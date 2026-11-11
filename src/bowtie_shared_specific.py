#!/usr/bin/env python3
"""Build shared-specific panels with Bowtie1 all-hit exact mapping (stdlib only)."""
import argparse
import csv
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
from collections import defaultdict


def table(path, required):
    with open(path) as handle:
        reader = csv.DictReader((s for s in handle if s.strip() and not s.startswith('#')), delimiter='\t')
        if not set(required).issubset(reader.fieldnames or []):
            raise ValueError('{} requires tab-separated columns: {}'.format(path, ', '.join(required)))
        rows = list(reader)
    if not rows:
        raise ValueError('Empty table: {}'.format(path))
    return rows


def identifier(value):
    if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*', value or ''):
        raise ValueError('Invalid identifier: {!r}'.format(value))
    return value


def resolve(manifest, value):
    path = Path(value)
    return path if path.is_absolute() else Path(manifest).resolve().parent / path


def candidates(manifest):
    records, targets, lengths = {}, set(), set()
    for row in table(manifest, ['target_id', 'shared_tsv']):
        target = identifier(row['target_id'])
        if target in targets:
            raise ValueError('Duplicate target: ' + target)
        targets.add(target)
        seen = set()
        with open(resolve(manifest, row['shared_tsv'])) as handle:
            for line in handle:
                if not line.strip() or line.startswith('#'):
                    continue
                fields = line.split()
                if fields[0].lower() == 'kmer':
                    continue
                seq = fields[0].upper()
                if len(fields) != 2 or not re.fullmatch('[ACGT]+', seq) or int(fields[1]) < 1:
                    raise ValueError('Expected KMC dump: ACGT kmer and positive integer count')
                canonical = min(seq, seq.translate(str.maketrans('ACGT', 'TGCA'))[::-1])
                if canonical in seen:
                    raise ValueError('Duplicate/reverse-complement kmer in ' + target)
                seen.add(canonical)
                lengths.add(len(seq))
                records['{}|kmer_{:06d}'.format(target, len(seen))] = (target, seq, fields[1])
    if len(lengths) > 1:
        raise ValueError('All candidates must have the same k-mer length')
    return records, targets


def load_regions(path, targets, genomes):
    regions = defaultdict(list)
    for row in table(path, ['target_id', 'genome_id', 'reference_name', 'start', 'end', 'coordinate_system']):
        target, genome = row['target_id'], row['genome_id']
        if target not in targets or genome not in genomes:
            raise ValueError('Unknown target/genome in region table: {}/{}'.format(target, genome))
        start, end = int(row['start']), int(row['end'])
        system = row['coordinate_system']
        if system == '1-based-closed':
            start -= 1
        elif system != '0-based-half-open':
            raise ValueError('Unsupported coordinate system: ' + system)
        if start < 0 or end <= start or not row['reference_name']:
            raise ValueError('Invalid region: ' + str(row))
        regions[genome, target, row['reference_name']].append((start, end))
    for genome in genomes:
        for target in targets:
            if not any(g == genome and t == target for g, t, _ in regions):
                raise ValueError('Missing target regions: {}/{}'.format(genome, target))
    return regions


def classify_sam(path, genome, records, regions, hits_path):
    counts = {q: [0, 0] for q in records}
    seen, unmapped, contigs = set(), set(), {}
    with open(path) as handle, open(hits_path, 'w') as output:
        writer = csv.writer(output, delimiter='\t', lineterminator='\n')
        writer.writerow(['target_id', 'kmer_id', 'kmer', 'genome_id', 'reference_name', 'start_0', 'end_0', 'strand', 'hit_class'])
        for line in handle:
            if line.startswith('@'):
                if line.startswith('@SQ\t'):
                    tags = dict(f.split(':', 1) for f in line.rstrip().split('\t')[1:])
                    contigs[tags['SN']] = int(tags['LN'])
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 11 or f[0] not in records:
                raise ValueError('Malformed SAM or unknown query: ' + line[:150])
            q, flag = f[0], int(f[1])
            seen.add(q)
            if flag & 4:
                unmapped.add(q)
                continue
            target, seq, _ = records[q]
            if f[5] != '{}M'.format(len(seq)) or 'NM:i:0' not in f[11:]:
                raise ValueError('Expected full-length exact Bowtie1 alignment: ' + q)
            start, end = int(f[3]) - 1, int(f[3]) - 1 + len(seq)
            if f[2] not in contigs or start < 0 or end > contigs[f[2]]:
                raise ValueError('Alignment outside reference bounds: ' + q)
            inside = any(a <= start and end <= b for a, b in regions.get((genome, target, f[2]), []))
            counts[q][0 if inside else 1] += 1
            writer.writerow([target, q, seq, genome, f[2], start, end, '-' if flag & 16 else '+', 'target' if inside else 'non_target'])
    if seen != set(records) or any(sum(counts[q]) for q in unmapped):
        raise ValueError('Incomplete or inconsistent SAM query records')
    for (g, _, chrom), intervals in regions.items():
        if g == genome and (chrom not in contigs or any(b > contigs[chrom] for _, b in intervals)):
            raise ValueError('Target region absent from/outside SAM reference: ' + chrom)
    return counts


def index_exists(prefix):
    return any(all(Path(str(prefix) + suffix + ext).is_file() for suffix in ('.1', '.2', '.3', '.4', '.rev.1', '.rev.2')) for ext in ('.ebwt', '.ebwtl'))


def run(args):
    records, targets = candidates(args.targets)
    refs = table(args.references, ['genome_id', 'fasta', 'index_prefix'])
    genomes = [identifier(r['genome_id']) for r in refs]
    if len(genomes) != len(set(genomes)):
        raise ValueError('Duplicate genome_id')
    regions = load_regions(args.regions, targets, set(genomes))
    if args.threads < 1:
        raise ValueError('threads must be positive')
    if not shutil.which(args.bowtie):
        raise ValueError('Bowtie1 not found; set --bowtie or BOWTIE_BIN')
    version = subprocess.check_output([args.bowtie, '--version'], universal_newlines=True)
    if not re.search(r'\bversion\s+1\.', version):
        raise ValueError('Bowtie1 is required (not Bowtie2)')
    out = Path(args.out_dir).resolve()
    out.mkdir(parents=True, exist_ok=False)
    query = out / 'all_targets.shared.fa'
    with query.open('w') as handle:
        for q, (_, seq, _) in records.items():
            handle.write('>{}\n{}\n'.format(q, seq))
    retained = set(records)
    with (out / 'classification.tsv').open('w') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        writer.writerow(['target_id', 'kmer_id', 'kmer', 'genome_id', 'target_hits', 'non_target_hits', 'status', 'retained'])
        for ref in refs:
            genome = ref['genome_id']
            prefix = resolve(args.references, ref['index_prefix']) if ref['index_prefix'] else out / 'index' / genome
            if not index_exists(prefix):
                fasta = resolve(args.references, ref['fasta'])
                if not ref['fasta'] or not fasta.is_file():
                    raise ValueError('Missing index and reference FASTA for ' + genome)
                prefix.parent.mkdir(parents=True, exist_ok=True)
                with (out / (genome + '.build.log')).open('w') as log:
                    subprocess.run([args.bowtie_build, str(fasta), str(prefix)], stdout=log, stderr=log, check=True)
                if not index_exists(prefix):
                    raise ValueError('Incomplete Bowtie1 index: ' + str(prefix))
            if not records:
                continue
            sam = out / (genome + '.sam')
            with (out / (genome + '.bowtie.log')).open('w') as log:
                subprocess.run([args.bowtie, '-f', '-S', '-v', '0', '-a', '-p', str(args.threads), str(prefix), str(query), str(sam)], stdout=log, stderr=log, check=True)
            counts = classify_sam(sam, genome, records, regions, out / (genome + '.hits.tsv'))
            for q, (on, off) in counts.items():
                status = 'mixed' if on and off else 'non_target_only' if off else 'target_only' if on else 'unmapped'
                keep = bool(on and not off)
                target, seq, _ = records[q]
                writer.writerow([target, q, seq, genome, on, off, status, int(keep)])
                if not keep:
                    retained.discard(q)
            if not args.keep_sam:
                sam.unlink()
    with (out / 'summary.tsv').open('w') as summary:
        summary.write('target_id\tshared_kmers\tshared_specific_kmers\n')
        for target in sorted(targets):
            total = sum(t == target for t, _, _ in records.values())
            kept = [q for q in records if records[q][0] == target and q in retained]
            with (out / (target + '.shared-specific.tsv')).open('w') as tsv, (out / (target + '.shared-specific.fa')).open('w') as fasta:
                for q in kept:
                    _, seq, count = records[q]
                    tsv.write('{}\t{}\n'.format(seq, count))
                    fasta.write('>{}\n{}\n'.format(q, seq))
            summary.write('{}\t{}\t{}\n'.format(target, total, len(kept)))
    (out / 'SUCCESS').write_text('Completed all references.\n')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--targets', required=True, help='TSV: target_id, shared_tsv (KMC dump)')
    p.add_argument('--references', required=True, help='TSV: genome_id, fasta, index_prefix; paths relative to this table; blank index builds locally')
    p.add_argument('--regions', required=True, help='TSV: target_id, genome_id, reference_name, start, end, coordinate_system')
    p.add_argument('--out-dir', required=True, help='New output directory; existing directories are never overwritten')
    p.add_argument('-t', '--threads', type=int, default=8)
    p.add_argument('--bowtie', default=os.environ.get('BOWTIE_BIN', 'bowtie'))
    p.add_argument('--bowtie-build', default=os.environ.get('BOWTIE_BUILD_BIN'), help='Defaults to bowtie-build alongside --bowtie, then PATH')
    p.add_argument('--keep-sam', action='store_true')
    try:
        args = p.parse_args()
        if not args.bowtie_build:
            executable = shutil.which(args.bowtie)
            sibling = Path(executable).parent / 'bowtie-build' if executable else None
            args.bowtie_build = str(sibling) if sibling and sibling.is_file() else 'bowtie-build'
        run(args)
    except (ValueError, OSError, subprocess.CalledProcessError) as exc:
        sys.exit('ERROR: {}'.format(exc))


if __name__ == '__main__':
    main()
