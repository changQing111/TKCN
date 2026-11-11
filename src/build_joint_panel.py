#!/usr/bin/env python3
"""Merge preselected k-mer FASTAs; retain every gene membership and source record."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import re


def fasta(path):
    name, parts = None, []
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if name is not None:
                    yield name, ''.join(parts).upper()
                name, parts = line[1:], []
                if not name:
                    raise ValueError('Empty FASTA header: ' + str(path))
            else:
                if name is None:
                    raise ValueError('Sequence before header: ' + str(path))
                parts.append(line)
    if name is not None:
        yield name, ''.join(parts).upper()


def build(manifest, output, k):
    manifest = Path(manifest).resolve()
    with manifest.open() as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    if not rows or not {'gene_id', 'role', 'fasta'}.issubset(rows[0]):
        raise ValueError('Manifest requires gene_id, role, fasta columns')
    genes, members, sources, stats, provenance = set(), {}, [], [], []
    for row in rows:
        gene, role = row['gene_id'], row['role']
        if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*', gene) or gene in genes:
            raise ValueError('Invalid or duplicate gene_id: ' + gene)
        if role not in ('single_copy', 'target'):
            raise ValueError('role must be single_copy or target')
        genes.add(gene)
        path = (manifest.parent / row['fasta']).resolve()
        unique, total, reversed_count = set(), 0, 0
        for header, seq in fasta(path):
            if len(seq) != k or not re.fullmatch('[ACGT]+', seq):
                raise ValueError('{}: {} must contain exactly {} ACGT bases'.format(path, header, k))
            rc = seq.translate(str.maketrans('ACGT', 'TGCA'))[::-1]
            canonical = min(seq, rc)
            reversed_count += canonical != seq
            total += 1
            unique.add(canonical)
            members.setdefault(canonical, set()).add((gene, role))
            sources.append((canonical, gene, header, '-' if canonical != seq else '+'))
        if not total:
            raise ValueError('Empty FASTA: ' + str(path))
        stats.append((gene, role, total, len(unique), total-len(unique), reversed_count))
        provenance.append(dict(gene_id=gene, role=role, fasta=str(path), sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
    output = Path(output)
    output.mkdir(parents=True, exist_ok=False)
    (output / 'by_gene').mkdir()
    ids = {seq: 'panel_{:06d}'.format(i) for i, seq in enumerate(sorted(members), 1)}
    def write_table(name, header, data):
        with (output / name).open('w') as handle:
            writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
            writer.writerow(header)
            writer.writerows(data)
    with (output / 'all_panels.k{}.fa'.format(k)).open('w') as handle:
        for seq in sorted(members):
            handle.write('>{}\n{}\n'.format(ids[seq], seq))
    for gene in sorted(genes):
        with (output / 'by_gene' / (gene + '.fa')).open('w') as handle:
            for seq in sorted(members):
                if any(g == gene for g, _ in members[seq]):
                    handle.write('>{}|{}\n{}\n'.format(gene, ids[seq], seq))
    write_table('membership.tsv', ['panel_id', 'kmer', 'gene_id', 'role'],
                ((ids[s], s, g, r) for s in sorted(members) for g, r in sorted(members[s])))
    write_table('source_records.tsv', ['panel_id', 'gene_id', 'original_header', 'canonical_orientation'],
                ((ids[s], g, h, strand) for s, g, h, strand in sources))
    write_table('summary.tsv', ['gene_id', 'role', 'input_records', 'unique_canonical_kmers', 'duplicate_records', 'reverse_complemented_records'], stats)
    write_table('cross_gene_overlap.tsv', ['panel_id', 'kmer', 'gene_ids'],
                ((ids[s], s, ','.join(sorted(g for g, _ in members[s]))) for s in sorted(members) if len(members[s]) > 1))
    metadata = dict(k=k, unique_kmers=len(members), cross_gene_shared_kmers=sum(len(v)>1 for v in members.values()), inputs=provenance,
                    policy='No additional sequence-quality filtering; canonicalize and deduplicate only. Use reads database as left input with intersect -ocleft.')
    (output / 'provenance.json').write_text(json.dumps(metadata, indent=2) + '\n')
    print(json.dumps({key: metadata[key] for key in ('k', 'unique_kmers', 'cross_gene_shared_kmers')}))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--manifest', required=True, help='TSV: gene_id, role, fasta; relative paths are relative to the manifest')
    p.add_argument('--out-dir', required=True, help='New output directory')
    p.add_argument('-k', type=int, default=31)
    args = p.parse_args()
    build(args.manifest, args.out_dir, args.k)
