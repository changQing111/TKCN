#!/usr/bin/env python3
"""Merge a single-copy membership TSV and a target panel TSV for joint counting."""
import argparse
from collections import defaultdict, Counter
import csv
import hashlib
import json
from pathlib import Path


def read_table(path, required):
    with open(path) as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        fields = reader.fieldnames
        if not set(required).issubset(fields or []):
            raise ValueError('Missing columns in {}: {}'.format(path, required))
        rows = list(reader)
    if not rows:
        raise ValueError('Empty panel: ' + str(path))
    return fields, rows


def canonical(seq):
    seq = seq.upper()
    if len(seq) != 31 or set(seq) - set('ACGT'):
        raise ValueError('Expected an ACGT 31-mer: ' + seq)
    return min(seq, seq.translate(str.maketrans('ACGT', 'TGCA'))[::-1])


def merge(single_path, target_path, out):
    single_fields, singles = read_table(single_path, ['panel_id', 'kmer', 'gene_id', 'role'])
    target_fields, targets = read_table(target_path, ['qname', 'sequence', 'gene_id'])
    if any(r['role'] != 'single_copy' for r in singles):
        raise ValueError('Single-copy input contains other roles')
    if {r['gene_id'] for r in singles} & {r['gene_id'] for r in targets}:
        raise ValueError('Gene identifiers collide between roles')
    members = defaultdict(set)
    gene_seqs = defaultdict(set)
    totals = Counter()
    for role, rows, col in [('single_copy', singles, 'kmer'), ('target', targets, 'sequence')]:
        for row in rows:
            if not row['gene_id'] or any(c.isspace() for c in row['gene_id']):
                raise ValueError('Invalid gene_id')
            seq = canonical(row[col])
            key = (row['gene_id'], role)
            members[seq].add(key)
            gene_seqs[key].add(seq)
            totals[key] += 1
    sequences = sorted(members)
    ids = {s: 'joint_{:07d}'.format(i) for i, s in enumerate(sequences, 1)}
    out = Path(out)
    out.mkdir(parents=True, exist_ok=False)

    def write(name, fields, rows):
        with (out/name).open('w') as h:
            w = csv.writer(h, delimiter='\t', lineterminator='\n')
            w.writerow(fields)
            w.writerows(rows)

    with (out/'all_panels.k31.fa').open('w') as h:
        for s in sequences:
            h.write('>{}\n{}\n'.format(ids[s], s))
    write('membership.tsv', ['panel_id', 'kmer', 'gene_id', 'role'],
          ((ids[s], s, g, r) for s in sequences for g, r in sorted(members[s])))
    write('gene_summary.tsv', ['gene_id', 'role', 'input_rows', 'unique_canonical_kmers', 'duplicate_rows', 'cross_gene_kmers'],
          ((g, r, totals[g,r], len(seqs), totals[g,r]-len(seqs), sum(len(members[s])>1 for s in seqs)) for (g,r), seqs in sorted(gene_seqs.items())))
    write('cross_gene_overlap.tsv', ['panel_id', 'kmer', 'gene_count', 'genes', 'roles'],
          ((ids[s], s, len(members[s]), ','.join(g for g,r in sorted(members[s])), ','.join(sorted({r for g,r in members[s]}))) for s in sequences if len(members[s])>1))
    # Preserve every original annotation and source identifier, including duplicate rows.
    for name, fields, rows, col in [('target_source_annotations.tsv', target_fields, targets, 'sequence'), ('single_copy_source_annotations.tsv', single_fields, singles, 'kmer')]:
        write(name, ['joint_panel_id', 'canonical_kmer'] + fields,
              ([ids[canonical(row[col])], canonical(row[col])] + [row[f] for f in fields] for row in rows))
    info = dict(k=31, unique_kmers=len(sequences), membership_rows=sum(len(v) for v in members.values()),
                single_copy_genes=sum(r=='single_copy' for g,r in gene_seqs), target_genes=sum(r=='target' for g,r in gene_seqs),
                single_copy_input_rows=len(singles), target_input_rows=len(targets),
                single_copy_unique_kmers=len({canonical(r['kmer']) for r in singles}),
                target_unique_kmers=len({canonical(r['sequence']) for r in targets}),
                cross_gene_kmers=sum(len(v)>1 for v in members.values()),
                cross_role_kmers=sum(len({r for g,r in v})>1 for v in members.values()),
                inputs=[dict(path=str(Path(p).resolve()), sha256=hashlib.sha256(Path(p).read_bytes()).hexdigest()) for p in (single_path, target_path)],
                policy='Canonicalize and deduplicate sequences; preserve all memberships and annotations. No additional quality filtering.')
    (out/'provenance.json').write_text(json.dumps(info, indent=2)+'\n')
    print(json.dumps(info, indent=2))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--single-membership', required=True)
    p.add_argument('--target-panel', required=True)
    p.add_argument('--out-dir', required=True)
    a = p.parse_args()
    merge(a.single_membership, a.target_panel, a.out_dir)
