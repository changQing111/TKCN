#!/usr/bin/env python3
"""Intersect canonical k-mer sets across individual copy FASTAs, without KMC."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import resource
import time


def read_copy(path):
    records, parts, header = [], [], None
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if header is not None:
                    records.append(''.join(parts).upper())
                header, parts = line[1:], []
                if not header:
                    raise ValueError('Empty header: ' + str(path))
            else:
                if header is None:
                    raise ValueError('Sequence before header: ' + str(path))
                parts.append(line)
    if header is not None:
        records.append(''.join(parts).upper())
    if len(records) != 1 or not records[0]:
        raise ValueError('Expected exactly one nonempty copy per FASTA: ' + str(path))
    return records[0]


def kmers(seq, k):
    complement = str.maketrans('ACGT', 'TGCA')
    result, skipped = set(), 0
    for i in range(len(seq)-k+1):
        word = seq[i:i+k]
        if set(word) - set('ACGT'):
            skipped += 1
        else:
            result.add(min(word, word.translate(complement)[::-1]))
    return result, skipped


def run(args):
    start = time.perf_counter()
    paths = sorted(Path(args.input_dir).glob(args.pattern))
    if not paths or args.k < 1:
        raise ValueError('No input FASTAs or invalid k')
    out = Path(args.out_dir)
    if out.exists():
        raise ValueError('Output directory already exists')
    shared, stats, provenance = None, [], []
    # All inputs are validated even if the running intersection becomes empty.
    for path in paths:
        seq = read_copy(path)
        if len(seq) < args.k:
            raise ValueError('Copy shorter than k: ' + str(path))
        current, skipped = kmers(seq, args.k)
        shared = current if shared is None else shared.intersection(current)
        stats.append((path.name, len(seq), len(current), skipped, len(shared)))
        provenance.append(dict(path=str(path.resolve()), sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
    compute_seconds = time.perf_counter()-start
    out.mkdir(parents=True)
    prefix = args.target_id + '_shared.k' + str(args.k)
    with (out/(prefix+'.tsv')).open('w') as tsv, (out/(prefix+'.fa')).open('w') as fasta:
        for i, word in enumerate(sorted(shared), 1):
            tsv.write(word+'\t1\n')
            fasta.write('>{}_k{}\n{}\n'.format(args.target_id, i, word))
    with (out/'per_copy.tsv').open('w') as h:
        w=csv.writer(h,delimiter='\t',lineterminator='\n')
        w.writerow(['copy_fasta','length','unique_kmers','ambiguous_windows_skipped','running_intersection'])
        w.writerows(stats)
    report=dict(target_id=args.target_id,k=args.k,copy_count=len(paths),shared_kmers=len(shared),
                input_bytes=sum(p.stat().st_size for p in paths),compute_seconds=compute_seconds,
                inputs=provenance,output_count_policy='presence only; count column set to 1')
    report['total_seconds']=time.perf_counter()-start
    report['peak_rss_kib_linux']=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='inputs'},indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input-dir',required=True)
    p.add_argument('--pattern',default='*.fasta',help='Nonrecursive glob; one record/copy per file')
    p.add_argument('--out-dir',required=True)
    p.add_argument('--target-id',required=True)
    p.add_argument('-k',type=int,default=31)
    run(p.parse_args())
