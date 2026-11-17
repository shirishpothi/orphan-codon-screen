"""Per-genome CGG decoding vs CGG skipping, for a VALIDATE_ALL run.

usage: tanb77.py SAMPLE.tsv MERGED_DIR
Writes MERGED_DIR/cgg_by_genome.tsv and prints totals by family and by CGG-tRNA status.
"""
import sys, os, json, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from orphans import reads

SAMPLE, M = sys.argv[1:3]
tax = {l.split('\t')[0]: l.split('\t')[1] for l in open(SAMPLE)}
scan = {d['acc']: d for d in map(json.loads, open(f'{M}/scan.jsonl')) if 'codons' in d}
rows = []
for v in map(json.loads, open(f'{M}/val.jsonl')):
    d = scan[v['acc']]
    readers = [k for k in d['trnas'] if 'CGG' in reads(k.split(':')[1])]
    o = v['orph']['CGG']
    t = tax[v['acc']].split(';')
    rows.append(dict(acc=v['acc'], family=t[4][3:], genus=t[5][3:], ntrna=sum(d['trnas'].values()),
                     cgg_reader=','.join(sorted(readers)) or '-', cgg_uses=o['uses'], cgg_pm=o['pm'],
                     n_plus=v['n_plus'], cgg_fs=o['fs'], cgg_gly=o['gly'],
                     other_fs=v['n_plus'] - o['fs'], other_gly=sum(v['plus_gly'].values()) - o['gly']))
cols = list(rows[0])
with open(f'{M}/cgg_by_genome.tsv', 'w') as w:
    w.write('\t'.join(cols) + '\n')
    for r in sorted(rows, key=lambda r: (r['family'], r['genus'], r['acc'])):
        w.write('\t'.join(str(r[c]) for c in cols) + '\n')


def agg(key):
    a = collections.defaultdict(lambda: collections.Counter())
    for r in rows:
        if r['ntrna'] < 25:
            continue
        k = key(r); a[k]['genomes'] += 1
        for f in ('cgg_fs', 'n_plus', 'cgg_gly', 'other_fs', 'other_gly'):
            a[k][f] += r[f]
        a[k]['skippers'] += r['cgg_fs'] >= 3 and r['cgg_fs'] >= 0.2 * r['n_plus']
    print(f"{'group':40} genomes skippers  CGG_fs/+1_fs  Gly@CGG  Gly@other")
    for k, c in sorted(a.items()):
        print(f"{str(k):40} {c['genomes']:7} {c['skippers']:8}  {c['cgg_fs']:5}/{c['n_plus']:<6}  "
              f"{c['cgg_gly']:4}/{c['cgg_fs']:<4} {c['other_gly']:4}/{c['other_fs']}")
    print()


print(f'{len(rows)} genomes aligned; totals use genomes with >= 25 tRNAs\n')
agg(lambda r: ('has CGG reader' if r['cgg_reader'] != '-' else 'no CGG reader'))
agg(lambda r: (r['family'], r['cgg_reader'] != '-'))
