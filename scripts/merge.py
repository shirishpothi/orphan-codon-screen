"""Merge shard outputs and summarise.

usage: merge.py SAMPLE.tsv RESULTS_DIR OUTDIR
Writes scan.jsonl, val.jsonl, hits.tsv, hotspots.tsv, lineage.tsv and summary.md to OUTDIR.
"""
import sys, os, glob, json, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from orphans import orphans

SAMPLE, RES, OUT = sys.argv[1:4]
os.makedirs(OUT, exist_ok=True)
RRN = {'ATT', 'ATC', 'GCA', 'GCG', 'GCT', 'GCC'}
tax = {l.split('\t')[0]: l.rstrip('\n').split('\t') for l in open(SAMPLE)}


def cat(pattern, dest):
    rows = []
    with open(dest, 'w') as w:
        for fn in sorted(glob.glob(os.path.join(RES, '**', pattern), recursive=True)):
            for l in open(fn):
                w.write(l)
                rows.append(json.loads(l))
    return rows


scans = cat('scan_*.jsonl', f'{OUT}/scan.jsonl')
vals = cat('val_*.jsonl', f'{OUT}/val.jsonl')
ok = [d for d in scans if 'codons' in d]
err = collections.Counter(d.get('error', '')[:40] for d in scans if 'codons' not in d)


def short(acc, lv=(1, 6)):
    return ';'.join(tax[acc][1].split(';')[lv[0]:lv[1]])


# 1. orphan codon is a +1 frameshift hotspot (the TANB77/CGG signature)
hits = []
for v in vals:
    for c, o in v['orph'].items():
        if o['fs'] >= 5 and o['fs'] >= 0.2 * max(1, v['n_plus']):
            hits.append((c, v['acc'], short(v['acc']), o['uses'], o['pm'], o['fs'], v['n_plus'], o['gly'], o['aa']))
hits.sort(key=lambda h: (h[0], h[2]))
with open(f'{OUT}/hits.tsv', 'w') as w:
    w.write('codon\tacc\ttaxonomy\tuses\tper_mille\tfs_at_codon\tn_plus\tgly_at_codon\thomolog_aa_at_inframe_codon\n')
    for h in hits:
        w.write('\t'.join(map(str, h)) + '\n')

# 2. any codon (orphan or not) that concentrates frameshifts in either direction
hot = []
for v in vals:
    for kind, key in (('+1', 'plus_codons'), ('-1', 'minus_codons')):
        n = v['n_plus'] if kind == '+1' else v['n_minus']
        for c, k in v[key].items():
            if k >= 15 and k >= 0.3 * n:
                hot.append((kind, c, v['acc'], short(v['acc']), k, n, c in v['orph']))
with open(f'{OUT}/hotspots.tsv', 'w') as w:
    w.write('direction\tcodon\tacc\ttaxonomy\tat_codon\ttotal\tis_orphan\n')
    for h in sorted(hot):
        w.write('\t'.join(map(str, h)) + '\n')

# 3. lineage-consistent orphan codons (order level)
tot, orph, use = collections.Counter(), collections.Counter(), collections.defaultdict(list)
for d in ok:
    if sum(d['trnas'].values()) < 25 or d['acc'] not in tax:
        continue
    g = ';'.join(tax[d['acc']][1].split(';')[:4])
    tot[g] += 1
    for c, n, pm in orphans(d, min_count=10):
        if c not in RRN and pm <= 5:
            orph[(g, c)] += 1
            use[(g, c)].append(pm)
lin = sorted(((c, v, tot[g], g) for (g, c), v in orph.items() if v >= 3 and v / tot[g] >= 0.5), key=lambda r: -r[1])
with open(f'{OUT}/lineage.tsv', 'w') as w:
    w.write('codon\tgenera_orphan\tgenera_scanned\torder\tmedian_per_mille\n')
    for c, v, t, g in lin:
        u = sorted(use[(g, c)])
        w.write(f'{c}\t{v}\t{t}\t{g}\t{u[len(u) // 2]:.2f}\n')

# summary
hit_orders = collections.Counter((h[0], ';'.join(h[2].split(';')[:3])) for h in hits)
with open(f'{OUT}/summary.md', 'w') as w:
    w.write(f'# Orphan-codon screen\n\nGenomes scanned: {len(ok)} of {len(scans)} (errors: {dict(err)})\n\n')
    w.write(f'Genomes with an orphan codon sent to frameshift alignment: {len(vals)}\n\n')
    w.write(f'## Orphan codon = +1 frameshift hotspot ({len(hits)} genomes)\n\n| codon | phylum;class;order | genomes |\n|---|---|---|\n')
    for (c, o), n in hit_orders.most_common():
        w.write(f'| {c} | {o} | {n} |\n')
    w.write(f'\n## Any-codon frameshift hotspots ({len(hot)})\n\n')
    for h in sorted(hot)[:60]:
        w.write(f'- {h[0]} {h[1]} {h[3]}: {h[4]}/{h[5]} (orphan={h[6]})\n')
    w.write(f'\n## Lineage-consistent orphan codons ({len(lin)})\n\n')
    for c, v, t, g in lin[:60]:
        w.write(f'- {c} {v}/{t} genera {g}\n')
print(open(f'{OUT}/summary.md').read())
