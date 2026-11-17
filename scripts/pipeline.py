"""One shard of the orphan-codon screen.

For each genome in the shard:
  1. download it from NCBI,
  2. call genes (pyrodigal) and tRNAs (pyaragorn), count codons and anticodons,
  3. if some sense codon is used but no tRNA in the genome can read it ("orphan"),
     align the genome to Swiss-Prot with DIAMOND in frameshift mode and record
     where the +1/-1 frameshifts fall and what the homologs have at the orphan codon.

usage: pipeline.py SAMPLE.tsv SHARD NSHARDS OUTDIR DIAMOND_DB
Set VALIDATE_ALL=1 to align every genome (not only those with an orphan codon);
CGG is then always reported in 'orph', whether or not a tRNA reads it.
"""
import sys, os, re, json, gzip, time, tempfile, subprocess, collections, urllib.request
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pyrodigal, pyaragorn
from orphans import orphans
from fscount import sites
from codonaa import table

RRN = {'ATT', 'ATC', 'GCA', 'GCG', 'GCT', 'GCC'}  # Ile-GAU / Ala-UGC sit in rRNA operons, often lost from MAGs
SAMPLE, SHARD, NSHARDS, OUT, DB = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4], sys.argv[5]
TMP = tempfile.mkdtemp()
VALIDATE_ALL = os.environ.get('VALIDATE_ALL') == '1'


def url_for(acc, name):
    acc = acc.split('_', 1)[1]
    n = re.sub(r'[^A-Za-z0-9._-]', '_', name)
    p = f"{acc[:3]}/{acc[4:7]}/{acc[7:10]}/{acc[10:13]}"
    return f"https://ftp.ncbi.nlm.nih.gov/genomes/all/{p}/{acc}_{n}/{acc}_{n}_genomic.fna.gz"


def fetch(row):
    acc, name = row[0], row[4]
    out = os.path.join(TMP, acc + '.fna.gz')
    u = url_for(acc, name)
    for t in range(6):
        if t == 1:  # name mangling differs: look the directory up
            try:
                a = acc.split('_', 1)[1]
                base = u.rsplit('/', 2)[0] + '/'
                html = urllib.request.urlopen(base, timeout=60).read().decode()
                m = re.search(r'href="(' + re.escape(a) + r'_[^"/]+)/"', html)
                if not m:
                    return None
                u = base + m.group(1) + '/' + m.group(1) + '_genomic.fna.gz'
            except Exception:
                time.sleep(3)
                continue
        try:
            d = urllib.request.urlopen(u, timeout=120).read()
            open(out, 'wb').write(d)
            return out
        except Exception as e:
            if '404' in str(e) and t >= 1:
                return None
            time.sleep(2 + 3 * t)
    return None


def readfa(fn):
    seqs, cur = [], []
    for line in gzip.open(fn, 'rt'):
        if line.startswith('>'):
            if cur:
                seqs.append(''.join(cur))
            cur = []
        else:
            cur.append(line.strip().upper())
    if cur:
        seqs.append(''.join(cur))
    return seqs


def scan(acc, fn):
    seqs = readfa(fn)
    gf = pyrodigal.GeneFinder(meta=False, mask=True)
    gf.train(*[s.encode() for s in seqs])
    cod = collections.Counter()
    ngenes = 0
    for s in seqs:
        for g in gf.find_genes(s.encode()):
            nt = g.sequence()
            ngenes += 1
            cod.update(nt[i:i + 3] for i in range(3, len(nt) - 3, 3))
    rf = pyaragorn.RNAFinder(translation_table=11, tmrna=False, trna=True, linear=False)
    trnas = collections.Counter()
    for s in seqs:
        for g in rf.find_rna(s):
            if isinstance(g, pyaragorn.TRNAGene):
                trnas[f"{g.amino_acid}:{g.anticodon.upper()}"] += 1
    return dict(acc=acc, ngenes=ngenes, table=gf.training_info.translation_table, codons=dict(cod), trnas=dict(trnas))


def validate(d, fn, sig):
    out = fn + '.dmnd.tsv'
    subprocess.run(['diamond', 'blastx', '-q', fn, '-d', DB, '-F', '15', '--fast', '--max-target-seqs', '1',
                    '--range-culling', '-e', '1e-5', '--threads', '1', '--query-gencode', str(d['table']),
                    '-f', '6', 'qseqid', 'sseqid', 'pident', 'qstart', 'qend', 'sstart', 'send', 'evalue',
                    'bitscore', 'qframe', 'qseq_translated', 'full_sseq', 'sseq', 'qseq', 'btop',
                    '-o', out, '--quiet'], check=True)
    s = sites(out)
    plus = [x for x in s if x['kind'] == '+1']
    minus = [x for x in s if x['kind'] == '-1']
    t = table(out)
    n_align = sum(1 for _ in open(out))
    os.remove(out)
    res = dict(acc=d['acc'], n_align=n_align, n_plus=len(plus), n_minus=len(minus),
               plus_codons=dict(collections.Counter(x['c0'] for x in plus)),
               minus_codons=dict(collections.Counter(x['c0'] for x in minus)),
               plus_gly=dict(collections.Counter(x['c0'] for x in plus if x['snext'] == 'G')),
               orph={})
    for cod, n, pm in sig:
        res['orph'][cod] = dict(uses=n, pm=round(pm, 3), fs=sum(1 for x in plus if x['c0'] == cod),
                                gly=sum(1 for x in plus if x['c0'] == cod and x['snext'] == 'G'),
                                aa=t[cod].most_common(4))
    return res


def work(args):
    acc, fn = args
    try:
        d = scan(acc, fn)
        v = None
        if VALIDATE_ALL:
            tot = sum(d['codons'].values())
            sig = [x for x in orphans(d, min_count=0) if x[0] not in RRN and x[0] != 'CGG']
            n = d['codons'].get('CGG', 0)
            v = validate(d, fn, sig + [('CGG', n, 1000 * n / max(1, tot))])
        elif sum(d['trnas'].values()) >= 25:
            sig = [x for x in orphans(d, min_count=20) if x[0] not in RRN and x[2] <= 5]
            if sig and len(sig) <= 4:
                v = validate(d, fn, sig)
        return json.dumps(d), (json.dumps(v) if v else None)
    except Exception as e:
        return json.dumps({'acc': acc, 'error': str(e)[:300]}), None
    finally:
        try:
            os.remove(fn)
        except OSError:
            pass


if __name__ == '__main__':
    rows = [l.rstrip('\n').split('\t') for i, l in enumerate(open(SAMPLE)) if i % NSHARDS == SHARD]
    os.makedirs(OUT, exist_ok=True)
    so = open(f'{OUT}/scan_{SHARD:03d}.jsonl', 'w')
    vo = open(f'{OUT}/val_{SHARD:03d}.jsonl', 'w')
    ncpu = os.cpu_count() or 4
    t0 = time.time()
    with ThreadPoolExecutor(8) as tp, ProcessPoolExecutor(ncpu) as pp:
        futs = []
        def drain(k):
            while len(futs) > k:
                s, v = futs.pop(0).result()
                so.write(s + '\n'); so.flush()
                if v:
                    vo.write(v + '\n'); vo.flush()
        for r, fn in zip(rows, tp.map(fetch, rows)):
            if fn is None:
                so.write(json.dumps({'acc': r[0], 'error': 'download'}) + '\n')
                continue
            futs.append(pp.submit(work, (r[0], fn)))
            drain(3 * ncpu)
        drain(0)
    print(f'shard {SHARD}: {len(rows)} genomes in {time.time() - t0:.0f}s')
