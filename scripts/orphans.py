"""Call 'orphan' sense codons: used in genes, but no tRNA anticodon in the genome can read them (generous wobble rules)."""
import sys, json, collections
B='TCAG'; COMP={'A':'T','T':'A','G':'C','C':'G'}
CODE=dict(zip([a+b+c for a in B for b in B for c in B],'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'))
FOUR={k[:2] for k in CODE if len({CODE[k[:2]+x] for x in B})==1}
def reads(ac):
    """codons read by anticodon ac (5'->3', positions 34-36), generous rules"""
    if len(ac)!=3 or any(c not in COMP for c in ac): return set()
    p1,p2=COMP[ac[2]],COMP[ac[1]]; w=ac[0]; box=p1+p2
    if w=='C': th={'G'}
    elif w=='G': th={'C','T'}
    elif w=='A': th={'T','C','A'}          # inosine
    else: th={'A','G','T','C'} if box in FOUR else {'A','G'}   # U34: superwobble in 4-fold boxes
    r={box+t for t in th}
    if ac=='CAT': r.add('ATA')               # lysidine/agmatidine Ile2 (can't tell from anticodon)
    return r
def load(fn, taxfile):
    tax={l.split('\t')[0]:l.split('\t') for l in open(taxfile)}
    for l in open(fn):
        d=json.loads(l)
        if 'codons' not in d: continue
        d['meta']=tax.get(d['acc']); yield d
def orphans(d, min_count=30):
    readable=set()
    for k in d['trnas']: readable|=reads(k.split(':')[1])
    tot=sum(d['codons'].values()); out=[]
    for c,n in d['codons'].items():
        if c not in CODE or CODE[c]=='*' or c in ('ATG','ATA'): continue
        if c not in readable and n>=min_count: out.append((c,n,1000*n/tot))
    return out
if __name__=='__main__':
    fn,taxfile=sys.argv[1],sys.argv[2]
    n=0; byc=collections.Counter(); rows=[]
    for d in load(fn,taxfile):
        if d['table']!=11: pass
        n+=1; o=orphans(d)
        for c,k,pm in o:
            byc[c]+=1; m=d['meta']
            rows.append((c,k,round(pm,2),d['acc'],m[1],m[2],m[3],m[6],sum(d['trnas'].values()),d['table']))
    print('genomes',n); print(byc.most_common(30))
    with open(sys.argv[3],'w') as w:
        for r in sorted(rows,key=lambda r:(r[0],r[4])): w.write('\t'.join(map(str,r))+'\n')
