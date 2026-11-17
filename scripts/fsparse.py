import sys, re, collections
CODE=dict(zip([a+b+c for a in 'TCAG' for b in 'TCAG' for c in 'TCAG'],'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'))
def tr(s): return ''.join(CODE.get(s[i:i+3],'X') for i in range(0,len(s)-2,3))
def rc(s): return s[::-1].translate(str.maketrans('ACGTN','TGCAN'))
TOK=re.compile(r'(\d+)|([A-Z*\\/-])([A-Z*\\/-])')
def events(qseq, btop):
    """walk alignment; yield (kind, nt_index_in_qseq, subject_aa_after) for each frameshift."""
    i=0; out=[]; toks=TOK.findall(btop)
    for k,(num,a,b) in enumerate(toks):
        if num: i+=3*int(num); continue
        if a=='-': continue           # gap in query
        if a=='\\': out.append(('+1',i)); i+=1; continue
        if a=='/': out.append(('-1',i)); i+=2; continue
        i+=3
    return out,i
if __name__=='__main__':
    n=collections.Counter(); bad=0
    for l in open(sys.argv[1]):
        f=l.rstrip('\n').split('\t'); q=f[13]; bt=f[14]
        ev,i=events(q,bt)
        if i!=len(q): bad+=1
    print('length mismatches',bad)
    f=open(sys.argv[1]).readline().split('\t'); q=f[13]
    print(tr(q)[:30], tr(rc(q))[:30], f[10][:30])
