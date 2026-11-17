"""For aligned in-frame codons, tabulate subject (homolog) residue per codon. Like Codetta, from DIAMOND alignments."""
import sys, collections
sys.path.insert(0, __file__.rsplit('/',1)[0])
from fsparse import TOK, rc, CODE
def table(fn):
    t=collections.defaultdict(collections.Counter)
    for l in open(fn):
        f=l.rstrip('\n').split('\t'); q=f[13]; bt=f[14]; s=f[12]
        if int(f[9])<0: q=rc(q)
        i=0
        for num,a,b in TOK.findall(bt):
            if num:
                for k in range(int(num)):
                    c=q[i:i+3]; t[c][CODE.get(c,'X')]+=1; i+=3
                continue
            if a=='-': continue
            if a=='\\': i+=1; continue
            if a=='/': i-=1; continue
            c=q[i:i+3]
            if b!='-': t[c][b]+=1
            i+=3
    return t
if __name__=='__main__':
    t=table(sys.argv[1])
    for c in sys.argv[2:]:
        n=sum(t[c].values()); print(c, n, [(a,round(v/n,2)) for a,v in t[c].most_common(5)] if n else '')
