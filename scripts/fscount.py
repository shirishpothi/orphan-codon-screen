"""Per-genome frameshift site table from DIAMOND -F output: codon at the shifted base, plus subject aa there."""
import sys, re, collections, json
sys.path.insert(0, __file__.rsplit('/',1)[0])
from fsparse import TOK, rc, tr
def sites(fn):
    res=[]
    for l in open(fn):
        f=l.rstrip('\n').split('\t'); q=f[13]; bt=f[14]
        if int(f[9])<0: q=rc(q)
        i=0; s_aa=[]  # walk and keep subject aa near shift
        toks=TOK.findall(bt)
        for k,(num,a,b) in enumerate(toks):
            if num: i+=3*int(num); continue
            if a=='-': continue
            if a in '\\/':
                # subject residue aligned to the next query codon
                nxt=None
                for num2,a2,b2 in toks[k+1:]:
                    if num2: nxt='match'; break
                    if a2 not in '\\/': nxt=b2; break
                q_next=tr(q[i+1:i+4]) if a=='\\' else tr(q[i-1:i+2])
                res.append(dict(kind='+1' if a=='\\' else '-1', c0=q[i:i+3], cm1=q[i-1:i+2], cm2=q[i-2:i+1], ctx=q[max(0,i-6):i+7], qnext=q_next, snext=nxt if nxt!='match' else q_next))
                i+= 1 if a=='\\' else -1; continue
            i+=3
    return res
if __name__=='__main__':
    r=sites(sys.argv[1]); p=[x for x in r if x['kind']=='+1']
    print('total',len(r),'+1',len(p))
    for key in ('c0','cm1','cm2'):
        print(key, collections.Counter(x[key] for x in p).most_common(5))
    print('subject aa at +1', collections.Counter(x['snext'] for x in p).most_common(5))
