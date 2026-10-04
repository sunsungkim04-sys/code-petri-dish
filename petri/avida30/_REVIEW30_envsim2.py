import sys
M=(1<<64)-1
class R:
    def __init__(s,seed): s.s=seed&M
    def nxt(s):
        s.s=(s.s+0x9E3779B97F4A7C15)&M; z=s.s
        z=((z^(z>>30))*0xBF58476D1CE4E5B9)&M
        z=((z^(z>>27))*0x94D049BB133111EB)&M
        return z^(z>>31)
    def u(s): return (s.nxt()>>11)*(1.0/9007199254740992.0)
    def below(s,n): return int(s.u()*n)
W=60;C=W*W;K=26
nb=[]
for c in range(C):
    x,y=c%W,c//W
    nb.append([((y-1)%W)*W+x, y*W+(x+1)%W, ((y+1)%W)*W+x, y*W+(x-1)%W])
for seed in range(1,9):
    base=(seed*2654435761+0x5EED)&M
    mix=R(base); mix.nxt(); mix.nxt(); env=R(mix.nxt())
    pool=[0]*(C*K)
    for t in range(C*200):
        c=env.below(C); k=env.below(K); pool[c*K+k]+=1
    for k in [ord(ch)-97 for ch in 'wzcagcccczvfcaxgab']:
        while True:
            c=env.below(C)
            if pool[c*K+k]>0: pool[c*K+k]-=1; break
    first=None; vals=[]
    for upd in range(0,150):
        if upd==0: env.below(1200)   # lifespan draw at first check
        d=env.u(); vals.append(d)
        if d<0.0003 and first is None: first=upd
        for q in range(int(C*0.15)):
            i=env.below(C); k=env.below(K)
            if pool[i*K+k]>0:
                t=nb[i][env.below(4)]; pool[i*K+k]-=1; pool[t*K+k]+=1
    print(seed,'first death-draw<3e-4 at update',first, 'min draw', min(vals))
