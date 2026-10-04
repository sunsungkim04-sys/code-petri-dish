# independent recount (reviewer) — own loop parser, own bootstrap RNG
import json, os, random, statistics as st
B = os.path.expanduser("~/petri")
def L_of(code):
    # loop = from after the last 's' before first 'l' through that 'l'; needs a 'c' inside; count non-'x'
    if 'l' not in code: return None
    li = code.find('l'); pre = code[:li]
    start = pre.rfind('s') + 1
    seg = code[start:li+1]
    if 'c' not in seg: return None
    return len(seg) - seg.count('x')
cells = {}
spec = {"DF": ("evo17", "", [0.001,0.003,0.005,0.01]), "DO": ("evo29/do", "_do", [0.001,0.003,0.005,0.01,0.03,0.05]), "FF": ("evo29/ff", "_ff", [0.001,0.003,0.005,0.01])}
ext = {}
for r,(d,suf,mus) in spec.items():
    for mu in mus:
        for s in range(1,21):
            z = json.load(open(f"{B}/{d}/mat_d8_mu{str(mu).replace('.','p')}_s{s:05d}{suf}.json"))
            assert abs(z["opts"]["mu"]-mu)<1e-15
            if z["extinct_at"] >= 0: ext[(r,mu,s)] = (z["extinct_at"], z["checksum"]); continue
            num = den = 0; Le = [0,0]
            for code, n in z["final_counts"]:
                L = L_of(code)
                if L is None: continue
                num += L*n; den += n
                if 'e' not in code: Le[0]+=L*n; Le[1]+=n
            cells[(r,mu,s)] = (num/den, Le[0]/Le[1])
print("extinct", len(ext)); [print("  ",k,v) for k,v in sorted(ext.items())]
def m(r,mu,seeds,i=0):
    v=[cells[(r,mu,s)][i] for s in seeds if (r,mu,s) in cells]; return sum(v)/len(v) if v else float('nan')
def q(seeds,i=0):
    rDF=m("DF",.01,seeds,i)-m("DF",.001,seeds,i); rDO=m("DO",.01,seeds,i)-m("DO",.001,seeds,i); rFF=m("FF",.01,seeds,i)-m("FF",.001,seeds,i)
    return dict(rDF=rDF,rDO=rDO,I=rDF-rDO,R=(rDF-rDO)/rDF,D3=m("DF",.01,seeds,i)-m("DO",.01,seeds,i),D4=m("DF",.01,seeds,i)-m("DO",.05,seeds,i),rFF=rFF,I_FF=rDF-rFF,R_FF=(rDF-rFF)/rDF)
S=list(range(1,21))
for r,(d,suf,mus) in spec.items():
    print(r, " ".join("%g:%.3f(n%d)"%(mu*100, m(r,mu,S), sum((r,mu,s) in cells for s in S)) for mu in mus))
def show(tag,seeds,i=0,boot=True):
    p=q(seeds,i); out={}
    if boot:
        rng=random.Random(12345); bs={k:[] for k in p}; nan=0
        for _ in range(4000):
            smp=[rng.choice(seeds) for _ in seeds]; o=q(smp,i)
            for k in o:
                if o[k]==o[k]: bs[k].append(o[k])
                else: nan+=1
        for k in p:
            v=sorted(bs[k]); out[k]="%.3f[%.3f,%.3f]"%(p[k],v[int(.025*len(v))],v[int(.975*len(v))-1])
        print(tag, out, "nan", nan)
    else: print(tag, {k:round(v,3) for k,v in p.items()})
show("ALL seeds (prereg def)", S)
show("ALL noe", S, 1)
common=[s for s in S if all(k in cells for k in [("DF",.001,s),("DF",.01,s),("DO",.001,s),("DO",.01,s)])]
print("common DF/DO 0.1&1%:", len(common), sorted(set(S)-set(common)))
show("COMMON set", common)
common5=[s for s in common if ("DO",.05,s) in cells]
print("common incl DO5%:", len(common5)); show("COMMON+DO5", common5)
commonFF=[s for s in common if ("FF",.001,s) in cells and ("FF",.01,s) in cells]; show("COMMON+FF", commonFF)
