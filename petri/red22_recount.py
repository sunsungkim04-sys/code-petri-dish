#!/usr/bin/env python3
"""해부22 독립 재계수 — 사전등록 §4 정의만으로 작성 (판정기 red22_analyze.py 는 열지 않음).
정의(원래 판 _RESULT_reduced_d.json 에서 0.633/0.875/0.578/+0.261/+0.232/0.639/1.079/1.151 재현으로 확인):
  R̄(code) = 시드 평균 R (rule=orig, mu=0)
  비 = (R̄-1)*L / S_dish ; S_dish = _RESULT_pairs14.json['code|<fam>|<code>']['S']
  M = 13 코드 비 중앙값 ; F = mean(racld 6) - mean(rascld 7)
  A1 = ln(class-mean R̄ L=4) - ln(L=5) [racld] ; ln(L=2) - ln(L=3) [rascld]
  D4 = OLS slope ln Ū(orig, mu=0.003) vs ln R̄(orig, mu=0), U>0 코드만
  D1 = 가족별 등급 평균 S=(R̄_cls-1)L 의 max/min
  부트스트랩: 시드 2,000회 재추출, 백분위 95%.
"""
import json, math, sys, random, statistics as st, hashlib, os

root = sys.argv[1]
petri = os.path.dirname(os.path.abspath(root.rstrip('/')))
seeds = list(range(2201, 2221))
P = json.load(open(os.path.join(petri, '_RESULT_pairs14.json')))
data = {}
cells = {}
codes = None
probs = []
for s in seeds:
    d = json.load(open(os.path.join(root, 's%d.json' % s)))
    a = d['args']
    exp = dict(cells=1000, density=8, grow=6000, ticks=3000, seeds=1, seed0=s, codes_of='racld,rascld',
               mus='0,0.003', rules='orig,rem', resident='racld', frac=0.1)
    for k, v in exp.items():
        if a.get(k) != v: probs.append('seed %d arg %s=%r' % (s, k, a.get(k)))
    if len(d['arm']) != 52: probs.append('seed %d arms %d' % (s, len(d['arm'])))
    if codes is None: codes = d['codes']
    elif codes != d['codes']: probs.append('codes differ %d' % s)
    for r in d['arm']:
        if r['seed'] != s: probs.append('seed field %d' % s)
        if not r['m0'] > 0: probs.append('m0 %d' % s)
        if not (r['R'] == r['R'] and r['R'] > 0): probs.append('R nan %d %s' % (s, r['code']))
    data[s] = d['arm']
    cells[s] = d['cells'][0]
    if cells[s]['pop'] <= 0: probs.append('pop %d' % s)
# identical arm lists across seeds
sig = {}
for s in seeds:
    h = hashlib.sha256(json.dumps([{k: v for k, v in r.items() if k != 'seed'} for r in data[s]], sort_keys=True).encode()).hexdigest()
    if h in sig: probs.append('identical arms %d %d' % (sig[h], s))
    sig[h] = s
impl = open(os.path.join(root, 'impl.sha256')).read().split()[0]
if impl != 'd9473aa0c1faaecaeb36df11d6cb3e6a62194f86e65c34ab0ec5a258da70ea43': probs.append('impl')
print('G0-like problems:', len(probs), probs[:5])
print('codes', len(codes))

fam = lambda c: c.replace('n', '', 1)
CODES = list(codes)

def get(s, c, rule, mu, k):
    for r in data[s]:
        if r['code'] == c and r['rule'] == rule and r['mu'] == mu: return r[k]
    raise KeyError

# pre-index per seed
R0 = {s: {c: get(s, c, 'orig', 0.0, 'R') for c in CODES} for s in seeds}
U3 = {s: {c: get(s, c, 'orig', 0.003, 'U') for c in CODES} for s in seeds}
Rrem0 = {s: {c: get(s, c, 'rem', 0.0, 'R') for c in CODES} for s in seeds}

def slope(x, y):
    mx = sum(x) / len(x); my = sum(y) / len(y)
    return sum((a - mx) * (b - my) for a, b in zip(x, y)) / sum((a - mx) ** 2 for a in x)

def stats(sample):
    Rb = {c: sum(R0[s][c] for s in sample) / len(sample) for c in CODES}
    Ub = {c: sum(U3[s][c] for s in sample) / len(sample) for c in CODES}
    rat = {c: (Rb[c] - 1) * codes[c] / P['code|%s|%s' % (fam(c), c)]['S'] for c in CODES}
    M = st.median(rat.values())
    fr = st.mean([rat[c] for c in CODES if fam(c) == 'racld'])
    fs = st.mean([rat[c] for c in CODES if fam(c) == 'rascld'])
    def cls(f, L): cs = [c for c in CODES if fam(c) == f and codes[c] == L]; return sum(Rb[c] for c in cs) / len(cs)
    a1r = math.log(cls('racld', 4)) - math.log(cls('racld', 5))
    a1s = math.log(cls('rascld', 2)) - math.log(cls('rascld', 3))
    d1r = max((cls('racld', L) - 1) * L for L in (4, 5)) / min((cls('racld', L) - 1) * L for L in (4, 5))
    d1s = max((cls('rascld', L) - 1) * L for L in (2, 3)) / min((cls('rascld', L) - 1) * L for L in (2, 3))
    cs = [c for c in CODES if Ub[c] > 0]
    d4 = slope([math.log(Rb[c]) for c in cs], [math.log(Ub[c]) for c in cs])
    out = dict(M=M, F=fr - fs, fam_racld=fr, fam_rascld=fs, A1_racld=a1r, A1_rascld=a1s, D1_racld=d1r, D1_rascld=d1s,
               D4=d4, D4_n=len(cs), rmin=min(rat.values()), rmax=max(rat.values()))
    for c in CODES: out['ratio_' + c] = rat[c]
    out['Mcode'] = min(CODES, key=lambda c: abs(rat[c] - M))
    return out

pt = stats(seeds)
rng = random.Random(22)
B = [stats([rng.choice(seeds) for _ in seeds]) for _ in range(2000)]
def ci(k):
    v = sorted(b[k] for b in B if b[k] == b[k])
    return v[int(0.025 * len(v))], v[int(0.975 * len(v)) - 1], len(B) - len(v)
for k in ['M', 'F', 'fam_racld', 'fam_rascld', 'A1_racld', 'A1_rascld', 'D4', 'D1_racld', 'D1_rascld', 'rmin', 'rmax'] + ['ratio_' + c for c in CODES]:
    lo, hi, nd = ci(k)
    print('%-16s %.4f [%.4f, %.4f] nan_drop=%d' % (k, pt[k], lo, hi, nd))
print('M code', pt['Mcode'], 'D4 n', pt['D4_n'], 'boot D4_n range', min(b['D4_n'] for b in B), max(b['D4_n'] for b in B))
print('rem rule mean R per code (descr):', {c: round(sum(Rrem0[s][c] for s in seeds) / 20, 3) for c in CODES})
pops = [cells[s]['pop'] for s in seeds]; fr = [cells[s]['free_used_frac'] for s in seeds]
print('pop %d~%d mean %.1f ; free %.4f~%.4f mean %.4f' % (min(pops), max(pops), st.mean(pops), min(fr), max(fr), st.mean(fr)))
# 판정
Mlo, Mhi, _ = ci('M'); Flo, _, _ = ci('F'); a, _, _ = ci('A1_racld'); b, _, _ = ci('A1_rascld'); dlo, dhi, _ = ci('D4')
print('P1', 0.5 < Mlo and Mhi < 1.0, 'P2', Flo > 0, 'P3', a > 0 and b > 0, 'P4', 0.484 <= dlo and dhi <= 0.884, '0.684 in D4 CI', dlo <= 0.684 <= dhi)
