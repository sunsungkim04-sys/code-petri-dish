#!/usr/bin/env python3
"""해부27 독립 재계수 — mini27_analyze.py · _RESULT_mini27* 를 열지 않고
prereg27_frozen.md §1-§5 (와 그것이 가리키는 해부26 §1 지표 정의)만으로 다시 센다.

series 열 = [t, n, m, mf, nf] (mini27.py census:523-533).
s = ln((m_T+0.5)/m_0) - ln((n_T-m_T+0.5)/(n_0-m_0)) · Δ = s(mut) - mean s(ref0..4).
R = draws / positions (계통·고리 키 '1_2' 끼운 · '0_4' 상주). I = err_realized/positions/μ.
재추출 = 배경 2,000 · 95% 백분위(G1 99%) · 이 스크립트만의 난수(sha256('recount27|칸')).
의존: numpy. 쓰기: _RECOUNT27_independent.txt 만.
usage: /usr/bin/python3 recount27_independent.py  (petri/ 에서)
"""
import json, os, sys, hashlib, math
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "mini27_out")
REPRO = os.path.join(HERE, "mini27_repro")
IDENT = os.path.join(HERE, "ident27_lab101")
M26 = os.path.join(HERE, "mini26_out")
NB = 2000
NBG = 20
MAX_DROP = 2
EXT = [0.0, 0.001625, 0.00325, 0.004875, 0.0065, 0.008125, 0.0121875, 0.01625]
ORIG = [0.0, 0.001, 0.002, 0.003, 0.004, 0.005, 0.0075, 0.01]
lines = []
problems = []


def P(s=""):
    lines.append(s)


def rng(name):
    h = hashlib.sha256(("recount27|" + name).encode()).digest()
    return np.random.default_rng(int.from_bytes(h[:8], "little"))


def load(d, rho, seed):
    p = os.path.join(d, "mini27_r%d_s%d.json" % (rho, seed))
    return json.load(open(p)) if os.path.exists(p) else None


def sval(series):
    t0, tT = series[0], series[-1]
    n0, m0, nT, mT = t0[1], t0[2], tT[1], tT[2]
    return math.log((mT + 0.5) / m0) - math.log((nT - mT + 0.5) / (n0 - m0))


def sval_f(series):  # 가임만(mf, nf)
    t0, tT = series[0], series[-1]
    n0, m0, nT, mT = t0[4], t0[3], tT[4], tT[3]
    return math.log((mT + 0.5) / m0) - math.log((nT - mT + 0.5) / (n0 - m0))


def mukey(mu):
    return round(mu * 1e7)


# ---------- 배경 고르기 + G0 ----------
BG = {}
DATA = {}
for rho in (8, 32):
    cands = sorted(int(f.split("_s")[1].split(".")[0]) for f in os.listdir(OUT)
                   if f.startswith("mini27_r%d_s" % rho))
    sel = []
    for sd in cands:
        if not (27001 <= sd <= 27030):
            problems.append("시드 대역 밖 %d" % sd)
        d = load(OUT, rho, sd)
        q = all(d["by_geom"][g]["grow"]["qualified"] for g in ("glob", "loc", "place"))
        if q and len(sel) < NBG:
            sel.append(sd)
            DATA[(rho, sd)] = d
    # 끊김 없이 시드 순 앞 20 인지
    first = []
    for sd in cands:
        d = load(OUT, rho, sd)
        if all(d["by_geom"][g]["grow"]["qualified"] for g in ("glob", "loc", "place")):
            first.append(sd)
    if first[:NBG] != sel or len(sel) < NBG:
        problems.append("배경 선택 ρ%d" % rho)
    BG[rho] = sel

cfg_exp = dict(K=1000, grow=6000, assay=3000, every=250, n0=100, frac=0.1, nref=5, min_pop=50)
narms = dict(glob=180, loc=180, place=54)
init_skip = {}
for (rho, sd), d in DATA.items():
    if d["version"] != "mini27 v1.0.0":
        problems.append("버전 %s" % d["version"])
    c = d["config"]
    for k, v in cfg_exp.items():
        if c[k] != v:
            problems.append("설정 %s %s" % (k, c[k]))
    if abs(c["death"] - 1 / 450) > 1e-15 or str(c["mat_radius"]) != "1" or d["grid"] != [40, 25]:
        problems.append("설정 d/반경/격자")
    if len(d["mus"]) != 15 or d["rules"] != ["roll", "find"]:
        problems.append("μ/규칙")
    for g in ("glob", "loc", "place"):
        G = d["by_geom"][g]
        gr = G["grow"]
        if not gr["cons_ok"] or gr["cons_checks"] < 26 or gr["cons_bad"]:
            problems.append("보존 성장 %d %d %s" % (rho, sd, g))
        init_skip[(rho, sd, g)] = gr.get("init_skipped", 0)
        arms = G["arms"]
        if len(arms) != narms[g]:
            problems.append("팔 수 %d %d %s %d" % (rho, sd, g, len(arms)))
        inj = set((a["n_inject"], a["skipped"]) for a in arms)
        if len(inj) != 1:
            problems.append("끼우기 동일성 %d %d %s" % (rho, sd, g))
        for a in arms:
            if not a["cons_ok"] or a["cons_checks"] < 14 or a["cons_bad"]:
                problems.append("보존 팔 %d %d %s" % (rho, sd, g))
        # μ0 규칙 동일성
        A0 = {(a["rule"], a["kind"]): a for a in arms if a["mu"] == 0.0}
        for kind in ["mut"] + ["ref%d" % i for i in range(5)]:
            r, f = A0.get(("roll", kind)), A0.get(("find", kind))
            if r is None or f is None:
                problems.append("μ0 팔 없음 %s %s" % (g, kind)); continue
            if r["series"] != f["series"] or set(r["counters"]) != set(f["counters"]) or \
               {k: v["positions"] for k, v in r["counters"].items()} != {k: v["positions"] for k, v in f["counters"].items()}:
                problems.append("μ0 규칙 동일성 %d %d %s %s" % (rho, sd, g, kind))

# 재현
for rho in (8, 32):
    a, b = load(OUT, rho, 27001), load(REPRO, rho, 27001)
    def strip(x):
        if isinstance(x, dict):
            return {k: strip(v) for k, v in x.items() if k != "seconds"}
        if isinstance(x, list):
            return [strip(v) for v in x]
        return x
    if b is None or json.dumps(strip(a), sort_keys=True) != json.dumps(strip(b), sort_keys=True):
        problems.append("재현 ρ%d" % rho)

# G2a: ident27_lab101 glob vs mini26_out (성장 · 팔)
g2a_same = 0
for rho in (8, 32):
    a = load(IDENT, rho, 26001)
    b = json.load(open(os.path.join(M26, "mini26_r%d_s26001.json" % rho)))
    G = a["by_geom"]["glob"]
    ga = {k: v for k, v in G["grow"].items()}
    gb = {k: v for k, v in b["grow"].items()}
    if json.dumps(ga, sort_keys=True) != json.dumps(gb, sort_keys=True):
        problems.append("G2a 성장 ρ%d" % rho)
    if len(G["arms"]) != len(b["arms"]):
        problems.append("G2a 팔 수 ρ%d" % rho)
    for x, y in zip(G["arms"], b["arms"]):
        if json.dumps(x, sort_keys=True) == json.dumps(y, sort_keys=True):
            g2a_same += 1
        else:
            problems.append("G2a 팔 ρ%d" % rho)

# 해부26 결과 해시
h26 = hashlib.sha256(open(os.path.join(HERE, "_RESULT_mini26.json"), "rb").read()).hexdigest()
if not h26.startswith("e863ec85"):
    problems.append("해부26 해시")
ref26 = json.load(open(os.path.join(HERE, "_RESULT_mini26.json")))["P4"]["lnS"]

P("배경 ρ8 %s~%s (%d) · ρ32 %s~%s (%d)" % (BG[8][0], BG[8][-1], len(BG[8]), BG[32][0], BG[32][-1], len(BG[32])))
P("G0+G2a 문제 %d %s · G2a 바이트 같은 팔 %d · 해부26 해시 %s…" % (len(problems), problems[:5], g2a_same, h26[:8]))
P("처음 배치 건너뜀(서술): " + " ".join("ρ%d·%s=%d" % (r, g, sum(init_skip[(r, s, g)] for s in BG[r]))
                                      for r in (8, 32) for g in ("glob", "loc", "place")))


# ---------- 칸 행렬 ----------
def arm(rho, sd, g, rule, mu, kind):
    for a in DATA[(rho, sd)]["by_geom"][g]["arms"]:
        if a["rule"] == rule and mukey(a["mu"]) == mukey(mu) and a["kind"] == kind:
            return a
    return None


def dead(a):
    return a["series"][-1][1] == 0 or a["extinct_tick"] != -1


def cell(rho, g, rule, mu, fertile=False):
    """배경마다 Δ (뺀 배경은 nan), ref-Δ(G1), 뺀 수, 끼운 계통 소멸 수"""
    f = sval_f if fertile else sval
    D, G1, drop, lost = [], [], 0, 0
    for sd in BG[rho]:
        m = arm(rho, sd, g, rule, mu, "mut")
        refs = [arm(rho, sd, g, rule, mu, "ref%d" % i) for i in range(5)]
        if m is None or any(r is None for r in refs):
            D.append(np.nan); G1.append(np.nan); drop += 1; continue
        if dead(m) or any(dead(r) for r in refs):
            D.append(np.nan); G1.append(np.nan); drop += 1; continue
        if m["series"][-1][2] == 0:
            lost += 1
        sr = [f(r["series"]) for r in refs]
        D.append(f(m["series"]) - np.mean(sr))
        G1.append(sr[0] - np.mean(sr[1:]))
    return np.array(D), np.array(G1), drop, lost


def boot_mean(x, name, q=(2.5, 97.5)):
    x = x[~np.isnan(x)]
    r = rng(name)
    idx = r.integers(0, len(x), size=(NB, len(x)))
    bm = x[idx].mean(1)
    return x.mean(), np.percentile(bm, q[0]), np.percentile(bm, q[1])


def cross(curve, ladder):
    c = np.asarray(curve, float)
    if np.any(np.isnan(c)):
        return np.nan
    if c[0] <= 0:
        return 0.0
    for i in range(1, len(c)):
        if c[i] <= 0:
            x0, x1, y0, y1 = ladder[i - 1], ladder[i], c[i - 1], c[i]
            return x0 + (x1 - x0) * y0 / (y0 - y1)
    return math.inf


def curve_mat(rho, g, rule, ladder):
    cols = []
    drops = []
    for mu in ladder:
        D, _, dr, _ = cell(rho, g, rule, mu)
        cols.append(D); drops.append(dr)
    return np.array(cols).T, drops  # 배경 × μ


def boot_cross(Mx, ladder, name, idx=None):
    nb = Mx.shape[0]
    if idx is None:
        idx = rng(name).integers(0, nb, size=(NB, nb))
    out = np.array([cross(np.nanmean(Mx[i], 0), ladder) for i in idx])
    return out


# ---------- G1 ----------
P()
for mu in (0.0, 0.01625):
    _, G1, dr, _ = cell(8, "loc", "roll", mu)
    m, lo, hi = boot_mean(G1, "G1|%s" % mu, (0.5, 99.5))
    P("G1 loc·roll·ρ8 μ%.4f%%: %+.3f [99%% %+.3f, %+.3f] 뺀 %d → %s" %
      (mu * 100, m, lo, hi, dr, "0 품음" if lo <= 0 <= hi else "0 안 품음"))

# ---------- Δ 표 ----------
P()
drop_tot = {}
lost_tot = {}
conds = [("loc", "roll"), ("loc", "find"), ("glob", "roll"), ("glob", "find"), ("place", "roll")]
DT = {}
for rho in (8, 32):
    for g, rule in conds:
        lad = EXT if g == "place" and rule == "roll" else sorted(set(EXT) | set(ORIG))
        if g == "place" and rule == "find":
            lad = [0.0]
        row = []
        for mu in lad:
            D, _, dr, lo_ = cell(rho, g, rule, mu)
            drop_tot[(rho, g, rule)] = drop_tot.get((rho, g, rule), 0) + dr
            lost_tot[(rho, g, rule)] = lost_tot.get((rho, g, rule), 0) + lo_
            m, lo, hi = boot_mean(D, "D|%d|%s|%s|%s" % (rho, g, rule, mu))
            DT[(rho, g, rule, mukey(mu))] = (m, lo, hi)
            tag = ("E" if mu in EXT else "") + ("O" if mu in ORIG else "")
            row.append("%s%.4g:%+.3f[%+.3f,%+.3f]" % (tag, mu * 100, m, lo, hi))
        P("Δ ρ%d %s·%s: " % (rho, g, rule) + " ".join(row))
P("뺀 배경 합 %s · 끼운 계통 소멸 합 %s" % (sum(drop_tot.values()), sum(lost_tot.values())))


# ---------- P1 · K1 · K2 ----------
def pform(g, tagname):
    Mx, drops = curve_mat(8, g, "roll", EXT)
    d0 = DT[(8, g, "roll", 0)]
    dE = DT[(8, g, "roll", mukey(0.01625))]
    c = cross(np.nanmean(Mx, 0), EXT)
    bc = boot_cross(Mx, EXT, "cross|%s" % g)
    finf = np.mean(np.isinf(bc))
    und = drops[0] > MAX_DROP or drops[-1] > MAX_DROP
    return d0, dE, c, finf, und, bc, Mx


P()
res = {}
for g, nm in (("loc", "P1"), ("glob", "K1"), ("place", "K2")):
    d0, dE, c, finf, und, bc, Mx = pform(g, nm)
    rev = d0[1] > 0 and dE[2] < 0
    if nm == "P1":
        v = "판정 불가" if und else ("✅" if rev else ("❌ (가)" if d0[1] <= 0 else "❌ (나)"))
    else:
        v = "판정 불가" if und else ("반전" if rev else ("반전 없음" if (dE[1] > 0 or finf >= 0.975) else "불확정"))
    res[nm] = (v, d0, dE, c, finf, bc, Mx)
    P("%s %s·roll·ρ8: Δ(0) %+.3f [%+.3f, %+.3f] · Δ(1.625%%) %+.3f [%+.3f, %+.3f] · 교차 %s · 재추출 inf 비율 %.3f · 교차 재추출 유한 중앙 %s → %s" %
      (nm, g, *d0, *dE, c, finf, (np.median(bc[np.isfinite(bc)]) if np.isfinite(bc).any() else "—"), v))

# ---------- G2b ----------
Mo, _ = curve_mat(8, "glob", "roll", ORIG)
c_o = cross(np.nanmean(Mo, 0), ORIG)


def lnRS(rho, g, mu=0.0, rule="roll"):
    lnR, lnS, keep = [], [], 0
    R2s, R4s = [], []
    for sd in BG[rho]:
        a = arm(rho, sd, g, rule, mu, "mut")
        c = a["counters"]
        if "1_2" not in c or "0_4" not in c:
            continue
        p2, p4 = c["1_2"]["positions"], c["0_4"]["positions"]
        if p2 < 100 or p4 < 100:
            continue
        R2, R4 = c["1_2"]["draws"] / p2, c["0_4"]["draws"] / p4
        if R2 <= 1 or R4 <= 1:
            continue
        keep += 1
        R2s.append(R2); R4s.append(R4)
        lnR.append(math.log(R2) - math.log(R4))
        lnS.append(math.log((R2 - 1) * 2) - math.log((R4 - 1) * 4))
    return np.array(lnR), np.array(lnS), keep, np.median(R2s), np.median(R4s)


gR, gS, gk, gR2, gR4 = lnRS(8, "glob")
mS, loS, hiS = boot_mean(gS, "G2b|lnS")
ov = not (hiS < ref26[1] or loS > ref26[2])
g2b = "✅" if (math.isinf(c_o) and c_o > 0 and ov) else "🟡"
P()
P("G2b glob·roll·ρ8 원 사다리 교차 %s · ln S비 %+.3f [%+.3f, %+.3f] = 비 %.3f [%.3f, %.3f] · 해부26 [%.3f, %.3f] 겹침 %s → %s" %
  (c_o, mS, loS, hiS, math.exp(mS), math.exp(loS), math.exp(hiS), math.exp(ref26[1]), math.exp(ref26[2]), ov, g2b))

# ---------- P3 (배경 짝) ----------
Mr, dr_r = curve_mat(8, "loc", "roll", EXT)
Mf, dr_f = curve_mat(8, "loc", "find", EXT)
idx = rng("P3|pair").integers(0, NBG, size=(NB, NBG))
cr = boot_cross(Mr, EXT, None, idx)
cf = boot_cross(Mf, EXT, None, idx)
valid = np.isfinite(cr) & (cr > 0)
push = valid & ((np.isinf(cf)) | (cf >= 2 * cr))
stay = valid & np.isfinite(cf) & (cf < 2 * cr)
p_push, p_stay, p_valid = push.mean(), stay.mean(), valid.mean()
p3raw = "✅" if p_push >= 0.975 else ("❌" if p_stay >= 0.975 else "🟡")
P1v = res["P1"][0]
p3 = ("판정 보류(P1 ❌) — 서술: %s" % p3raw) if P1v.startswith("❌") else p3raw
c_rf = cross(np.nanmean(Mr, 0), EXT); c_ff = cross(np.nanmean(Mf, 0), EXT)
P("P3 loc: c_roll %s · c_find %s · P(밀림) %.3f · P(남음) %.3f · P(c_roll 유효) %.3f → %s" %
  (c_rf, c_ff, p_push, p_stay, p_valid, p3))

# ---------- P4 ----------
lR, lS, lk, lR2, lR4 = lnRS(8, "loc")
pR, pS, pk, pR2, pR4 = lnRS(8, "place")
r1, r2 = rng("P4a|loc"), rng("P4a|glob")
bl = lR[r1.integers(0, len(lR), size=(NB, len(lR)))].mean(1)
bg_ = gR[r2.integers(0, len(gR), size=(NB, len(gR)))].mean(1)
diff = bl - bg_
dpt = lR.mean() - gR.mean()
dlo, dhi = np.percentile(diff, [2.5, 97.5])
p4a = "판정 불가" if (lk < 15 or gk < 15) else ("✅" if dlo > 0 else "❌")
m4, lo4, hi4 = boot_mean(lS, "P4b|lnS")
a, b = math.log(0.8), math.log(1.25)
p4b = "판정 불가" if lk < 15 else ("✅" if (lo4 >= a and hi4 <= b) else ("❌" if (hi4 < a or lo4 > b) else "🟡"))
p4 = "✅" if (p4a == "✅" and p4b == "✅") else ("❌" if "❌" in (p4a, p4b) else "🟡")
P()
P("P4 자격 loc %d · glob %d · place %d" % (lk, gk, pk))
P("R 중앙 μ0 roll ρ8: loc L2 %.2f L4 %.2f · glob %.2f %.2f · place %.2f %.2f" % (lR2, lR4, gR2, gR4, pR2, pR4))
P("P4a loc lnR차 %+.3f · glob %+.3f · loc−glob %+.3f [%+.3f, %+.3f] → %s" % (lR.mean(), gR.mean(), dpt, dlo, dhi, p4a))
P("P4b loc ln S비 %+.5f [%+.5f, %+.5f] = 비 %.3f [%.3f, %.3f] · ln0.8 %+.5f · 상한−ln0.8 %+.5f → %s" %
  (m4, lo4, hi4, math.exp(m4), math.exp(lo4), math.exp(hi4), a, hi4 - a, p4b))
# P4b 경계 안정성: 다른 난수 흐름 10 개
alt = []
for k in range(10):
    _, l_, h_ = boot_mean(lS, "P4b|alt%d" % k)
    alt.append(h_)
alt = np.array(alt)
P("P4b 상한 흐름 10개: 최소 %+.5f 최대 %+.5f · ln0.8 아래(❌) %d/10" % (alt.min(), alt.max(), int((alt < a).sum())))
P("P4 → %s" % p4)

# ---------- 배경 조작 확인 · I ----------
P()
for rho in (8, 32):
    s = []
    for g in ("glob", "loc", "place"):
        Rs = []
        for sd in BG[rho]:
            c = DATA[(rho, sd)]["by_geom"][g]["grow"]["counters"]["0_4"]
            Rs.append(c["draws"] / c["positions"])
        s.append("%s %.2f" % (g, np.median(Rs)))
    P("배경 R 중앙(상주) ρ%d: %s" % (rho, " · ".join(s)))


def Iratio(rho, g, mu, key):
    v = []
    for sd in BG[rho]:
        c = arm(rho, sd, g, "roll", mu, "mut")["counters"].get(key)
        if c and c["positions"] > 0:
            v.append(c["err_realized"] / c["positions"] / mu)
    return np.median(v)


P("I loc·roll·ρ8 μ1.625%%: L2 %.2f · L4 %.2f · 비 %.2f" %
  (Iratio(8, "loc", 0.01625, "1_2"), Iratio(8, "loc", 0.01625, "0_4"),
   Iratio(8, "loc", 0.01625, "1_2") / Iratio(8, "loc", 0.01625, "0_4")))

# ---------- 결과 · 해석 ----------
K1v, K2v = res["K1"][0], res["K2"][0]
if P1v == "판정 불가":
    br = "해석하지 않음"
elif P1v == "✅":
    if K1v == "반전":
        br = "B"
    elif K1v == "판정 불가" or K2v == "판정 불가":
        br = "'필요한 성분' 안 씀"
    elif K1v == "반전 없음" and K2v == "반전 없음":
        br = "A"
    elif K2v == "반전":
        br = "C"
    else:
        br = "국소에서만 구간상 반전 · 대조 불확정"
else:
    if K1v == "반전":
        br = "E"
    elif K1v == "판정 불가":
        br = "E/D 구분 불가"
    else:
        br = "D(가)" if "가" in P1v else "D(나)"
P()
P("결과: 국소 결핍 모형에서 P1 %s · G2b %s · P3 %s · P4 %s (P4a %s · P4b %s)" % (P1v, g2b, p3, p4, p4a, p4b))
P("K1 %s · K2 %s → 갈래 %s" % (K1v, K2v, br))

txt = "\n".join(lines) + "\n"
open(os.path.join(HERE, "_RECOUNT27_independent.txt"), "w").write(txt)
sys.stdout.write(txt)
