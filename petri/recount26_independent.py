#!/usr/bin/env python3
"""해부26 독립 재계수 (recount26_independent.py)

판정기 mini26_analyze.py 와 _RESULT_mini26* 를 열지 않고, 사전등록 노트
(사전등록/해부26-사전등록-2026-10-02.md §1 · §4 · §5)의 정의 · 판정선만으로
mini26_out/ · mini26_repro/ 원자료 JSON 에서 G0(일부) · G1 · P1~P4 와 서술 칸을 다시 계산한다.

정의 (노트에서 옮김):
- s = ln((m_T+0.5)/m_0) - ln((n_T-m_T+0.5)/(n_0-m_0)) · series 행 = [t, n, m, mf, nf] (mini26.py census)
- Δ = s(mut) - mean s(ref0..4) · 배경마다
- 교차 = 평균 Δ 곡선(μ %)의 첫 + → ≤0 부호 변화 선형 보간 · Δ(0) 평균 ≤ 0 → 0 · 1% 까지 + → +inf
- 재추출 2,000 · 95% 백분위 (G1 99%) · 난수는 이 스크립트 고유 (판정기와 다름)
- 칸에서 팔 하나라도 개체 전멸(extinct_tick >= 0)이면 그 배경을 그 칸에서 뺀다 · > 2 이면 판정 불가
- P4: roll · ρ8 · μ0 · 끼운 팔 · R = draws/positions · 계통 1_2 (L2) · 0_4 (L4) · 자격 positions ≥ 100 · R > 1 (두 계통) · k ≥ 15
"""
import hashlib, json, math, os, sys
import numpy as np

ROOT = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(ROOT, "mini26_out")
REPRO = os.path.join(ROOT, "mini26_repro")
MUS = [0.0, 0.001, 0.002, 0.003, 0.004, 0.005, 0.0075, 0.01]
MUP = [m * 100 for m in MUS]
RHOS = [8, 32]
RULES = ["roll", "find"]
NB = 2000
NBG = 20
MAXDROP = 2
LO, HI = math.log(0.8), math.log(1.25)
rng_master = np.random.default_rng(int(hashlib.sha256(b"recount26-independent").hexdigest()[:12], 16))


def load(path):
    with open(path) as f:
        return json.load(f)


# ---------------- 적재 · 배경 선택 ----------------
data = {}
for rho in RHOS:
    seeds = sorted(int(fn.split("_s")[1].split(".")[0]) for fn in os.listdir(OUT) if fn.startswith("mini26_r%d_" % rho))
    qual = []
    for s in seeds:
        d = load(os.path.join(OUT, "mini26_r%d_s%d.json" % (rho, s)))
        g = d["grow"]
        if g["extinct_tick"] < 0 and g["n_end"] >= 50 and g["qualified"]:
            qual.append((s, d))
        if len(qual) == NBG:
            break
    data[rho] = qual
    print("배경 ρ%d: %d · 시드 %d~%d" % (rho, len(qual), qual[0][0], qual[-1][0]))

# ---------------- G0 (재계수자가 볼 수 있는 부분) ----------------
problems = []
for rho in RHOS:
    for s, d in data[rho]:
        c = d["config"]
        if d["version"] != "mini26 v1.0.0": problems.append((rho, s, "version"))
        exp = dict(K=1000, grow=6000, assay=3000, every=250, n0=100, frac=0.1, nref=5, min_pop=50)
        for k, v in exp.items():
            if c[k] != v: problems.append((rho, s, "config " + k))
        if abs(c["death"] - 1 / 450) > 1e-15: problems.append((rho, s, "death"))
        if d["mus"] != MUS or d["rules"] != RULES: problems.append((rho, s, "mus/rules"))
        g = d["grow"]
        if not g["cons_ok"] or g["cons_bad"] or g["cons_checks"] < 26: problems.append((rho, s, "grow cons"))
        if len(d["arms"]) != 96: problems.append((rho, s, "arms"))
        inj = set()
        for a in d["arms"]:
            if not a["cons_ok"] or a["cons_bad"] or a["cons_checks"] < 14: problems.append((rho, s, "arm cons"))
            inj.add((a["n_inject"], a["skipped"]))
        if len(inj) != 1: problems.append((rho, s, "inject"))
        # μ0 규칙 동일성: 궤적 · 완료 위치 · 키 집합
        for kind in ["mut"] + ["ref%d" % i for i in range(5)]:
            ar = [a for a in d["arms"] if a["mu"] == 0.0 and a["kind"] == kind]
            r = {a["rule"]: a for a in ar}
            if r["roll"]["series"] != r["find"]["series"]: problems.append((rho, s, "mu0 series " + kind))
            if set(r["roll"]["counters"]) != set(r["find"]["counters"]): problems.append((rho, s, "mu0 keys " + kind))
            for k in r["roll"]["counters"]:
                if r["roll"]["counters"][k]["positions"] != r["find"]["counters"][k]["positions"]:
                    problems.append((rho, s, "mu0 pos " + kind))
# 재현 (seconds 제외 동일)
for rho in RHOS:
    fn = "mini26_r%d_s26001.json" % rho
    a, b = load(os.path.join(OUT, fn)), load(os.path.join(REPRO, fn))
    a.pop("seconds", None); b.pop("seconds", None)
    if json.dumps(a, sort_keys=True) != json.dumps(b, sort_keys=True): problems.append((rho, 26001, "repro"))
print("G0(재계수 범위) 문제:", len(problems), problems[:5])


# ---------------- 지표 ----------------
def s_of(series, fert=False):
    r0, rT = series[0], series[-1]
    if fert:
        n0, m0, nT, mT = r0[4], r0[3], rT[4], rT[3]
    else:
        n0, m0, nT, mT = r0[1], r0[2], rT[1], rT[2]
    return math.log((mT + 0.5) / m0) - math.log((nT - mT + 0.5) / (n0 - m0))


def arm(d, rule, mu, kind):
    for a in d["arms"]:
        if a["rule"] == rule and a["mu"] == mu and a["kind"] == kind:
            return a
    raise KeyError


# D[rho][rule] : (NBG, 8) Δ · NaN 은 뺀 배경
D, DF, G1D = {}, {}, {}
drops, lost = {}, {}
for rho in RHOS:
    for rule in RULES:
        M = np.full((NBG, 8), np.nan); MF = np.full((NBG, 8), np.nan); G = np.full((NBG, 8), np.nan)
        for i, (s, d) in enumerate(data[rho]):
            for j, mu in enumerate(MUS):
                arms = [arm(d, rule, mu, "mut")] + [arm(d, rule, mu, "ref%d" % k) for k in range(5)]
                if any(a["extinct_tick"] >= 0 for a in arms):
                    drops[(rho, rule, j)] = drops.get((rho, rule, j), 0) + 1
                    continue
                ss = [s_of(a["series"]) for a in arms]
                M[i, j] = ss[0] - np.mean(ss[1:])
                G[i, j] = ss[1] - np.mean(ss[2:])
                try:
                    sf = [s_of(a["series"], True) for a in arms]
                    MF[i, j] = sf[0] - np.mean(sf[1:])
                except (ValueError, ZeroDivisionError):
                    pass
                if arms[0]["series"][-1][2] == 0:
                    lost[(rho, rule, j)] = lost.get((rho, rule, j), 0) + 1
        D[(rho, rule)], DF[(rho, rule)], G1D[(rho, rule)] = M, MF, G
print("뺀 배경 합:", sum(drops.values()), drops, "· 끼운 계통 소멸 합:", sum(lost.values()), lost)


def cens_flip(curve):
    curve = np.asarray(curve, float)
    if not np.all(np.isfinite(curve)):
        return np.nan
    if curve[0] <= 0:
        return 0.0
    for j in range(1, len(curve)):
        if curve[j] <= 0:
            x0, x1, y0, y1 = MUP[j - 1], MUP[j], curve[j - 1], curve[j]
            return x0 + (x1 - x0) * y0 / (y0 - y1)
    return math.inf


def ci(v, lvl=95):
    v = np.asarray(v, float); v = v[np.isfinite(v)] if lvl else v
    a = (100 - lvl) / 2
    return np.percentile(v, a), np.percentile(v, 100 - a)


def boot_idx(n, rng):
    return rng.integers(0, n, size=(NB, n))


def fmt(x):
    return "%+.3f" % x


# ---------------- G1 ----------------
rng = np.random.default_rng(rng_master.integers(2**63))
G = G1D[(8, "roll")]
for j in (0, 7):
    v = G[:, j]; v = v[np.isfinite(v)]
    bi = boot_idx(len(v), rng)
    b = v[bi].mean(1)
    lo, hi = np.percentile(b, 0.5), np.percentile(b, 99.5)
    print("G1 μ%s%%: %s [99%% %s, %s] → %s" % (MUP[j], fmt(v.mean()), fmt(lo), fmt(hi), "0 품음" if lo <= 0 <= hi else "0 안 품음"))

# ---------------- Δ 표 · 교차 ----------------
print("\nΔ 표 (평균 [95%])")
boots = {}
for rho in RHOS:
    for rule in RULES:
        M = D[(rho, rule)]
        rng = np.random.default_rng(rng_master.integers(2**63))
        bi = boot_idx(NB, rng)[:1] if False else rng.integers(0, NBG, size=(NB, NBG))
        boots[(rho, rule)] = bi
        B = np.nanmean(M[bi], axis=1)  # (NB, 8)
        row = []
        for j in range(8):
            lo, hi = np.percentile(B[:, j], [2.5, 97.5])
            row.append("%s[%s,%s]" % (fmt(np.nanmean(M[:, j])), fmt(lo), fmt(hi)))
        c = cens_flip(np.nanmean(M, axis=0))
        cb = np.array([cens_flip(B[k]) for k in range(NB)])
        print("ρ%d %s: " % (rho, rule) + " | ".join(row))
        fin = cb[np.isfinite(cb)]
        print("   교차 %s · 재추출 inf 비율 %.3f · 유한 %d%s" % (c, np.mean(np.isinf(cb)), len(fin),
              (" · 유한 구간 [%.3f, %.3f]" % tuple(np.percentile(fin, [2.5, 97.5]))) if len(fin) else ""))
        MFm = DF[(rho, rule)]
        print("   (S3 가임만) " + " ".join(fmt(x) for x in np.nanmean(MFm, axis=0)) + " · 교차 %s" % cens_flip(np.nanmean(MFm, axis=0)))

# ---------------- P1 ----------------
M = D[(8, "roll")]
B = np.nanmean(M[boots[(8, "roll")]], axis=1)
l0 = np.percentile(B[:, 0], 2.5); h1 = np.percentile(B[:, 7], 97.5)
P1 = l0 > 0 and h1 < 0
print("\nP1: Δ(0) 하한 %s · Δ(1%%) 상한 %s → %s" % (fmt(l0), fmt(h1), "✅" if P1 else "❌"))


# ---------------- P2 · P3 (서술: P1 ❌ 이면 보류) ----------------
def events(crf, cff):
    ok = np.isfinite(crf) & (crf > 0)
    push = ok & (np.isinf(cff) | (np.isfinite(cff) & (cff >= 2 * crf)))
    stay = ok & np.isfinite(cff) & (cff < 2 * crf)
    return push.mean(), stay.mean(), ok.mean()


rng = np.random.default_rng(rng_master.integers(2**63))
bi = rng.integers(0, NBG, size=(NB, NBG))
crf = np.array([cens_flip(np.nanmean(D[(8, "roll")][b], 0)) for b in bi])
cff = np.array([cens_flip(np.nanmean(D[(8, "find")][b], 0)) for b in bi])
print("P2(짝): c_rf %s · c_ff %s · P(밀림) %.3f · P(남음) %.3f · P(c_rf 유효) %.3f" % (
    cens_flip(np.nanmean(D[(8, "roll")], 0)), cens_flip(np.nanmean(D[(8, "find")], 0)), *events(crf, cff)))
bj = rng.integers(0, NBG, size=(NB, NBG))
c32 = np.array([cens_flip(np.nanmean(D[(32, "roll")][b], 0)) for b in bj])
print("P3(독립): c_rf8 %s · c_rf32 %s · P(밀림) %.3f · P(남음) %.3f · P(c_rf 유효) %.3f · c32 재추출 inf 비율 %.3f" % (
    cens_flip(np.nanmean(D[(8, "roll")], 0)), cens_flip(np.nanmean(D[(32, "roll")], 0)), *events(crf, c32), np.mean(np.isinf(c32))))

# ---------------- P4 ----------------
lr, ls, R2s, R4s = [], [], [], []
for s, d in data[8]:
    c = arm(d, "roll", 0.0, "mut")["counters"]
    a2, a4 = c["1_2"], c["0_4"]
    R2, R4 = a2["draws"] / a2["positions"], a4["draws"] / a4["positions"]
    if a2["positions"] < 100 or a4["positions"] < 100 or R2 <= 1 or R4 <= 1:
        continue
    R2s.append(R2); R4s.append(R4)
    lr.append(math.log(R2) - math.log(R4))
    ls.append(math.log(((R2 - 1) * 2) / ((R4 - 1) * 4)))
lr, ls = np.array(lr), np.array(ls)
k = len(lr)
rng = np.random.default_rng(rng_master.integers(2**63))
bi = rng.integers(0, k, size=(NB, k))
a_lo, a_hi = np.percentile(lr[bi].mean(1), [2.5, 97.5])
b_lo, b_hi = np.percentile(ls[bi].mean(1), [2.5, 97.5])
P4a = a_lo > 0
if LO <= b_lo and b_hi <= HI: P4b = "✅"
elif b_hi < LO or b_lo > HI: P4b = "❌"
else: P4b = "🟡"
print("\nP4 자격 %d/20" % k)
print("P4a: %s [%s, %s] → %s · R 중앙 L2 %.2f · L4 %.2f" % (fmt(lr.mean()), fmt(a_lo), fmt(a_hi), "✅" if P4a else "❌", np.median(R2s), np.median(R4s)))
print("P4b: ln 비 %s [%s, %s] · 비 %.3f [%.3f, %.3f] → %s" % (fmt(ls.mean()), fmt(b_lo), fmt(b_hi), math.exp(ls.mean()), math.exp(b_lo), math.exp(b_hi), P4b))
P4 = "✅" if (P4a and P4b == "✅") else ("❌" if (not P4a or P4b == "❌") else "🟡")
print("P4 →", P4)

# ---------------- 서술: 배경 조작 · R/T/I 사다리 · S1 ③ ----------------
for rho in RHOS:
    Rg = [d["grow"]["counters"]["0_4"]["draws"] / d["grow"]["counters"]["0_4"]["positions"] for s, d in data[rho]]
    print("배경 R(상주) ρ%d 중앙 %.2f" % (rho, np.median(Rg)))
print("\nR 사다리(끼운 팔 · 계통별 중앙) · T/R − L 최소/최대")
trl = []
for rho in RHOS:
    for rule in RULES:
        row = []
        for mu in MUS:
            r2 = np.median([arm(d, rule, mu, "mut")["counters"]["1_2"]["draws"] / arm(d, rule, mu, "mut")["counters"]["1_2"]["positions"] for s, d in data[rho]])
            r4 = np.median([arm(d, rule, mu, "mut")["counters"]["0_4"]["draws"] / arm(d, rule, mu, "mut")["counters"]["0_4"]["positions"] for s, d in data[rho]])
            row.append("%.2f/%.2f" % (r2, r4))
            for s, d in data[rho]:
                for key, L in (("1_2", 2), ("0_4", 4)):
                    cc = arm(d, rule, mu, "mut")["counters"][key]
                    trl.append(cc["copy_ticks"] / cc["attempts"] - L)
        print("ρ%d %s: %s" % (rho, rule, " ".join(row)))
print("T/R − L: 최소 %.3g 최대 %.3g" % (min(trl), max(trl)))

# S1 ③
mu = 0.003
X, Y, bg = [], [], []
for i, (s, d) in enumerate(data[8]):
    c = arm(d, "roll", mu, "mut")["counters"]
    for key in ("1_2", "0_4"):
        cc = c[key]
        tot = cc["kids_ok"] + cc["kids_sterile"]
        if cc["kids_sterile"] > 0 and tot > 0 and cc["positions"] > 0:
            X.append(math.log(cc["draws"] / cc["positions"])); Y.append(math.log(cc["kids_sterile"] / tot)); bg.append(i)
X, Y, bg = map(np.array, (X, Y, bg))
slope = np.polyfit(X, Y, 1)[0]
rng = np.random.default_rng(rng_master.integers(2**63))
sl = []
for _ in range(NB):
    pick = rng.integers(0, NBG, NBG)
    idx = np.concatenate([np.where(bg == p)[0] for p in pick])
    sl.append(np.polyfit(X[idx], Y[idx], 1)[0])
lo, hi = np.percentile(sl, [2.5, 97.5])
print("\nS1 ③ 기울기 %s [%s, %s] · 점 %d · 점추정 띠[0.48,0.88] 안: %s · 구간 전체 띠 안: %s" % (
    fmt(slope), fmt(lo), fmt(hi), len(X), 0.48 <= slope <= 0.88, 0.48 <= lo and hi <= 0.88))

# 끼운 개체 수
for rho in RHOS:
    ni = [d["arms"][0]["n_inject"] for s, d in data[rho]]
    sk = sum(d["arms"][0]["skipped"] for s, d in data[rho])
    print("끼운 개체 ρ%d %d~%d · 건너뜀 %d" % (rho, min(ni), max(ni), sk))
