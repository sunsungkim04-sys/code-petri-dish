"""해부 3A 판정 — 에너지 빈도 해부. 사전등록 해부3-사전등록-2026-09-16.md §2 · §4 · §5. 결과 전에 썼다.

  Δ 는 해부 1 과 같다 — 같은 배경의 WT 팔 11개(ref + neutral 10) s 평균을 뺀다. 가짜 돌연변이는 자기를 뺀 10개 평균.
  자료: frac 0.01 · 0.5 는 이번 발사(dms_freq_*) · frac 0.1 은 해부 1(dms_main · μ 0)과 해부 2A(dms_mu1 · μ 1%)의 같은 코드 팔.
  판정(μ 마다 · 초점 코드마다): 짝 차이 Δ(50%) − Δ(1%) 의 배경 재추출 부트스트랩 95%
     상한 < 0 → 드물 때 더 낫다 · 하한 > 0 → 흔할 때 더 낫다 · 0 품음 → 빈도 의존 검출 안 됨
  덧붙임: Δ(50%) 구간이 0 아래거나 0 을 품으면 '흔해지면 이점이 사라진다'
  점검 ①②③ 은 배경을 빼고 보고 · ④ 가짜 돌연변이 10개 중 3개 이상이면 그 조합 보류 · ⑤ 양성 대조 reascl 이 해로움이어야 함
실행: python3 freq_analyze.py [petri 폴더, 기본 ~/petri]  →  _RESULT_freq.json
"""
import glob
import json
import math
import os
import random
import sys

BASE = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1 else "~/petri")
COND, WT = "energy", "reascld"
FOCAL = ["reashcld", "reaschld", "reasccld"]
POS = "reascl"
N_BG, N_BOOT = 60, 2000
RNG = random.Random(20260920)
SRC = [
    (0.01, 0.0, "dms_freq_f0p01_mu0"),
    (0.01, 0.01, "dms_freq_f0p01_mu0p01"),
    (0.1, 0.0, "dms_main"),
    (0.1, 0.01, "dms_mu1"),
    (0.5, 0.0, "dms_freq_f0p5_mu0"),
    (0.5, 0.01, "dms_freq_f0p5_mu0p01"),
]


def pct(xs, p):
    xs = sorted(xs)
    k = (len(xs) - 1) * p
    f = math.floor(k)
    c = min(f + 1, len(xs) - 1)
    return xs[f] + (xs[c] - xs[f]) * (k - f)


def s_of(arm):
    r = arm["series"]
    n0, m0 = r[0][1], r[0][2]
    nT, mT = r[-1][1], r[-1][2]
    return math.log((mT + 0.5) / m0) - math.log((nT - mT + 0.5) / (n0 - m0))


def load(folder, frac, mu):
    out, probs = {}, []
    for f in sorted(glob.glob(os.path.join(BASE, folder, f"dms_{COND}_{WT}_s*.json"))):
        z = json.load(open(f))
        why = []
        if not z["checksum_match"]:
            why.append("② 체크섬")
        if not z["fidelity"]["same"]:
            why.append("① 복제 충실도")
        if abs(z.get("mu_assay", 0.0) - mu) > 1e-12:
            why.append(f"μ {z.get('mu_assay')}")
        if abs(z["frac"] - frac) > 1e-12:
            why.append(f"frac {z['frac']}")
        if not z["grow_top"] or z["grow_top"][0][0] != WT:
            why.append("배경 규칙")
        if any(not a["cons_ok"] for a in z["arms"]):
            why.append("③ 재료 보존")
        n_expect = 11 + (z["n_mutants_total"] if z["arms_mode"] == "all" else len(z["codes"]))
        if len(z["arms"]) != n_expect:
            why.append(f"팔 수 {len(z['arms'])}≠{n_expect}")
        if why:
            probs.append((z["seed"], why))
            continue
        wt_arms = [a for a in z["arms"] if a["kind"] in ("ref", "neutral")]
        s_wt = {a["labels"][0]: s_of(a) for a in wt_arms}
        base = sum(s_wt.values()) / len(s_wt)
        d = {a["key"]: s_of(a) - base for a in z["arms"] if a["kind"] == "mut"}
        for lab, v in s_wt.items():
            if lab == "wt":
                continue
            others = [x for k, x in s_wt.items() if k != lab]
            d[lab] = v - sum(others) / len(others)
        out[z["seed"]] = d
    return out, probs


data, halt = {}, False
print(f"자료 — 초점 {', '.join(FOCAL)} · 양성 대조 {POS}")
for frac, mu, folder in SRC:
    d, probs = load(folder, frac, mu)
    data[(frac, mu)] = d
    print(f"  frac {frac:<5} μ {mu:<5} {folder:<22}: 쓴 배경 {len(d)}/{N_BG}" + ("" if not probs else " · 뺀 배경 " + " · ".join(f"시드 {s}:{'/'.join(w)}" for s, w in probs)))
    if len(d) != N_BG:
        halt = True
common = sorted(set.intersection(*[set(v) for v in data.values()])) if data else []
print(f"  짝지은 배경 {len(common)}")
if halt or len(common) < N_BG:
    print("🚨 계획한 배경이 다 모이지 않았다 — 판정하지 않는다")
    sys.exit(1)
idxs = [[RNG.randrange(len(common)) for _ in range(len(common))] for _ in range(N_BOOT)]


def ci(vals):
    n = len(vals)
    bm = [sum(vals[i] for i in ix) / n for ix in idxs]
    return sum(vals) / n, pct(bm, .025), pct(bm, .975)


def cls(lo, hi):
    return "해로움" if hi < 0 else ("이로움" if lo > 0 else "구별 안 됨")


out = {"n_bg": len(common), "checks": {}, "focal": {}}
print("\n== 점검 ④ 가짜 돌연변이 10개 중 효과 있음 · ⑤ 양성 대조")
hold = {}
for frac, mu, folder in SRC:
    d = data[(frac, mu)]
    n_false = 0
    for k in range(1, 11):
        m, lo, hi = ci([d[s][f"neutral:{k}"] for s in common])
        n_false += cls(lo, hi) != "구별 안 됨"
    m, lo, hi = ci([d[s][POS] for s in common])
    pos_cls = cls(lo, hi)
    h = n_false >= 3 or pos_cls != "해로움"
    hold[(frac, mu)] = h
    out["checks"][f"{frac}_{mu}"] = {"n_false": n_false, "pos": [m, lo, hi, pos_cls], "hold": h}
    print(f"  frac {frac:<5} μ {mu:<5}: ④ {n_false}/10 {'🚨 보류' if n_false >= 3 else '통과'} · ⑤ {POS} Δ {m:+.2f} [{lo:+.2f}, {hi:+.2f}] {pos_cls}{' 🚨' if pos_cls != '해로움' else ''}")

print("\n== 판정 — 초점 코드의 Δ 와 빈도 의존")
for mu in (0.0, 0.01):
    print(f"  μ {mu}")
    for code in FOCAL:
        row = {}
        for frac in (0.01, 0.1, 0.5):
            row[frac] = ci([data[(frac, mu)][s][code] for s in common])
        diff = ci([data[(0.5, mu)][s][code] - data[(0.01, mu)][s][code] for s in common])
        if hold.get((0.5, mu)) or hold.get((0.01, mu)):
            v = "판정 보류 — 점검 ④/⑤ 발동"
        elif diff[2] < 0:
            v = "드물 때 더 낫다"
        elif diff[1] > 0:
            v = "흔할 때 더 낫다"
        else:
            v = "빈도 의존 검출 안 됨"
        gone = "흔해지면 이점이 사라진다" if row[0.5][2] < 0 or (row[0.5][1] <= 0 <= row[0.5][2]) else "50% 에서도 이득이 남는다"
        out["focal"][f"{code}_{mu}"] = {"by_frac": {str(k): v2 for k, v2 in row.items()}, "diff": diff, "verdict": v, "at50": gone}
        print(f"    {code}: 1% {row[0.01][0]:+.2f} [{row[0.01][1]:+.2f}, {row[0.01][2]:+.2f}] · 10% {row[0.1][0]:+.2f} [{row[0.1][1]:+.2f}, {row[0.1][2]:+.2f}] · "
              f"50% {row[0.5][0]:+.2f} [{row[0.5][1]:+.2f}, {row[0.5][2]:+.2f}]")
        print(f"       짝 차이 50%−1% {diff[0]:+.2f} [{diff[1]:+.2f}, {diff[2]:+.2f}] → {v} · {gone}")
json.dump(out, open(os.path.join(BASE, "_RESULT_freq.json"), "w"), ensure_ascii=False, indent=1)
print("\n저장: _RESULT_freq.json")
