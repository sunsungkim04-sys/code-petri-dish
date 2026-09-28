"""해부 2A 판정 — 복제 오류를 켠 해부를 해부 1(μ = 0)과 같은 배경끼리 짝지어 비교. 사전등록 해부2-사전등록-2026-09-15.md §2 · §4. 결과 전에 썼다.

  Δ 는 dms_analyze.py 와 같다(같은 배경 WT 팔 11개 s 평균을 뺀다). 배경 시드가 같은 것끼리 짝.
  초점 1 — mat 의 씨앗 rascld(racld 에 s 끼어듦): Δ(μ0) · Δ(μ1) · 짝 차이 Δ(μ1) − Δ(μ0), 배경 재추출 부트스트랩 2,000회 95%
    Δ(μ1) 상한 < 0                         → '복제 오류를 켜면 racld 가 씨앗보다 낫다(뒤집힌다)'
    아니고 짝 차이 상한 < 0                 → '복제 오류가 씨앗의 이점을 줄이지만 뒤집지는 않는다'
    짝 차이 하한 > 0                        → '복제 오류가 오히려 씨앗 쪽으로 민다'
    짝 차이가 0 을 품음                     → '복제 오류로는 설명 안 됨'
  초점 2 — energy 에서 μ0 이로움 7개가 μ1 에서도 이로움인가(개수) · 코드별 짝 차이
  전체 — WT 마다 분류 개수 μ0 대 μ1 · 분류가 바뀐 코드
  보류 — μ1 판정(_RESULT_dms_mu1.json)의 점검 ④ 발동 · ⑤ 실패면 그 WT 초점 판정을 보류한다
실행: python3 dms2_compare.py [petri 폴더, 기본 ~/petri]
"""
import glob
import json
import math
import os
import random
import sys
from collections import Counter

BASE = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1 else "~/petri")
WTS = [("space", "reasccld"), ("mat", "racld"), ("energy", "reascld")]
N_BOOT = 2000
RNG = random.Random(20260917)


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


def per_bg(folder, cond, wt):
    out = {}
    for f in sorted(glob.glob(os.path.join(BASE, folder, f"dms_{cond}_{wt}_s*_all.json"))):
        z = json.load(open(f))
        wt_arms = [a for a in z["arms"] if a["kind"] in ("ref", "neutral")]
        base = sum(s_of(a) for a in wt_arms) / len(wt_arms)
        out[z["seed"]] = {a["key"]: s_of(a) - base for a in z["arms"] if a["kind"] == "mut"}
    return out


def ci(vals, idxs):
    n = len(vals)
    bm = [sum(vals[i] for i in idx) / n for idx in idxs]
    return sum(vals) / n, pct(bm, .025), pct(bm, .975)


R0 = json.load(open(os.path.join(BASE, "_RESULT_dms_main.json")))
R1 = json.load(open(os.path.join(BASE, "_RESULT_dms_mu1.json")))
out = {}
for cond, wt in WTS:
    W0, W1 = R0["wts"][wt], R1["wts"][wt]
    hold = bool(W1["hold"]) or not R1["positive_all_harmful"]
    d0, d1 = per_bg("dms_main", cond, wt), per_bg("dms_mu1", cond, wt)
    seeds = sorted(set(d0) & set(d1))
    n = len(seeds)
    idxs = [[RNG.randrange(n) for _ in range(n)] for _ in range(N_BOOT)]
    c0 = Counter(r["cls"] for r in W0["mutants"])
    c1 = Counter(r["cls"] for r in W1["mutants"])
    l0 = sum(r["lethal"] for r in W0["mutants"])
    l1 = sum(r["lethal"] for r in W1["mutants"])
    print(f"\n==================== {cond} · WT {wt} · 짝지은 배경 {n}{' · 🚨 μ1 판정 보류(점검 ④/⑤) — 초점 판정 안 함' if hold else ''}")
    print(f"  분류 μ0: 해로움 {c0['해로움']}(치명 {l0}) · 구별 안 됨 {c0['구별 안 됨']} · 이로움 {c0['이로움']}")
    print(f"  분류 μ1: 해로움 {c1['해로움']}(치명 {l1}) · 구별 안 됨 {c1['구별 안 됨']} · 이로움 {c1['이로움']}")
    cls0 = {r["code"]: r for r in W0["mutants"]}
    cls1 = {r["code"]: r for r in W1["mutants"]}
    changed = [(k, cls0[k]["cls"], cls1[k]["cls"]) for k in cls0 if k in cls1 and cls0[k]["cls"] != cls1[k]["cls"]]
    print(f"  분류가 바뀐 코드 {len(changed)}: " + (" · ".join(f"{k} {a}→{b} ({cls0[k]['mean']:+.2f}→{cls1[k]['mean']:+.2f})" for k, a, b in changed[:20]) if changed else "없음"))
    res = {"n_bg": n, "hold": hold, "counts_mu0": dict(c0), "counts_mu1": dict(c1), "lethal": [l0, l1], "changed": changed}
    if cond == "mat" and not hold:
        x0 = [d0[s]["rascld"] for s in seeds]
        x1 = [d1[s]["rascld"] for s in seeds]
        dd = [b - a for a, b in zip(x0, x1)]
        m0, lo0, hi0 = ci(x0, idxs)
        m1, lo1, hi1 = ci(x1, idxs)
        md, lod, hid = ci(dd, idxs)
        if hi1 < 0:
            v = "복제 오류를 켜면 racld 가 씨앗보다 낫다(뒤집힌다)"
        elif hid < 0:
            v = "복제 오류가 씨앗의 이점을 줄이지만 뒤집지는 않는다"
        elif lod > 0:
            v = "복제 오류가 오히려 씨앗 쪽으로 민다"
        else:
            v = "복제 오류로는 설명 안 됨"
        print(f"  초점 1 씨앗 rascld: Δ μ0 {m0:+.3f} [{lo0:+.3f}, {hi0:+.3f}] · Δ μ1 {m1:+.3f} [{lo1:+.3f}, {hi1:+.3f}] · 짝 차이 {md:+.3f} [{lod:+.3f}, {hid:+.3f}] → {v}")
        res["focal_rascld"] = {"mu0": [m0, lo0, hi0], "mu1": [m1, lo1, hi1], "diff": [md, lod, hid], "verdict": v}
    if cond == "energy" and not hold:
        good0 = [r["code"] for r in W0["mutants"] if r["cls"] == "이로움"]
        stay = [k for k in good0 if cls1[k]["cls"] == "이로움"]
        print(f"  초점 2 μ0 이로움 {len(good0)}개 중 μ1 에서도 이로움 {len(stay)}")
        rows = []
        for k in good0:
            dd = [d1[s][k] - d0[s][k] for s in seeds]
            md, lod, hid = ci(dd, idxs)
            rows.append({"code": k, "mu0": cls0[k]["mean"], "mu1": cls1[k]["mean"], "cls1": cls1[k]["cls"], "diff": [md, lod, hid]})
            print(f"    {k}: μ0 {cls0[k]['mean']:+.2f} → μ1 {cls1[k]['mean']:+.2f} {cls1[k]['cls']} · 짝 차이 {md:+.2f} [{lod:+.2f}, {hid:+.2f}]")
        res["focal_energy"] = {"n_good0": len(good0), "n_stay": len(stay), "rows": rows}
    new_good = sorted((r for r in W1["mutants"] if r["cls"] == "이로움"), key=lambda r: -r["mean"])
    print("  μ1 이로움 전부: " + (" · ".join(f"{r['code']}({'/'.join(r['labels'])}) {r['mean']:+.2f}[{r['lo']:+.2f},{r['hi']:+.2f}]" for r in new_good) if new_good else "없음"))
    out[wt] = res
json.dump(out, open(os.path.join(BASE, "_RESULT_dms2_compare.json"), "w"), ensure_ascii=False, indent=1)
print("\n저장: _RESULT_dms2_compare.json")
