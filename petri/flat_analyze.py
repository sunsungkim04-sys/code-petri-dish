"""해부 3B 판정 — 평평함 직접 재기. 사전등록 해부3-사전등록-2026-09-16.md §3 · §4 · §5. 결과 전에 썼다.

  같은 mat 배경 20접시 · 같은 설정(끼우기 10% · μ 0 · 3,000틱)에서 잰 두 WT 의 한 글자 이웃을 비교한다.
    racld  = 해부 1 자료(dms_main)      ·  rascld = 이번 발사(dms_flat)
  평평함(주 판정) = 그 WT 의 이웃 전부의 Δ 평균(배경마다 하나) → 배경 재추출 부트스트랩 2,000회
     (씨앗 평평함 − racld 평평함) 상한 < 0 → 'racld 가 더 평평하다' · 하한 > 0 → '씨앗이 더 평평하다' · 0 품음 → 차이 검출 안 됨
  덧붙임: 치명 이웃 비율(Wilson) · 분류 개수 · 일관성 확인(씨앗의 이웃 racld 의 Δ 가 해부 1 의 +0.39 와 부호 반대 · 비슷한 크기인가)
  점검 ①②③ 배경 제외 · ④ 가짜 돌연변이 10개 중 3개 이상이면 보류 · ⑤ 양성 대조(rascl = d 빠뜨림) 해로움
실행: python3 flat_analyze.py [petri 폴더, 기본 ~/petri]  →  _RESULT_flat.json
"""
import glob
import json
import math
import os
import random
import sys
from collections import Counter

BASE = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1 else "~/petri")
COND = "mat"
SETS = {"racld": "dms_main", "rascld": "dms_flat"}
N_BG, N_BOOT = 20, 2000
RNG = random.Random(20260921)


def pct(xs, p):
    xs = sorted(xs)
    k = (len(xs) - 1) * p
    f = math.floor(k)
    c = min(f + 1, len(xs) - 1)
    return xs[f] + (xs[c] - xs[f]) * (k - f)


def wilson(k, n, z=1.96):
    if n == 0:
        return (float("nan"), float("nan"))
    p = k / n
    c = (p + z * z / (2 * n)) / (1 + z * z / n)
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n)
    return (max(0.0, c - h), min(1.0, c + h))


def s_of(arm):
    r = arm["series"]
    n0, m0 = r[0][1], r[0][2]
    nT, mT = r[-1][1], r[-1][2]
    return math.log((mT + 0.5) / m0) - math.log((nT - mT + 0.5) / (n0 - m0))


def load(wt, folder):
    out, probs, ext = {}, [], Counter()
    for f in sorted(glob.glob(os.path.join(BASE, folder, f"dms_{COND}_{wt}_s*_all.json"))):
        z = json.load(open(f))
        why = []
        if not z["checksum_match"]:
            why.append("② 체크섬")
        if not z["fidelity"]["same"]:
            why.append("① 복제 충실도")
        if abs(z.get("mu_assay", 0.0)) > 1e-12 or abs(z["frac"] - 0.1) > 1e-12 or z["ticks"] != 3000:
            why.append("설정 다름")
        if any(not a["cons_ok"] for a in z["arms"]):
            why.append("③ 재료 보존")
        if len(z["arms"]) != 11 + z["n_mutants_total"]:
            why.append(f"팔 수 {len(z['arms'])}")
        if why:
            probs.append((z["seed"], why))
            continue
        wt_arms = [a for a in z["arms"] if a["kind"] in ("ref", "neutral")]
        s_wt = {a["labels"][0]: s_of(a) for a in wt_arms}
        base = sum(s_wt.values()) / len(s_wt)
        d = {}
        for a in z["arms"]:
            if a["kind"] == "mut":
                d[a["key"]] = s_of(a) - base
                if a["series"][-1][2] == 0:
                    ext[a["key"]] += 1
        for lab, v in s_wt.items():
            if lab == "wt":
                continue
            others = [x for k, x in s_wt.items() if k != lab]
            d[lab] = v - sum(others) / len(others)
        out[z["seed"]] = d
    return out, probs, ext


D, E, halt = {}, {}, False
for wt, folder in SETS.items():
    d, probs, ext = load(wt, folder)
    D[wt], E[wt] = d, ext
    print(f"{wt:>7} ({folder}): 쓴 배경 {len(d)}/{N_BG}" + ("" if not probs else " · 뺀 배경 " + " · ".join(f"시드 {s}:{'/'.join(w)}" for s, w in probs)))
    if len(d) != N_BG:
        halt = True
common = sorted(set(D["racld"]) & set(D["rascld"]))
print(f"짝지은 배경 {len(common)} · 시드 {common[:5]}…")
if halt or len(common) != N_BG:
    print("🚨 두 WT 의 배경이 같지 않거나 다 모이지 않았다 — 판정하지 않는다")
    sys.exit(1)
idxs = [[RNG.randrange(len(common)) for _ in range(len(common))] for _ in range(N_BOOT)]


def ci(vals):
    n = len(vals)
    bm = [sum(vals[i] for i in ix) / n for ix in idxs]
    return sum(vals) / n, pct(bm, .025), pct(bm, .975)


def cls(lo, hi):
    return "해로움" if hi < 0 else ("이로움" if lo > 0 else "구별 안 됨")


print("\n== 점검 ④⑤")
hold = False
for wt in SETS:
    d = D[wt]
    n_false = sum(cls(*ci([d[s][f"neutral:{k}"] for s in common])[1:]) != "구별 안 됨" for k in range(1, 11))
    pos = "rascl" if wt == "rascld" else "racl"
    m, lo, hi = ci([d[s][pos] for s in common])
    pc = cls(lo, hi)
    h = n_false >= 3 or pc != "해로움"
    hold = hold or h
    print(f"  {wt:>7}: ④ {n_false}/10 {'🚨 보류' if n_false >= 3 else '통과'} · ⑤ {pos} Δ {m:+.2f} [{lo:+.2f}, {hi:+.2f}] {pc}{' 🚨' if pc != '해로움' else ''}")

print("\n== 평평함 — 이웃 전부의 Δ 평균")
flat, res = {}, {}
for wt in SETS:
    keys = [k for k in D[wt][common[0]] if not k.startswith("neutral:")]
    per_bg = [sum(D[wt][s][k] for k in keys) / len(keys) for s in common]
    flat[wt] = per_bg
    m, lo, hi = ci(per_bg)
    n_leth = sum(E[wt][k] == len(common) for k in keys)
    wl = wilson(n_leth, len(keys))
    classes = Counter(cls(*ci([D[wt][s][k] for s in common])[1:]) for k in keys)
    res[wt] = {"n_neighbors": len(keys), "flatness": [m, lo, hi], "lethal": [n_leth, len(keys), wl[0], wl[1]], "classes": dict(classes),
               "quartiles": [pct([sum(D[wt][s][k] for s in common) / len(common) for k in keys], q) for q in (.25, .5, .75)]}
    print(f"  {wt:>7}: 이웃 {len(keys)} · 평평함 {m:+.3f} [{lo:+.3f}, {hi:+.3f}] · 치명 {n_leth}/{len(keys)} = {n_leth / len(keys):.0%} [{wl[0]:.0%}, {wl[1]:.0%}] · "
          f"이웃 Δ 사분위 {res[wt]['quartiles'][0]:+.2f} / {res[wt]['quartiles'][1]:+.2f} / {res[wt]['quartiles'][2]:+.2f} · 분류 {dict(classes)}")
diff = ci([a - b for a, b in zip(flat["rascld"], flat["racld"])])
v = "판정 보류 — 점검 ④/⑤ 발동" if hold else ("racld 가 더 평평하다" if diff[2] < 0 else ("씨앗이 더 평평하다" if diff[1] > 0 else "차이 검출 안 됨"))
print(f"  씨앗 − racld: {diff[0]:+.3f} [{diff[1]:+.3f}, {diff[2]:+.3f}] → {v}")

m, lo, hi = ci([D["rascld"][s]["racld"] for s in common])
print(f"\n== 일관성 확인: 씨앗의 이웃 racld 의 Δ {m:+.3f} [{lo:+.3f}, {hi:+.3f}] (해부 1 에서 racld 의 이웃 rascld 는 +0.39)")
json.dump({"n_bg": len(common), "hold": hold, "wt": res, "diff": diff, "verdict": v, "mirror": [m, lo, hi]},
          open(os.path.join(BASE, "_RESULT_flat.json"), "w"), ensure_ascii=False, indent=1)
print("저장: _RESULT_flat.json")
