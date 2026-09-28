"""해부 2 설계용 탐색 — '씨앗 시절' 배경을 몇 틱으로 잡을지. 본실험 기록(시드 1~100)만 읽는다 · 돌연변이 결과는 보지 않는다.

해부 1 과 같은 배경 접시(끝 우세 코드가 WT 인 접시를 시드 순으로 space 20 · mat 20 · energy 60)에서:
  - 틱별 개체 수(중앙 · 5%) — 10% 를 끼워도 표지 계통이 충분히 큰가
  - 틱별 씨앗 rascld · WT 의 개체 비율(중앙) — WT 가 아직 퍼지기 전인가
  - WT 가 처음 상위 1위에 오른 틱(500틱 기록, 중앙 · 5%)
실행: python3 dms2_early_explore.py [petri 폴더, 기본 ~/petri]
"""
import json
import math
import os
import sys

BASE = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1 else "~/petri")
WTS = [("space", "reasccld", 20), ("mat", "racld", 20), ("energy", "reascld", 60)]
TICKS = [500, 1000, 1500, 2000, 2500, 3000, 4000, 5000]


def pct(xs, p):
    xs = sorted(x for x in xs if x == x)
    if not xs:
        return float("nan")
    k = (len(xs) - 1) * p
    f = math.floor(k)
    c = min(f + 1, len(xs) - 1)
    return xs[f] + (xs[c] - xs[f]) * (k - f)


for cond, wt, nbg in WTS:
    zs = []
    for s in range(1, 101):
        z = json.load(open(os.path.join(BASE, "main", f"{cond}_d8_mu0p01_s{s:05d}.json")))
        if z["extinct_at"] < 0 and z["final_counts"] and z["final_counts"][0][0] == wt:
            zs.append(z)
        if len(zs) == nbg:
            break
    print(f"\n== {cond} · WT {wt} · 배경 {len(zs)} (시드 {zs[0]['opts']['seed']}~{zs[-1]['opts']['seed']})")
    first_top = []
    for z in zs:
        t1 = next((t for t, top in z["tops"] if top and top[0][0] == wt), None)
        first_top.append(t1 if t1 is not None else float("nan"))
    print(f"   WT 가 처음 1위인 틱: 중앙 {pct(first_top, .5):,.0f} · 5% {pct(first_top, .05):,.0f} · 95% {pct(first_top, .95):,.0f}")
    for t in TICKS:
        pops, seed_sh, wt_sh = [], [], []
        for z in zs:
            pop = next(r[1] for r in z["samples"] if r[0] == t)
            top = dict(next(tp for tt, tp in z["tops"] if tt == t))
            pops.append(pop)
            seed_sh.append(top.get("rascld", 0) / pop if pop else float("nan"))
            wt_sh.append(top.get(wt, 0) / pop if pop else float("nan"))
        print(f"   틱 {t:>5,}: 개체 중앙 {pct(pops, .5):>6,.0f} (5% {pct(pops, .05):>5,.0f}) · 씨앗 비율 중앙 {pct(seed_sh, .5):.0%} · WT 비율 중앙 {pct(wt_sh, .5):.1%}")
