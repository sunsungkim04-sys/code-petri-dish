"""이긴 코드 해부 사후 진단 — 설계 노트 밖, 판정에 쓰지 않는다. 판정 결과를 믿어도 되는지와 읽는 법을 보려고 본 것.

  (1) 소멸 시각 — '치명' 은 *못 번식해서* 사라진 것인가 *밀려서* 사라진 것인가.
      양성 대조(d 빠뜨림 = 자식을 떼지 못함)는 수명(300~600틱) 안에 사라져야 한다. 치명 돌연변이가 그보다 한참 늦게 사라지면
      번식은 하는데 경쟁에서 밀린 것이다. 비-WT 코드만 망가뜨리는 계측기 결함이면 양성 대조와 같은 시각에 사라질 것이다.
  (2) 가장 덜 해로운 돌연변이 10개 · 무늬 x 끼어듦 전부 · 이로움 전부
실행: python3 dms_diagnostics.py [petri 폴더, 기본 ~/petri]  (먼저 dms_analyze.py 를 돌려 _RESULT_dms_main.json 이 있어야 한다)
"""
import glob
import json
import math
import os
import sys
from collections import defaultdict

BASE = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1 else "~/petri")
RES = json.load(open(os.path.join(BASE, "_RESULT_dms_main.json")))
WTS = [("space", "reasccld"), ("mat", "racld"), ("energy", "reascld")]


def pct(xs, p):
    xs = sorted(xs)
    if not xs:
        return float("nan")
    k = (len(xs) - 1) * p
    f = math.floor(k)
    c = min(f + 1, len(xs) - 1)
    return xs[f] + (xs[c] - xs[f]) * (k - f)


for cond, wt in WTS:
    W = RES["wts"][wt]
    muts = {r["code"]: r for r in W["mutants"]}
    stops = defaultdict(list)
    for f in sorted(glob.glob(os.path.join(BASE, "dms_main", f"dms_{cond}_{wt}_s*_all.json"))):
        z = json.load(open(f))
        for a in z["arms"]:
            if a["kind"] == "mut":
                stops[a["key"]].append(a["stopped"])
    pos = W["positive"]["code"]
    pos_stop = [s for s in stops[pos] if s >= 0]
    print(f"\n==================== {cond} · WT {wt} · 배경 {W['n_bg']}")
    print(f"(1) 양성 대조 {pos} 소멸 틱: 중앙 {pct(pos_stop, .5):,.0f} · 범위 {min(pos_stop):,}~{max(pos_stop):,} ({len(pos_stop)}/{len(stops[pos])})")
    lethal = [r for r in muts.values() if r["lethal"] and r["code"] != pos]
    med = [(pct([s for s in stops[r["code"]] if s >= 0], .5), r) for r in lethal]
    late = [m for m, _ in med if m > max(pos_stop)]
    print(f"    치명 돌연변이 {len(lethal)}개(양성 대조 제외)의 소멸 틱 중앙값 분포: 5% {pct([m for m, _ in med], .05):,.0f} · 중앙 {pct([m for m, _ in med], .5):,.0f} · 95% {pct([m for m, _ in med], .95):,.0f} · "
          f"양성 대조의 가장 늦은 소멸보다 늦게 사라진 것 {len(late)}/{len(lethal)}")
    fast = sorted(med, key=lambda x: x[0])[:6]
    print("    가장 빨리 사라진 치명: " + " · ".join(f"{r['code']}({'/'.join(r['labels'])}) {m:,.0f}" for m, r in fast))
    slow = sorted(med, key=lambda x: -x[0])[:6]
    print("    가장 늦게 사라진 치명: " + " · ".join(f"{r['code']}({'/'.join(r['labels'])}) {m:,.0f}" for m, r in slow))
    harm = sorted((r for r in muts.values() if r["cls"] == "해로움"), key=lambda r: -r["mean"])[:10]
    print("(2) 가장 덜 해로운 10: " + " · ".join(f"{r['code']}({'/'.join(r['labels'])}) {r['mean']:+.2f} 소멸 {r['extinct']}/{r['n']}" for r in harm))
    xins = sorted((r for r in muts.values() if any(l.startswith("ins:") and l.endswith(":x") for l in r["labels"])), key=lambda r: min(int(l.split(":")[1]) for l in r["labels"]))
    print("    무늬 x 끼어듦: " + " · ".join(f"{r['code']} {r['mean']:+.2f} 소멸 {r['extinct']}/{r['n']} 틱 {pct([s for s in stops[r['code']] if s >= 0], .5):,.0f}" for r in xins))
    good = sorted((r for r in muts.values() if r["cls"] == "이로움"), key=lambda r: -r["mean"])
    print("    이로움: " + (" · ".join(f"{r['code']}({'/'.join(r['labels'])}) {r['mean']:+.2f}[{r['lo']:+.2f},{r['hi']:+.2f}]" for r in good) if good else "없음"))
    undecided = [r for r in muts.values() if r["cls"] == "구별 안 됨"]
    if undecided:
        print("    구별 안 됨: " + " · ".join(f"{r['code']}({'/'.join(r['labels'])}) {r['mean']:+.2f}[{r['lo']:+.2f},{r['hi']:+.2f}]" for r in undecided))
