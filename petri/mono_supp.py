"""해부 4B 보강 — 판정 아님. 사전등록 규칙이 '보류'(시드에 따라 갈린 이웃 14개 > 5)를 냈을 때의 보강 측정.

  같은 이웃들을 시드 10개로 늘려 다시 키우고(μ 0 · 4,000틱), 이진 분류 대신 **중앙값**으로 본다.
  (1) 중앙 개체 ≥ 10 으로 분류한 U_dead — 이진 분류의 흔들림을 줄인 판
  (2) 애매한 이웃 수(중앙 개체가 5~20 사이)
  (3) 등급 부하 U_graded = 1 − Σ w·(그 이웃 중앙 개체 / WT 중앙 개체, 1 로 자름) — 죽고 사는 이분법 대신 연속값
  무게 w 는 sim.js 오류 규칙 그대로: 바뀜 1/11 · 빠뜨림 1/3 · 끼어듦 1/33
실행: python3 mono_supp.py [petri 폴더]
"""
import glob
import json
import math
import os
import sys
from collections import defaultdict

BASE = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1 else "~/petri")
LET = "nsracldjehx"
W = {"sub": 1 / 11, "del": 1 / 3, "ins": 1 / 33}
CODES = ["racld", "rascld"]


def med(xs):
    xs = sorted(xs)
    n = len(xs)
    return xs[n // 2] if n % 2 else (xs[n // 2 - 1] + xs[n // 2]) / 2


pops = defaultdict(list)
for f in sorted(glob.glob(os.path.join(BASE, "mono_nb", "mono_mat_*.json"))):
    z = json.load(open(f))
    pops[z["code"]].append(z["samples"][-1][1])
print(f"코드 {len(pops)} · 시드 수 중앙 {med([len(v) for v in pops.values()]):.0f}")


def neighbours(wt):
    out = {}
    def add(key, label):
        if key == wt or not (2 <= len(key) <= 48):
            return
        out.setdefault(key, []).append(label)
    for p in range(len(wt)):
        for ch in LET:
            if ch != wt[p]:
                add(wt[:p] + ch + wt[p + 1:], "sub")
    for p in range(len(wt)):
        add(wt[:p] + wt[p + 1:], "del")
    for p in range(len(wt) + 1):
        for ch in LET:
            add(wt[:p] + ch + wt[p:], "ins")
    return out


out = {}
for wt in CODES:
    nb = neighbours(wt)
    wt_med = med(pops[wt])
    tot = sum(sum(W[l] for l in labs) for labs in nb.values())
    dead = graded = 0.0
    n_dead = ambiguous = 0
    for k, labs in nb.items():
        w = sum(W[l] for l in labs)
        m = med(pops[k]) if k in pops else 0
        if m < 10:
            dead += w
            n_dead += 1
        if 5 <= m <= 20:
            ambiguous += 1
        graded += w * min(1.0, m / wt_med)
    out[wt] = {"wt_median_pop": wt_med, "n": len(nb), "n_dead": n_dead, "U_dead_median": dead / tot,
               "U_graded": 1 - graded / tot, "ambiguous": ambiguous}
    print(f"  {wt:>7}: WT 중앙 개체 {wt_med:.0f} · 이웃 {len(nb)} · 중앙 <10 인 이웃 {n_dead} ({n_dead / len(nb):.0%}) · "
          f"U_dead(중앙) {dead / tot:.3f} · U_graded {1 - graded / tot:.3f} · 애매(5~20) {ambiguous}")
a, b = out["racld"], out["rascld"]
print(f"  → U_dead 차이(racld − 씨앗) {a['U_dead_median'] - b['U_dead_median']:+.3f} · U_graded 차이 {a['U_graded'] - b['U_graded']:+.3f}")
print("  (🟡 보강 측정 — 사전등록 판정은 '보류' 그대로다)")
json.dump(out, open(os.path.join(BASE, "_RESULT_mono_supp.json"), "w"), ensure_ascii=False, indent=1)
