#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""🟡 탐색 (판정 아님) — 해부 12 P 가 이미 잰 이웃 528개의 Δ 를 **복사 고리 틱 수**로 설명할 수 있나.

해부 11 C 의 고리 사다리는 `racld` 한 가족 안에서 네 칸(2 · 3 · 4 · 5)이었다. 여기서는 **같은 자료**를
다시 읽어, 고리 길이를 코드마다 계산하고 Δ 와의 관계를 본다. 새 실행 없음 · 판정 아님 —
해부 13(둘째 코드 가족)의 예측을 세우기 위한 것이고, 사후 분석임을 사전등록에 적는다.

고리 틱 수 정의 (`sim.js` step() 의 의미 그대로)
  l 은 **바로 앞의 s 다음 자리**로 돌아간다. s 가 없으면 코드 맨 앞.
  고리 = [되돌아갈 자리 … l] 구간. 그 안에서 **x 는 틱을 안 쓴다**(step 이 건너뜀).
  고리 안에 c 가 없으면 복사 고리가 아니다 → 따로 센다.
실행: python3 loop13_explore.py [base]
"""
import glob
import json
import math
import os
import statistics
import sys

BASE = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1 else "~/petri")
ROOT = os.environ.get("PETRI_ROOT", "inv12")
NCODES = ["racld", "rascld", "rsacld", "racldx"]


def loop_ticks(code):
    """(고리 틱 수, 고리가 c 를 품나, l 이 있나). 없으면 (None, …)."""
    if "l" not in code or "c" not in code:
        return None, False, "l" in code
    li = code.index("l")
    si = code.rfind("s", 0, li)
    start = si + 1 if si >= 0 else 0
    if start > li:                      # l 이 s 바로 뒤 — 제자리 고리
        return None, False, True
    seg = code[start:li + 1]
    return sum(1 for ch in seg if ch != "x"), ("c" in seg), True


def s_of(arm):
    r = arm["series"]
    n0, m0 = r[0][1], r[0][2]
    nT, mT = r[-1][1], r[-1][2]
    return math.log((mT + 0.5) / m0) - math.log((nT - mT + 0.5) / (n0 - m0))


D = {}
for code in NCODES:
    for f in sorted(glob.glob(os.path.join(BASE, "%snbr_%s" % (ROOT, code), "mu0", "dms_mat_racld_s*_nbr_%s.json" % code))):
        z = json.load(open(f))
        wt = [a for a in z["arms"] if a["kind"] in ("ref", "neutral")]
        base = sum(s_of(a) for a in wt) / len(wt)
        for a in z["arms"]:
            if a["kind"] == "mut":
                D.setdefault(a["key"], []).append((s_of(a) - base, a["stopped"] >= 0 and a["series"][-1][2] == 0))

print("== 🟡 탐색 — 이웃 Δ 를 복사 고리 틱 수로 설명할 수 있나 (해부 12 P 자료 재분석 · 판정 아님)")
print("  서로 다른 이웃 코드 %d · 배경당 반복 %d" % (len(D), len(next(iter(D.values())))))

rows = []
for key, v in D.items():
    d = statistics.median(x[0] for x in v)
    dead = sum(1 for x in v if x[1]) * 2 > len(v)
    L, has_c, has_l = loop_ticks(key)
    rows.append({"key": key, "d": d, "dead": dead, "L": L, "c": has_c, "l": has_l, "len": len(key)})

nol = [r for r in rows if not r["l"]]
noc = [r for r in rows if r["l"] and (r["L"] is None or not r["c"])]
ok = [r for r in rows if r["L"] is not None and r["c"]]
print("  고리 계산: 복사 고리 있음 %d · l 이 없음 %d · 고리에 c 없음/제자리 %d" % (len(ok), len(nol), len(noc)))
print("  l 없는 이웃 Δ 중앙값 %+.2f (불임 %d/%d) · 고리에 c 없는 이웃 %+.2f (불임 %d/%d)"
      % (statistics.median([r["d"] for r in nol]) if nol else float("nan"), sum(r["dead"] for r in nol), len(nol),
         statistics.median([r["d"] for r in noc]) if noc else float("nan"), sum(r["dead"] for r in noc), len(noc)))

print("\n  고리 틱 | 코드 수 | 불임 | Δ 중앙값 | Δ 사분위 | 살아남은 것의 Δ 중앙값")
by = {}
for r in ok:
    by.setdefault(r["L"], []).append(r)
for L in sorted(by):
    v = by[L]
    ds = sorted(x["d"] for x in v)
    alive = [x["d"] for x in v if not x["dead"]]
    q1 = ds[len(ds) // 4]; q3 = ds[3 * len(ds) // 4]
    print("  %7d | %7d | %4d | %+8.2f | %+.2f / %+.2f | %s"
          % (L, len(v), sum(x["dead"] for x in v), statistics.median(ds), q1, q3,
             ("%+.2f (n %d)" % (statistics.median(alive), len(alive))) if alive else "—"))

alive = [r for r in ok if not r["dead"]]
if len(alive) > 3:
    xs = [r["L"] for r in alive]; ys = [r["d"] for r in alive]
    mx = sum(xs) / len(xs); my = sum(ys) / len(ys)
    sxx = sum((x - mx) ** 2 for x in xs)
    b = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sxx
    res = [y - (my + b * (x - mx)) for x, y in zip(xs, ys)]
    se = math.sqrt(sum(r * r for r in res) / (len(xs) - 2) / sxx)
    print("\n  살아남은 이웃 %d개에서 Δ ~ 고리 틱: 기울기 %+.3f ± %.3f (SE) · 고리 범위 %d~%d"
          % (len(alive), b, se, min(xs), max(xs)))
    print("  🟡 사후 분석이다. 해부 11 C 의 사다리(틱당 −0.510)와 부호를 견주기만 한다.")
    lens = sorted({r["len"] for r in alive})
    for L0 in lens:
        g = [r for r in alive if r["len"] == L0]
        if len(g) > 3:
            gx = [r["L"] for r in g]; gy = [r["d"] for r in g]
            if len(set(gx)) > 1:
                gmx = sum(gx) / len(gx); gmy = sum(gy) / len(gy)
                gb = sum((x - gmx) * (y - gmy) for x, y in zip(gx, gy)) / sum((x - gmx) ** 2 for x in gx)
                print("    길이 %d 만: n %d · 기울기 %+.3f · 고리 %s" % (L0, len(g), gb, sorted(set(gx))))

print("\n  🚩 읽기 주의: 이웃은 고리 길이 말고도 다르다(글자 종류 · 자리). 이 기울기는 **교란된 값**이고,")
print("     해부 13 은 이것을 예측의 출처로만 쓰고 새 자료에서 판정한다.")
