# -*- coding: utf-8 -*-
"""해부 6 설계용 탐색 — 판정 아님. 해부6 노트 §1c 의 재생성기.
   mat racld 이웃의 침입 Δ 가 μ 0(해부 1) → μ 1%(해부 2A) 로 얼마나 움직였나 — 씨앗만의 일인가, 살아 있는 이웃 전부의 일인가.
   입력: 볼트의 _RESULT_dms_main.json · _RESULT_dms_mu1.json
실행: python3 explore6_shift.py  →  _RESULT_explore6_shift.txt 로 저장해 둔다
"""
import json
import statistics as st

a = {m["code"]: m for m in json.load(open("_RESULT_dms_main.json"))["wts"]["racld"]["mutants"]}
b = {m["code"]: m for m in json.load(open("_RESULT_dms_mu1.json"))["wts"]["racld"]["mutants"]}
rows = sorted(((a[c]["mean"], b[c]["mean"], c) for c in a
               if c in b and not a[c]["lethal"] and not b[c]["lethal"]), reverse=True)
print("두 μ 모두 치명 아닌 이웃 %d개 (전체 %d)" % (len(rows), len(a)))
print("%-10s %8s %8s %8s" % ("코드", "Δ μ0", "Δ μ1", "이동"))
for x, y, c in rows[:15]:
    print("%-10s %+8.2f %+8.2f %+8.2f" % (c, x, y, y - x))
sh = sorted(y - x for x, y, c in rows)
print("이동 중앙값 %+.2f · 사분위 %+.2f ~ %+.2f · 음수 %d/%d"
      % (st.median(sh), sh[len(sh) // 4], sh[3 * len(sh) // 4], sum(1 for v in sh if v < 0), len(sh)))
near = [y - x for x, y, c in rows if x > -2]
print("μ0 에서 Δ > −2 인 이웃 %d개의 이동 중앙값 %+.2f" % (len(near), st.median(near)))
na = json.load(open("_RESULT_dms_main.json"))["wts"]["racld"]["neutral"]
nb = json.load(open("_RESULT_dms_mu1.json"))["wts"]["racld"]["neutral"]
print("가짜 돌연변이(WT) 이동:", " ".join("%+.2f" % (q["mean"] - p["mean"]) for p, q in zip(na, nb)))
