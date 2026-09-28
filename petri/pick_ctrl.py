# -*- coding: utf-8 -*-
"""해부 5 설계용 — 대조 짝으로 쓸 복제자 고르기. 규칙: mono_nb(μ 0 · 시드 1~10 · 4,000틱)에서
   끝 개체 수 중앙값 ≥ 10 인 코드를 길이별로 중앙값 큰 순서로 세우고 위에서 고른다(두 WT 는 뺀다)."""
import glob
import json
import os
import statistics as st
from collections import defaultdict

rows = defaultdict(list)
for f in glob.glob(os.path.expanduser("~/petri/mono_nb/*.json")):
    z = json.load(open(f))
    rows[z["code"]].append(z["samples"][-1][1] if z["extinct_at"] < 0 else 0)

by_len = defaultdict(list)
for code, pops in rows.items():
    if len(pops) < 3:
        continue
    m = st.median(pops)
    if m >= 10:
        by_len[len(code)].append((m, code, len(pops)))

print("복제하는 코드 수 (길이별) — mono_nb %d 코드 중" % len(rows))
for L in sorted(by_len):
    print("  길이 %d: %3d개 · 상위 6 = %s" % (L, len(by_len[L]),
          ", ".join("%s(%.0f)" % (c, m) for m, c, _ in sorted(by_len[L], reverse=True)[:6])))
WT = {"racld", "rascld"}
print()
for L in (5, 6):
    top = [(m, c) for m, c in sorted(((m, c) for m, c, _ in by_len.get(L, [])), reverse=True) if c not in WT]
    print("  길이 %d 상위(WT 제외): %s" % (L, ", ".join("%s(%.0f)" % (c, m) for m, c in top[:5])))
