# -*- coding: utf-8 -*-
"""W1 재계수 — 🟡 사후(판정 아님). 설계검토-2026-09-17.md W1.
해부 5 §6a 의 '재료 총량은 어느 칸에서도 제약이 아니다(밀도 4 도 89% 남음)' 는 열한 글자 전부를 분모에 넣은 값이다.
코드가 **쓰는 글자**만 세면 어떤가 — mono.js v1.1.0 기록의 free_by_letter_end(끝 시점 글자별 자유 재료)로 다시 센다.
  자유 몫(글자 k) = free_by_letter_end[k] ÷ (mat_total0 / 11)   (처음 배치는 글자마다 같은 기대 수)
  코드 안 글자 = 코드에 든 글자 종류 · 밖 글자 = 나머지 · 전체 = free_end / mat_total0 (노트 값)
  자료: dens/mat_d{04,08,14,20,32}/mono_mat_<code>_mu0_s*.json (해부 5A · μ 0 · 살아남은 접시)
실행: python3 w1_recount.py [petri 폴더]  →  표준 출력
"""
import glob
import json
import os
import sys
from collections import defaultdict

BASE = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1 else "~/petri")
LET = "nsracldjehx"
acc = defaultdict(list)
for f in glob.glob(os.path.join(BASE, "dens", "mat_d*", "mono_mat_*_mu0_s*.json")):
    z = json.load(open(f))
    if z["extinct_at"] >= 0:
        continue
    per = z["mat_total0"] / 11
    fr = [x / per for x in z["free_by_letter_end"]]
    code = z["code"]
    inside = [fr[LET.index(ch)] for ch in set(code)]
    outside = [fr[i] for i, ch in enumerate(LET) if ch not in code]
    acc[(z["density"], code)].append((sum(inside) / len(inside), min(inside), sum(outside) / len(outside), z["free_end"] / z["mat_total0"], z["ticks_run"]))
print("== W1 재계수 (🟡 사후 · 판정 아님) — 끝 시점 자유 재료 몫: 코드가 쓰는 글자만 대 전체")
print("   %-4s %-8s %3s | %-16s | %-6s | %-12s" % ("밀도", "코드", "n", "안 글자 평균(최소)", "밖 글자", "전체(노트 값)"))
for k in sorted(acc):
    v = acc[k]
    n = len(v)
    m = lambda j: sum(x[j] for x in v) / n
    print("   d%-3d %-8s %3d | %.2f (%.2f)        | %.2f   | %.2f" % (k[0], k[1], n, m(0), m(1), m(2), m(3)))
print("\n읽기(🟡): 전체 몫은 코드가 안 쓰는 글자 5~7종이 부풀린 값이다. 코드가 쓰는 글자로 보면 밀도 8 은 85% 가 몸에 묶여 있고(공급 천장),")
print("          밀도 4 는 개체가 적어 오히려 76% 가 남는다. '재료 총량은 어느 칸에서도 제약이 아니다' 는 글자별로는 성립하지 않는다.")
