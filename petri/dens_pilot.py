# -*- coding: utf-8 -*-
"""해부 5 설계용 탐색 — 판정 아님.
   ① 해부 4A 기록(μ 0 · 밀도 8)에서 개체당 출산율·점유율을 본다 → 무엇이 제약인가
   ② 밀도 4·8·14·20·32 와 space 에서 접시가 사는지, 재료가 바닥나는지 짧게 본다(시드 1~3)"""
import glob
import json
import os
import statistics as st

DISH = 10732

print("== ① 해부 4A 기록 다시 보기 (μ 0 · 밀도 8 · 시드 1~30 · 20,000틱)")
print("   %-8s %8s %8s %9s %10s %12s" % ("코드", "끝 개체", "점유%", "자유재료%", "개체당재료", "개체당출산/천틱"))
for code in ("racld", "rascld"):
    pops, occ, birth = [], [], []
    for s in range(1, 31):
        f = "mono_curve/mono_mat_%s_mu0_s%05d.json" % (code, s)
        if not os.path.exists(f):
            continue
        z = json.load(open(f))
        sm = z["samples"]
        last = sm[-1]
        pops.append(last[1]); occ.append(100.0 * last[1] / DISH)
        # 정상 상태(1만~2만틱)의 개체당 출산율
        a = [r for r in sm if r[0] == 10000][0]
        db, dt = last[4] - a[4], last[0] - a[0]
        mid = st.mean(r[1] for r in sm if 10000 <= r[0] <= 20000)
        birth.append(1000.0 * db / (dt * mid))
    print("   %-8s %8.0f %8.1f %9s %10s %12.2f"
          % (code, st.median(pops), st.median(occ), "-", "-", st.median(birth)))

print("\n   (자유재료 · 개체당재료는 v1.1.0 으로 다시 돌린 시드 1~3 에서)")
for code in ("racld", "rascld"):
    fr, per = [], []
    for s in (1, 2, 3):
        f = "mono_regress5/mono_mat_%s_mu0_s%05d.json" % (code, s)
        z = json.load(open(f))
        fr.append(100.0 * z["free_end"] / z["mat_total0"])
        per.append(z["tied_up_end"] / z["samples"][-1][1])
    print("   %-8s 자유재료 %.1f%% · 개체당 묶인 글자 %.2f (코드 길이 %d)" % (code, st.mean(fr), st.mean(per), len(code)))
