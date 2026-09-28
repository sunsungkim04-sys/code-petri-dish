# -*- coding: utf-8 -*-
"""해부 5 사후 진단 — 사전등록 밖. **판정 아님.**
  (mono.js 의 파일명은 μ=0 을 "mu0" 으로 쓴다 — "mu0p0" 이 아니다)
  ① 천장 점검: 칸마다 P(1만)/P(2만) — 1 에 가까우면 1만 틱에 이미 천장이라 D 가 0 으로 눌린다
  ② 오류 아래 순도: μ 마다 끝에 '출발 코드 그대로' 인 개체 몫 — 강건성의 직접 흔적
실행: python3 dens_diag.py  →  화면
"""
import glob
import json
import math
import os
import statistics as st

BASE = os.path.expanduser("~/petri")
RUNGS = [("mat_d04", "밀도 4"), ("mat_d08", "밀도 8"), ("mat_d14", "밀도 14"),
         ("mat_d20", "밀도 20"), ("mat_d32", "밀도 32"), ("space", "space(재료없음)")]
MUS = [0.0, 0.001, 0.002, 0.003, 0.004, 0.005, 0.0075, 0.01]


def at(z, t):
    for r in z["samples"]:
        if r[0] == t:
            return r
    return None


print("== ① 천장 점검 — 1만 틱에 이미 다 찼나 (판정 아님)")
print("   %-16s %-8s %9s %9s %8s  %s" % ("칸", "코드", "P(1만)", "P(2만)", "1만/2만", "읽기"))
for folder, name in RUNGS:
    for code in ("racld", "rascld"):
        rs = [json.load(open(f)) for f in sorted(glob.glob(os.path.join(BASE, "dens", folder, "mono_*_%s_mu0_s*.json" % code)))]
        p10 = st.median(at(z, 10000)[1] for z in rs if at(z, 10000))
        p20 = st.median(z["samples"][-1][1] for z in rs)
        r = p10 / p20 if p20 else float("nan")
        note = "이미 천장 — D 가 0 으로 눌린다" if r > 0.97 else ("아직 자라는 중" if r < 0.93 else "거의 천장")
        print("   %-16s %-8s %9.0f %9.0f %8.3f  %s" % (name, code, p10, p20, r, note))

print("\n== ② 오류 아래 순도 — 끝에 출발 코드 그대로인 개체 몫 (판정 아님)")
print("   %6s | %-26s | %-26s" % ("μ", "racld  (몫 · 개체수)", "씨앗 rascld"))
for mu in MUS:
    row = []
    for code in ("racld", "rascld"):
        fr, wt = [], []
        for f in sorted(glob.glob(os.path.join(BASE, "mufine", "mono_mat_%s_mu%s_s*.json"
                                               % (code, ("0" if mu == 0 else str(mu).replace(".", "p")))))):
            z = json.load(open(f))
            if z["extinct_at"] >= 0:
                continue
            last = z["samples"][-1]
            if last[1]:
                fr.append(100.0 * last[3] / last[1]); wt.append(last[3])
        row.append("%6.1f%% · %6.0f" % (st.median(fr), st.median(wt)) if fr else "  (없음)")
    print("   %5.2f%% | %-26s | %-26s" % (mu * 100, row[0], row[1]))
