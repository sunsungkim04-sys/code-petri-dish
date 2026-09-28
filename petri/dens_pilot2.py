# -*- coding: utf-8 -*-
"""해부 5 설계용 탐색 — 밀도 사다리 파일럿 요약. 판정 아님(시드 3개)."""
import glob
import json
import os
import statistics as st

DISH = 10732
RUNGS = [("mat_d04", "밀도 4"), ("mat_d08", "밀도 8"), ("mat_d14", "밀도 14"),
         ("mat_d20", "밀도 20"), ("mat_d32", "밀도 32"), ("space", "space(재료없음)")]


def load(folder, code):
    out = []
    for f in sorted(glob.glob(os.path.expanduser("~/petri/densp/%s/mono_*_%s_mu0_s*.json" % (folder, code)))):
        out.append(json.load(open(f)))
    return out


def at(z, t):
    for r in z["samples"]:
        if r[0] == t:
            return r
    return None


print("== 밀도 사다리 파일럿 (μ 0 · 시드 1~3 · 20,000틱) — 판정 아님")
print("%-16s %-7s %8s %8s %9s %8s %9s %6s" % ("칸", "코드", "P(1만)", "P(2만)", "점유%", "자유재료%", "개체당글자", "멸종"))
for folder, name in RUNGS:
    for code in ("racld", "rascld"):
        zs = load(folder, code)
        if not zs:
            print("%-16s %-7s  (없음)" % (name, code)); continue
        p10 = st.median(at(z, 10000)[1] for z in zs)
        p20 = st.median(z["samples"][-1][1] for z in zs)
        occ = 100.0 * p20 / DISH
        fr = "-" if zs[0]["free_end"] is None else "%.1f" % st.mean(100.0 * z["free_end"] / z["mat_total0"] for z in zs)
        per = "-" if zs[0]["tied_up_end"] is None else "%.2f" % st.mean(z["tied_up_end"] / z["samples"][-1][1] for z in zs)
        ext = sum(1 for z in zs if z["extinct_at"] >= 0)
        print("%-16s %-7s %8.0f %8.0f %9.1f %8s %9s %6d" % (name, code, p10, p20, occ, fr, per, ext))
    a = load(folder, "racld"); b = load(folder, "rascld")
    if a and b:
        import math
        d10 = st.mean(math.log(at(x, 10000)[1] + .5) - math.log(at(y, 10000)[1] + .5) for x, y in zip(a, b))
        d20 = st.mean(math.log(x["samples"][-1][1] + .5) - math.log(y["samples"][-1][1] + .5) for x, y in zip(a, b))
        print("%-16s %-7s   D(1만) %+.3f · D(2만) %+.3f" % ("", "→", d10, d20))
