#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""해부 14 후속 ③ — 헛손질을 요청 글자 종류별로. 🟡 탐색 · 판정 아님(해부 14 노트 §10).
  회귀: 팔 기록 = inv14ib 의 같은 팔 · 칸마다 stall_by 합 = copy_stall
  코드마다: 글자당(쓴 글자 1개당) 헛손질 수를 글자 종류별로(배경 20 합) · N3 실패 아홉 쌍에서 어느 글자의 헛손질이 다른가(ln 비)
실행: python3 stallby_analyze.py [petri]  (환경 PETRI_ROOT=inv14bl PETRI_REF=inv14)
"""
import glob, json, math, os, sys
BASE = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1 else "~/petri")
ROOT = os.environ.get("PETRI_ROOT", "inv14bl"); REF = os.environ.get("PETRI_REF", "inv14")
LET = "nsracldjehx"
NBASES = "racld rascld rsacld racldx acld rascled".split()
N3 = [("rascld", "rascldn", "rnascld"), ("rascld", "rasclnd", "rnascld"), ("acld", "acnld", "ancld"), ("acld", "acnld", "nacld"),
      ("rascled", "nrascled", "rasclend"), ("rascled", "ranscled", "rasclend"), ("rascled", "rascledn", "rasclend"), ("rascled", "rasclend", "rasclned"), ("rascled", "rasclend", "rnascled")]
def loop_ticks(code):
    if "l" not in code or "c" not in code: return None
    li = code.index("l"); si = code.rfind("s", 0, li); start = si + 1 if si >= 0 else 0
    if start > li: return None
    seg = code[start:li + 1]
    return sum(1 for ch in seg if ch != "x") if "c" in seg else None
def inserts(base, ch):
    out = {}
    for p in range(len(base) + 1):
        k = base[:p] + ch + base[p:]; L = loop_ticks(k)
        if L is not None: out.setdefault(k, L)
    return out
bad = []; g0 = [0, 0]; AGG = {}
for b in NBASES:
    codes = inserts(b, "n")
    for f in sorted(glob.glob(os.path.join(BASE, "%s_%s" % (ROOT, b), "mu0", "ib_inv_s*_mu0.json"))):
        z = json.load(open(f)); s = z["seed"]
        if not z.get("by_letter") or list(z["codes"]) != list(codes): bad.append("설정 %s %d" % (b, s)); continue
        rf = os.path.join(BASE, "%sib_%s" % (REF, b), "mu0", "ib_inv_s%05d_mu0.json" % s)
        ref = json.load(open(rf)) if os.path.exists(rf) else None
        for a in z["arms"]:
            if ref is not None:
                ra = [x for x in ref["arms"] if x["key"] == a["key"] and x["kind"] == a["kind"]]
                if ra:
                    g0[0] += 1
                    if ra[0]["series"] != a["series"]: g0[1] += 1
            for bn in a["bins"]:
                if sum(bn["stall_by"][1]) != bn["copy_stall"][1]: bad.append("합 %s %d %s" % (b, s, a["key"]))
            if a["kind"] != "mut": continue
            if a["stopped"] >= 0 and a["series"][-1][2] == 0: continue
            d = AGG.setdefault((b, a["key"]), {"ok": 0, "stall": 0, "by": [0] * 11, "n": 0})
            d["ok"] += sum(bn["copy_ok"][1] for bn in a["bins"]); d["stall"] += sum(bn["copy_stall"][1] for bn in a["bins"]); d["n"] += 1
            for bn in a["bins"]:
                for k in range(11): d["by"][k] += bn["stall_by"][1][k]
print("== 해부 14 후속 ③ — 요청 글자별 헛손질 · 루트 %s · 회귀 대조 %d 팔 · 불일치 %d · 점검 문제 %d" % (ROOT, g0[0], g0[1], len(bad)))
for x in bad[:5]: print("   ", x)
print("\n== 코드마다 쓴 글자 1개당 헛손질 (글자 종류별 · 배경 합)")
for b in NBASES:
    print("  밑 코드 `%s`" % b)
    for (bb, c), d in sorted(AGG.items()):
        if bb != b: continue
        per = {LET[k]: d["by"][k] / d["ok"] for k in range(11) if d["by"][k]}
        print("     %-9s L %d · 글자당 헛손질 %.2f · %s" % (c, loop_ticks(c), d["stall"] / d["ok"], " ".join("%s %.2f" % (kv) for kv in sorted(per.items(), key=lambda kv: -kv[1]))))
print("\n== N3 실패 쌍: 글자별 헛손질(글자당)의 ln 비 (앞/뒤) — 어느 글자가 다른가")
for b, x, y in N3:
    dx, dy = AGG.get((b, x)), AGG.get((b, y))
    if not dx or not dy: continue
    row = []
    for k in range(11):
        px, py = dx["by"][k] / dx["ok"], dy["by"][k] / dy["ok"]
        if px > 0 and py > 0: row.append((LET[k], math.log(px / py), px - py))
    tot = math.log((dx["stall"] / dx["ok"]) / (dy["stall"] / dy["ok"]))
    row.sort(key=lambda t: -abs(t[2]))
    print("   %-9s / %-9s 전체 ln 비 %+.3f · 차(글자당) 큰 순: %s" % (x, y, tot, " · ".join("%s %+.3f(ln %+.2f)" % (l, dd, lr) for l, lr, dd in row[:5])))
json.dump({"%s|%s" % k: v for k, v in AGG.items()}, open("_RESULT_stallby.json", "w"), indent=1)
print("기록 → _RESULT_stallby.json")
