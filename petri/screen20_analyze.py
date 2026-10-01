#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""해부 20 선별 판정 — μ 0 에서 계통이 무너지는 밑 코드는 짝 계측(assay)이 못 미덥다 → 짝 단계에서 뺀다
사전등록: 해부20-사전등록-2026-10-01.md §3 (규칙은 새 코드의 결과를 보기 전에 고정)

  자료: inv20S_<밑>/mu0/dms_mat_<밑>_s*_list.json  (선별 배경 10 · 밑 코드 = 후보 + 기존 여섯)
  무너짐(collapse): 한 팔(ref · neutral · mut — 전부 표지 계통)의 끝 개체 수 / 처음 개체 수 < COLLAPSE(0.10)  (멸종 = 0 포함)
  K_B = 선별 배경 중 무너진 팔이 하나라도 있는 배경 수
  규칙  K_B ≥ K_MAX(1) 이면 그 밑 코드는 **짝 단계에서 뺀다**. 아니면 통과.
  규칙 대조(서술 · 규칙 자체는 안 바뀐다):
     S+  `acld` 는 빠져야 한다(기존 자료 40 배경 중 40 에서 무너짐)
     S−  기존 다섯(racld rascld rsacld racldx rascled)은 통과해야 한다(기존 자료 40 배경 중 0 · 최소 비 0.43)
  이 스크립트는 **s(로그비) · Δ · ΔΔ 를 계산하지도 인쇄하지도 않는다.** 밑 코드마다 K 와 최소 비만 인쇄한다.
실행: python3 screen20_analyze.py [petri 폴더]   (환경: PETRI_ROOT · PETRI_NBG · PETRI_PILOT · PETRI_SEEDS)
  → _RESULT_screen20.json (파일럿이면 _RESULT_pilot20_screen.json) — inv20_launch.sh 의 P 단계가 이 파일의 pass 를 읽는다
"""
import glob
import json
import os
import sys

BASE = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1 else "~/petri")
ROOT = os.environ.get("PETRI_ROOT", "inv20")
PILOT = os.environ.get("PETRI_PILOT") == "1"
NBG = int(os.environ.get("PETRI_NBG", "2" if PILOT else "10"))
OLD6 = "racld rascld rsacld racldx acld rascled".split()
RELIABLE5 = [b for b in OLD6 if b != "acld"]
# ---- 규칙 (사전등록 §3 · 결과 전 고정) ----
COLLAPSE = 0.10
K_MAX = 1


def loop_ticks(code):
    if "l" not in code or "c" not in code:
        return None
    li = code.index("l"); si = code.rfind("s", 0, li); start = si + 1 if si >= 0 else 0
    if start > li:
        return None
    seg = code[start:li + 1]
    return sum(1 for ch in seg if ch != "x") if "c" in seg else None


def inserts(base, ch):
    out = {}
    for p in range(len(base) + 1):
        k = base[:p] + ch + base[p:]; L = loop_ticks(k)
        if L is not None:
            out.setdefault(k, L)
    return out


def expected_seeds():
    # 검토 개정 R3(10-01): 판정 모드에서는 배경 목록 · 배경 수를 환경으로 바꿀 수 없다(판정선을 움직이는 손잡이를 닫는다)
    if not PILOT and (os.environ.get("PETRI_SEEDS") or os.environ.get("PETRI_NBG")):
        print("🚨 판정 모드에서 PETRI_SEEDS · PETRI_NBG 를 쓸 수 없다 — 파일럿이면 PETRI_PILOT=1"); sys.exit(1)
    if os.environ.get("PETRI_SEEDS"):
        return sorted(int(x) for x in os.environ["PETRI_SEEDS"].split())
    picked = []
    for s in range(100, 52, -1):
        z = json.load(open(os.path.join(BASE, "main", "mat_d8_mu0p01_s%05d.json" % s)))
        if z["extinct_at"] < 0 and z["final_counts"] and z["final_counts"][0][0] == "racld":
            picked.append(s)
    return sorted(picked[20:30])


cand = json.load(open(os.path.join(BASE, "_RESULT_cand20.json")))["picked"]
BASES = cand + OLD6
SEEDS = expected_seeds()
bad, OUT = [], {"root": ROOT, "pilot": PILOT, "collapse": COLLAPSE, "k_max": K_MAX, "candidates": cand, "seeds": SEEDS}
print("== 해부 20 선별 %s· 루트 %sS · 후보 %s · 대조 %s · 선별 배경 %s" % ("파일럿 " if PILOT else "", ROOT, " ".join(cand), " ".join(OLD6), SEEDS))
if len(SEEDS) != NBG:
    bad.append("기대 배경 수 %d ≠ %d" % (len(SEEDS), NBG))
res = {}
for b in BASES:
    codes = list(inserts(b, "n"))
    fs = sorted(glob.glob(os.path.join(BASE, "%sS_%s" % (ROOT, b), "mu0", "dms_mat_%s_s*_list.json" % b)))
    seen = []
    K, minr, n_arms, n_ext = 0, float("inf"), 0, 0
    for f in fs:
        z = json.load(open(f)); s = z["seed"]; seen.append(s)
        if abs(z["mu_assay"]) > 1e-12 or z["wt"] != b or z["arms_mode"] != "list" or z["cond"] != "mat": bad.append("설정 %s %d" % (b, s))
        if z.get("find_first") or z.get("remember_die") or z.get("density") or z.get("grow_mu") is not None or z.get("age0") is not None: bad.append("규칙 · 세계 옵션 %s %d" % (b, s))
        if not z["checksum_match"]: bad.append("배경 체크섬 %s %d" % (b, s))
        if not z["fidelity"]["same"]: bad.append("충실도 %s %d" % (b, s))
        if z["grow_extinct"] >= 0 or z["grow_top"][0][0] != "racld": bad.append("배경 %s %d" % (b, s))
        if not all(a["cons_ok"] for a in z["arms"]): bad.append("재료 보존 %s %d" % (b, s))
        muts = [a for a in z["arms"] if a["kind"] == "mut"]
        if [a["key"] for a in muts] != codes or len(z["arms"]) != 11 + len(codes): bad.append("팔 목록 %s %d" % (b, s))
        col = False
        for a in z["arms"]:
            m0, mT = a["series"][0][2], a["series"][-1][2]
            if m0 <= 0: bad.append("처음 개체 0 %s %d %s" % (b, s, a["key"])); continue
            r = mT / m0; minr = min(minr, r); n_arms += 1
            if mT == 0: n_ext += 1
            if r < COLLAPSE: col = True
        K += col
    if sorted(seen) != SEEDS: bad.append("배경 목록 %s %s ≠ %s" % (b, sorted(seen), SEEDS))
    res[b] = {"K": K, "n_bg": len(fs), "min_ratio": minr if n_arms else None, "n_arms": n_arms, "n_extinct_arms": n_ext,
              "pass": K < K_MAX, "role": "candidate" if b in cand else ("S+" if b == "acld" else "S-"), "codes": codes}

print("  점검 문제 %d" % len(bad))
for x in bad[:15]: print("    ", x)
complete = not bad
if bad and not PILOT:
    print("🚨 점검 실패 — 선별하지 않는다(짝 단계 발사 불가)")
print("\n  밑 코드     역할        K(무너진 배경)  최소 끝/처음 비  멸종 팔   규칙")
for b in BASES:
    r = res[b]
    print("  %-10s %-10s %3d / %-3d      %s     %4d / %-4d  %s" % (b, r["role"], r["K"], r["n_bg"], ("%.3f" % r["min_ratio"]) if r["min_ratio"] is not None else "  -  ",
                                                       r["n_extinct_arms"], r["n_arms"], "통과" if r["pass"] else "제외"))
s_plus = not res["acld"]["pass"]
s_minus = all(res[b]["pass"] for b in RELIABLE5)
print("\n  규칙 대조 S+ (`acld` 제외됨): %s · S− (기존 다섯 통과): %s" % ("✅" if s_plus else "❌", "✅" if s_minus else "❌ " + " ".join(b for b in RELIABLE5 if not res[b]["pass"])))
passed = [b for b in BASES if res[b]["pass"]]
cand_pass = [b for b in cand if res[b]["pass"]]
print("  짝 단계로: 후보 %s · 대조 %s" % (" ".join(cand_pass) or "(없음)", " ".join(b for b in passed if b not in cand)))
if not cand_pass:
    print("  ⚠️ 통과한 후보가 없다 — 사전등록 §5 의 '선별 전멸' 해석을 적용한다(짝 단계는 대조만 돌릴지 사람이 정한다)")
OUT.update({"per_base": res, "S_plus": s_plus, "S_minus": s_minus, "pass": passed if complete else [], "candidates_pass": cand_pass,
            "complete": complete, "bad": bad})
fn = "_RESULT_pilot20_screen.json" if PILOT else "_RESULT_screen20.json"
json.dump(OUT, open(os.path.join(BASE, fn), "w"), ensure_ascii=False, indent=1)
print("기록 →", fn)
if bad and not PILOT:
    sys.exit(1)
