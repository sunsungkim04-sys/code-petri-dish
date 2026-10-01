#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""해부 20 판정 — 가족 밖 밑 코드에서 짝지은 고리 효과(고리 한 틱당 Δ)가 같은 부호인가
사전등록: 해부20-사전등록-2026-10-01.md §4 (판정선은 결과 전 고정 · 해부 13 · 15 E 와 같은 통계)

  자료: inv20P_<밑>/mu0/dms_mat_<밑>_s*_list.json  (짝 배경 20 · 밑 코드 = _RESULT_screen20.json 의 pass)
  s  = ln((m_T + 0.5)/m_0) − ln((n_T − m_T + 0.5)/(n_0 − m_0))     (해부 13 과 같다)
  ΔΔ = s(고리 L+1 코드) − s(고리 L 코드) · 배경 짝 · 배경 재추출 2,000 · 95% 구간
  짝 = 이웃 고리 등급(틱 차 1)의 모든 (짧은, 긴) 쌍 · 음성 쌍 = 고리 같은 모든 쌍 (둘 다 `n` 삽입만)

  W4 (서술)           고리당 복사가 2 인 밑 코드는 W1 × 2 = '글자당 틱' 단위 값도 적는다
  W1 (주 · 후보마다)  합친 틱당 ΔΔ: 배경마다 그 밑 코드의 모든 짝 ΔΔ/ΔL 평균 → 재추출 — 상한 < 0
  N1 (후보마다)       모든 짝의 ΔΔ 상한 < 0                               (해부 13 N1)
  N2 (서술)           틱당 점추정이 [−1.0, −0.2] 안인 짝 수                  (해부 13 N2 · 여기선 판정 아님)
  N3 (후보마다)       모든 음성 쌍 |ΔΔ| ≤ 0.15                              (해부 13 N3)
  PC (양성 대조)      짝 단계에 들어온 기존 밑 코드마다 W1 상한 < 0 — 하나라도 ❌ 면 후보 판정을 해석하지 않는다(보류)
  읽기  후보가 W1 · N1 · N3 전부 ✅ → '분리된 재현' · W1 ✅ 이지만 N1 또는 N3 ❌ → '부호만' · W1 ❌ → '재현 안 됨'
실행: python3 pairs20_analyze.py [petri 폴더]   (환경: PETRI_ROOT · PETRI_NBG · PETRI_PILOT · PETRI_SEEDS · PETRI_SCREEN)
  → _RESULT_pairs20.json (파일럿이면 _RESULT_pilot20.json)
"""
import glob
import json
import math
import os
import random
import sys

BASE = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1 else "~/petri")
ROOT = os.environ.get("PETRI_ROOT", "inv20")
PILOT = os.environ.get("PETRI_PILOT") == "1"
NBG = int(os.environ.get("PETRI_NBG", "2" if PILOT else "20"))
SCREEN = os.environ.get("PETRI_SCREEN", "_RESULT_pilot20_screen.json" if PILOT else "_RESULT_screen20.json")
# ---- 판정선 (사전등록 §4 · 결과 전 고정) ----
TOL_NEG, PT_LO, PT_HI = 0.15, -1.0, -0.2        # 해부 13 그대로
N_BOOT = 2000
RNG_SEED = 20261001     # 재추출 난수는 (시드, 이름)으로 키 — 통과 밑 코드 집합이 바뀌어도 한 칸의 구간이 안 바뀐다
bad, notes = [], []


def loop_ticks(code):
    if "l" not in code or "c" not in code:
        return None
    li = code.index("l"); si = code.rfind("s", 0, li); start = si + 1 if si >= 0 else 0
    if start > li:
        return None
    seg = code[start:li + 1]
    return sum(1 for ch in seg if ch != "x") if "c" in seg else None


def copies_per_loop(code):
    li = code.index("l"); si = code.rfind("s", 0, li); start = si + 1 if si >= 0 else 0
    return code[start:li].count("c")


def inserts(base, ch):
    out = {}
    for p in range(len(base) + 1):
        k = base[:p] + ch + base[p:]; L = loop_ticks(k)
        if L is not None:
            out.setdefault(k, L)
    return out


def s_of(arm):
    r = arm["series"]; n0, m0 = r[0][1], r[0][2]; nT, mT = r[-1][1], r[-1][2]
    return math.log((mT + 0.5) / m0) - math.log((nT - mT + 0.5) / (n0 - m0))


def pct(xs, p):
    xs = sorted(xs); k = (len(xs) - 1) * p; f = math.floor(k); c = min(f + 1, len(xs) - 1)
    return xs[f] + (xs[c] - xs[f]) * (k - f)


def boot_mean(v, name):
    rng = random.Random("%d|%s" % (RNG_SEED, name))
    m = sum(v) / len(v)
    bs = [sum(v[rng.randrange(len(v))] for _ in v) / len(v) for _ in range(N_BOOT)]
    return m, pct(bs, 0.025), pct(bs, 0.975)


def fmt(m, lo, hi):
    return "%+.3f [%+.3f, %+.3f]" % (m, lo, hi)


def pairs_of(codes):
    out = []
    Ls = sorted(set(codes.values()))
    for L1, L2 in zip(Ls, Ls[1:]):
        if L2 - L1 != 1: continue
        for a in codes:
            for b in codes:
                if codes[a] == L1 and codes[b] == L2: out.append((a, b))
    return out


def same_loop_pairs(codes):
    cl = sorted(codes); out = []
    for i, a in enumerate(cl):
        for b in cl[i + 1:]:
            if codes[a] == codes[b]: out.append((a, b))
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
    return sorted(picked[:20])


scr = json.load(open(os.path.join(BASE, SCREEN)))
if not scr.get("complete"):
    print("🚨 선별 판정(%s)이 완결되지 않았다 — 판정하지 않는다" % SCREEN); sys.exit(1)
CAND = scr["candidates"]
PASS = scr["pass"]
SEEDS = expected_seeds()
if set(SEEDS) & set(scr["seeds"]) and not PILOT:
    bad.append("짝 배경과 선별 배경이 겹친다 %s" % sorted(set(SEEDS) & set(scr["seeds"])))
# 검토 개정 R4(10-01): 판정 모드에서는 선별 JSON 이 본 선별(파일럿 아님 · 선별 배경 10 = 100 에서 내려가며 21~30 번째)이어야 한다
if not PILOT:
    _all = []
    for _s in range(100, 52, -1):
        _z = json.load(open(os.path.join(BASE, "main", "mat_d8_mu0p01_s%05d.json" % _s)))
        if _z["extinct_at"] < 0 and _z["final_counts"] and _z["final_counts"][0][0] == "racld":
            _all.append(_s)
    if scr.get("pilot") or sorted(scr["seeds"]) != sorted(_all[20:30]) or scr.get("root") != ROOT:
        bad.append("선별 JSON 이 본 선별이 아니다(pilot %s · 배경 %s · 루트 %s)" % (scr.get("pilot"), scr["seeds"], scr.get("root")))
OUT = {"root": ROOT, "pilot": PILOT, "screen": SCREEN, "pass": PASS, "candidates": CAND, "seeds": SEEDS,
       "tol_neg": TOL_NEG, "per_tick_band": [PT_LO, PT_HI]}
print("== 해부 20 %s— 루트 %sP · 짝 단계 밑 코드 %s (후보 %s) · 배경 %s" % ("파일럿 " if PILOT else "판정 ", ROOT, " ".join(PASS), " ".join(b for b in PASS if b in CAND) or "없음", SEEDS))

S = {}          # b -> {seed: {code: s}}
COLL = {}       # b -> 무너진 팔 수(서술)
for b in PASS:
    codes = list(inserts(b, "n"))
    fs = sorted(glob.glob(os.path.join(BASE, "%sP_%s" % (ROOT, b), "mu0", "dms_mat_%s_s*_list.json" % b)))
    d, seen, nc = {}, [], 0
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
        d[s] = {a["key"]: s_of(a) for a in muts}
        nc += sum(1 for a in muts if a["series"][-1][2] < 0.10 * a["series"][0][2])
    if sorted(seen) != SEEDS: bad.append("배경 목록 %s (%d/%d)" % (b, len(seen), len(SEEDS)))
    S[b] = d; COLL[b] = nc

# 같은 배경의 밑 코드 파일은 배경이 같아야 한다(grow_checksum) — 짝 단계 전 밑 코드 사이
gcs = {}
for b in PASS:
    for f in glob.glob(os.path.join(BASE, "%sP_%s" % (ROOT, b), "mu0", "dms_mat_%s_s*_list.json" % b)):
        z = json.load(open(f)); gcs.setdefault(z["seed"], set()).add(z["grow_checksum"])
ndiff = sum(1 for v in gcs.values() if len(v) > 1)
print("  배경 동일성: 배경 %d 중 밑 코드 사이 grow_checksum 이 다른 배경 %d" % (len(gcs), ndiff))
if ndiff: bad.append("밑 코드 사이 배경 체크섬 불일치 %d" % ndiff)
print("  점검 문제 %d" % len(bad))
for x in bad[:15]: print("    ", x)
if bad and not PILOT:
    print("🚨 점검 실패 — 판정하지 않는다"); sys.exit(1)

res = {}
for b in PASS:
    codes = inserts(b, "n"); role = "후보" if b in CAND else "대조"
    k = copies_per_loop(b)
    print("\n-- 밑 코드 `%s` (%s · 고리 %d틱 · 고리당 복사 %d · 등급 %s) · 짝 배경 %d · 짝 단계에서 무너진 팔 %d(서술)" % (b, role, loop_ticks(b), k, sorted(set(codes.values())), len(S[b]), COLL[b]))
    n1, n2, n3 = [], [], []
    per_bg = {}
    for a, c in pairs_of(codes):
        dL = codes[c] - codes[a]
        v = {s: S[b][s][c] - S[b][s][a] for s in S[b] if a in S[b][s] and c in S[b][s]}
        if len(v) < 2: continue
        m, lo, hi = boot_mean([v[s] for s in sorted(v)], "pair|%s|%s|%s" % (b, a, c))
        per = m / dL; ok1 = hi < 0; ok2 = PT_LO <= per <= PT_HI
        n1.append(ok1); n2.append(ok2)
        for s in v: per_bg.setdefault(s, []).append(v[s] / dL)
        print("     ΔΔ %-10s(%d) − %-10s(%d) = %s · 배경 %d · 틱당 %+.3f %s%s" % (c, codes[c], a, codes[a], fmt(m, lo, hi), len(v), per,
                                                                           "✅" if ok1 else "❌ 상한≥0", "" if ok2 else " (틱당 N2 밖)"))
        OUT["pair|%s|%s|%s" % (b, a, c)] = {"dd": m, "ci": [lo, hi], "per_tick": per, "n_bg": len(v), "N1": ok1, "N2": ok2}
    negs = []
    for a, c in same_loop_pairs(codes):
        v = [S[b][s][c] - S[b][s][a] for s in sorted(S[b]) if a in S[b][s] and c in S[b][s]]
        if len(v) < 2: continue
        m, lo, hi = boot_mean(v, "neg|%s|%s|%s" % (b, a, c)); ok3 = abs(m) <= TOL_NEG; n3.append(ok3); negs.append(abs(m))
        print("     [음성] %-10s vs %-10s (%d틱) = %s %s" % (a, c, codes[a], fmt(m, lo, hi), "✅" if ok3 else "❌ 0.15 밖"))
        OUT["neg|%s|%s|%s" % (b, a, c)] = {"dd": m, "ci": [lo, hi], "N3": ok3}
    if len(per_bg) >= 2:
        w = [sum(x) / len(x) for s, x in sorted(per_bg.items())]
        wm, wlo, whi = boot_mean(w, "W1|%s" % b); w1 = whi < 0
    else:
        wm = wlo = whi = float("nan"); w1 = False
    N1 = bool(n1) and all(n1); N3 = bool(n3) and all(n3)
    read = ("분리된 재현" if (w1 and N1 and N3) else ("부호만" if w1 else "재현 안 됨"))
    print("   W1 합친 틱당 ΔΔ %s %s · N1 %s(%d/%d 짝) · N2 서술 %d/%d · N3 %s(%d/%d 쌍 · 최대 |ΔΔ| %.3f) → %s"
          % (fmt(wm, wlo, whi), "✅" if w1 else "❌", "✅" if N1 else "❌", sum(n1), len(n1), sum(n2), len(n2),
             "✅" if N3 else "❌", sum(n3), len(n3), max(negs) if negs else float("nan"), read if role == "후보" else "(대조)"))
    # W4 서술: 고리당 복사 k 가 1 이 아니면 '글자당 틱' 으로 환산한 값도 적는다(한 틱 더 = 글자당 1/k 틱 더)
    if k != 1 and wm == wm:
        print("   W4 서술: 고리당 복사 %d — 글자당 틱으로 환산하면 %+.3f [%+.3f, %+.3f] (고리 틱당 × %d)" % (k, wm * k, wlo * k, whi * k, k))
    res[b] = {"role": role, "W1": {"m": wm, "ci": [wlo, whi], "pass": w1}, "N1": N1, "n1": [sum(n1), len(n1)], "N2": [sum(n2), len(n2)],
              "N3": N3, "n3": [sum(n3), len(n3)], "neg_max": max(negs) if negs else None, "read": read, "collapsed_arms": COLL[b],
              "loop": loop_ticks(b), "copies_per_loop": k}
OUT["per_base"] = res

print("\n== 종합")
ctrl = [b for b in PASS if b not in CAND]
pc = bool(ctrl) and all(res[b]["W1"]["pass"] for b in ctrl)
print("  PC 양성 대조(기존 밑 코드 %s 의 W1): %s" % (" ".join(ctrl), "✅" if pc else "❌ — 후보 판정 **보류**"))
OUT["PC"] = pc
cands = [b for b in PASS if b in CAND]
if not cands:
    print("  선별을 통과한 후보가 없다 — 사전등록 §5 '선별 전멸'")
for b in cands:
    r = res[b]
    print("  후보 `%s`: W1 %s · N1 %s · N3 %s → %s%s" % (b, "✅" if r["W1"]["pass"] else "❌", "✅" if r["N1"] else "❌", "✅" if r["N3"] else "❌", r["read"],
                                                  "" if pc else " (PC ❌ — 해석 보류)"))
n_sep = sum(1 for b in cands if res[b]["read"] == "분리된 재현"); n_sign = sum(1 for b in cands if res[b]["W1"]["pass"])
print("  후보 %d 중 W1(같은 부호) %d · 분리된 재현(W1 · N1 · N3) %d · 선별 제외 후보 %s" % (len(cands), n_sign, n_sep, " ".join(b for b in CAND if b not in PASS) or "없음"))
OUT["summary"] = {"candidates_measured": cands, "W1_pass": n_sign, "separated": n_sep, "excluded": [b for b in CAND if b not in PASS], "PC": pc}
if PILOT:
    print("  🟡 파일럿 — 판정 아님(배경 %d)" % len(SEEDS))
OUT["bad"] = bad; OUT["notes"] = notes
fn = "_RESULT_pilot20.json" if PILOT else "_RESULT_pairs20.json"
json.dump(OUT, open(os.path.join(BASE, fn), "w"), ensure_ascii=False, indent=1)
print("기록 →", fn)
