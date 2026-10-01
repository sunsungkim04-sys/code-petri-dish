#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""해부 24 판정 — 직접 계수(invbud R · T)를 둘째 · 셋째 배경 세트에서 반복한다. 배경 세트마다 따로 판정한다.
사전등록: 사전등록/해부24-사전등록-2026-10-01.md §3 · §4(🔒 판정선). 판정선은 아래 상수 — 결과를 보기 전에 잠갔다.

  세트 1  inv14ib_*  배경 3~22(오름차순 racld 정상 20)       해부 14 판정 그대로 → 여기서는 **양성 대조**(점추정이 _RESULT_pairs14.json 과 같아야 한다 · 판정 밖)
  세트 2  inv15Eib_* 배경 26~50 중 racld 정상 20(해부 15 E)  **사전등록된 재분석**(새 실행 없음) · G0a 기준 = inv15E_* (dms 목록 모드 μ 0)
  세트 3  inv24ib_*  배경 = _RESULT_pairs20.json seeds(100…76) 새 실행 · 밑 코드 = 해부 20 짝 단계 8 · G0a 기준 = inv20P_*
  파일럿  --pilot --root R --seeds a,b --bases …  (R ib_* · 기준 R D_*)

  정의(해부 14 · pairs14_analyze.py 에서 복사 — 아래 함수마다 출처 줄):
     R = (copy_ok + copy_stall + copy_del) / copy_ok     T = copy_phase / copy_ok      (끼운 계통 · 팔 3,000틱 합 · 배경마다 하나)
     L = loop_ticks(code) · k = copies_per_loop(code)(pairs20_analyze.py 50행 · 해부 14 의 밑 코드는 전부 k = 1)
  자기 검산  모든 코드 pooled T/R − L/k ∈ [−0.5, +1.5]  (k = 1 이면 해부 14 와 같은 식) — 벗어나면 멈춤
  G0   파일 점검(버전 · μ 0 · 기본 설정 · 체크섬 · 코드 목록 · 배경 목록 정확 일치) + G0a: 끼운 팔의 [t, n, m] · stopped = 같은 배경 · 코드의 dms 목록 모드 팔 · grow_checksum 같음
  B1 (= 해부 14 A1)   밑 코드마다 이웃 고리 등급: 배경별 등급 평균 ln R 차(짧은 − 긴) 하한 > 0
  B2 (= 해부 14 A2a)  모든 같은-고리 쌍 |배경 평균 ln(R_i/R_j)| ≤ TAU 0.10
  B3 (= 해부 14 A2b)  모든 같은-고리 쌍 |pooled T/R_i − T/R_j| ≤ TR_PAIR 0.10 틱
  B4 🆕 핵심 1       밑 코드마다 Q = max 같은-고리 |ln R 비| ÷ min 이웃 등급 효과 — 배경 재추출(공통 배경 한꺼번에) 상한 < Q_MAX 0.5
  B5 🆕 핵심 2       그 세트의 N3 실패 쌍마다 배경별 (ln T 비 − ln R 비) = ln(T/R)_a − ln(T/R)_c 의 95% 구간 ⊂ [−D_MAX, +D_MAX], D_MAX 0.04
  서술  B6 N3 실패 쌍의 부호 정합(R 이 큰 쪽 = 침입 적합도가 낮은 쪽인가) · A1b · A2c 점추정 · A3 꼴 인구학 ln 비
실행: python3 pairs24_analyze.py --set 2|3 [--base ~/petri]         → _RESULT_pairs24_set2.json · _RESULT_pairs24_set3.json (txt 는 stdout 을 받는다)
      python3 pairs24_analyze.py --set 1                            → _RESULT_pairs24_set1.json (양성 대조 · 판정 밖)
      python3 pairs24_analyze.py --pilot --root pilot24 --seeds 53,54 --bases racld,rascled,reasccld → _RESULT_pilot24.json
      python3 pairs24_analyze.py --reg pilot24reg                   → G0r 회귀만(지금 invbud.js 출력 = 옛 출력 · 판정 통계량 없음)
종료 코드: 0 정상 · 2 점검(G0 · 자기 검산) 실패로 판정하지 않음
"""
import argparse
import glob
import json
import math
import os
import random
import sys

ap = argparse.ArgumentParser()
ap.add_argument("--set", type=int, choices=[1, 2, 3])
ap.add_argument("--pilot", action="store_true")
ap.add_argument("--root")
ap.add_argument("--seeds")
ap.add_argument("--bases")
ap.add_argument("--reg")
ap.add_argument("--base", default="~/petri")
ap.add_argument("--out")
A = ap.parse_args()
BASE = os.path.expanduser(A.base)

# ---------------- 🔒 판정선 (사전등록 §4 · 10-01 · 둘째 · 셋째 세트의 판정 통계량을 보기 전) ----------------
TAU = 0.10                 # B2 — pairs14_analyze.py 44행 그대로
TR_PAIR = 0.10             # B3 — pairs14_analyze.py 45행 그대로
TR_LO, TR_HI = -0.5, 1.5   # 자기 검산 — pairs14_analyze.py 46행 그대로
Q_MAX = 0.5                # B4 — _RESULT_pairs14.txt 132~137행 A2c 점추정 0.13~0.27 · 해부 14 τ 의 근거('고리 한 틱 효과의 절반', 14 사전등록 §4b)와 같은 뜻
D_MAX = 0.04               # B5 — _RESULT_pairs14.txt A3 아홉 쌍 'T-R' 최대 |+0.020|(145행 ranscled/rasclend) 의 2배
N_BOOT = 2000
RNG_SEED = 20261001
OLD6 = ["racld", "rascld", "rsacld", "racldx", "acld", "rascled"]
FAM1 = {"racld", "rascld", "rsacld", "racldx"}
# 해부 14 N3 실패 아홉 쌍 — pairs14_analyze.py 50~52행 그대로(세트 1 양성 대조용)
N3_FAIL14 = {("acnld", "ancld"), ("acnld", "nacld"), ("nrascled", "rasclend"), ("ranscled", "rasclend"), ("rascledn", "rasclend"),
             ("rasclend", "rasclned"), ("rasclend", "rnascled"), ("rascldn", "rnascld"), ("rasclnd", "rnascld")}
KEYS_SUM = ["copy_ok", "copy_stall", "copy_del", "copy_phase", "births", "deaths", "death_age", "death_prestep", "death_other",
            "stepped", "not_stepped", "unknown"]          # pairs14_analyze.py 와 같은 칸
bad, notes = [], []


# ---------------- 함수 (출처 표기) ----------------
def loop_ticks(code):                       # pairs14_analyze.py 57~66행 그대로
    if "l" not in code or "c" not in code:
        return None
    li = code.index("l")
    si = code.rfind("s", 0, li)
    start = si + 1 if si >= 0 else 0
    if start > li:
        return None
    seg = code[start:li + 1]
    return sum(1 for ch in seg if ch != "x") if "c" in seg else None


def copies_per_loop(code):                  # pairs20_analyze.py 50~52행 그대로
    li = code.index("l"); si = code.rfind("s", 0, li); start = si + 1 if si >= 0 else 0
    return code[start:li].count("c")


def inserts(base, ch):                      # pairs14_analyze.py 69~76행 그대로
    out = {}
    for p in range(len(base) + 1):
        k = base[:p] + ch + base[p:]
        L = loop_ticks(k)
        if L is not None:
            out.setdefault(k, L)
    return out


def pct(xs, p):                             # pairs14_analyze.py 79~81행 그대로
    xs = sorted(xs); k = (len(xs) - 1) * p; f = math.floor(k); c = min(f + 1, len(xs) - 1)
    return xs[f] + (xs[c] - xs[f]) * (k - f)


def boot_mean(v, name):                     # pairs14 84~87행의 정의 · 난수는 이름으로 키(pairs20_analyze.py 74행 꼴 — 공용 흐름 결합 방지)
    rng = random.Random("%d|%s" % (RNG_SEED, name))
    m = sum(v) / len(v)
    bs = [sum(v[rng.randrange(len(v))] for _ in v) / len(v) for _ in range(N_BOOT)]
    return m, pct(bs, 0.025), pct(bs, 0.975)


def fmt(m, lo, hi):
    return "%+.3f [%+.3f, %+.3f]" % (m, lo, hi)


def same_loop_pairs(codes):                 # pairs15/20 의 same_loop_pairs 와 같음(정렬된 a < c)
    cl = sorted(codes); out = []
    for i, a in enumerate(cl):
        for c in cl[i + 1:]:
            if codes[a] == codes[c]: out.append((a, c))
    return out


def pick_asc(first, n):                     # inv14/15_launch.sh 의 배경 선택과 같음
    out = []
    for s in range(first, 101):
        z = json.load(open(os.path.join(BASE, "main", "mat_d8_mu0p01_s%05d.json" % s)))
        if z["extinct_at"] < 0 and z["final_counts"] and z["final_counts"][0][0] == "racld":
            out.append(s)
        if len(out) == n:
            break
    return out


def arm_cmp(a, b):                          # G0r — elapsed_ms 만 빼고 팔 전체 비교
    ka = {k: v for k, v in a.items() if k != "elapsed_ms"}; kb = {k: v for k, v in b.items() if k != "elapsed_ms"}
    return json.dumps(ka, sort_keys=True) == json.dumps(kb, sort_keys=True)


# ---------------- G0r 회귀 모드 ----------------
if A.reg:
    print("== 해부 24 G0r — 지금 invbud.js 의 출력이 옛 출력(해부 14 · 15 를 만든 c1a7b76b…)과 같은가")
    pairs = [("%s/s3_racld/ib_inv_s00003_mu0.json" % A.reg, "inv14ib_racld/mu0/ib_inv_s00003_mu0.json"),
             ("%s/s26_rascled/ib_inv_s00026_mu0.json" % A.reg, "inv15Eib_rascled/mu0/ib_inv_s00026_mu0.json")]
    ok_all = True
    for new, old in pairs:
        zn = json.load(open(os.path.join(BASE, new))); zo = json.load(open(os.path.join(BASE, old)))
        top = [k for k in set(zn) | set(zo) if k not in ("arms", "elapsed_ms") and zn.get(k) != zo.get(k)]
        na = len(zn["arms"]); same = sum(1 for x, y in zip(zn["arms"], zo["arms"]) if arm_cmp(x, y))
        ok = (not top) and na == len(zo["arms"]) and same == na
        ok_all &= ok
        print("  %s vs %s: 머리 칸 다른 것 %s · 팔 %d/%d 같음 %s" % (new, old, top or "없음", same, na, "✅" if ok else "❌"))
    print("  G0r 종합: %s" % ("✅" if ok_all else "❌ — 지금 invbud.js 로는 세트 1 · 2 와 같은 계측기라고 할 수 없다 · 멈춤"))
    sys.exit(0 if ok_all else 2)

# ---------------- 세트 설정 ----------------
if A.pilot:
    SET = "pilot"; ROOT = A.root
    SEEDS = sorted(int(x) for x in A.seeds.split(","))
    BASES = A.bases.split(",")
    IB_DIR = ROOT + "ib_%s/mu0"; REF_DIR = ROOT + "D_%s/mu0"
elif A.set == 1:
    SET = 1; SEEDS = pick_asc(1, 20); BASES = OLD6; IB_DIR = "inv14ib_%s/mu0"; REF_DIR = None
elif A.set == 2:
    SET = 2; SEEDS = pick_asc(26, 20); BASES = OLD6; IB_DIR = "inv15Eib_%s/mu0"; REF_DIR = "inv15E_%s/mu0"
elif A.set == 3:
    SET = 3
    p20 = json.load(open(os.path.join(BASE, "_RESULT_pairs20.json")))
    if p20.get("pilot") or len(p20["seeds"]) != 20: bad.append("_RESULT_pairs20.json 이 본 판정이 아니다")
    SEEDS = sorted(p20["seeds"]); BASES = list(p20["pass"])
    IB_DIR = "inv24ib_%s/mu0"; REF_DIR = "inv20P_%s/mu0"
else:
    print("--set 1|2|3 · --pilot · --reg 중 하나"); sys.exit(1)

# N3 실패 쌍과 그 쌍의 ΔΔ(= s(c) − s(a) · a < c 정렬) — 세트마다 **그 세트의 Δ 판정 파일**에서
N3 = {}   # (b, a, c) -> dd
if SET == 1:
    for fn, fam in (("_RESULT_pairs13.json", None), ("_RESULT_pairs13_racldfam.json", FAM1)):
        z = json.load(open(os.path.join(BASE, fn)))
        for k, v in z.items():
            if not k.startswith("neg|"): continue
            _, b, ch, a, c = k.split("|")
            if ch != "n" or (fam is not None) != (b in FAM1): continue
            if (a, c) in N3_FAIL14 or (c, a) in N3_FAIL14:
                N3[(b, a, c)] = -v["dd"]          # pairs13 의 dd = s(a) − s(c) (pairs13_analyze.py 131행) → 부호 뒤집어 s(c) − s(a)
elif SET == 2:
    z = json.load(open(os.path.join(BASE, "_RESULT_pairs15.json")))
    for k, v in z.items():
        if k.startswith("E2|") and v["fail"]:
            _, b, a, c = k.split("|"); N3[(b, a, c)] = v["dd"]   # pairs15_analyze.py 328행: d[s][c] − d[s][a]
elif SET == 3:
    for k, v in p20.items():
        if k.startswith("neg|") and not v["N3"]:
            _, b, a, c = k.split("|"); N3[(b, a, c)] = v["dd"]   # pairs20_analyze.py 193행: S[c] − S[a]
else:   # 파일럿: Δ 없음 — 모든 같은-고리 쌍을 B5 경로 연습에 쓴다(판정 아님)
    for b in BASES:
        for a, c in same_loop_pairs(inserts(b, "n")):
            N3[(b, a, c)] = None
OUT = {"set": SET, "seeds": SEEDS, "bases": BASES, "ib_dir": IB_DIR, "ref_dir": REF_DIR,
       "lines": {"TAU": TAU, "TR_PAIR": TR_PAIR, "TR": [TR_LO, TR_HI], "Q_MAX": Q_MAX, "D_MAX": D_MAX},
       "n3_pairs": ["%s|%s|%s" % k for k in sorted(N3)]}

# ---------------- 읽기 · G0 ----------------
IB = {}; DEAD = {}
g0a = {"checked": 0, "mismatch": [], "no_ref": 0, "gcs_diff": 0}
for b in BASES:
    codes = inserts(b, "n")
    IB[b] = {c: {} for c in codes}; DEAD[b] = {c: [] for c in codes}
    fs = sorted(glob.glob(os.path.join(BASE, IB_DIR % b, "ib_inv_s*_mu0.json")))
    got = sorted(json.load(open(f))["seed"] for f in fs)
    if got != SEEDS: bad.append("IB 배경 목록 %s: %s ≠ %s" % (b, got, SEEDS))
    for f in fs:
        z = json.load(open(f)); s = z["seed"]
        if z.get("invbud_version") != "1.0.0" or z.get("sim_version") != "0.3.5": bad.append("IB 버전 %s %d" % (b, s))
        if z["mode"] != "inv" or abs(z["mu"]) > 1e-12 or z.get("wt") != "racld": bad.append("IB 설정 %s %d" % (b, s))
        if any(z.get(k) is not None for k in ("by_letter", "density", "grow_mu", "ancestor")): bad.append("IB 기본값 아님 %s %d" % (b, s))
        if not z["checksum_match"]: bad.append("IB 배경 체크섬 %s %d" % (b, s))
        if list(z["codes"]) != list(codes): bad.append("IB 코드 목록 %s %d" % (b, s))
        if len(z["arms"]) != 1 + len(codes) or z["arms"][0]["kind"] != "ref": bad.append("IB 팔 수 %s %d" % (b, s))
        ref = None
        if REF_DIR:
            rp = os.path.join(BASE, REF_DIR % b, "dms_mat_%s_s%05d_list.json" % (b, s))
            if os.path.exists(rp):
                ref = json.load(open(rp))
                if ref.get("grow_checksum") != z.get("grow_checksum"): g0a["gcs_diff"] += 1
            else:
                g0a["no_ref"] += 1
        for a in z["arms"]:
            if a["kind"] != "mut": continue
            c = a["key"]
            if ref is not None:
                cand = [x for x in ref["arms"] if x["kind"] == "mut" and x["key"] == c]
                if not cand:
                    g0a["mismatch"].append("%s s%d %s: 기준 팔 없음" % (b, s, c))
                else:
                    g0a["checked"] += 1
                    if [r[:3] for r in cand[0]["series"]] != [r[:3] for r in a["series"]] or cand[0]["stopped"] != a["stopped"]:
                        g0a["mismatch"].append("%s s%d %s" % (b, s, c))
            d = {k: sum(bn[k][1] for bn in a["bins"]) for k in KEYS_SUM}      # 표지 1(끼운 계통) — pairs14 143행
            if a["stopped"] >= 0 and a["series"][-1][2] == 0:
                DEAD[b][c].append(s); continue
            if d["copy_del"] != 0: bad.append("IB μ0 인데 빠뜨림 %s %d %s" % (b, s, c))
            expo = d["stepped"] + d["not_stepped"]
            if expo and d["unknown"] / expo > 0.01: bad.append("IB 모름 몫 %s %d %s %.3f" % (b, s, c, d["unknown"] / expo))
            if d["copy_ok"] == 0: bad.append("IB copy_ok 0 %s %d %s" % (b, s, c)); continue
            d["R"] = (d["copy_ok"] + d["copy_stall"] + d["copy_del"]) / d["copy_ok"]
            d["T"] = d["copy_phase"] / d["copy_ok"]
            d["expo"] = expo
            IB[b][c][s] = d

for (b, a, c) in N3:      # N3 쌍의 코드가 이 세트에 없으면 판정 전에 멈춘다(조용히 빠지지 않게)
    if b not in IB or a not in IB[b] or c not in IB[b]:
        bad.append("N3 쌍 %s %s/%s 의 코드가 이 세트에 없다" % (b, a, c))
print("== 해부 24 %s — 세트 %s · 밑 코드 %s · 배경 %d %s" % ("파일럿" if A.pilot else ("양성 대조(판정 밖)" if SET == 1 else "판정"), SET, " ".join(BASES), len(SEEDS), SEEDS))
print("  N3 실패 쌍(그 세트의 Δ 판정에서): %d %s" % (len(N3), "· 파일럿은 모든 같은-고리 쌍(연습)" if A.pilot else ""))
if REF_DIR:
    print("  G0a 팔 기록 대조: %d 팔 · 불일치 %d · 기준 파일 없음 %d · grow_checksum 다른 배경 %d" % (g0a["checked"], len(g0a["mismatch"]), g0a["no_ref"], g0a["gcs_diff"]))
    for x in g0a["mismatch"][:5]: print("     ✗", x)
    if g0a["mismatch"] or g0a["no_ref"] or g0a["gcs_diff"]:
        bad.append("G0a 불일치 %d · 기준 없음 %d · 체크섬 %d" % (len(g0a["mismatch"]), g0a["no_ref"], g0a["gcs_diff"]))
else:
    print("  G0a 없음(세트 1 은 해부 14 G0a 에서 920 팔 대조됨 — 여기선 점추정 재현을 본다)")
ndead = sum(len(v) for b in DEAD for v in DEAD[b].values())
print("  끼운 계통이 사라진 팔 %d (그 칸만 뺀다 · 서술)" % ndead)
for b in DEAD:
    for c, v in DEAD[b].items():
        if v: print("     %s %s: 배경 %s" % (b, c, v))
print("  점검 문제 %d" % len(bad))
for x in bad[:15]: print("    ", x)
OUT.update({"G0a": {k: (len(v) if isinstance(v, list) else v) for k, v in g0a.items()}, "dead": {"%s|%s" % (b, c): v for b in DEAD for c, v in DEAD[b].items() if v}})
if bad and not A.pilot:
    print("🚨 점검 실패 — 판정하지 않는다"); OUT["bad"] = bad
    json.dump(OUT, open(A.out or "_RESULT_pairs24_set%s.json" % SET, "w"), indent=1, ensure_ascii=False); sys.exit(2)


def pooled(b, c):                           # pairs14_analyze.py 222~225행 그대로
    ds = IB[b][c].values()
    ok = sum(d["copy_ok"] for d in ds); st = sum(d["copy_stall"] for d in ds); de = sum(d["copy_del"] for d in ds); ph = sum(d["copy_phase"] for d in ds)
    return ((ok + st + de) / ok if ok else float("nan"), ph / ok if ok else float("nan"))


# ---------------- 자기 검산 ----------------
print("\n== 자기 검산 — pooled T/R − L/k ∈ [%.1f, %.1f] (k = 고리당 복사 · 해부 14 밑 코드는 k = 1)" % (TR_LO, TR_HI))
tr_fail = []
for b in BASES:
    codes = inserts(b, "n")
    for c, L in codes.items():
        if not IB[b][c]: continue
        k = copies_per_loop(c); R, T = pooled(b, c); dev = T / R - L / k
        ok = TR_LO <= dev <= TR_HI
        if not ok: tr_fail.append((b, c, dev))
        print("   %-10s L %d k %d · R %.3f · T %.2f · T/R %.3f · 차 %+.3f %s" % (c, L, k, R, T, T / R, dev, "✅" if ok else "❌"))
        OUT["code|%s|%s" % (b, c)] = {"L": L, "k": k, "R": R, "T": T, "TR": T / R, "S": (R - 1) * L / k, "n_bg": len(IB[b][c])}
OUT["tr_fail"] = tr_fail
if tr_fail and not A.pilot:
    print("🚨 자기 검산 실패 %d 코드 — 고리 배정 의심 · 판정하지 않는다" % len(tr_fail))
    json.dump(OUT, open(A.out or "_RESULT_pairs24_set%s.json" % SET, "w"), indent=1, ensure_ascii=False); sys.exit(2)


def common_seeds(b):
    codes = inserts(b, "n")
    return sorted(set.intersection(*[set(IB[b][c]) for c in codes if IB[b][c]])) if any(IB[b][c] for c in codes) else []


def grade_effects(b, seeds):                # 배경별 등급 평균 ln R 차 — pairs14 A1(246~263행)과 같은 계산
    codes = inserts(b, "n"); Ls = sorted(set(codes.values())); out = {}
    for L1, L2 in zip(Ls, Ls[1:]):
        c1 = [c for c in codes if codes[c] == L1]; c2 = [c for c in codes if codes[c] == L2]
        out[(L1, L2)] = [sum(math.log(IB[b][c][s]["R"]) for c in c1) / len(c1) - sum(math.log(IB[b][c][s]["R"]) for c in c2) / len(c2) for s in seeds]
    return out


# ---------------- B1 ----------------
print("\n== B1 (= 해부 14 A1) — 이웃 고리 등급의 배경별 등급 평균 ln R 차(짧은 − 긴) · 하한 > 0")
b1_all = True
for b in BASES:
    seeds = common_seeds(b); codes = inserts(b, "n")
    print("  밑 코드 `%s` · 등급 %s · k %d · 공통 배경 %d" % (b, sorted(set(codes.values())), copies_per_loop(b), len(seeds)))
    for (L1, L2), v in grade_effects(b, seeds).items():
        if not v: continue
        m, lo, hi = boot_mean(v, "B1|%s|%d|%d" % (b, L1, L2)); ok = lo > 0; b1_all &= ok
        print("     ln R(%d) − ln R(%d) = %s %s" % (L1, L2, fmt(m, lo, hi), "✅" if ok else "❌"))
        OUT["B1|%s|%d|%d" % (b, L1, L2)] = {"d": m, "ci": [lo, hi], "pass": ok, "n_bg": len(v)}
    Sc = {}
    for L in sorted(set(codes.values())):
        vals = [OUT["code|%s|%s" % (b, c)]["S"] for c in codes if codes[c] == L and ("code|%s|%s" % (b, c)) in OUT]
        if vals: Sc[L] = sum(vals) / len(vals)
    if Sc:
        print("     A1b 서술: 등급 평균 S = %s · max/min %.3f" % (" · ".join("L%d %.2f" % (L, v) for L, v in sorted(Sc.items())), max(Sc.values()) / min(Sc.values())))
        OUT["A1b|%s" % b] = {"S": {str(k): v for k, v in Sc.items()}, "ratio": max(Sc.values()) / min(Sc.values())}
OUT["B1"] = b1_all
print("  B1 종합: %s" % ("✅" if b1_all else "❌"))

# ---------------- B2 · B3 ----------------
print("\n== B2 (= A2a) |ln(R_i/R_j)| ≤ %.2f · B3 (= A2b) |T/R_i − T/R_j| ≤ %.2f 틱 — 모든 같은-고리 쌍 (★ = 그 세트의 N3 실패 쌍)" % (TAU, TR_PAIR))
b2_all, b3_all, n_pairs, absd = True, True, 0, []
for b in BASES:
    codes = inserts(b, "n")
    for a, c in same_loop_pairs(codes):
        seeds = sorted(set(IB[b][a]) & set(IB[b][c]))
        if not seeds: continue
        v = [math.log(IB[b][a][s]["R"]) - math.log(IB[b][c][s]["R"]) for s in seeds]
        m, lo, hi = boot_mean(v, "B2|%s|%s|%s" % (b, a, c))
        ka, kc = OUT["code|%s|%s" % (b, a)], OUT["code|%s|%s" % (b, c)]
        trd = abs(ka["TR"] - kc["TR"])
        ok2 = abs(m) <= TAU; ok3 = trd <= TR_PAIR; b2_all &= ok2; b3_all &= ok3; n_pairs += 1; absd.append(abs(m))
        star = "★" if (b, a, c) in N3 and not A.pilot else " "
        print("   %s %-10s vs %-10s (L %d) ln R 비 = %s · 배경 %d %s · T/R 차 %.3f %s" % (star, a, c, codes[a], fmt(m, lo, hi), len(v), "✅" if ok2 else "❌", trd, "✅" if ok3 else "❌"))
        OUT["B2|%s|%s|%s" % (b, a, c)] = {"d": m, "ci": [lo, hi], "n_bg": len(v), "pass": ok2, "tr_diff": trd, "pass_b3": ok3}
if absd:
    print("  같은-고리 쌍 %d · |ln 비| 중앙값 %.4f · 90분위 %.4f · 최대 %.4f" % (n_pairs, pct(absd, 0.5), pct(absd, 0.9), max(absd)))
OUT.update({"B2": b2_all, "B3": b3_all, "n_pairs": n_pairs})
print("  B2 종합: %s · B3 종합: %s" % ("✅" if b2_all else "❌", "✅" if b3_all else "❌"))

# ---------------- B4 ----------------
print("\n== B4 🆕 — 밑 코드마다 Q = max 같은-고리 |ln R 비| ÷ min 이웃 등급 효과 · 공통 배경 재추출 상한 < %.2f" % Q_MAX)
b4_all = True
for b in BASES:
    codes = inserts(b, "n"); seeds = common_seeds(b); sp = same_loop_pairs(codes)
    if not sp or len(set(codes.values())) < 2 or len(seeds) < 2:
        print("   `%s`: 같은-고리 쌍 또는 이웃 등급 없음 — 해당 없음" % b); continue
    pv = {(a, c): [math.log(IB[b][a][s]["R"]) - math.log(IB[b][c][s]["R"]) for s in seeds] for a, c in sp}
    ge = grade_effects(b, seeds)

    def qstat(ix):
        mx = max(abs(sum(v[i] for i in ix) / len(ix)) for v in pv.values())
        mn = min(sum(v[i] for i in ix) / len(ix) for v in ge.values())
        return mx / mn if mn > 0 else float("inf")
    n = len(seeds); q = qstat(list(range(n)))
    rng = random.Random("%d|B4|%s" % (RNG_SEED, b))
    bs = [qstat([rng.randrange(n) for _ in range(n)]) for _ in range(N_BOOT)]
    lo, hi = pct(bs, 0.025), pct(bs, 0.975); ok = hi < Q_MAX; b4_all &= ok
    print("   `%s` Q = %.3f [%.3f, %.3f] · 공통 배경 %d · 같은-고리 쌍 %d %s" % (b, q, lo, hi, n, len(sp), "✅" if ok else "❌"))
    OUT["B4|%s" % b] = {"Q": q, "ci": [lo, hi], "pass": ok, "n_bg": n}
OUT["B4"] = b4_all
print("  B4 종합: %s" % ("✅" if b4_all else "❌"))

# ---------------- B5 · B6 · A3 꼴 ----------------
print("\n== B5 🆕 — N3 실패 쌍마다 배경별 ln(T/R)_a − ln(T/R)_c (= ln T 비 − ln R 비) 의 구간 ⊂ [−%.2f, +%.2f]  ·  B6 서술 부호 정합" % (D_MAX, D_MAX))
b5_all, b5_n, coh, coh_res, n_res = True, 0, 0, 0, 0
for (b, a, c), dd in sorted(N3.items()):
    if b not in IB or a not in IB[b] or c not in IB[b]:
        continue
    seeds = sorted(set(IB[b][a]) & set(IB[b][c]))
    if len(seeds) < 2: continue
    v5 = [math.log(IB[b][a][s]["T"] / IB[b][a][s]["R"]) - math.log(IB[b][c][s]["T"] / IB[b][c][s]["R"]) for s in seeds]
    m5, lo5, hi5 = boot_mean(v5, "B5|%s|%s|%s" % (b, a, c)); ok5 = -D_MAX <= lo5 and hi5 <= D_MAX
    b5_all &= ok5; b5_n += 1
    row = {}
    for name, f in [("R", lambda d: d["R"]), ("T", lambda d: d["T"]), ("births", lambda d: d["births"] / d["expo"]), ("deaths", lambda d: d["deaths"] / d["expo"]),
                    ("prestep", lambda d: d["death_prestep"] / d["expo"])]:
        v = []
        for s in seeds:
            x, y = f(IB[b][a][s]), f(IB[b][c][s])
            if x > 0 and y > 0: v.append(math.log(x) - math.log(y))
        row[name] = boot_mean(v, "A3|%s|%s|%s|%s" % (name, b, a, c)) if v else (float("nan"),) * 3
    rec = {"B5": {"d": m5, "ci": [lo5, hi5], "pass": ok5}, "dd": dd, **{k: {"d": x[0], "ci": [x[1], x[2]]} for k, x in row.items()}}
    line = "   %-8s %-10s / %-10s B5 %s %s · ln R 비 %s · ln T 비 %s · 출생 %s · 사망 %s" % (b, a, c, fmt(m5, lo5, hi5), "✅" if ok5 else "❌",
                                                                                 fmt(*row["R"]), fmt(*row["T"]), fmt(*row["births"]), fmt(*row["deaths"]))
    if dd is not None:
        # B6: R 이 큰 쪽(헛손질이 많은 쪽)이 침입 적합도가 낮은가 — ln R(a/c) > 0 이면 s(c) − s(a) > 0 이어야 정합
        cohere = (row["R"][0] > 0) == (dd > 0); resolved = row["R"][1] > 0 or row["R"][2] < 0
        coh += cohere; n_res += resolved; coh_res += (cohere and resolved)
        rec.update({"coherent": cohere, "R_resolved": resolved, "dd_over_lnR": dd / row["R"][0] if row["R"][0] else None})
        line += " · ΔΔ s(c)−s(a) %+.3f · 부호 %s%s" % (dd, "정합" if cohere else "어긋남", "" if resolved else "(ln R 구간이 0 을 품음)")
    print(line)
    OUT["N3|%s|%s|%s" % (b, a, c)] = rec
OUT["B5"] = b5_all if b5_n else None
print("  B5 종합: %s" % (("✅" if b5_all else "❌") + " (%d 쌍)" % b5_n if b5_n else "해당 없음(N3 실패 쌍 0)"))
if not A.pilot and N3:
    print("  B6 서술: 부호 정합 %d/%d · ln R 구간이 0 을 안 품는 쌍 %d 중 정합 %d" % (coh, b5_n, n_res, coh_res))
    OUT["B6"] = {"coherent": coh, "n": b5_n, "resolved": n_res, "coherent_resolved": coh_res}

# ---------------- A2c 점추정 서술 ----------------
for b in BASES:
    eff = [v["d"] for k, v in OUT.items() if k.startswith("B1|%s|" % b)]
    same = [abs(v["d"]) for k, v in OUT.items() if k.startswith("B2|%s|" % b)]
    if eff and same and min(eff) > 0:
        OUT["A2c|%s" % b] = {"max_same": max(same), "min_eff": min(eff), "ratio": max(same) / min(eff)}

# ---------------- 세트 1 양성 대조: 해부 14 점추정 재현 ----------------
if SET == 1:
    z14 = json.load(open(os.path.join(BASE, "_RESULT_pairs14.json")))
    diffs = []
    for k, v in OUT.items():
        if k.startswith("code|"):
            r = z14.get(k)
            if r is None or abs(r["R"] - v["R"]) > 1e-9 or abs(r["T"] - v["T"]) > 1e-9: diffs.append(k)
        elif k.startswith("B1|"):
            r = z14.get("A1|" + k[3:])
            if r is None or abs(r["d"] - v["d"]) > 1e-9: diffs.append(k)
        elif k.startswith("B2|"):
            r = z14.get("A2|" + k[3:])
            if r is None or abs(r["d"] - v["d"]) > 1e-9 or abs(r["tr_diff"] - v["tr_diff"]) > 1e-9: diffs.append(k)
    print("\n== 양성 대조 — 해부 14 점추정(R · T · A1 · A2a · A2b) 재현: 다른 항목 %d %s" % (len(diffs), "✅" if not diffs else "❌ " + " ".join(diffs[:5])))
    OUT["repro14"] = {"diffs": diffs}

# ---------------- 종합 ----------------
print("\n== 종합 (세트 %s)" % SET)
if A.pilot:
    print("  파일럿 — 판정 없음. 위 값은 경로 연습 · 판정선을 바꾸는 근거로 쓰지 않는다(사전등록 §8).")
else:
    core = OUT["B2"] and OUT["B4"] and (OUT["B5"] is not False)
    print("  G0 ✅ · 자기 검산 ✅ · B1 %s · B2 %s · B3 %s · B4 %s · B5 %s" % tuple(("✅" if x else ("해당 없음" if x is None else "❌")) for x in (OUT["B1"], OUT["B2"], OUT["B3"], OUT["B4"], OUT["B5"])))
    print("  핵심(B2 · B4 · B5) %s%s" % ("✅" if core else "❌", " — 세트 1 은 양성 대조(판정 밖)" if SET == 1 else ""))
    OUT["core"] = core
OUT["bad"] = bad; OUT["notes"] = notes
fn = A.out or ("_RESULT_pilot24.json" if A.pilot else "_RESULT_pairs24_set%s.json" % SET)
json.dump(OUT, open(fn, "w"), indent=1, ensure_ascii=False)
print("기록 →", fn)
