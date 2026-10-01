#!/usr/bin/env python3
"""해부 23 독립 재계수 (판정기 noise23_analyze.py 를 열지 않고 사전등록 노트 §2 · §4 정의만으로 새로 씀).

실행: cd ~/petri && python3 /path/recount23_independent.py [BASE] > recount23.txt   (JSON 은 stdout 끝 'JSON:' 뒤)
- 원자료: inv15B_<b>/mu0 · inv15E_<b>/mu0 · inv20P_<b>/mu0 (*_list.json) · inv12nbr_racld/mu0 · main/
- 발표 N3 값(G0-c 대조용): _RESULT_pairs13.json · _RESULT_pairs13_racldfam.json · _RESULT_pairs15.json(E2|) · _RESULT_pairs20.json(neg|)
- 판정 출력(_RESULT_noise23*)은 읽지 않는다.
구현 메모(판정기와 일부러 다를 수 있는 곳):
- s 는 노트 §2 식 그대로 새로 씀(series[0]/[-1] 의 열 1 = n, 열 2 = m).
- 같은 고리 판정: '고리' = 첫 l 앞 마지막 s 다음부터 l 까지의 x 아닌 글자 수(c 가 있어야) — 해부 13·15·20 정의를 내 손으로 다시 씀.
- 재추출: random.Random("20261023|floor|<세트>|<밑>") 로 배경 20 개를 복원추출 2,000 번, 매번 σ̄² 재계산, 선형보간 분위수.
  난수 소비 순서가 판정기와 다를 수 있으므로 q95_hi 는 소수 셋째 자리 정도의 차는 '구현 차' 로 본다.
- R_var = Var_배경(d, ddof 1) / (2σ̄²).
"""
import glob, json, math, os, random, statistics, sys

BASE = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1 else "~/petri")
SETS = {
    "S1": ("inv15B", list(range(3, 23)), "racld rascld rsacld racldx acld rascled".split()),
    "S2": ("inv15E", [26, 27, 28, 30, 31, 32, 33, 35, 36, 38, 39, 40, 41, 42, 43, 44, 45, 47, 48, 50],
           "racld rascld rsacld racldx acld rascled".split()),
    "S3": ("inv20P", [76, 77, 79, 80, 81, 83, 84, 85, 86, 87, 88, 89, 92, 93, 94, 95, 96, 98, 99, 100],
           "rcald crald reasccld racld rascld rsacld racldx rascled".split()),
}
TOL = 0.15
Z = 1.959964
NB = 2000
problems = []


def loop_len(code):
    # 첫 'l' 까지, 그 앞 마지막 's' 다음부터의 구간에서 x 아닌 글자 수 · 구간에 c 가 있어야
    if "l" not in code or "c" not in code:
        return None
    li = code.index("l")
    st = code.rfind("s", 0, li) + 1
    seg = code[st:li + 1]
    if "c" not in seg:
        return None
    return len([ch for ch in seg if ch != "x"])


def n_inserts(b):
    out = []
    for p in range(len(b) + 1):
        k = b[:p] + "n" + b[p:]
        if loop_len(k) is not None and k not in out:
            out.append(k)
    return out


def s_arm(a):
    sr = a["series"]
    n0, m0 = sr[0][1], sr[0][2]
    nT, mT = sr[-1][1], sr[-1][2]
    return math.log((mT + 0.5) / m0) - math.log((nT - mT + 0.5) / (n0 - m0))


def q(xs, p):
    xs = sorted(xs)
    k = (len(xs) - 1) * p
    f = int(math.floor(k))
    c = min(f + 1, len(xs) - 1)
    return xs[f] + (xs[c] - xs[f]) * (k - f)


def mean(v):
    return sum(v) / len(v)


def var1(v):
    m = mean(v)
    return sum((x - m) ** 2 for x in v) / (len(v) - 1)


def Phi(x):
    return 0.5 * (1 + math.erf(x / math.sqrt(2)))


# ---- 배경 목록 규칙 점검(main/) ----
def rule_ok(s):
    z = json.load(open(os.path.join(BASE, "main", "mat_d8_mu0p01_s%05d.json" % s)))
    return z["extinct_at"] < 0 and bool(z["final_counts"]) and z["final_counts"][0][0] == "racld"


qual = [s for s in range(1, 101) if rule_ok(s)]
rule_lists = {
    "S1(1부터 오름 · 앞 20)": [s for s in qual if s >= 1][:20],
    "S2(26부터 오름 · 앞 20)": [s for s in qual if s >= 26][:20],
    "S3(100부터 내림 · 앞 20)": sorted(sorted([s for s in qual if s >= 53], reverse=True)[:20]),
}
for (nm, lst), st in zip(rule_lists.items(), ["S1", "S2", "S3"]):
    ok = lst == SETS[st][1]
    print("배경 규칙 %s = 하드코딩 %s : %s" % (nm, st, "일치" if ok else "불일치 %s" % lst))
    if not ok:
        problems.append("배경 규칙 %s" % st)
allbg = [s for st in SETS for s in SETS[st][1]]
if len(allbg) != len(set(allbg)):
    problems.append("세트끼리 배경 겹침")
print("세트끼리 배경 겹침: %s" % ("없음" if len(allbg) == len(set(allbg)) else "있음"))

# ---- 원자료 읽기 ----
CELL = {}
nfiles = 0
gck = {}
for st, (root, seeds, bases) in SETS.items():
    for b in bases:
        fs = sorted(glob.glob(os.path.join(BASE, "%s_%s" % (root, b), "mu0", "dms_mat_%s_s*_list.json" % b)))
        codes = n_inserts(b)
        wt = {}
        mut = {}
        K = 0
        crn_raw = {}
        seen = []
        for f in fs:
            z = json.load(open(f))
            nfiles += 1
            s = z["seed"]
            seen.append(s)
            if z["wt"] != b or z["arms_mode"] != "list" or abs(z["mu_assay"]) > 0 or z["sim_version"] != "0.3.5":
                problems.append("설정 %s %s %d" % (st, b, s))
            if not z["checksum_match"] or not z["fidelity"]["same"] or z["grow_extinct"] >= 0 or z["grow_top"][0][0] != "racld":
                problems.append("배경 %s %s %d" % (st, b, s))
            if not all(a["cons_ok"] for a in z["arms"]):
                problems.append("보존 %s %s %d" % (st, b, s))
            gck.setdefault((st, s), set()).add(z["grow_checksum"])
            ref = [a for a in z["arms"] if a["kind"] == "ref"]
            neu = [a for a in z["arms"] if a["kind"] == "neutral"]
            mu = [a for a in z["arms"] if a["kind"] == "mut"]
            if len(ref) != 1 or len(neu) != 10 or ref[0]["key"] != b or any(a["key"] != b for a in neu):
                problems.append("WT 팔 %s %s %d" % (st, b, s))
            if [a["labels"] for a in neu] != [["neutral:%d" % i] for i in range(1, 11)]:
                problems.append("neutral 이름 %s %s %d" % (st, b, s))
            ps = [a["post_seed"] for a in ref + neu]
            if len(set(ps)) != 11 or any(a["post_seed"] != ref[0]["post_seed"] for a in mu):
                problems.append("뒤 난수 %s %s %d" % (st, b, s))
            if len(set(json.dumps(a["inj"], sort_keys=True) for a in z["arms"])) != 1:
                problems.append("inj %s %s %d" % (st, b, s))
            if [a["key"] for a in mu] != codes:
                problems.append("mut 목록 %s %s %d" % (st, b, s))
            col = any(a["series"][-1][2] < 0.10 * a["series"][0][2] for a in z["arms"])
            K += col
            wt[s] = [s_arm(a) for a in ref + neu]
            mut[s] = {a["key"]: s_arm(a) for a in mu}
            crn_raw[s] = (s_arm(ref[0]), [s_arm(a) for a in neu])
        if sorted(seen) != seeds:
            problems.append("파일 시드 %s %s" % (st, b))
        CELL[(st, b)] = {"wt": wt, "mut": mut, "K": K, "codes": codes, "crn_raw": crn_raw}
ngc = sum(1 for v in gck.values() if len(v) > 1)
print("읽은 원자료 파일 %d · 같은 세트 같은 배경 grow_checksum 불일치 %d" % (nfiles, ngc))
if ngc:
    problems.append("grow_checksum")

# ---- G0-b ----
gb_ok = 0
gb_n = 0
for f in sorted(glob.glob(os.path.join(BASE, "inv12nbr_racld", "mu0", "*.json"))):
    z = json.load(open(f))
    if z["seed"] not in SETS["S1"][1]:
        continue
    gb_n += 1
    ref = [a for a in z["arms"] if a["kind"] == "ref"][0]
    slf = [a for a in z["arms"] if a["kind"] == "mut" and a["key"] == "racld"]
    if slf and slf[0]["post_seed"] == ref["post_seed"] and slf[0]["series"] == ref["series"]:
        gb_ok += 1
print("G0-b 같은 코드 같은 뒤 난수 동일 궤적 %d/%d" % (gb_ok, gb_n))
if gb_ok != 20 or gb_n != 20:
    problems.append("G0-b")


def same_loop(codes):
    L = {c: loop_len(c) for c in codes}
    cs = sorted(codes)
    return [(a, c) for i, a in enumerate(cs) for c in cs[i + 1:] if L[a] == L[c]]


# ---- 셈 ----
RES = {}
for (st, b), C in CELL.items():
    seeds = SETS[st][1]
    wt = C["wt"]
    sig2 = mean([var1(wt[s]) for s in seeds])
    sd0 = math.sqrt(2 * sig2 / len(seeds))
    q95 = Z * sd0
    rng = random.Random("20261023|floor|%s|%s" % (st, b))
    bq = []
    for _ in range(NB):
        pick = [seeds[rng.randrange(len(seeds))] for _ in seeds]
        sg = mean([var1(wt[s]) for s in pick])
        bq.append(Z * math.sqrt(2 * sg / len(seeds)))
    lo, hi = q(bq, 0.025), q(bq, 0.975)
    # 경험 55 쌍
    fl = []
    for i in range(11):
        for j in range(i + 1, 11):
            fl.append(mean([wt[s][j] - wt[s][i] for s in seeds]))
    q_emp = q([abs(x) for x in fl], 0.95)
    exc = sum(1 for x in fl if abs(x) > q95) / len(fl)
    ffr = 2 * (1 - Phi(TOL / sd0))
    ffr_emp = sum(1 for x in fl if abs(x) > TOL) / len(fl)
    pairs = {}
    for a, c in same_loop(C["codes"]):
        d = [C["mut"][s][c] - C["mut"][s][a] for s in seeds]
        dd = mean(d)
        x = abs(dd)
        rv = var1(d) / (2 * sig2) if sig2 > 0 else float("nan")
        cls = "N" if x <= q95 else ("A" if x <= hi else "P")
        rob = cls == "P" and rv == rv and x > hi * max(1.0, math.sqrt(rv))
        pairs[(a, c)] = {"dd": dd, "x": x, "Rvar": rv, "cls": cls, "rob": rob, "fail": x > TOL}
    # CRN 비
    crn = []
    for code in C["codes"]:
        vr = var1([C["mut"][s][code] - C["crn_raw"][s][0] for s in seeds])
        vn = mean([var1([C["mut"][s][code] - C["crn_raw"][s][1][k] for s in seeds]) for k in range(10)])
        crn.append(vr / vn)
    RES[(st, b)] = {"K": C["K"], "judged": C["K"] == 0, "sig2": sig2, "sd0": sd0, "q95": q95, "lo": lo, "hi": hi,
                    "q_emp": q_emp, "exc": exc, "ffr": ffr, "ffr_emp": ffr_emp, "pairs": pairs,
                    "crn_med": statistics.median(crn), "crn_lt1": sum(1 for r in crn if r < 1) / len(crn),
                    "rv_nonfail_med": (statistics.median([p["Rvar"] for p in pairs.values() if not p["fail"]])
                                       if any(not p["fail"] for p in pairs.values()) else float("nan"))}

# ---- G0-c 발표 N3 대조 ----
def P(f):
    return json.load(open(os.path.join(BASE, f)))


pub = {}
for f in ["_RESULT_pairs13.json", "_RESULT_pairs13_racldfam.json"]:
    for k, v in P(f).items():
        p = k.split("|")
        if p[0] == "neg" and len(p) == 5 and p[2] == "n":
            a, c = sorted([p[3], p[4]])
            sign = -1 if (p[3], p[4]) == (a, c) else 1  # 해부 13 = s(앞키) − s(뒤키)
            pub[("S1", p[1], a, c)] = sign * v["dd"]
for k, v in P("_RESULT_pairs15.json").items():
    p = k.split("|")
    if p[0] == "E2" and len(p) == 4:
        a, c = sorted([p[2], p[3]])
        pub[("S2", p[1], a, c)] = v["dd"] * (1 if (p[2], p[3]) == (a, c) else -1)
for k, v in P("_RESULT_pairs20.json").items():
    p = k.split("|")
    if p[0] == "neg" and len(p) == 4:
        a, c = sorted([p[2], p[3]])
        pub[("S3", p[1], a, c)] = v["dd"] * (1 if (p[2], p[3]) == (a, c) else -1)
maxdiff = 0.0
cnt_ok = 0
for (st, b), R in RES.items():
    mine = {(a, c): R["pairs"][(a, c)]["dd"] for (a, c) in R["pairs"]}
    theirs = {(a, c): v for (s2, b2, a, c), v in pub.items() if s2 == st and b2 == b}
    if set(mine) == set(theirs):
        cnt_ok += 1
    else:
        problems.append("G0-c 쌍 집합 %s %s (내 %d · 발표 %d)" % (st, b, len(mine), len(theirs)))
    for k in mine:
        if k in theirs:
            maxdiff = max(maxdiff, abs(mine[k] - theirs[k]))
print("G0-c 칸 쌍 집합 일치 %d/%d · 최대 차 %.3g" % (cnt_ok, len(RES), maxdiff))
if maxdiff >= 1e-9:
    problems.append("G0-c 값")
print("점검 문제 %d %s" % (len(problems), problems[:10]))

# ---- 인쇄 ----
print("\n세트 밑 K 판정 | σ̄² SD0 q95 [lo, hi] q95_emp 초과율 거짓실패(정규 · 경험) | 쌍 N/A/P/Prob | 실패 N/A/P/Prob | CRN중앙 <1몫 | 실패아닌 R_var중앙")
for (st, b), R in RES.items():
    ps = R["pairs"].values()
    fs = [p for p in ps if p["fail"]]
    cn = lambda L: "%d/%d/%d/%d" % (sum(p["cls"] == "N" for p in L), sum(p["cls"] == "A" for p in L),
                                   sum(p["cls"] == "P" for p in L), sum(p["rob"] for p in L))
    print("%s %-9s K%2d %s | %.5f %.4f %.4f [%.4f, %.4f] %.4f %.3f %.2e %.3f | %d %s | %d %s | %.3f %.2f | %.3f" % (
        st, b, R["K"], "판정" if R["judged"] else "제외", R["sig2"], R["sd0"], R["q95"], R["lo"], R["hi"], R["q_emp"],
        R["exc"], R["ffr"], R["ffr_emp"], len(R["pairs"]), cn(list(ps)), len(fs), cn(fs), R["crn_med"], R["crn_lt1"],
        R["rv_nonfail_med"]))

print("\nF (판정 칸 실패 쌍):")
F = []
for (st, b), R in RES.items():
    if not R["judged"]:
        continue
    for (a, c), p in sorted(R["pairs"].items()):
        if p["fail"]:
            F.append((st, b, a, c))
            print("  %s %-9s %-10s/%-10s ΔΔ %+.4f  %s  R_var %.3f  rob-문턱 %.4f  %s" % (
                st, b, a, c, p["dd"], p["cls"], p["Rvar"], R["hi"] * max(1, math.sqrt(p["Rvar"])), "P_rob" if p["rob"] else "-"))
print("F = %d (S1 %d · S2 %d · S3 %d)" % (len(F), sum(f[0] == "S1" for f in F), sum(f[0] == "S2" for f in F), sum(f[0] == "S3" for f in F)))

ALL3 = []
for b in SETS["S1"][2]:
    if all((st, b) in RES and RES[(st, b)]["judged"] for st in SETS):
        for k in RES[("S1", b)]["pairs"]:
            if all(k in RES[(st, b)]["pairs"] for st in SETS):
                ALL3.append((b,) + k)
vplus, vb, weak = [], [], []
for b, a, c in ALL3:
    pp = [RES[(st, b)]["pairs"][(a, c)] for st in SETS]
    sg = set(1 if p["dd"] > 0 else -1 for p in pp)
    if all(p["rob"] for p in pp) and len(sg) == 1:
        vplus.append((b, a, c))
    elif all(p["cls"] == "P" for p in pp) and len(sg) == 1:
        vb.append((b, a, c))
    nP = [st for st, p in zip(SETS, pp) if p["cls"] == "P"]
    if len(nP) == 2 and len(set(1 if RES[(st, b)]["pairs"][(a, c)]["dd"] > 0 else -1 for st in nP)) == 1:
        weak.append((b, a, c, nP, [st for st in nP if RES[(st, b)]["pairs"][(a, c)]["fail"]]))
print("\nALL3 배열 %d" % len(ALL3))
for b, a, c in vplus:
    print("  V+ 배열 %-8s %-10s/%-10s  %s" % (b, a, c, " · ".join("%s %+.3f(Rv %.2f)" % (st, RES[(st, b)]["pairs"][(a, c)]["dd"], RES[(st, b)]["pairs"][(a, c)]["Rvar"]) for st in SETS)))
print("세 세트 P_rob 같은 부호 %d · 세 세트 P 같은 부호이나 P_rob 아님 %d" % (len(vplus), len(vb)))
for w in weak:
    print("  'V+ 약' %s %s/%s P 세트 %s · F 에 넣은 세트 %s" % w)
fF = [RES[(f[0], f[1])]["pairs"][(f[2], f[3])] for f in F]
if vplus:
    V = "V+"
elif vb:
    V = "V?b"
elif not F:
    V = "V?f"
elif not any(p["cls"] == "P" for p in fF) and sum(p["cls"] == "N" for p in fF) * 2 >= len(fF):
    V = "V0"
else:
    V = "V?"
print("\n종합 = %s" % V)
crnflag = [(st, b, R["crn_med"]) for (st, b), R in RES.items() if R["judged"] and R["crn_med"] < 0.8 and R["crn_lt1"] >= 0.75]
print("CRN 결합 기준(중앙 < 0.8 · 1미만 몫 ≥ 0.75) 판정 칸: %s" % crnflag)
judged = [R for R in RES.values() if R["judged"]]
print("판정 칸 q95 범위 %.4f ~ %.4f · q95_hi 최대 %.4f · 거짓 실패율 최대 %.3g · 경험 거짓 실패 >0 칸 %d" % (
    min(r["q95"] for r in judged), max(r["q95"] for r in judged), max(r["hi"] for r in judged),
    max(r["ffr"] for r in judged), sum(1 for r in judged if r["ffr_emp"] > 0)))
out = {"V": V, "problems": problems, "nfiles": nfiles, "G0b": [gb_ok, gb_n], "G0c_maxdiff": maxdiff,
       "cells": {"%s|%s" % k: {kk: vv for kk, vv in R.items() if kk != "pairs"} for k, R in RES.items()},
       "pairs": {"%s|%s|%s|%s" % (k[0], k[1], a, c): p for k, R in RES.items() for (a, c), p in R["pairs"].items()},
       "F": ["|".join(f) for f in F], "ALL3": len(ALL3), "vplus": ["|".join(v) for v in vplus],
       "vb": ["|".join(v) for v in vb], "weak": [list(w[:3]) + [w[3], w[4]] for w in weak], "crnflag": crnflag}
print("JSON:" + json.dumps(out))
