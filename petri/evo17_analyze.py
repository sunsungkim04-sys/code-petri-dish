#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""해부 17 판정 — 진화 궤적: 오류율이 인구의 복사 고리 길이를 어디로 데려가나. 사전등록 해부17-사전등록-2026-09-28.md §4.
  G0   μ 1% 시드 1~20 의 체크섬 = main/ 같은 시드 체크섬 (20/20) · loops 마지막 항목 = final_counts 로 센 히스토그램
  V1   끝점 L̄(50k) 의 μ 4 층 중앙값 단조 증가 · Spearman(μ, L̄) 부트스트랩 하한 > 0 (접시 80)
  V2   μ 0.1% 끝점 최빈 등급 = 2 ≥ 15/20 · μ 1% 최빈 등급 ≥ 4 ≥ 15/20
  V3   서술: 0.3% · 0.5% 최빈 등급 분포 · 교대 시각 t4(등급 4 몫이 등급 2 몫을 처음 넘는 틱) · L̄(t) 중앙값 곡선
  V4   서술: 끝점 e 포함 몫 · 고리 없는 개체 몫
  V5   서술: μ 1% 시드 1~20 끝점 히스토그램 대 main 100 접시 final_counts 히스토그램
실행: python3 evo17_analyze.py [petri 폴더]  (환경: PETRI_ROOT PETRI_SEEDS PETRI_MUS PETRI_PILOT)
"""
import glob
import json
import math
import os
import random
import sys

BASE = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1 else "~/petri")
ROOT = os.environ.get("PETRI_ROOT", "evo17")
SEEDS = [int(x) for x in os.environ["PETRI_SEEDS"].split()] if os.environ.get("PETRI_SEEDS") else list(range(1, 21))
MUS = [float(x) for x in os.environ["PETRI_MUS"].split()] if os.environ.get("PETRI_MUS") else [0.001, 0.003, 0.005, 0.01]
PILOT = os.environ.get("PETRI_PILOT") == "1"
V2_MIN = 15
N_BOOT = 2000
RNG = random.Random(20261001)
bad, notes = [], []
OUT = {"root": ROOT, "pilot": PILOT, "mus": MUS, "seeds": SEEDS}


def loop_ticks(code):
    if "l" not in code or "c" not in code: return None
    li = code.index("l"); si = code.rfind("s", 0, li); start = si + 1 if si >= 0 else 0
    if start > li: return None
    seg = code[start:li + 1]
    return sum(1 for ch in seg if ch != "x") if "c" in seg else None


def mtag(mu): return "mu" + str(mu).replace(".", "p")


def hist_of(counts):
    h, none = {}, 0
    for code, n in counts:
        L = loop_ticks(code)
        if L is None: none += n
        else: h[L] = h.get(L, 0) + n
    return h, none


def lbar(h):
    tot = sum(h.values())
    return sum(int(L) * n for L, n in h.items()) / tot if tot else float("nan")


def shares(h):
    tot = sum(h.values()) or 1
    return {int(L): n / tot for L, n in h.items()}


def pct(xs, p):
    xs = sorted(xs); k = (len(xs) - 1) * p; f = math.floor(k); c = min(f + 1, len(xs) - 1)
    return xs[f] + (xs[c] - xs[f]) * (k - f)


def spearman(x, y):
    def rank(a):
        o = sorted(range(len(a)), key=lambda i: a[i]); r = [0.0] * len(a); i = 0
        while i < len(o):
            j = i
            while j + 1 < len(o) and a[o[j + 1]] == a[o[i]]: j += 1
            for k in range(i, j + 1): r[o[k]] = (i + j) / 2
            i = j + 1
        return r
    rx, ry = rank(x), rank(y); n = len(x); mx = sum(rx) / n; my = sum(ry) / n
    sxy = sum((a - mx) * (b - my) for a, b in zip(rx, ry)); sxx = sum((a - mx) ** 2 for a in rx); syy = sum((b - my) ** 2 for b in ry)
    return sxy / math.sqrt(sxx * syy) if sxx and syy else float("nan")


# ---------------- 읽기 ----------------
D = {}   # (mu, seed) -> dict
EXTINCT = []
g0 = {"checked": 0, "mismatch": [], "self": 0, "self_bad": []}
for mu in MUS:
    for s in SEEDS:
        f = os.path.join(BASE, ROOT, "mat_d8_%s_s%05d.json" % (mtag(mu), s))
        if not os.path.exists(f): bad.append("없음 %s" % f); continue
        z = json.load(open(f))
        # 09-28 판정 전 개정 1회: 사전등록이 초기 멸종을 안 다뤘다 — 멸종 접시(시조가 자리 잡기 전 사라짐)는 점검 실패가 아니라 **제외**하고 개수를 적는다.
        #   G0 체크섬 대조는 멸종 접시에도 한다(main 도 같은 틱에 멸종해야 한다). V2 의 절대 개수(≥ 15/20)는 그대로라 멸종은 통과에 불리하다.
        if z["extinct_at"] >= 0:
            EXTINCT.append((mu, s, z["extinct_at"]))
            if abs(mu - 0.01) < 1e-12:
                mf = os.path.join(BASE, "main", "mat_d8_mu0p01_s%05d.json" % s)
                if os.path.exists(mf):
                    g0["checked"] += 1
                    if json.load(open(mf))["checksum"] != z["checksum"]: g0["mismatch"].append(s)
            continue
        if not z.get("loops_on") or abs(z["opts"]["mu"] - mu) > 1e-12 or z["ticks_run"] != 50000: bad.append("설정 %s %d" % (mu, s)); continue
        loops = z["loops"]
        hf, nf = hist_of(z["final_counts"])
        last = loops[-1]
        g0["self"] += 1
        if last[0] != z["ticks_run"] or {int(k): v for k, v in last[1].items()} != hf or last[2] != nf: g0["self_bad"].append("%s %d" % (mu, s))
        if abs(mu - 0.01) < 1e-12:
            mf = os.path.join(BASE, "main", "mat_d8_mu0p01_s%05d.json" % s)
            if os.path.exists(mf):
                g0["checked"] += 1
                if json.load(open(mf))["checksum"] != z["checksum"]: g0["mismatch"].append(s)
        traj = [(t, {int(k): v for k, v in h.items()}, none) for t, h, none in loops]
        pop_e = z["samples"][-1][6] / z["samples"][-1][1] if z["samples"][-1][1] else float("nan")
        # 교대 시각
        t4 = None
        for t, h, none in traj:
            sh = shares(h)
            if sh.get(4, 0) + sh.get(5, 0) > sh.get(2, 0): t4 = t; break
        end_sh = shares(hf); mode = max(end_sh, key=end_sh.get) if end_sh else None
        D[(mu, s)] = {"lbar_end": lbar(hf), "shares_end": end_sh, "mode": mode, "t4": t4, "none_end": nf / sum(hf.values()) if hf else float("nan"),
                      "e_end": pop_e, "traj": [(t, lbar(h)) for t, h, none in traj], "traj_sh": [(t, shares(h)) for t, h, none in traj]}
print("== 해부 17 %s— 진화 궤적 · 루트 %s · μ %s · 시드 %d" % ("파일럿 " if PILOT else "판정 ", ROOT, MUS, len(SEEDS)))
print("  G0 체크섬 대조(main μ 1%%): %d · 불일치 %s · 자기 검산(loops 끝 = final_counts) %d · 불일치 %s" % (g0["checked"], g0["mismatch"], g0["self"], g0["self_bad"]))
if g0["mismatch"] or g0["self_bad"]: bad.append("G0")
print("  멸종 접시(제외) %d: %s" % (len(EXTINCT), ["μ %g 시드 %d @%d" % e for e in EXTINCT]))
OUT["extinct"] = EXTINCT
print("  점검 문제 %d" % len(bad))
for x in bad[:10]: print("    ", x)
if bad and not PILOT: print("🚨 점검 실패 — 판정하지 않는다"); sys.exit(1)

# ---------------- V1 ----------------
print("\n== V1 — 끝점 L̄(50k) 대 μ")
med = {}
for mu in MUS:
    v = [D[(mu, s)]["lbar_end"] for s in SEEDS if (mu, s) in D]
    if not v: continue
    med[mu] = pct(v, .5)
    print("   μ %.1f%%: L̄ 중앙값 %.3f · 사분위 [%.3f, %.3f] · 접시 %d" % (mu * 100, med[mu], pct(v, .25), pct(v, .75), len(v)))
    OUT["V1|%s" % mu] = {"median": med[mu], "q1": pct(v, .25), "q3": pct(v, .75), "n": len(v), "values": v}
mono = all(med[a] < med[b] for a, b in zip(MUS, MUS[1:]) if a in med and b in med)
xs = [mu for mu in MUS for s in SEEDS if (mu, s) in D]; ys = [D[(mu, s)]["lbar_end"] for mu in MUS for s in SEEDS if (mu, s) in D]
rho = spearman(xs, ys) if len(xs) > 3 else float("nan")
bs = []
for _ in range(N_BOOT):
    idx = [RNG.randrange(len(xs)) for _ in xs]; bs.append(spearman([xs[i] for i in idx], [ys[i] for i in idx]))
lo, hi = pct(bs, .025), pct(bs, .975)
v1 = mono and lo > 0
print("   단조 %s · Spearman %.3f [%.3f, %.3f] · 접시 %d → V1 %s" % ("✅" if mono else "❌", rho, lo, hi, len(xs), "✅" if v1 else "❌"))
OUT["V1"] = {"mono": mono, "rho": rho, "ci": [lo, hi], "pass": v1}

# ---------------- V2 · V3 ----------------
print("\n== V2 — 끝점 최빈 등급 (μ 0.1%% → 2 · μ 1%% → ≥ 4 · 각 ≥ %d/20) · V3 서술" % V2_MIN)
v2 = True
for mu in MUS:
    modes = [D[(mu, s)]["mode"] for s in SEEDS if (mu, s) in D]
    cnt = {}
    for m in modes: cnt[m] = cnt.get(m, 0) + 1
    t4s = [D[(mu, s)]["t4"] for s in SEEDS if (mu, s) in D]
    t4v = [t for t in t4s if t is not None]
    line = "   μ %.1f%%: 최빈 등급 %s · 교대 t4 %s" % (mu * 100, dict(sorted(cnt.items(), key=lambda kv: -kv[1])),
                                              ("중앙 %d · 범위 %d~%d · 없음 %d" % (pct(t4v, .5), min(t4v), max(t4v), len(t4s) - len(t4v))) if t4v else "없음 %d/%d" % (len(t4s), len(t4s)))
    if abs(mu - 0.001) < 1e-12:
        ok = cnt.get(2, 0) >= V2_MIN; v2 &= ok; line += " → 등급 2 %d/%d %s" % (cnt.get(2, 0), len(modes), "✅" if ok else "❌")
    if abs(mu - 0.01) < 1e-12:
        n4 = sum(n for m, n in cnt.items() if m is not None and m >= 4); ok = n4 >= V2_MIN; v2 &= ok; line += " → 등급 ≥4 %d/%d %s" % (n4, len(modes), "✅" if ok else "❌")
    print(line)
    OUT["V2|%s" % mu] = {"modes": {str(k): v for k, v in cnt.items()}, "t4": t4s}
print("  V2 종합: %s" % ("✅" if v2 else "❌")); OUT["V2"] = v2

# ---------------- L̄(t) 중앙값 곡선 · V4 · V5 ----------------
print("\n== L̄(t) 중앙값 (2,500틱마다) · V4 끝점 e 몫 · 고리 없는 몫")
for mu in MUS:
    ks = [s for s in SEEDS if (mu, s) in D]
    if not ks: continue
    T = [t for t, _ in D[(mu, ks[0])]["traj"]]
    curve = [pct([dict(D[(mu, s)]["traj"])[t] for s in ks], .5) for t in T]
    OUT["curve|%s" % mu] = {"t": T, "median": curve, "q1": [pct([dict(D[(mu, s)]["traj"])[t] for s in ks], .25) for t in T], "q3": [pct([dict(D[(mu, s)]["traj"])[t] for s in ks], .75) for t in T]}
    print("   μ %.1f%%: %s" % (mu * 100, " ".join("%.2f" % curve[i] for i in range(0, len(T), 5))))
    e = [D[(mu, s)]["e_end"] for s in ks]; nn = [D[(mu, s)]["none_end"] for s in ks]
    print("           V4: e 포함 몫 중앙 %.3f [%.3f, %.3f] · 고리 없는 몫 중앙 %.3f" % (pct(e, .5), pct(e, .25), pct(e, .75), pct(nn, .5)))
    OUT["V4|%s" % mu] = {"e_end_median": pct(e, .5), "none_end_median": pct(nn, .5)}
mf = sorted(glob.glob(os.path.join(BASE, "main", "mat_d8_mu0p01_s*.json")))
if mf:
    v100 = []
    for f in mf:
        z = json.load(open(f))
        if z["extinct_at"] >= 0: continue
        h, _ = hist_of(z["final_counts"]); v100.append(lbar(h))
    v20 = [D[(0.01, s)]["lbar_end"] for s in SEEDS if (0.01, s) in D]
    if v20:
        print("\n== V5 서술 — μ 1%% 끝점 L̄: 시드 1~20 중앙 %.3f 대 main %d 접시 중앙 %.3f [%.3f, %.3f]" % (pct(v20, .5), len(v100), pct(v100, .5), pct(v100, .25), pct(v100, .75)))
        OUT["V5"] = {"n20_median": pct(v20, .5), "n100": len(v100), "n100_median": pct(v100, .5), "n100_q1": pct(v100, .25), "n100_q3": pct(v100, .75)}

print("\n== 종합: G0 %s · V1 %s · V2 %s" % ("✅" if not (g0["mismatch"] or g0["self_bad"]) else "❌", "✅" if v1 else "❌", "✅" if v2 else "❌"))
OUT["bad"] = bad
fn = "_RESULT_pilot17.json" if PILOT else "_RESULT_evo17.json"
json.dump(OUT, open(fn, "w"), indent=1, ensure_ascii=False); print("기록 →", fn)
