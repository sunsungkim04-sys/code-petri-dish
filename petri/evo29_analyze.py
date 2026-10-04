#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""해부 29 판정 — 진화 궤적의 고리 길어짐이 한 번 굴림 규칙에서 사라지나 (M3). 사전등록 해부29-사전등록-2026-10-03.md §4.

규칙: DF = draw-first(원래 · evo17/ 기록) · DO = draw-once(--remember-die 1 · evo29/do/) · FF = find-first(--find-first 1 · evo29/ff/)
양: 접시마다 끝점 L̄ = 고리 있는 개체의 개체 가중 평균 고리 길이(final_counts · evo17 과 같은 loop_ticks).
칸 평균 m(규칙, μ) = 그 칸의 멸종 안 한 접시 L̄ 의 산술 평균. 오름 rise(규칙) = m(규칙, 1%) − m(규칙, 0.1%).
재추출: 시드 1~20 을 복원 추출(2,000회) — 한 재추출 안에서 모든 규칙 · μ 가 **같은 시드 목록**을 쓴다(시드 짝). 칸 평균은 그 목록에서 멸종 안 한 접시로.
  난수 random.Random("haebu29-drawonce-evo") · 구간은 2.5 · 97.5 백분위.

판정선(결과 전 고정)
  G0   회귀: G0a 깃발 없음(sim 0.3.6) ≡ evo17 (μ 1% 시드 1 · 2 · μ 0.1% 시드 1 — 기록 전부) · G0b --find-first ≡ rewind8ff (시드 101 · 102)
       · G0c --remember-die 두 번 같다 · G0d dms.js 식 켜기 체크섬 = G0c · G0e 판정기가 evo17 칸 중앙값을 다시 계산해 _RESULT_evo17.json 과 1e-12 안
       · 점검: 설정(규칙 표시 · μ · 5만 틱 · sim 0.3.6) · 재료 보존 위반 0 · 자기 검산(loops 끝 = final_counts) — 하나라도 어긋나면 판정하지 않는다
  D1   (주) 상호작용 I = rise(DF) − rise(DO) · 95% 하한 > 0
  D2   (주) 지운 몫 R = I / rise(DF) · 종합 ✅ = D1 ∧ R 95% 하한 ≥ 0.5 · 🟡 = D1 만 · ❌ = D1 실패
       등급(서술): R 하한 ≥ 0.75 '대부분 사라진다' · 0.5~0.75 '절반 넘게' · 그 밖 '일부'
  D3   (부) μ 1% 끝점 차 m(DF,1%) − m(DO,1%) · 95% 하한 > 0
  D4   (부 · 실현 오류 맞춤) m(DF,1%) − m(DO,5%) · 95% 하한 > 0 — DO 명목 5% 는 DF 명목 1% 와 실현 오류가 비슷하다(§3)
       서술: rise_hi(DO) = m(DO,5%) − m(DO,1%) · m(DO,3%)
  D5   (부 · 눈 가리지 못함) I_FF = rise(DF) − rise(FF) · 95% 하한 > 0 · R_FF
  서술: 칸마다 중앙값 [사분위] · 최빈 등급 · 멸종 · Spearman(μ, L̄) 규칙별(4 μ · 접시 재추출) · e 포함 몫 · e 없는 개체만의 L̄ 과 그 오름 ·
        고리 없는 몫 · 교대 t4 · L̄(t) 중앙값 곡선
실행: python3 evo29_analyze.py [petri 폴더]  (환경: PETRI_ROOT PETRI_SEEDS PETRI_MUS_DO PETRI_MUS_FF PETRI_DF_ROOT PETRI_PILOT)
"""
import glob
import json
import math
import os
import random
import sys

BASE = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1 else "~/petri")
ROOT = os.environ.get("PETRI_ROOT", "evo29")
REG = ROOT + "reg"
DF_ROOT = os.environ.get("PETRI_DF_ROOT", "evo17")
SEEDS = [int(x) for x in os.environ["PETRI_SEEDS"].split()] if os.environ.get("PETRI_SEEDS") else list(range(1, 21))
MUS_DF = [0.001, 0.003, 0.005, 0.01]
MUS_DO = [float(x) for x in os.environ["PETRI_MUS_DO"].split()] if os.environ.get("PETRI_MUS_DO") else [0.001, 0.003, 0.005, 0.01, 0.03, 0.05]
MUS_FF = [float(x) for x in os.environ["PETRI_MUS_FF"].split()] if os.environ.get("PETRI_MUS_FF") else [0.001, 0.003, 0.005, 0.01]
PILOT = os.environ.get("PETRI_PILOT") == "1"
if PILOT and os.environ.get("PETRI_MUS_DF"): MUS_DF = [float(x) for x in os.environ["PETRI_MUS_DF"].split()]
LO, HI = 0.001, 0.01          # 오름의 양 끝
HIGH = 0.05                   # D4 실현 오류 맞춤 칸
R_PASS = 0.5                  # D2 종합 문턱(R 하한)
N_BOOT = 2000
RNG = random.Random("haebu29-drawonce-evo")
bad, notes = [], []
OUT = {"root": ROOT, "pilot": PILOT, "seeds": SEEDS, "mus": {"DF": MUS_DF, "DO": MUS_DO, "FF": MUS_FF}}


def loop_ticks(code):
    if "l" not in code or "c" not in code: return None
    li = code.index("l"); si = code.rfind("s", 0, li); start = si + 1 if si >= 0 else 0
    if start > li: return None
    seg = code[start:li + 1]
    return sum(1 for ch in seg if ch != "x") if "c" in seg else None


def mtag(mu): return "mu" + str(mu).replace(".", "p")


def hist_of(counts, skip_e=False):
    h, none = {}, 0
    for code, n in counts:
        if skip_e and "e" in code: continue
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
    xs = sorted(x for x in xs if x == x)
    if not xs: return float("nan")
    k = (len(xs) - 1) * p; f = math.floor(k); c = min(f + 1, len(xs) - 1)
    return xs[f] + (xs[c] - xs[f]) * (k - f)


def mean(xs):
    xs = [x for x in xs if x == x]
    return sum(xs) / len(xs) if xs else float("nan")


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


def ck(name, ok):
    print(("  ✅ " if ok else "  🚨 ") + name)
    if not ok: bad.append(name)


RECORD_KEYS = ("checksum", "samples", "tops", "specials", "snap", "final_counts", "births", "deaths", "ticks_run", "extinct_at")

# ---------------- G0 회귀 ----------------
print("== 해부 29 %s— 진화 궤적 · 한 번 굴림 규칙 · 루트 %s · 시드 %d · DF 기준 %s" % ("파일럿 " if PILOT else "판정 ", ROOT, len(SEEDS), DF_ROOT))
print("\n== G0 회귀")
for mu, s in ((0.01, 1), (0.01, 2), (0.001, 1)):
    a = os.path.join(BASE, REG, "df", "mat_d8_%s_s%05d.json" % (mtag(mu), s)); b = os.path.join(BASE, "evo17", "mat_d8_%s_s%05d.json" % (mtag(mu), s))
    if not (os.path.exists(a) and os.path.exists(b)): ck("G0a 파일 없음 μ %g 시드 %d" % (mu, s), False); continue
    za, zb = json.load(open(a)), json.load(open(b))
    same = all(za[k] == zb[k] for k in RECORD_KEYS + ("loops",)) and za["opts"] == zb["opts"]
    ck("G0a 깃발 없음 sim %s ≡ evo17 sim %s · μ %g 시드 %d · 체크섬 %s/%s · 기록 전부 %s" % (za["sim_version"], zb["sim_version"], mu, s, za["checksum"], zb["checksum"], same),
       same and za["sim_version"] == "0.3.6" and "remember_die" not in za and "find_first" not in za)
for s in (101, 102):
    a = os.path.join(BASE, REG, "ff", "mat_d8_mu0p01_s%05d_ff.json" % s); b = os.path.join(BASE, "rewind8ff", "mat_d8_mu0p01_s%05d_ff.json" % s)
    if not (os.path.exists(a) and os.path.exists(b)): ck("G0b 파일 없음 시드 %d" % s, False); continue
    za, zb = json.load(open(a)), json.load(open(b))
    same = all(za[k] == zb[k] for k in RECORD_KEYS)
    ck("G0b --find-first sim %s ≡ rewind8ff sim %s · 시드 %d · 체크섬 %s/%s · 기록 전부 %s" % (za["sim_version"], zb["sim_version"], s, za["checksum"], zb["checksum"], same),
       same and za.get("find_first") is True)
a = os.path.join(BASE, REG, "do_a", "mat_d8_mu0p01_s00901_do.json"); b = os.path.join(BASE, REG, "do_b", "mat_d8_mu0p01_s00901_do.json")
ca = None
if os.path.exists(a) and os.path.exists(b):
    za, zb = json.load(open(a)), json.load(open(b)); ca = za["checksum"]
    ck("G0c --remember-die 두 번 같다 · 체크섬 %s/%s · remember_die %s · opts.rememberDie %s" % (za["checksum"], zb["checksum"], za.get("remember_die"), za["opts"].get("rememberDie")),
       all(za[k] == zb[k] for k in RECORD_KEYS) and za.get("remember_die") is True and za["opts"].get("rememberDie") is True and za["ticks_run"] == 5000)
else: ck("G0c 파일 없음", False)
f = os.path.join(BASE, REG, "dms_style.txt")
if os.path.exists(f):
    parts = open(f).read().split()
    ck("G0d dms.js 식 켜기 (%s) 체크섬 = G0c (%s)" % (" ".join(parts), ca), len(parts) == 3 and parts[0] == "0.3.6" and parts[1] == "5000" and parts[2] == ca)
else: ck("G0d 파일 없음", False)
# G0e 계측기 자기 검산: evo17 칸 중앙값을 이 판정기의 정의로 다시 계산
ev = os.path.join(BASE, "_RESULT_evo17.json")
if os.path.exists(ev) and not PILOT:
    E = json.load(open(ev)); dev = 0.0
    for mu in MUS_DF:
        v = []
        for s in range(1, 21):
            z = json.load(open(os.path.join(BASE, "evo17", "mat_d8_%s_s%05d.json" % (mtag(mu), s))))
            if z["extinct_at"] >= 0: continue
            v.append(lbar(hist_of(z["final_counts"])[0]))
        dev = max(dev, abs(pct(v, .5) - E["V1|%s" % mu]["median"]))
    ck("G0e evo17 칸 중앙값 재계산 최대 차 %.2e (≤ 1e-12)" % dev, dev <= 1e-12)
elif PILOT: print("  (파일럿: G0e 는 판정 시드 기록으로만 한다 — 건너뜀)")

# ---------------- 읽기 ----------------
DIRS = {"DF": (DF_ROOT, ""), "DO": (os.path.join(ROOT, "do"), "_do"), "FF": (os.path.join(ROOT, "ff"), "_ff")}
MUS = {"DF": MUS_DF, "DO": MUS_DO, "FF": MUS_FF}
D, EXT = {}, []
selfchk = {"n": 0, "bad": []}
for rule in ("DF", "DO", "FF"):
    d, suf = DIRS[rule]
    for mu in MUS[rule]:
        for s in SEEDS:
            f = os.path.join(BASE, d, "mat_d8_%s_s%05d%s.json" % (mtag(mu), s, suf))
            if not os.path.exists(f): bad.append("없음 %s" % f); continue
            z = json.load(open(f))
            flag_ok = (z.get("remember_die") is True) == (rule == "DO") and (z.get("find_first") is True) == (rule == "FF")
            ver_ok = z["sim_version"] == ("0.3.5" if (rule == "DF" and DF_ROOT in ("evo17", "pilot17b")) else "0.3.6")
            if not flag_ok or not ver_ok or abs(z["opts"]["mu"] - mu) > 1e-12 or z["ticks_planned"] != 50000 or not z.get("loops_on") or z["conservation_violations"]:
                bad.append("설정 %s μ %g 시드 %d" % (rule, mu, s)); continue
            if z["extinct_at"] >= 0: EXT.append((rule, mu, s, z["extinct_at"])); continue
            if z["ticks_run"] != 50000: bad.append("틱 %s μ %g 시드 %d" % (rule, mu, s)); continue
            hf, nf = hist_of(z["final_counts"]); last = z["loops"][-1]
            selfchk["n"] += 1
            if last[0] != z["ticks_run"] or {int(k): v for k, v in last[1].items()} != hf or last[2] != nf: selfchk["bad"].append("%s %g %d" % (rule, mu, s))
            he, _ = hist_of(z["final_counts"], skip_e=True)
            traj = [(t, {int(k): v for k, v in h.items()}) for t, h, _n in z["loops"]]
            t4 = None
            for t, h in traj:
                sh = shares(h)
                if sh.get(4, 0) + sh.get(5, 0) > sh.get(2, 0): t4 = t; break
            esh = shares(hf); mode = max(esh, key=esh.get) if esh else None
            last_s = z["samples"][-1]
            D[(rule, mu, s)] = {"L": lbar(hf), "L_noe": lbar(he), "mode": mode, "t4": t4, "e": last_s[6] / last_s[1] if last_s[1] else float("nan"),
                                "none": nf / (sum(hf.values()) + nf) if hf else float("nan"), "traj": [(t, lbar(h)) for t, h in traj]}
ck("자기 검산(loops 끝 = final_counts) %d · 불일치 %s" % (selfchk["n"], selfchk["bad"]), not selfchk["bad"])
print("  멸종(제외) %d: %s" % (len(EXT), ["%s μ %g 시드 %d @%d" % e for e in EXT]))
print("  점검 문제 %d" % len(bad))
for x in bad[:12]: print("    ", x)
OUT["extinct"] = EXT; OUT["bad"] = bad
if bad and not PILOT: print("🚨 점검 실패 — 판정하지 않는다"); json.dump(OUT, open("_RESULT_evo29_halted.json", "w"), indent=1, ensure_ascii=False); sys.exit(1)


# ---------------- 재추출 틀 ----------------
def cell_mean(rule, mu, seeds, key="L"):
    return mean([D[(rule, mu, s)][key] for s in seeds if (rule, mu, s) in D])


def stats(seeds, key="L"):
    m = lambda r, mu: cell_mean(r, mu, seeds, key)
    o = {}
    o["rise_DF"] = m("DF", HI) - m("DF", LO)
    if HI in MUS_DO and LO in MUS_DO:
        o["rise_DO"] = m("DO", HI) - m("DO", LO); o["I"] = o["rise_DF"] - o["rise_DO"]; o["R"] = o["I"] / o["rise_DF"] if o["rise_DF"] else float("nan")
        o["D3"] = m("DF", HI) - m("DO", HI)
    if HIGH in MUS_DO and HI in MUS_DO:
        o["D4"] = m("DF", HI) - m("DO", HIGH); o["rise_hi_DO"] = m("DO", HIGH) - m("DO", HI)
    if HI in MUS_FF and LO in MUS_FF:
        o["rise_FF"] = m("FF", HI) - m("FF", LO); o["I_FF"] = o["rise_DF"] - o["rise_FF"]; o["R_FF"] = o["I_FF"] / o["rise_DF"] if o["rise_DF"] else float("nan")
    return o


def boot(key="L"):
    pt = stats(SEEDS, key); bs = {k: [] for k in pt}
    rng = random.Random(RNG.random())
    for _ in range(N_BOOT):
        smp = [SEEDS[rng.randrange(len(SEEDS))] for _ in SEEDS]
        o = stats(smp, key)
        for k in bs: bs[k].append(o.get(k, float("nan")))
    return {k: (pt[k], pct(bs[k], .025), pct(bs[k], .975)) for k in pt}


# ---------------- 칸 표 (서술) ----------------
print("\n== 칸 — 끝점 L̄ (평균 · 중앙값 [사분위]) · 최빈 등급 · e 포함 몫 · e 없는 개체 L̄ · 고리 없는 몫 · 교대 t4")
for rule in ("DF", "DO", "FF"):
    for mu in MUS[rule]:
        ks = [s for s in SEEDS if (rule, mu, s) in D]
        if not ks: print("   %s μ %.1f%%: 접시 0" % (rule, mu * 100)); continue
        v = [D[(rule, mu, s)]["L"] for s in ks]
        cnt = {}
        for s in ks: cnt[D[(rule, mu, s)]["mode"]] = cnt.get(D[(rule, mu, s)]["mode"], 0) + 1
        t4 = [D[(rule, mu, s)]["t4"] for s in ks]; t4v = [t for t in t4 if t is not None]
        row = {"n": len(ks), "mean": mean(v), "median": pct(v, .5), "q1": pct(v, .25), "q3": pct(v, .75), "modes": {str(k): n for k, n in cnt.items()},
               "e": pct([D[(rule, mu, s)]["e"] for s in ks], .5), "L_noe": pct([D[(rule, mu, s)]["L_noe"] for s in ks], .5),
               "none": pct([D[(rule, mu, s)]["none"] for s in ks], .5), "t4_none": len(t4) - len(t4v), "t4_median": pct(t4v, .5) if t4v else None, "values": v}
        OUT["cell|%s|%s" % (rule, mu)] = row
        print("   %s μ %.1f%%: n %2d · 평균 %.3f · 중앙 %.3f [%.3f, %.3f] · 최빈 %s · e %.3f · e없는 L̄ %.3f · 고리없음 %.3f · t4 %s" % (
            rule, mu * 100, row["n"], row["mean"], row["median"], row["q1"], row["q3"], dict(sorted(cnt.items(), key=lambda kv: -kv[1])),
            row["e"], row["L_noe"], row["none"], ("중앙 %d · 없음 %d" % (row["t4_median"], row["t4_none"])) if t4v else "없음 %d/%d" % (len(t4), len(t4))))
        T = [t for t, _ in D[(rule, mu, ks[0])]["traj"]]
        OUT["curve|%s|%s" % (rule, mu)] = {"t": T, "median": [pct([dict(D[(rule, mu, s)]["traj"]).get(t, float("nan")) for s in ks], .5) for t in T]}

# ---------------- 판정 ----------------
B = boot("L"); OUT["boot"] = {k: list(v) for k, v in B.items()}
f3 = lambda k: "%+.3f [%+.3f, %+.3f]" % B[k]
print("\n== 판정 (시드 재추출 %d · 칸 평균)" % N_BOOT)
print("   rise(DF) = m(DF,1%%) − m(DF,0.1%%) = %s" % f3("rise_DF"))
verdict = None
if "I" in B:
    print("   rise(DO) = %s" % f3("rise_DO"))
    d1 = B["I"][1] > 0
    print("   D1 (주) I = rise(DF) − rise(DO) = %s → %s" % (f3("I"), "✅" if d1 else "❌"))
    r_lo = B["R"][1]
    grade = "대부분 사라진다" if r_lo >= 0.75 else ("절반 넘게 사라진다" if r_lo >= 0.5 else "일부만 사라진다")
    print("   D2 (주) R = I / rise(DF) = %.3f [%.3f, %.3f] → 등급 '%s'" % (B["R"][0], B["R"][1], B["R"][2], grade))
    verdict = "✅" if d1 and r_lo >= R_PASS else ("🟡" if d1 else "❌")
    print("   종합: %s  (✅ = D1 ∧ R 하한 ≥ %.2f · 🟡 = D1 만 · ❌ = D1 실패)" % (verdict, R_PASS))
    print("   D3 (부) m(DF,1%%) − m(DO,1%%) = %s → %s" % (f3("D3"), "✅" if B["D3"][1] > 0 else "❌"))
    OUT["verdict"] = {"D1": d1, "R_grade": grade, "overall": verdict, "D3": B["D3"][1] > 0}
if "D4" in B:
    print("   D4 (부 · 실현 오류 맞춤) m(DF,1%%) − m(DO,5%%) = %s → %s · 서술 rise_hi(DO) = m(DO,5%%) − m(DO,1%%) = %s" % (f3("D4"), "✅" if B["D4"][1] > 0 else "❌", f3("rise_hi_DO")))
    OUT.setdefault("verdict", {})["D4"] = B["D4"][1] > 0
if "I_FF" in B:
    print("   D5 (부 · 눈 가리지 못함) rise(FF) = %s · I_FF = %s → %s · R_FF %.3f [%.3f, %.3f]" % (f3("rise_FF"), f3("I_FF"), "✅" if B["I_FF"][1] > 0 else "❌", *B["R_FF"]))
    OUT.setdefault("verdict", {})["D5"] = B["I_FF"][1] > 0

# ---------------- 서술: Spearman · e 없는 개체 ----------------
print("\n== 서술 — Spearman(μ, L̄) 규칙별 (명목 0.1~1% 네 층 · 접시 재추출)")
for rule in ("DF", "DO", "FF"):
    pairs = [(mu, D[(rule, mu, s)]["L"]) for mu in MUS[rule] if mu <= HI + 1e-12 for s in SEEDS if (rule, mu, s) in D]
    if len(pairs) < 4 or len(set(p[0] for p in pairs)) < 2: continue
    rho = spearman([p[0] for p in pairs], [p[1] for p in pairs]); bs = []
    rng = random.Random("haebu29-spearman-" + rule)
    for _ in range(N_BOOT):
        smp = [pairs[rng.randrange(len(pairs))] for _ in pairs]; bs.append(spearman([p[0] for p in smp], [p[1] for p in smp]))
    print("   %s: %.3f [%.3f, %.3f] · 접시 %d" % (rule, rho, pct(bs, .025), pct(bs, .975), len(pairs)))
    OUT["spearman|%s" % rule] = [rho, pct(bs, .025), pct(bs, .975), len(pairs)]
Be = boot("L_noe"); OUT["boot_noe"] = {k: list(v) for k, v in Be.items()}
print("\n== 서술 — e 없는 개체만의 L̄ 로 같은 양 (e 혼입 점검)")
for k in ("rise_DF", "rise_DO", "I", "R", "rise_FF", "I_FF"):
    if k in Be: print("   %-8s %.3f [%.3f, %.3f]" % ((k,) + Be[k]))

print("\n== 종합: G0 %s · 종합 %s" % ("✅" if not [b for b in bad if b.startswith("G0")] else "❌", verdict))
fn = "_RESULT_pilot29.json" if PILOT else "_RESULT_evo29.json"
json.dump(OUT, open(os.path.join(BASE, fn), "w"), indent=1, ensure_ascii=False); print("기록 →", fn)
