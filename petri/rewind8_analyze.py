# -*- coding: utf-8 -*-
"""해부 8 A 판정 — 재료 먼저 세계(sim v0.3.3 findFirst)에서 처음부터 되감으면 mat 의 도착점이 바뀌나.
사전등록 해부8-사전등록-2026-09-17.md §2 · §5. 결과 전에 썼다(파일럿 시드 9101~9110 은 봤다 — 노트 §5).

  재료 먼저: rewind8ff/mat_d8_mu0p01_s{101..200}_ff.json   (mat 밀도 8 · μ 1% · 50,000틱)
  원래 규칙: rewind2/mat_d8_mu0p01_s{101..200}.json         (되감기 2 기록 · sim 0.3.1 — v0.3.3 기본값과 궤적 같음 3/3 확인)
  같은 시드 = 같은 첫 재료 배치(첫 복사 전까지 난수가 같다). 그 뒤 궤적은 갈라지므로 두 묶음은 **짝이 아니라 독립 표본**으로 다룬다.
  정의는 되감기 2(rewind2_analyze.py)와 같다:
    끝 우세 = 마지막 틱의 가장 흔한 코드 · 안정성 = 40,000틱 이후 상위 기록(500틱마다) 중 1위가 끝 우세인 몫
    길 = 1위 코드가 기록 두 번 연속 1위가 된 순서(씨앗 제외)

  A1 도착점   : p = 살아남은 접시 중 끝 우세가 racld 인 몫. 차(재료 먼저 − 원래) — 두 묶음 각각 접시 재추출 2,000회 95%
                상한 < 0 → 'racld 도착점이 줄어든다' · 하한 > 0 → '늘어난다' · 그 밖 '차이 검출 안 됨'
                덧붙임: p_재료먼저 Wilson 상한 < 0.5 → 'racld 는 더 이상 흔한 도착점이 아니다'
  A2 무엇이   : 재료 먼저 묶음에 되감기 2 Q2' 규칙 그대로 — 안정성 중앙 < 0.90 이면 '떠돈다',
                아니면 가장 흔한 끝 우세 공유 몫 Wilson 하한 ≥ 0.90 '한 곳' · 상한 ≤ 0.50 '갈라짐' · 그 밖 '대체로 한 곳'
  A3 씨앗 자리 : 길이 비어 있는(씨앗이 한 번도 1위를 내주지 않은) 접시 몫 — 두 묶음 차 (A1 과 같은 재추출) · 서술 + 구간
  A4 서술     : 멸종 · 끝 racld 몫(개체 중) · 끝 씨앗 몫 · e(먹기) 포함 몫 · 평균 길이 · 코드 종류 · 처음 racld 가 상위 3 에 든 틱
  완결성: 200 파일 · 설정 · 재료 보존 위반 0 — 아니면 판정하지 않는다
실행: python3 rewind8_analyze.py [petri 폴더]  →  _RESULT_rewind8.json
마른 실행(판정 시드를 보지 않게 파일럿 시드로): PETRI_ROOT=pilotA PETRI_SEEDS=9101-9103 PETRI_ORIG=pilotA_orig PETRI_ORIG_V=0.3.3 python3 rewind8_analyze.py
"""
import json
import math
import os
import random
import sys

BASE = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1 else "~/petri")
ROOT = os.environ.get("PETRI_ROOT", "rewind8ff")
ORIG = os.environ.get("PETRI_ORIG", "rewind2")
ORIG_V = set(os.environ.get("PETRI_ORIG_V", "0.3.1").split(","))
_s = os.environ.get("PETRI_SEEDS", "101-200").split("-")
SEEDS = list(range(int(_s[0]), int(_s[1]) + 1))
TICKS, MU, DENS = 50000, 0.01, 8
LAG_FROM = 40000
STABLE, ONE_PLACE, SPLIT = 0.90, 0.90, 0.50
N_BOOT = 2000
RNG = random.Random(20260917)
SEED_CODE, TARGET = "rascld", "racld"


def pct(xs, p):
    xs = sorted(xs)
    k = (len(xs) - 1) * p
    f = math.floor(k)
    c = min(f + 1, len(xs) - 1)
    return xs[f] + (xs[c] - xs[f]) * (k - f)


def wilson(k, n, z=1.96):
    if n == 0:
        return (float("nan"), float("nan"))
    p = k / n
    c = (p + z * z / (2 * n)) / (1 + z * z / n)
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n)
    return (max(0.0, c - h), min(1.0, c + h))


def dominant(z):
    return z["final_counts"][0][0]


def stability(z):
    fin = dominant(z)
    recs = [top[0][0] for t, top in z["tops"] if t >= LAG_FROM and top]
    return sum(k == fin for k in recs) / len(recs) if recs else float("nan")


def path(z):
    order, seen, prev, run = [], set(), None, 0
    for t, top in z["tops"]:
        k = top[0][0] if top else None
        if k is None:
            prev, run = None, 0
            continue
        run = run + 1 if k == prev else 1
        if run == 2 and k != SEED_CODE and k not in seen:
            seen.add(k)
            order.append(k)
        prev = k
    return tuple(order)


def load(folder, suffix, want_ff, versions):
    ds, probs = [], []
    for s in SEEDS:
        f = os.path.join(BASE, folder, "mat_d8_mu0p01_s%05d%s.json" % (s, suffix))
        if not os.path.exists(f):
            probs.append("없음 " + os.path.basename(f))
            continue
        z = json.load(open(f))
        o = z["opts"]
        if (z["sim_version"] not in versions or z["ticks_planned"] != TICKS or o["mu"] != MU or o["density"] != DENS
                or z["cond"] != "mat" or o["seed"] != s or bool(o.get("findFirst", False)) != want_ff
                or bool(z.get("find_first", False)) != want_ff):
            probs.append("설정 다름 " + os.path.basename(f))
        if z["conservation_violations"]:
            probs.append("재료 보존 위반 " + os.path.basename(f))
        ds.append(z)
    return ds, probs


ff, p1 = load(ROOT, "_ff", True, {"0.3.3", "0.3.4", "0.3.5"})
orig, p2 = load(ORIG, "", False, ORIG_V)
print("== 해부 8 A 재료 먼저 세계에서 처음부터 되감기 — 판정")
print("  재료 먼저 %s %d/%d · 원래 규칙 %s %d/%d" % (ROOT, len(ff), len(SEEDS), ORIG, len(orig), len(SEEDS)))
if p1 or p2:
    print("🚨 완결성 가드 %d건 — 판정하지 않는다" % (len(p1) + len(p2)))
    for p in (p1 + p2)[:20]:
        print("   ", p)
    sys.exit(1)

sets = {"orig": [z for z in orig if z["extinct_at"] < 0], "ff": [z for z in ff if z["extinct_at"] < 0]}
out = {"n": {k: len(v) for k, v in sets.items()}, "extinct": {"orig": len(orig) - len(sets["orig"]), "ff": len(ff) - len(sets["ff"])}}
print("  생존: 원래 %d/%d · 재료 먼저 %d/%d" % (len(sets["orig"]), len(orig), len(sets["ff"]), len(ff)))

ind = {k: [(dominant(z) == TARGET, len(path(z)) == 0) for z in v] for k, v in sets.items()}
idx = {k: [[RNG.randrange(len(v)) for _ in v] for _ in range(N_BOOT)] for k, v in ind.items()}


def share(k, j, ix=None):
    v = ind[k]
    rows = [v[i] for i in ix] if ix is not None else v
    return sum(r[j] for r in rows) / len(rows)


for j, name in ((0, "A1 끝 우세 racld 몫"), (1, "A3 씨앗이 1위를 한 번도 안 내준 몫")):
    po, pf = share("orig", j), share("ff", j)
    d = [share("ff", j, idx["ff"][b]) - share("orig", j, idx["orig"][b]) for b in range(N_BOOT)]
    lo, hi = pct(d, .025), pct(d, .975)
    kf = sum(r[j] for r in ind["ff"])
    wl, wh = wilson(kf, len(ind["ff"]))
    if j == 0:
        verdict = ("racld 도착점이 줄어든다" if hi < 0 else "racld 도착점이 늘어난다" if lo > 0 else "차이 검출 안 됨")
        extra = " · 재료 먼저 Wilson [%.2f, %.2f] → %s" % (wl, wh, "racld 는 더 이상 흔한 도착점이 아니다" if wh < 0.5 else "(상한 ≥ 0.5)")
    else:
        verdict = ("씨앗이 자리를 더 지킨다" if lo > 0 else "씨앗이 자리를 덜 지킨다" if hi < 0 else "차이 검출 안 됨") + " (서술)"
        extra = " · 재료 먼저 Wilson [%.2f, %.2f]" % (wl, wh)
    print("\n-- %s: 원래 %.3f (%d/%d) → 재료 먼저 %.3f (%d/%d) · 차 %+.3f [%+.3f, %+.3f] → **%s**%s"
          % (name, po, sum(r[j] for r in ind["orig"]), len(ind["orig"]), pf, kf, len(ind["ff"]), pf - po, lo, hi, verdict, extra))
    out["A1" if j == 0 else "A3"] = {"orig": po, "ff": pf, "diff": [pf - po, lo, hi], "ff_wilson": [wl, wh], "verdict": verdict}

# A2 — 되감기 2 Q2' 규칙
for k in ("ff", "orig"):
    v = sets[k]
    stab = [stability(z) for z in v]
    doms = {}
    for z in v:
        doms[dominant(z)] = doms.get(dominant(z), 0) + 1
    top = sorted(doms.items(), key=lambda x: (-x[1], x[0]))
    code, cnt = top[0]
    lo, hi = wilson(cnt, len(v))
    med = pct(stab, .5)
    if med < STABLE:
        verdict = "끝 우세 코드가 자리 잡지 않는다(떠돈다)"
    elif lo >= ONE_PLACE:
        verdict = "한 곳"
    elif hi <= SPLIT:
        verdict = "갈라진다"
    else:
        verdict = "대체로 한 곳"
    print("\n-- A2 %s: 안정성 중앙 %.3f (1 미만 %d/%d) · 가장 흔한 끝 우세 %s %d/%d = %.0f%% [%.0f%%, %.0f%%] → %s%s"
          % ("재료 먼저" if k == "ff" else "원래(참고)", med, sum(x < 1 for x in stab), len(v), code, cnt, len(v),
             100 * cnt / len(v), 100 * lo, 100 * hi, "**%s**" % verdict if k == "ff" else verdict, ""))
    print("        끝 우세 상위: " + " · ".join("%s %d" % t for t in top[:6]))
    out.setdefault("A2", {})[k] = {"stab_median": med, "top": top[:10], "wilson": [lo, hi], "verdict": verdict}

# A4 서술
print("\n-- A4 서술 (살아남은 접시 · 중앙 [사분위])")


def q(xs):
    return "%.3f [%.3f, %.3f]" % (pct(xs, .5), pct(xs, .25), pct(xs, .75))


for k in ("orig", "ff"):
    v = sets[k]
    fr = lambda z, pred: sum(c for key, c in z["final_counts"] if pred(key)) / max(1, sum(c for _, c in z["final_counts"]))
    rac = [fr(z, lambda s: s == TARGET) for z in v]
    sed = [fr(z, lambda s: s == SEED_CODE) for z in v]
    eat = [fr(z, lambda s: "e" in s) for z in v]
    ln = [z["samples"][-1][3] for z in v]
    kinds = [z["samples"][-1][2] for z in v]
    first = []
    for z in v:
        t0 = next((t for t, top in z["tops"] if any(key == TARGET for key, _ in top[:3])), None)
        if t0 is not None:
            first.append(t0)
    paths = {}
    for z in v:
        paths[path(z)] = paths.get(path(z), 0) + 1
    ptop = sorted(paths.items(), key=lambda x: -x[1])[:4]
    print("   %s: 멸종 %d · 끝 racld 몫 %s · 끝 씨앗 몫 %s · e 포함 몫 %s · 평균 길이 %s · 코드 종류 %s"
          % ("원래    " if k == "orig" else "재료 먼저", out["extinct"][k], q(rac), q(sed), q(eat), q(ln), q(kinds)))
    print("            racld 가 상위 3 에 처음 든 틱: %d/%d 접시 · 중앙 %s · 흔한 길 %s"
          % (len(first), len(v), ("%.0f" % pct(first, .5)) if first else "-",
             " · ".join("%s %d" % ("→".join(p) if p else "(씨앗 그대로)", c) for p, c in ptop)))
    out.setdefault("A4", {})[k] = {"racld_share": [pct(rac, .5), pct(rac, .25), pct(rac, .75)],
                                   "seed_share": [pct(sed, .5), pct(sed, .25), pct(sed, .75)],
                                   "eat_share": [pct(eat, .5), pct(eat, .25), pct(eat, .75)],
                                   "mean_len": pct(ln, .5), "kinds": pct(kinds, .5), "racld_top3": len(first),
                                   "paths": [["→".join(p), c] for p, c in ptop]}

json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "_RESULT_rewind8.json"), "w"), ensure_ascii=False, indent=1)
print("\n저장: _RESULT_rewind8.json")
