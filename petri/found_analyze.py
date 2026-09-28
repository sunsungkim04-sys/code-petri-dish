"""해부 2B 판정 — 씨앗 시절 창업 경쟁. 사전등록 해부2-사전등록-2026-09-15.md §3 · §4. 결과 전에 썼다.

  팔의 퍼짐 s = ln((m1_T + 0.5) / (m0_T + 0.5)) − ln(m1_0 / m0_0)    m1 = 짝 코드 계통(표지 1) · m0 = 씨앗 계통(표지 0)
  배경마다 · μ 마다 Δ = (시험 팔 5개 s 평균) − (중립 팔 5개 s 평균)  →  배경 평균 · 배경 재추출 부트스트랩 2,000회 95%
  판정: 짝이 낫다 = 하한 > 0 · 씨앗이 낫다 = 상한 < 0 · 구별 안 됨 = 0 을 품음
  μ 효과: 배경마다 Δ(μ 1%) − Δ(μ 0) 의 평균 · 95%
  점검 ① 복제 충실도 · ③ 재료 보존 — 실패한 배경은 빼고 보고 · 멸종해 팔이 없는 배경도 따로 센다
       ② 본실험 기록 대조(T0 개체 수 · 상위 20) — 하나라도 불일치면 전부 중단
       ④ 중립 팔 s 의 배경 평균 구간이 0 을 품어야 한다(μ 마다) — 아니면 반반 나누기가 치우친 것 → 그 조건 판정 보류
       ⑤ 양성 대조: space(짝 reascld) · energy(짝 rahcld) 가 μ 0 에서 '짝이 낫다' 여야 한다 — 아니면 mat 해석 보류
실행: python3 found_analyze.py [petri 폴더, 기본 ~/petri]
"""
import glob
import json
import math
import os
import random
import sys

BASE = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1 else "~/petri")
PAIRS = [("space", "reascld", 500, "양성 대조"), ("energy", "rahcld", 500, "양성 대조"), ("mat", "racld", 2000, "시험")]
MUS = [0.0, 0.01]
REPS = 5
N_BOOT = 2000
RNG = random.Random(20260918)


def pct(xs, p):
    xs = sorted(xs)
    k = (len(xs) - 1) * p
    f = math.floor(k)
    c = min(f + 1, len(xs) - 1)
    return xs[f] + (xs[c] - xs[f]) * (k - f)


def s_of(arm):
    r = arm["series"]
    _, _, m1_0, m0_0 = r[0]
    _, _, m1_T, m0_T = r[-1]
    return math.log((m1_T + 0.5) / (m0_T + 0.5)) - math.log(m1_0 / m0_0)


def boot_ci(vals, idxs):
    n = len(vals)
    bm = [sum(vals[i] for i in idx) / n for idx in idxs]
    return sum(vals) / n, pct(bm, .025), pct(bm, .975)


def verdict(lo, hi):
    return "짝이 낫다" if lo > 0 else ("씨앗이 낫다" if hi < 0 else "구별 안 됨")


out = {"pairs": {}}
halt = False
pos_ok = True
for cond, b, t0, role in PAIRS:
    files = sorted(glob.glob(os.path.join(BASE, "found", f"found_{cond}_{b}_t{t0}_s*.json")))
    print(f"\n==================== {cond} · 씨앗 rascld 대 짝 {b} · T0 {t0:,} · {role} · 파일 {len(files)}/100")
    if len(files) != 100:
        print("  🚨 계획한 100 배경이 다 모이지 않았다 — 판정하지 않는다")
        halt = True
        continue
    bgs, dropped, extinct = [], [], 0
    for f in files:
        z = json.load(open(f))
        # 결과 전 교정(09-15, 발사 직후 로그에서 발견): T0 전에 멸종한 접시는 본실험 기록이 멸종 틱에서 끝나 T0 기록이 없다
        # → found.js 는 불일치로 적는다. 이 경우는 두 기록의 멸종 틱이 같은지로 대조한다(같으면 궤적 일치 · 팔 없음으로 센다).
        if z["grow_extinct"] >= 0:
            main_ext = json.load(open(os.path.join(BASE, "main", f"{cond}_d8_mu0p01_s{z['seed']:05d}.json")))["extinct_at"]
            if main_ext != z["grow_extinct"]:
                print(f"  🚨 점검 ② 멸종 틱 불일치 — 시드 {z['seed']} 여기 {z['grow_extinct']} · 본실험 {main_ext} → 전부 중단")
                halt = True
            extinct += 1
            continue
        if not z["early"]["checked"] or not z["early"]["match"]:
            print(f"  🚨 점검 ② 본실험 기록과 불일치 — 시드 {z['seed']} {z['early']} → 전부 중단")
            halt = True
        if not z["arms"]:
            extinct += 1
            continue
        why = []
        if not z["fidelity"]["same"]:
            why.append("① 복제 충실도")
        if len(z["arms"]) != len(MUS) * REPS * 2:
            why.append(f"팔 수 {len(z['arms'])}")
        if any(not a["cons_ok"] for a in z["arms"]):
            why.append("③ 재료 보존")
        if any(a["conv"]["to_a"] - a["conv"]["skipped"] <= 0 or a["series"][0][2] == 0 or a["series"][0][3] == 0 for a in z["arms"]):
            why.append("바꾼 계통이 비었다")
        if why:
            dropped.append((z["seed"], why))
        else:
            bgs.append(z)
    if halt:
        continue
    nb = len(bgs)
    print(f"  쓴 배경 {nb} · T0 전 멸종/개체 부족 {extinct} · 뺀 배경 {len(dropped)} " + " · ".join(f"시드 {s}: {'/'.join(w)}" for s, w in dropped))
    pops = [z["grow_pop"] for z in bgs]
    skips = sum(a["conv"]["skipped"] for z in bgs for a in z["arms"])
    print(f"  T0 개체 수 중앙 {pct(pops, .5):,.0f} [5% {pct(pops, .05):,.0f}] · 재료 부족으로 못 바꾼 개체 누계 {skips}")
    idxs = [[RNG.randrange(nb) for _ in range(nb)] for _ in range(N_BOOT)]
    res = {"n_bg": nb, "extinct_before": extinct, "dropped": dropped, "mu": {}}
    deltas = {}
    for mu in MUS:
        neut = [sum(s_of(a) for a in z["arms"] if a["kind"] == "neutral" and a["mu"] == mu) / REPS for z in bgs]
        test = [sum(s_of(a) for a in z["arms"] if a["kind"] == "test" and a["mu"] == mu) / REPS for z in bgs]
        d = [t - n_ for t, n_ in zip(test, neut)]
        deltas[mu] = d
        nm, nlo, nhi = boot_ci(neut, idxs)
        dm, dlo, dhi = boot_ci(d, idxs)
        check4 = nlo <= 0 <= nhi
        fixed = sum(1 for z in bgs for a in z["arms"] if a["kind"] == "test" and a["mu"] == mu and a["series"][-1][3] == 0)
        lost = sum(1 for z in bgs for a in z["arms"] if a["kind"] == "test" and a["mu"] == mu and a["series"][-1][2] == 0)
        share = [a["series"][-1][2] / max(1, a["series"][-1][2] + a["series"][-1][3]) for z in bgs for a in z["arms"] if a["kind"] == "test" and a["mu"] == mu]
        v = verdict(dlo, dhi) if check4 else "판정 보류 — 점검 ④ 발동(중립 팔이 한쪽으로 치우침)"
        res["mu"][str(mu)] = {"neutral": [nm, nlo, nhi], "delta": [dm, dlo, dhi], "check4": check4, "fixed": fixed, "lost": lost,
                              "share_median": pct(share, .5), "verdict": v}
        print(f"  μ {mu:<4}: 중립 s {nm:+.3f} [{nlo:+.3f}, {nhi:+.3f}] {'✓' if check4 else '🚨'} · Δ {dm:+.3f} [{dlo:+.3f}, {dhi:+.3f}] → {v}")
        print(f"          시험 팔 끝 짝 계통 몫 중앙 {pct(share, .5):.0%} · 짝이 씨앗 계통을 없앤 팔 {fixed}/{nb * REPS} · 짝 계통이 사라진 팔 {lost}/{nb * REPS}")
    dd = [x1 - x0 for x0, x1 in zip(deltas[0.0], deltas[0.01])]
    mm, mlo, mhi = boot_ci(dd, idxs)
    res["mu_effect"] = [mm, mlo, mhi]
    print(f"  μ 효과 Δ(1%) − Δ(0): {mm:+.3f} [{mlo:+.3f}, {mhi:+.3f}] → " +
          ("복제 오류가 짝 쪽으로 민다" if mlo > 0 else ("복제 오류가 씨앗 쪽으로 민다" if mhi < 0 else "복제 오류 효과 검출 안 됨")))
    if role == "양성 대조" and not res["mu"]["0.0"]["verdict"] == "짝이 낫다":
        pos_ok = False
    out["pairs"][cond] = res

if not halt:
    print(f"\n점검 ⑤ 양성 대조(space · energy 가 μ 0 에서 짝이 낫다): {'통과' if pos_ok else '🚨 발동 — mat 해석 보류'}")
    out["positive_ok"] = pos_ok
    json.dump(out, open(os.path.join(BASE, "_RESULT_found.json"), "w"), ensure_ascii=False, indent=1)
else:
    print("\n🚨 판정 중단 사유가 있다 — 위를 볼 것")
    sys.exit(1)
