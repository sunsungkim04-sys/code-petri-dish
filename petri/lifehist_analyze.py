# -*- coding: utf-8 -*-
"""해부 8 B 판정 — 오류가 없을 때(μ 0) 씨앗의 침입 Δ ≈ +0.39 는 어떻게 쌓이나. 사전등록 해부8-사전등록-2026-09-17.md §3 · §5.
결과 전에 썼다(파일럿 = 되감기 2 배경 101~103 — 노트 §5).

  배경: 본실험 mat 접시 중 끝 우세가 racld 인 접시를 시드 순으로 늘어놓아 **21~40 번째** 20개(1~20 번째는 해부 1 · 6C · 7 이 썼다)
  팔: ref + neutral 10(WT = racld) · rascld(씨앗) · μ 0 · 6,000틱 · lifehist.js v1.0.0
  s(t) = ln((m_t + 0.5)/m_0) − ln((r_t + 0.5)/r_0)   (dms.js 와 같다 · m = 끼운 계통 · r = 나머지)
  Δ(t) = s_씨앗(t) − (WT 팔 11개 s(t) 평균) · 배경 재추출 2,000회 95%
  g1(t) = ln((m_t + 0.5)/m_0) · g0(t) = ln((r_t + 0.5)/r_0) → Δ = [g1 씨앗 − g1 WT] − [g0 씨앗 − g0 WT]

  B0 재현   : Δ(3000) 하한 > 0 → '새 배경에서도 씨앗이 이긴다' — 아니면 B1~B4 는 적되 '이득을 설명' 으로 읽지 않는다
  B1 복리   : 배경마다 Δ(t) 를 t = 500, 750, …, 6000 에 직선 적합 → 기울기 b. 꾸준한 몫 S = 3000·b̄ / Δ̄(3000)
              S 하한 ≥ 2/3 → '대부분 꾸준히 쌓인다(복리)' · S 상한 ≤ 1/3 → '대부분 끼운 직후에 생긴다' · 그 밖 '둘 다 보탠다'
              서술: Δ(6000)/Δ(3000) · 기울기 앞 절반[500, 3000] 대 뒤 절반[3000, 6000] 비
  B2 공동체 : 공동체 몫 Q0 = −(g0 씨앗 − g0 WT)(3000) / Δ(3000)
              하한 > 0 → '나머지 공동체 쪽도 보탠다' · 상한 < 0 → '공동체 쪽은 거꾸로 작용한다' · 그 밖 '검출 안 됨'
  B3 사망 쪽 : [0, 3000) 의 250틱 칸마다 끼운 계통(씨앗 팔) 대 끼운 계통(WT 팔 11개 합)의 조사망률 차를 나이 칸(50틱)으로
              Kitagawa 분해: 구성 몫 Σ(c₁−c₂)(m₁+m₂)/2 · 나이별 비율 몫 Σ(m₁−m₂)(c₁+c₂)/2 (한쪽 노출 0 인 나이 칸은 비율 차 0)
              사망을 원인 둘로 따로 분해한다(마른 실행 뒤 고침 — 파일럿에서 '그 밖 사망'의 나이별 비율이 씨앗 팔에서 낮았다):
    B3a 나이 사망 : 구성 몫 비 W_age = 구성 / (구성 + 비율)
              W_age 하한 ≥ 2/3 → '나이 사망 쪽 이득은 대부분 나이 구성(더 젊은 계통) 몫 — 출생 이득의 메아리'
              W_age 상한 ≤ 1/3 → '나이별 수명 사망률 자체가 다르다' · 그 밖 '둘 다'
              (규칙상 수명은 태어날 때 코드와 무관하게 뽑힌다 — W_age ≈ 1 이 아니면 이해나 계측을 의심한다)
    B3b 그 밖 사망(먹힘 · 무작위): 나이별 비율 몫(성장 단위) 하한 > 0 → '씨앗 계통은 같은 나이에서 덜 먹힌다(덜 죽는다)'
              상한 < 0 → '더 죽는다' · 그 밖 '검출 안 됨' · 구성 몫은 서술 · 같은 팔 나머지 계통의 그 밖 사망률 비와 먹기(e) 개체 수는 서술
  B4 출생 쪽 : 같은 분해를 출생률(부모의 나이)에. 나이별 비율 몫(성장 단위 = 250 × 칸 합)
              하한 > 0 → '씨앗은 같은 나이에서 더 자주 낳는다' · 상한 < 0 → '같은 나이에서는 덜 낳는다(출생 이득은 구성 탓)' · 그 밖 '검출 안 됨'
              구성 몫은 서술(더 젊은 계통은 아직 못 낳는 나이가 많아 음수일 수 있다)
  B5 서술   : 닫힘 — 출생 몫 + 사망 몫(계통 셈 적분) 대 g1 차(개체 수), 그리고 − g0 차를 더한 것 대 Δ(3000)
              나이별 출생 비율 몫 중 300틱 이상 나이의 몫(수명 끝에 한 번 더 낳는 몫) · [3000, 6000) 의 같은 분해 · 끼운 계통 몫(3000 · 6000)
  점검: 20 배경 · 체크섬 · 끝 우세 racld · 팔 12 · 재료 보존 · 부모 잃음 0 · 칸 셈(끝 = 시작 + 출생 − 사망) · 놓친 출생 < 1% · 멈춘 팔 0
실행: python3 lifehist_analyze.py [petri 폴더]  →  _RESULT_lifehist8.json
마른 실행: PETRI_ROOT=pilotB PETRI_BGS=101,102,103 PETRI_MAIN=rewind2 python3 lifehist_analyze.py
"""
import json
import math
import os
import random
import sys

BASE = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1 else "~/petri")
ROOT = os.environ.get("PETRI_ROOT", "lh8")
MAIN = os.environ.get("PETRI_MAIN", "main")
T, BIN = 6000, 250
N_BOOT = 2000
RNG = random.Random(20260919)
TS = list(range(500, T + 1, 250))


def pct(xs, p):
    xs = sorted(x for x in xs if x == x)
    k = (len(xs) - 1) * p
    f = math.floor(k)
    c = min(f + 1, len(xs) - 1)
    return xs[f] + (xs[c] - xs[f]) * (k - f)


if os.environ.get("PETRI_BGS"):
    BGS = [int(x) for x in os.environ["PETRI_BGS"].split(",")]
else:
    pick = []
    for s in range(1, 101):
        z = json.load(open(os.path.join(BASE, MAIN, "mat_d8_mu0p01_s%05d.json" % s)))
        if z["extinct_at"] < 0 and z["final_counts"] and z["final_counts"][0][0] == "racld":
            pick.append(s)
    BGS = pick[20:40]

Z, probs = [], []
for s in BGS:
    f = os.path.join(BASE, ROOT, "lh_mat_s%05d_t%d.json" % (s, T))
    if not os.path.exists(f):
        probs.append("없음 %s" % f)
        continue
    z = json.load(open(f))
    kinds = [(a["kind"], a["key"]) for a in z["arms"]]
    if (z["lifehist_version"] != "1.0.0" or z["ticks"] != T or z["bin"] != BIN or not z["checksum_match"]
            or z["grow_top"][0][0] != "racld" or kinds != [("ref", "racld")] + [("neutral", "racld")] * 10 + [("mut", "rascld")]):
        probs.append("설정 · 체크섬 · 팔 %d" % s)
    births = missed = 0
    for a in z["arms"]:
        if not a["cons_ok"] or a["lost_parent"] or a["stopped"] >= 0:
            probs.append("재료 · 부모 · 멈춤 %d %s" % (s, a["labels"]))
        ser = {r[0]: r for r in a["series"]}
        for bn in a["bins"]:
            births += sum(bn["b"][0]) + sum(bn["b"][1])
            t1 = bn["t"] + BIN
            for L in (0, 1):
                end = ser[t1][2] if L == 1 else ser[t1][1] - ser[t1][2]
                if bn["n0"][L] + sum(bn["b"][L]) - sum(bn["da"][L]) - sum(bn["do"][L]) != end:
                    probs.append("칸 셈 %d %s t%d" % (s, a["labels"], bn["t"]))
        missed += a["missed_births"]
    if missed / max(1, births + missed) >= 0.01:
        probs.append("놓친 출생 %d %.3f" % (s, missed / (births + missed)))
    Z.append(z)
print("== 해부 8 B 침입 Δ 가 쌓이는 경로 — 판정")
print("  배경 %s · 파일 %d/%d · 점검 문제 %d" % (",".join(map(str, BGS)), len(Z), len(BGS), len(probs)))
if probs or len(Z) != len(BGS):
    print("🚨 점검 실패 — 판정하지 않는다")
    for p in probs[:20]:
        print("   ", p)
    sys.exit(1)


def g(a, t):
    ser = {r[0]: r for r in a["series"]}
    n0, m0 = a["series"][0][1], a["series"][0][2]
    x = ser[t]
    return math.log((x[2] + 0.5) / m0), math.log((x[1] - x[2] + 0.5) / (n0 - m0))


def kit(bins1, bins2, key, t0, t1):
    """칸 [t0, t1) 의 Kitagawa — 돌려주는 값은 성장 단위(250 × 칸 합) 의 (구성, 비율, 나이 ≥300 비율 몫) — 계통 1 끼리."""
    st = rt = rt_old = 0.0
    for b1, b2 in zip(bins1, bins2):
        if not (t0 <= b1["t"] < t1):
            continue
        o1, o2 = b1["ot"][1], b2["ot"][1]
        e1, e2 = b1[key][1] if key != "d" else [x + y for x, y in zip(b1["da"][1], b1["do"][1])], \
            b2[key][1] if key != "d" else [x + y for x, y in zip(b2["da"][1], b2["do"][1])]
        T1, T2 = sum(o1), sum(o2)
        for a in range(len(o1)):
            c1, c2 = o1[a] / T1, o2[a] / T2
            m1 = e1[a] / o1[a] if o1[a] else None
            m2 = e2[a] / o2[a] if o2[a] else None
            if m1 is None and m2 is None:
                continue
            if m1 is None:
                m1 = m2
            if m2 is None:
                m2 = m1
            st += BIN * (c1 - c2) * (m1 + m2) / 2
            r = BIN * (m1 - m2) * (c1 + c2) / 2
            rt += r
            if a >= 6:
                rt_old += r
    return st, rt, rt_old


def pooled(arms):
    """WT 팔 여럿의 칸 기록을 합친다."""
    out = []
    for j in range(len(arms[0]["bins"])):
        bn = {"t": arms[0]["bins"][j]["t"]}
        for key in ("ot", "b", "da", "do"):
            bn[key] = [[sum(a["bins"][j][key][L][k] for a in arms) for k in range(12)] for L in (0, 1)]
        out.append(bn)
    return out


rows = []
for z in Z:
    wt = [a for a in z["arms"] if a["kind"] in ("ref", "neutral")]
    sd = z["arms"][-1]
    dl = {}
    for t in TS:
        g1s, g0s = g(sd, t)
        gw = [g(a, t) for a in wt]
        g1w = sum(x[0] for x in gw) / len(gw)
        g0w = sum(x[1] for x in gw) / len(gw)
        dl[t] = (g1s - g1w, g0s - g0w)
    delta = {t: v[0] - v[1] for t, v in dl.items()}
    xs = TS
    mx = sum(xs) / len(xs)
    my = sum(delta[t] for t in xs) / len(xs)
    slope = sum((t - mx) * (delta[t] - my) for t in xs) / sum((t - mx) ** 2 for t in xs)
    h1 = [t for t in xs if t <= 3000]
    h2 = [t for t in xs if t >= 3000]

    def sl(ts):
        mx_ = sum(ts) / len(ts)
        my_ = sum(delta[t] for t in ts) / len(ts)
        return sum((t - mx_) * (delta[t] - my_) for t in ts) / sum((t - mx_) ** 2 for t in ts)

    pw = pooled(wt)
    kd = kit(sd["bins"], pw, "d", 0, 3000)
    ka = kit(sd["bins"], pw, "da", 0, 3000)
    ko = kit(sd["bins"], pw, "do", 0, 3000)
    # 나머지 계통(0)의 그 밖 사망률 · 먹기 개체 수 — 씨앗 팔 대 WT 팔 평균, [0, 3000)
    rest = {}
    for nm, arms_ in (("s", [sd]), ("w", wt)):
        ot0 = sum(sum(bn["ot"][0]) for a in arms_ for bn in a["bins"] if bn["t"] < 3000)
        do0 = sum(sum(bn["do"][0]) for a in arms_ for bn in a["bins"] if bn["t"] < 3000)
        ne = [sum(a["bins"][j]["ne"][0] + a["bins"][j]["ne"][1] for a in arms_) / len(arms_) for j in (0, 12, 23)]
        rest["do0_" + nm] = do0 / ot0
        rest["ne_" + nm] = ne
    kb = kit(sd["bins"], pw, "b", 0, 3000)
    kd2 = kit(sd["bins"], pw, "d", 3000, 6000)
    kb2 = kit(sd["bins"], pw, "b", 3000, 6000)
    s3 = {r[0]: r for r in sd["series"]}
    rows.append({
        "seed": z["seed"], "d3000": delta[3000], "d6000": delta[6000], "slope": slope, "s1": sl(h1), "s2": sl(h2),
        "g1": dl[3000][0], "g0": dl[3000][1],
        # 사망 몫은 성장에 − 로 들어간다
        "dS": -kd[0], "dR": -kd[1], "bS": kb[0], "bR": kb[1], "bR_old": kb[2],
        "aS": -ka[0], "aR": -ka[1], "oS": -ko[0], "oR": -ko[1], **rest,
        "dS2": -kd2[0], "dR2": -kd2[1], "bS2": kb2[0], "bR2": kb2[1],
        "f3000": s3[3000][2] / s3[3000][1], "f6000": s3[6000][2] / s3[6000][1], "f0": s3[0][2] / s3[0][1],
    })

idx = [[RNG.randrange(len(rows)) for _ in rows] for _ in range(N_BOOT)]


def mean(key, ix=None):
    rr = [rows[i] for i in ix] if ix is not None else rows
    return sum(r[key] for r in rr) / len(rr)


def boot(f):
    v = [f(ix) for ix in idx]
    return f(None), pct(v, .025), pct(v, .975)


out = {"bgs": BGS, "rows": rows}
d3 = boot(lambda ix: mean("d3000", ix))
v0 = "새 배경에서도 씨앗이 이긴다" if d3[1] > 0 else "씨앗 이득이 검출되지 않는다 — B1~B4 를 이득의 설명으로 읽지 않는다"
print("\n-- B0 재현: Δ(3000) %+.3f [%+.3f, %+.3f] → **%s** · 배경 20 중 양수 %d" % (d3 + (v0, sum(r["d3000"] > 0 for r in rows))))
out["B0"] = {"d3000": d3, "verdict": v0}

S = boot(lambda ix: 3000 * mean("slope", ix) / mean("d3000", ix))
v1 = ("대부분 꾸준히 쌓인다(복리)" if S[1] >= 2 / 3 else "대부분 끼운 직후에 생긴다" if S[2] <= 1 / 3 else "꾸준히 쌓이는 몫과 끼운 직후 몫이 둘 다 보탠다")
r63 = boot(lambda ix: mean("d6000", ix) / mean("d3000", ix))
rh = boot(lambda ix: mean("s2", ix) / mean("s1", ix))
print("-- B1 복리: 기울기 %+.4f /1000틱 · 꾸준한 몫 S = %.2f [%.2f, %.2f] → **%s**" % (1000 * mean("slope"), S[0], S[1], S[2], v1))
print("   (서술) Δ(6000) %+.3f · Δ(6000)/Δ(3000) = %.2f [%.2f, %.2f] · 기울기 뒤/앞 절반 = %.2f [%.2f, %.2f] · 끼운 계통 몫 %.3f → %.3f → %.3f"
      % ((mean("d6000"),) + r63 + rh + (mean("f0"), mean("f3000"), mean("f6000"))))
out["B1"] = {"S": S, "verdict": v1, "ratio_6000_3000": r63, "slope_half_ratio": rh, "slope_per_1000": 1000 * mean("slope"),
             "tag1_share": [mean("f0"), mean("f3000"), mean("f6000")]}

Q0 = boot(lambda ix: -mean("g0", ix) / mean("d3000", ix))
v2 = "나머지 공동체 쪽도 보탠다" if Q0[1] > 0 else "공동체 쪽은 거꾸로 작용한다" if Q0[2] < 0 else "공동체 몫 검출 안 됨"
print("-- B2 공동체: g1 차 %+.3f · g0 차 %+.3f · 공동체 몫 Q0 = %.2f [%.2f, %.2f] → **%s**" % (mean("g1"), mean("g0"), Q0[0], Q0[1], Q0[2], v2))
out["B2"] = {"g1": mean("g1"), "g0": mean("g0"), "Q0": Q0, "verdict": v2}

dsum = boot(lambda ix: mean("dS", ix) + mean("dR", ix))
asum = boot(lambda ix: mean("aS", ix) + mean("aR", ix))
osum = boot(lambda ix: mean("oS", ix) + mean("oR", ix))
print("-- B3 사망 쪽: 사망 몫 %+.3f [%+.3f, %+.3f] = 나이 사망 %+.3f [%+.3f, %+.3f] + 그 밖 사망 %+.3f [%+.3f, %+.3f]"
      % (dsum + asum + osum))
Wa = boot(lambda ix: mean("aS", ix) / (mean("aS", ix) + mean("aR", ix)))
v3a = ("나이 사망 쪽 이득은 대부분 나이 구성(더 젊은 계통) 몫 — 출생 이득의 메아리" if Wa[1] >= 2 / 3
       else "나이별 수명 사망률 자체가 다르다" if Wa[2] <= 1 / 3 else "구성과 나이별 수명 사망률이 둘 다")
print("   B3a 나이 사망: 구성 %+.3f · 나이별 %+.3f · W_age = %.2f [%.2f, %.2f] → **%s**" % ((mean("aS"), mean("aR")) + Wa + (v3a,)))
oR = boot(lambda ix: mean("oR", ix))
v3b = ("씨앗 계통은 같은 나이에서 덜 먹힌다(덜 죽는다)" if oR[1] > 0 else "씨앗 계통은 같은 나이에서 더 죽는다" if oR[2] < 0 else "나이별 그 밖 사망률 차 검출 안 됨")
rr = boot(lambda ix: math.log(mean("do0_s", ix) / mean("do0_w", ix)))
print("   B3b 그 밖 사망: 구성 %+.3f · 나이별 %+.3f [%+.3f, %+.3f] → **%s**" % ((mean("oS"),) + oR + (v3b,)))
print("       (서술) 나머지 계통의 그 밖 사망률 ln(씨앗 팔/WT 팔) %+.3f [%+.3f, %+.3f] · 먹기(e) 개체 수 0 · 3000 · 5750틱: 씨앗 팔 %s · WT 팔 %s"
      % (rr + ("/".join("%.0f" % mean_i for mean_i in [sum(r["ne_s"][j] for r in rows) / len(rows) for j in range(3)]),
               "/".join("%.0f" % mean_i for mean_i in [sum(r["ne_w"][j] for r in rows) / len(rows) for j in range(3)]))))
out["B3"] = {"total": dsum, "age": asum, "other": osum, "W_age": Wa, "verdict_age": v3a, "other_rate": oR, "verdict_other": v3b,
             "rest_other_lnratio": rr, "struct": mean("dS"), "rate": mean("dR")}

bR = boot(lambda ix: mean("bR", ix))
bS = boot(lambda ix: mean("bS", ix))
v4 = ("씨앗은 같은 나이에서 더 자주 낳는다" if bR[1] > 0 else "같은 나이에서는 덜 낳는다(출생 이득은 구성 탓)" if bR[2] < 0 else "나이별 출생률 차 검출 안 됨")
bsum = boot(lambda ix: mean("bS", ix) + mean("bR", ix))
print("-- B4 출생 쪽: 출생 몫 %+.3f [%+.3f, %+.3f] = 구성 %+.3f [%+.3f, %+.3f] + 나이별 %+.3f [%+.3f, %+.3f] → **%s**"
      % (bsum + bS + bR + (v4,)))
out["B4"] = {"total": bsum, "struct": bS, "rate": bR, "verdict": v4}

old = boot(lambda ix: mean("bR_old", ix) / mean("bR", ix))
print("\n-- B5 서술")
print("   닫힘: 출생 몫 + 사망 몫 = %+.3f · g1 차(개체 수) %+.3f · 여기에 − g0 차 = %+.3f · Δ(3000) %+.3f"
      % (bsum[0] + dsum[0], mean("g1"), bsum[0] + dsum[0] - mean("g0"), d3[0]))
print("   나이별 출생 비율 몫 중 300틱 이상 나이의 몫 %.2f [%.2f, %.2f]" % old)
print("   [3000, 6000): 출생 몫 구성 %+.3f · 나이별 %+.3f · 사망 몫 구성 %+.3f · 나이별 %+.3f"
      % (mean("bS2"), mean("bR2"), mean("dS2"), mean("dR2")))
out["B5"] = {"closure": [bsum[0] + dsum[0], mean("g1"), bsum[0] + dsum[0] - mean("g0"), d3[0]], "old_share": old,
             "late": {"bS": mean("bS2"), "bR": mean("bR2"), "dS": mean("dS2"), "dR": mean("dR2")}}

json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "_RESULT_lifehist8.json"), "w"), ensure_ascii=False, indent=1)
print("\n저장: _RESULT_lifehist8.json")
