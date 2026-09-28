# -*- coding: utf-8 -*-
"""해부 6A 판정 — 계보 추적. 사전등록 해부6-사전등록-2026-09-16.md §2 · §5. 결과 전에 썼다.

  접시마다: u_obs = WT 부모의 출생 중 변종 몫 · u_nom = 1 − (1 − μ(1/3 + 1/3 + 10/11))^L (헛손질이 없을 때의 몫)
            부풀림 = u_obs / u_nom · r1 = k1 / kWT (한 걸음 변종의 자식 수 ÷ WT 의 자식 수, 창 [5000, T−700])
            순도 = 끝 개체 중 WT 몫
  두 코드 비교는 같은 시드끼리 로그 비 · 시드 재추출 부트스트랩 2,000회 95% · 멸종 접시는 그 짝을 뺀다(수 보고)

  A1 부풀림     : 밀도 8 에서 두 코드 모두 모든 μ 에서 부풀림 하한 > 3 → '헛손질이 복제 오류를 부풀린다'
                  밀도 32 에서 점추정 < 1.5 → '재료가 넉넉하면 부풀림이 없다'
  A2 다시 시도  : 밀도 8 · ln(u_obs(rascld) / u_obs(rsacld)) 하한 > 0 (같은 6글자 · 고리 2틱 대 3틱)
                  → '같은 길이에서 빨리 다시 시도하는 코드가 변종을 더 낳는다'
  A3 밀도       : S = u_obs(rascld)/u_obs(racld). ln S(밀도 8) − ln S(밀도 32) 하한 > 0 → '씨앗의 추가 변종 공급은 헛손질에서 온다'
                  초과분 e = ln S − ln S_len (S_len = u_nom(6)/u_nom(5), 길이만의 비). 밀도 32 에서 사라진 몫 1 − e32/e8 ≥ 2/3 이면
                  '넉넉하면 초과분이 대부분 사라진다' · 그 밖 '일부만 사라진다'
                  (마른 실행에서 고침 — 처음 판정선 'S(32) 구간이 S_len 을 품는다' 는 n 이 크면 근사 오차만으로 떨어진다)
  A4 몫 가르기  : 밀도 8 · |ln S| 와 |ln R| (R = r1(rascld)/r1(racld)) — 더 큰 쪽이 순도 차이의 주된 몫
  A5 순도       : 밀도 8 · 순도 rascld < racld (로그 비 상한 < 0)
  A6 먹혔나     : 밀도 8 씨앗 접시 끝에서 `racld` 가 차지한 몫(시드 평균) — 20% 이상이면 🚩 '씨앗 접시의 순도 손실 일부는 racld 에게 먹힌 것'
                  (마른 실행에서 더함 — 밀도 32 racld 접시가 변종에게 통째로 먹힌 것을 보고)
  점검: 파일 완결성 · 놓친 출생 몫 < 1% · μ 마다 u_obs > 0
실행: python3 lineage_analyze.py [petri 폴더]  →  _RESULT_lineage.json
"""
import glob
import json
import math
import os
import random
import sys

BASE = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1 else "~/petri")
ROOT = os.environ.get("PETRI_ROOT", "lin6")
SEEDS = [int(x) for x in os.environ["PETRI_SEEDS"].split(",")] if os.environ.get("PETRI_SEEDS") else list(range(601, 631))
MUS = [0.001, 0.005, 0.01]
PLAN = {"d08": ["racld", "rascld", "rsacld"], "d32": ["racld", "rascld"]}
N_BOOT = 2000
RNG = random.Random(20261001)


def pct(xs, p):
    xs = sorted(xs)
    k = (len(xs) - 1) * p
    f = math.floor(k)
    c = min(f + 1, len(xs) - 1)
    return xs[f] + (xs[c] - xs[f]) * (k - f)


def boot(vals):
    if not vals:
        return (float("nan"),) * 3
    m = sum(vals) / len(vals)
    bs = []
    for _ in range(N_BOOT):
        bs.append(sum(vals[RNG.randrange(len(vals))] for _ in vals) / len(vals))
    return m, pct(bs, 0.025), pct(bs, 0.975)


def unom(mu, L):
    return 1 - (1 - mu * (1 / 3 + 1 / 3 + 10 / 11)) ** L


def mtag(mu):
    return str(mu).replace(".", "p")


D, missing = {}, []
for dd, codes in PLAN.items():
    for code in codes:
        for mu in MUS:
            for s in SEEDS:
                f = os.path.join(BASE, ROOT, dd, "lin_%s_mu%s_s%05d.json" % (code, mtag(mu), s))
                if os.path.exists(f):
                    D[(dd, code, mu, s)] = json.load(open(f))
                else:
                    missing.append((dd, code, mu, s))
print("== 해부 6A 계보 추적 — 판정")
n_plan = sum(len(c) for c in PLAN.values()) * len(MUS) * len(SEEDS)
miss_frac = max((z["missed_births"] / max(z["births"], 1) for z in D.values()), default=0)
zero_u = [k for k, z in D.items() if z["extinct_at"] < 0 and z["wt_parent_births"] and z["wt_parent_mut_births"] == 0]
ext = sum(1 for z in D.values() if z["extinct_at"] >= 0)
print("  파일 %d/%d · 멸종 접시 %d · 놓친 출생 몫 최대 %.3f%% · u_obs = 0 인 접시 %d"
      % (len(D), n_plan, ext, 100 * miss_frac, len(zero_u)))
if missing or miss_frac >= 0.01 or zero_u:
    print("🚨 점검 실패 — 판정하지 않는다")
    for m in missing[:10]:
        print("    없음", m)
    sys.exit(1)


def stats(z):
    if z["extinct_at"] >= 0 or not z["wt_parent_births"]:
        return None
    c = z["classes"]
    kw = c["wt"][1] / c["wt"][0] if c["wt"][0] else float("nan")
    k1 = c["d1"][1] / c["d1"][0] if c["d1"][0] else float("nan")
    last = z["series"][-1]
    uo = z["wt_parent_mut_births"] / z["wt_parent_births"]
    top = dict(z["top_end"])
    return {"u": uo, "amp": uo / unom(z["mu"], len(z["code"])), "r1": k1 / kw if kw > 0 else float("nan"),
            "kw": kw, "k1": k1, "pure": last[2] / last[1] if last[1] else float("nan"), "n1": c["d1"][0],
            "racld_share": top.get("racld", 0) / last[1] if last[1] else float("nan"),
            "top_mut": next(((k, v / last[1]) for k, v in z["top_end"] if k != z["code"]), ("-", 0.0))}


S = {k: stats(z) for k, z in D.items()}
out = {"table": {}, "A1": {}, "A2": {}, "A3": {}, "A4": {}, "A5": {}}

print("\n-- 코드별 (시드 평균 · 95%)")
print("   %-4s %-7s %6s | %-24s | %-22s | %-22s | %-22s" % ("밀도", "코드", "μ", "부풀림 u_obs/u_nom", "u_obs", "r1 = k1/kWT", "순도(끝 WT 몫)"))
for dd, codes in PLAN.items():
    for code in codes:
        for mu in MUS:
            xs = [S[(dd, code, mu, s)] for s in SEEDS if S[(dd, code, mu, s)]]
            row = {}
            for key in ("amp", "u", "r1", "pure"):
                v = [x[key] for x in xs if not math.isnan(x[key])]
                row[key] = boot(v)
            out["table"]["%s|%s|%s" % (dd, code, mu)] = dict(row, n=len(xs))
            print("   %-4s %-7s %5.2f%% | %6.2f [%5.2f,%6.2f] | %.4f [%.4f,%.4f] | %.3f [%.3f,%.3f] | %.3f [%.3f,%.3f]  n=%d"
                  % (dd, code, mu * 100, *row["amp"], *row["u"], *row["r1"], *row["pure"], len(xs)))


def paired_by_seed(dd, a, b, key, mu):
    v = {}
    for s in SEEDS:
        x, y = S[(dd, a, mu, s)], S[(dd, b, mu, s)]
        if x and y and not (math.isnan(x[key]) or math.isnan(y[key])) and x[key] > 0 and y[key] > 0:
            v[s] = math.log(x[key] / y[key])
    return v


def paired(dd, a, b, key, mu):
    return list(paired_by_seed(dd, a, b, key, mu).values())


# A1
print("\n-- A1 부풀림")
lo8 = min(out["table"]["d08|%s|%s" % (c, mu)]["amp"][1] for c in PLAN["d08"] for mu in MUS)
hi32 = max(out["table"]["d32|%s|%s" % (c, mu)]["amp"][0] for c in PLAN["d32"] for mu in MUS)
a1a = lo8 > 3
a1b = hi32 < 1.5
print("   밀도 8 부풀림 하한의 최소 %.2f → %s" % (lo8, "**헛손질이 복제 오류를 부풀린다**" if a1a else "부풀림이 3배를 넘는다고 할 수 없다"))
print("   밀도 32 부풀림 점추정의 최대 %.2f → %s" % (hi32, "**재료가 넉넉하면 부풀림이 없다**" if a1b else "밀도 32 에서도 부풀림이 남는다"))
out["A1"] = {"min_lo_d8": lo8, "max_pt_d32": hi32, "amplified_d8": a1a, "none_d32": a1b}

# A2
print("\n-- A2 다시 시도 (밀도 8 · 같은 6글자: rascld 고리 2틱 대 rsacld 3틱)")
a2 = True
for mu in MUS:
    m, lo, hi = boot(paired("d08", "rascld", "rsacld", "u", mu))
    a2 = a2 and lo > 0
    out["A2"][str(mu)] = [m, lo, hi]
    print("   μ %.2f%%: ln(u_obs 씨앗/rsacld) %+.3f [%+.3f, %+.3f] (비 %.2f)" % (mu * 100, m, lo, hi, math.exp(m)))
print("   → %s" % ("**같은 길이에서 빨리 다시 시도하는 코드가 변종을 더 낳는다**(세 μ 모두)" if a2 else "세 μ 모두에서 그렇다고 할 수 없다"))
out["A2"]["verdict"] = a2

# A3
print("\n-- A3 밀도 (S = u_obs 씨앗/racld)")
slen = unom(0.01, 6) / unom(0.01, 5)
a3a, a3b = True, True
for mu in MUS:
    v8 = paired_by_seed("d08", "rascld", "racld", "u", mu)
    v32 = paired_by_seed("d32", "rascld", "racld", "u", mu)
    s8, s32 = boot(list(v8.values())), boot(list(v32.values()))
    dif = boot([v8[k] - v32[k] for k in v8 if k in v32])        # 같은 시드끼리
    sl = math.log(unom(mu, 6) / unom(mu, 5))
    a3a = a3a and dif[1] > 0
    e8, e32 = s8[0] - sl, s32[0] - sl
    gone = 1 - e32 / e8 if e8 > 0 else float("nan")
    a3b = a3b and gone >= 2 / 3
    out["A3"][str(mu)] = {"lnS_d8": s8, "lnS_d32": s32, "diff": dif, "ln_Slen": sl, "excess_gone": gone}
    print("   μ %.2f%%: S 밀도8 %.2f [%.2f,%.2f] · 밀도32 %.2f [%.2f,%.2f] · 길이만 %.2f · ln 차 %+.3f [%+.3f,%+.3f] · 초과분 사라진 몫 %.2f"
          % (mu * 100, math.exp(s8[0]), math.exp(s8[1]), math.exp(s8[2]), math.exp(s32[0]), math.exp(s32[1]), math.exp(s32[2]),
             math.exp(sl), *dif, gone))
print("   → %s" % ("**씨앗의 추가 변종 공급은 헛손질에서 온다**(세 μ 모두 밀도 8 > 밀도 32)" if a3a else "밀도 8 과 32 의 차이가 세 μ 모두에서 서지 않는다"))
print("   → %s" % ("**넉넉하면 초과분이 대부분(≥2/3) 사라진다**(세 μ 모두)" if a3b else "초과분이 일부만 사라지는 μ 가 있다"))
out["A3"]["verdict_diff"] = a3a
out["A3"]["verdict_len"] = a3b

# A4
print("\n-- A4 순도 차이의 몫 (밀도 8 · 씨앗 대 racld)")
for mu in MUS:
    ls = boot(paired("d08", "rascld", "racld", "u", mu))
    lr = boot(paired("d08", "rascld", "racld", "r1", mu))
    main = "공급(변종을 더 낳는다)" if abs(ls[0]) > abs(lr[0]) else "번식(변종이 더 잘 산다)"
    out["A4"][str(mu)] = {"lnS": ls, "lnR": lr, "main": main}
    print("   μ %.2f%%: ln S %+.3f [%+.3f,%+.3f] · ln R %+.3f [%+.3f,%+.3f] → 주된 몫 = **%s**" % (mu * 100, *ls, *lr, main))

# A5
print("\n-- A5 순도 (밀도 8 · 씨앗 대 racld · 로그 비)")
a5 = True
for mu in MUS:
    p = boot(paired("d08", "rascld", "racld", "pure", mu))
    a5 = a5 and p[2] < 0
    out["A5"][str(mu)] = p
    print("   μ %.2f%%: ln(순도 씨앗/racld) %+.3f [%+.3f, %+.3f] (비 %.3f)" % (mu * 100, *p, math.exp(p[0])))
print("   → %s" % ("**씨앗 접시가 덜 순수하다**(세 μ 모두)" if a5 else "세 μ 모두에서 그렇다고 할 수 없다"))
out["A5"]["verdict"] = a5

# A6
print("\n-- A6 씨앗 접시가 racld 에게 먹혔나 (밀도 8 · 끝 개체 중 racld 몫 · 가장 흔한 변종)")
out["A6"] = {}
flag = False
for mu in MUS:
    xs = [S[("d08", "rascld", mu, s)] for s in SEEDS if S[("d08", "rascld", mu, s)]]
    rs = boot([x["racld_share"] for x in xs])
    from collections import Counter
    tops = Counter(x["top_mut"][0] for x in xs).most_common(3)
    flag = flag or rs[0] >= 0.20
    out["A6"][str(mu)] = {"racld_share": rs, "top_mut_modes": tops}
    print("   μ %.2f%%: racld 몫 %.3f [%.3f, %.3f] · 접시마다 가장 흔한 변종(최빈) %s" % (mu * 100, *rs, tops))
print("   → %s" % ("🚩 **씨앗 접시의 순도 손실 일부는 racld 에게 먹힌 것이다**(20% 이상인 μ 가 있다)" if flag else "racld 에게 먹힌 것은 아니다(모든 μ 에서 20% 미만)"))
out["A6"]["flag"] = flag

json.dump(out, open(os.path.join(BASE, os.environ.get("PETRI_OUT", "_RESULT_lineage.json")), "w"), ensure_ascii=False, indent=1)
print("\n저장: " + os.environ.get("PETRI_OUT", "_RESULT_lineage.json"))
