# -*- coding: utf-8 -*-
"""해부 7 ② 공급 · ③a 교체 경로 판정 — 사전등록 해부7-사전등록-2026-09-16.md §2 · §3 · §4. 결과 전에 썼다.

  같은 시드(701~730)를 원래 규칙 · 재료 먼저 규칙으로 돌려 짝짓는다. 밀도 8 · μ 0.1 · 0.5 · 1%.
  u_obs = WT 부모 출생 중 변종 몫 · u_nom = 1 − (1 − μ(1/3 + 1/3 + 10/11))^L · 부풀림 = u_obs/u_nom
  순도 = (끝 WT 수 + 0.5)/(끝 개체 + 0.5)  (0 을 피하려는 반 개 — 미리 정함)
  시드 재추출 부트스트랩 2,000회 95% · 멸종 접시는 그 짝을 뺀다(수 보고)

  F1 조작 점검  : 재료 먼저 · 세 코드 · 세 μ 부풀림 점추정이 **모두 < 1.5** — 아니면 ② 공급 판정은 하지 않는다
  F2 공급       : 초과분 e = ln(u 씨앗/u racld) − ln S_len. 사라진 몫 1 − e_재료먼저/e_원래 ≥ 2/3 (세 μ 모두)
                  → '재료 먼저 세계에서 씨앗의 추가 변종 공급이 대부분 사라진다'
  F3 순도       : ln(순도 씨앗/racld) 의 (재료 먼저 − 원래) 하한 > 0 (세 μ 모두) → '격차가 준다'
  P1 경로       : 원래 규칙 씨앗 접시 · 새로 생긴 racld(부모가 racld 아님) 중 부모가 씨앗 그대로인 몫
                  (접시별 몫의 평균, 새로 생김이 있는 접시만) — **μ 마다** ≥ 0.5 '씨앗에서 글자 하나 빠져 바로 생긴다' · 아니면 '다른 변종을 거쳐 생긴다'
                  세 μ 가 같으면 그 읽기 · 다르면 'μ 에 따라 다르다' (마른 실행에서 고침 — 세 μ 평균이 0.73 · 0.25 · 0.02 를 가렸다)
  P2 퍼짐       : 원래 규칙 · racld 출생 중 부모도 racld 인 몫 — **μ 마다** ≥ 0.9 '생긴 뒤 스스로 퍼진다' · < 0.5 '계속 새로 생겨 쌓인다' · 사이 '섞임'
  P3 제자리     : 원래 규칙 · 같은 접시에서 ln(racld 로 태어난 개체의 **같은 코드 자식 수** / 씨앗으로 태어난 개체의 같은 코드 자식 수) 하한 > 0
                  (세 μ 모두) → '씨앗 접시 안에서 racld 가 제 코드를 더 많이 남긴다'
                  자식 총수의 비(번식력)는 서술로 함께 적는다 (마른 실행에서 고침 — 총수 비는 음수인데 racld 가 퍼졌다: 번식력과 충실도가 섞였다)
  P4 규칙 효과  : 끝 racld 몫(씨앗 접시) 재료 먼저 − 원래 — 서술(구간)
  (μ 는 따로 판정하고 '세 μ 모두' 가 붙은 줄은 셋 다 맞아야 한다)
실행: python3 ff_lin_analyze.py [petri 폴더]  →  _RESULT_ff_lin.json
"""
import json
import math
import os
import random
import sys

BASE = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1 else "~/petri")
ROOT = os.environ.get("PETRI_ROOT", "lin7")
SEEDS = [int(x) for x in os.environ["PETRI_SEEDS"].split(",")] if os.environ.get("PETRI_SEEDS") else list(range(701, 731))
MUS = [0.001, 0.005, 0.01]
CODES = {"racld": "rascld", "rascld": "racld", "rsacld": "racld"}
RULES = ["orig", "ff"]
N_BOOT = 2000
RNG = random.Random(20261006)


def pct(xs, p):
    xs = sorted(xs)
    k = (len(xs) - 1) * p
    f = math.floor(k)
    c = min(f + 1, len(xs) - 1)
    return xs[f] + (xs[c] - xs[f]) * (k - f)


def boot(v):
    if not v:
        return (float("nan"),) * 3
    m = sum(v) / len(v)
    bs = [sum(v[RNG.randrange(len(v))] for _ in v) / len(v) for _ in range(N_BOOT)]
    return m, pct(bs, 0.025), pct(bs, 0.975)


def unom(mu, L):
    return 1 - (1 - mu * (2 / 3 + 10 / 11)) ** L


def fname(rule, code, mu, s):
    return os.path.join(BASE, ROOT, rule, "lin_%s_mu%s_s%05d%s.json" % (code, str(mu).replace(".", "p"), s, "_ff" if rule == "ff" else ""))


D, missing, bad = {}, [], []
for rule in RULES:
    for code, watch in CODES.items():
        for mu in MUS:
            for s in SEEDS:
                f = fname(rule, code, mu, s)
                if not os.path.exists(f):
                    missing.append((rule, code, mu, s)); continue
                z = json.load(open(f))
                if z["sim_version"] not in ("0.3.3", "0.3.4", "0.3.5") or z["find_first"] != (rule == "ff") or z["watch"] != watch:
                    bad.append((rule, code, mu, s))
                D[(rule, code, mu, s)] = z
ext = sum(1 for z in D.values() if z["extinct_at"] >= 0)
miss = max((z["missed_births"] / max(z["births"], 1) for z in D.values()), default=1)
print("== 해부 7 ② 공급 · ③a 교체 경로 — 판정")
print("  파일 %d/%d · 설정 어긋남 %d · 멸종 %d · 놓친 출생 몫 최대 %.3f%%"
      % (len(D), len(RULES) * len(CODES) * len(MUS) * len(SEEDS), len(bad), ext, 100 * miss))
if missing or bad or miss >= 0.01:
    print("🚨 점검 실패 — 판정하지 않는다")
    for m in (missing + bad)[:10]:
        print("    ", m)
    sys.exit(1)


def st(z):
    if z["extinct_at"] >= 0 or not z["wt_parent_births"]:
        return None
    last = z["series"][-1]
    u = z["wt_parent_mut_births"] / z["wt_parent_births"]
    c = z["classes"]
    kw = c["wt"][1] / c["wt"][0] if c["wt"][0] else float("nan")
    wc = z["watch_class"]
    kx = wc[1] / wc[0] if wc[0] else float("nan")
    wp = z["watch_births_by_parent"]
    denovo = wp["wt"] + wp["other"]
    return {"u": u, "amp": u / unom(z["mu"], len(z["code"])), "pure": (last[2] + 0.5) / (last[1] + 0.5),
            "watch_share": last[6] / last[1] if last[1] else float("nan"),
            "from_wt": wp["wt"] / denovo if denovo else None,
            "inherit": wp["same"] / (wp["same"] + denovo) if (wp["same"] + denovo) else None,
            "k_ratio": math.log(kx / kw) if kx > 0 and kw > 0 else None, "denovo": denovo, "first": z["watch_first_denovo_tick"],
            "same_ratio": (math.log((z["same_kids_watch"] / wc[0]) / (z["same_kids_wt"] / c["wt"][0]))
                           if wc[0] and c["wt"][0] and z["same_kids_watch"] and z["same_kids_wt"] else None)}


S = {k: st(z) for k, z in D.items()}
out = {"table": {}, "F1": {}, "F2": {}, "F3": {}, "P": {}}

print("\n-- 부풀림 · 순도 (시드 평균)")
print("   %-5s %-7s %6s | %-22s | %-10s" % ("규칙", "코드", "μ", "부풀림 [95%]", "순도"))
for rule in RULES:
    for code in CODES:
        for mu in MUS:
            xs = [S[(rule, code, mu, s)] for s in SEEDS if S[(rule, code, mu, s)]]
            a = boot([x["amp"] for x in xs]); p = boot([x["pure"] for x in xs])
            out["table"]["%s|%s|%s" % (rule, code, mu)] = {"amp": a, "pure": p, "n": len(xs)}
            print("   %-5s %-7s %5.2f%% | %5.2f [%5.2f,%5.2f] | %.4f  n=%d" % (rule, code, mu * 100, *a, p[0], len(xs)))

f1 = max(out["table"]["ff|%s|%s" % (c, mu)]["amp"][0] for c in CODES for mu in MUS)
f1ok = f1 < 1.5
out["F1"] = {"max_amp_ff": f1, "passes": f1ok}
print("\n-- F1 조작 점검: 재료 먼저 부풀림 점추정 최대 %.2f → %s" % (f1, "**통과 — 재료 먼저 세계에서는 부풀림이 사라졌다**" if f1ok else "🚨 조작이 먹지 않았다 — ② 공급 판정은 하지 않는다"))


def pair(rule, a, b, key, mu):
    v = {}
    for s in SEEDS:
        x, y = S[(rule, a, mu, s)], S[(rule, b, mu, s)]
        if x and y and x[key] > 0 and y[key] > 0:
            v[s] = math.log(x[key] / y[key])
    return v


if f1ok:
    print("\n-- F2 씨앗의 추가 공급 (초과분 e = ln S − ln S_len)")
    f2 = True
    for mu in MUS:
        sl = math.log(unom(mu, 6) / unom(mu, 5))
        eo = {s: v - sl for s, v in pair("orig", "rascld", "racld", "u", mu).items()}
        ef = {s: v - sl for s, v in pair("ff", "rascld", "racld", "u", mu).items()}
        mo, mf = boot(list(eo.values())), boot(list(ef.values()))
        common = [s for s in eo if s in ef]
        gone_bs = []
        for _ in range(N_BOOT):
            ix = [common[RNG.randrange(len(common))] for _ in common]
            a = sum(eo[s] for s in ix) / len(ix); b = sum(ef[s] for s in ix) / len(ix)
            if a > 0:
                gone_bs.append(1 - b / a)
        gone = 1 - mf[0] / mo[0] if mo[0] > 0 else float("nan")
        f2 = f2 and gone >= 2 / 3
        out["F2"][str(mu)] = {"e_orig": mo, "e_ff": mf, "gone": gone, "gone_ci": [pct(gone_bs, 0.025), pct(gone_bs, 0.975)]}
        print("   μ %.2f%%: e 원래 %+.3f [%+.3f,%+.3f] · 재료먼저 %+.3f [%+.3f,%+.3f] · 사라진 몫 %.2f [%.2f, %.2f]"
              % (mu * 100, *mo, *mf, gone, pct(gone_bs, 0.025), pct(gone_bs, 0.975)))
    out["F2"]["verdict"] = f2
    print("   → %s" % ("**재료 먼저 세계에서 씨앗의 추가 변종 공급이 대부분 사라진다**(세 μ 모두 ≥ 2/3)" if f2 else "대부분 사라진다고 할 수 없는 μ 가 있다"))

    print("\n-- F3 순도 격차 ln(순도 씨앗/racld) — 재료 먼저 − 원래")
    f3 = True
    for mu in MUS:
        po, pf = pair("orig", "rascld", "racld", "pure", mu), pair("ff", "rascld", "racld", "pure", mu)
        d = boot([pf[s] - po[s] for s in po if s in pf])
        f3 = f3 and d[1] > 0
        out["F3"][str(mu)] = {"orig": boot(list(po.values())), "ff": boot(list(pf.values())), "diff": d}
        print("   μ %.2f%%: 원래 %+.2f · 재료먼저 %+.2f · 차 %+.2f [%+.2f, %+.2f]"
              % (mu * 100, out["F3"][str(mu)]["orig"][0], out["F3"][str(mu)]["ff"][0], *d))
    out["F3"]["verdict"] = f3
    print("   → %s" % ("**재료 먼저 세계에서 순도 격차가 준다**(세 μ 모두)" if f3 else "세 μ 모두에서 준다고 할 수 없다"))

print("\n-- ③a 씨앗 접시가 racld 로 바뀌는 경로 (원래 규칙 · 재료 먼저는 서술)")
for rule in RULES:
    for mu in MUS:
        xs = [S[(rule, "rascld", mu, s)] for s in SEEDS if S[(rule, "rascld", mu, s)]]
        fw = boot([x["from_wt"] for x in xs if x["from_wt"] is not None])
        ih = boot([x["inherit"] for x in xs if x["inherit"] is not None])
        kr = boot([x["k_ratio"] for x in xs if x["k_ratio"] is not None])
        sr = boot([x["same_ratio"] for x in xs if x["same_ratio"] is not None])
        sh = boot([x["watch_share"] for x in xs])
        fs = sorted(x["first"] for x in xs if x["first"] is not None and x["first"] >= 0)
        out["P"]["%s|%s" % (rule, mu)] = {"from_wt": fw, "inherit": ih, "k_ratio": kr, "same_ratio": sr, "share": sh,
                                          "first_median": fs[len(fs) // 2] if fs else None, "n_denovo_dishes": len(fs)}
        print("   %-5s μ %.2f%%: 새로 생김 중 씨앗 부모 %.2f [%.2f,%.2f] · 물려받음 %.3f [%.3f,%.3f] · ln 같은코드자식비 %+.2f [%+.2f,%+.2f] · (ln 자식총수비 %+.2f) · 끝 racld 몫 %.3f · 처음 생긴 틱(중앙) %s"
              % (rule, mu * 100, *fw, *ih, *sr, kr[0], sh[0], out["P"]["%s|%s" % (rule, mu)]["first_median"]))
po = [out["P"]["orig|%s" % mu] for mu in MUS]
l1 = ["씨앗에서 글자 하나 빠져 바로 생긴다" if p["from_wt"][0] >= 0.5 else "다른 변종을 거쳐 생긴다" for p in po]
l2 = ["생긴 뒤 스스로 퍼진다" if p["inherit"][0] >= 0.9 else ("계속 새로 생겨 쌓인다" if p["inherit"][0] < 0.5 else "섞임") for p in po]
p1 = l1[0] if len(set(l1)) == 1 else "μ 에 따라 다르다"
p2 = l2[0] if len(set(l2)) == 1 else "μ 에 따라 다르다"
p3 = all(p["same_ratio"][1] > 0 for p in po)
out["P"]["P1"] = {"per_mu": l1, "verdict": p1}
out["P"]["P2"] = {"per_mu": l2, "verdict": p2}
out["P"]["P3"] = {"verdict": p3}
print("   P1 새로 생긴 racld 의 부모 (μ 마다: %s) → **%s**" % (" · ".join(l1), p1))
print("   P2 racld 출생 중 물려받은 몫 (μ 마다: %s) → **%s**" % (" · ".join(l2), p2))
print("   P3 → %s" % ("**씨앗 접시 안에서 racld 가 제 코드를 더 많이 남긴다**(세 μ 모두)" if p3 else "세 μ 모두에서 제 코드를 더 많이 남긴다고 할 수 없다"))
print("   P4 끝 racld 몫 재료 먼저 − 원래:")
for mu in MUS:
    d = boot([S[("ff", "rascld", mu, s)]["watch_share"] - S[("orig", "rascld", mu, s)]["watch_share"]
              for s in SEEDS if S[("ff", "rascld", mu, s)] and S[("orig", "rascld", mu, s)]])
    out["P"]["P4|%s" % mu] = d
    print("      μ %.2f%%: %+.3f [%+.3f, %+.3f]" % (mu * 100, *d))

json.dump(out, open(os.path.join(BASE, os.environ.get("PETRI_OUT", "_RESULT_ff_lin.json")), "w"), ensure_ascii=False, indent=1)
print("\n저장: " + os.environ.get("PETRI_OUT", "_RESULT_ff_lin.json"))
