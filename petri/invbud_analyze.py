# -*- coding: utf-8 -*-
"""해부 7 ③b 판정 — 오류가 없을 때(μ 0) 씨앗이 찬 접시 침입에서 이기는 이유. 사전등록 해부7-사전등록-2026-09-16.md §3 · §4.
결과 전에 썼다.

  침입: 20 배경 · 팔 ref(racld) · rascld · rsacld · racldx · 3,000틱 · 표지 1(끼운 계통)만 본다 — 팔끼리는 **같은 개체를 같은 자리에
        끼운** 짝이다. 비교는 늘 'X 팔의 표지 1' 대 'ref 팔의 표지 1' (배경 안 짝 · 배경 재추출 2,000회 95%).
  단일: racld · rascld · rsacld × 시드 801~820 · 15,000틱 이후 칸 (같은 시드끼리 짝)

  G0 회귀   : 침입 팔 80개의 계통 개체 수 기록이 해부 1(dms_main)의 같은 팔과 같다 — 하나라도 다르면 중단
  G1 무엇이 : 개체-틱당 출생 · 사망을 시간으로 적분해 계통 성장을 가른다(두 팔의 공통 칸만)
              Δ출생 = ∫b(씨앗 팔) − ∫b(ref 팔) · Δ사망 = −(∫d(씨앗 팔) − ∫d(ref 팔))
              |Δ출생| > |Δ사망| → '주로 출생 쪽' · 아니면 '주로 사망 쪽' (점추정 · 구간 병기)
              사망 쪽이면 원인(나이 · 차례 전 · 그 밖)별 몫도 적는다
              서술: 같은 방식의 표지0(나머지 공동체) 쪽 차이 — s 는 두 계통의 차이라 표지1 쪽만으로는 Δ 가 다 설명되지 않는다
              (마른 실행에서 더함 — 표지1 쪽 합이 Δ 의 절반뿐이었다)
  G2 가로채기 : 글자당 시간 t = 복사 중인 틱 / 복사 성공
              섞인 비 = t(씨앗 팔 표지1) / t(ref 팔 표지1) · 단일 비 = t(단일 씨앗) / t(단일 racld)
              ln 섞인 비 − ln 단일 비 : 상한 < 0 → '섞이면 씨앗이 글자를 더 빨리 얻는다 — 가로채기'
                                     하한 > 0 → '섞이면 씨앗이 글자를 더 늦게 얻는다'
                                     그 밖   → '섞여도 글자 얻는 속도 비는 같다'
              (두 비는 서로 다른 접시라 따로 재추출한다) · rsacld 도 같은 방식으로 서술
  G3 서술   : 출생/틱 · 사망/틱 · 평균 나이 · 출생당 헛손질 · 막힘 · 글자당 시간 · 출생당 개체-틱
  점검      : 배경 체크섬 · 팔 4 · μ0 빠뜨림 0 · 모름 몫 < 1% · 단일 접시 완결성(멸종은 그 짝을 빼고 수 보고)
실행: python3 invbud_analyze.py [petri 폴더]  →  _RESULT_invbud.json
"""
import json
import math
import os
import random
import sys

BASE = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1 else "~/petri")
ROOT = os.environ.get("PETRI_ROOT", "ib7")
NBG = int(os.environ.get("PETRI_NBG", "20"))
MSEEDS = [int(x) for x in os.environ["PETRI_SEEDS"].split(",")] if os.environ.get("PETRI_SEEDS") else list(range(801, 821))
CODES = ["rascld", "rsacld", "racldx"]
MONO = ["racld", "rascld", "rsacld"]
N_BOOT = 2000
RNG = random.Random(20261008)


def pct(xs, p):
    xs = sorted(xs)
    k = (len(xs) - 1) * p
    f = math.floor(k)
    c = min(f + 1, len(xs) - 1)
    return xs[f] + (xs[c] - xs[f]) * (k - f)


def boot(v):
    m = sum(v) / len(v)
    bs = [sum(v[RNG.randrange(len(v))] for _ in v) / len(v) for _ in range(N_BOOT)]
    return m, pct(bs, 0.025), pct(bs, 0.975), bs


pick = []
for s in range(1, 101):
    z = json.load(open(os.path.join(BASE, "main", "mat_d8_mu0p01_s%05d.json" % s)))
    if z["extinct_at"] < 0 and z["final_counts"] and z["final_counts"][0][0] == "racld":
        pick.append(s)
    if len(pick) == NBG:
        break
INV, bad = {}, []
for s in pick:
    f = os.path.join(BASE, ROOT, "inv", "ib_inv_s%05d_mu0.json" % s)
    if not os.path.exists(f):
        bad.append("없음 %s" % f); continue
    z = json.load(open(f))
    INV[s] = z
    if not z["checksum_match"]: bad.append("배경 체크섬 %d" % s)
    if len(z["arms"]) != 4: bad.append("팔 수 %d" % s)
reg_n = reg_ok = 0
for s, z in INV.items():
    old = {(a["kind"], a["key"], a["post_seed"]): a["series"] for a in json.load(open(os.path.join(BASE, "dms_main", "dms_mat_racld_s%05d_all.json" % s)))["arms"]}
    for a in z["arms"]:
        reg_n += 1
        reg_ok += old.get((a["kind"], a["key"], a["post_seed"])) == a["series"]
MO, mono_missing = {}, []
for c in MONO:
    for s in MSEEDS:
        f = os.path.join(BASE, ROOT, "mono", "ib_mono_%s_s%05d.json" % (c, s))
        if os.path.exists(f):
            MO[(c, s)] = json.load(open(f))
        else:
            mono_missing.append((c, s))
dels = sum(sum(b["copy_del"]) for z in INV.values() for a in z["arms"] for b in a["bins"]) + \
       sum(sum(b["copy_del"]) for z in MO.values() for b in z["bins"])


def tot(bins, L, keys, lim=None):
    bb = bins if lim is None else bins[:lim]
    return {k: sum(b[k][L] for b in bb) for k in keys}


unk = 0
for z in list(INV.values()):
    for a in z["arms"]:
        t = tot(a["bins"], 1, ["unknown", "alloc_dead", "stepped"])
        unk = max(unk, (t["unknown"] + t["alloc_dead"]) / max(t["stepped"], 1))
mext = [k for k, z in MO.items() if z["extinct_at"] >= 0]
print("== 해부 7 ③b μ 0 침입에서 씨앗이 이기는 이유 — 판정")
print("  침입 파일 %d/%d · 단일 파일 %d/%d · 점검 문제 %d · G0 회귀 %d/%d 팔 · 빠뜨림 %d · 표지1 모름 몫 최대 %.2f%% · 단일 멸종 %d"
      % (len(INV), NBG, len(MO), len(MONO) * len(MSEEDS), len(bad), reg_ok, reg_n, dels, 100 * unk, len(mext)))
if bad or mono_missing or reg_ok != reg_n or reg_n == 0 or dels or unk >= 0.01:
    print("🚨 점검 실패 — 판정하지 않는다")
    for b in (bad + [str(m) for m in mono_missing])[:10]:
        print("    ", b)
    sys.exit(1)
bgs = sorted(INV)
KEYS = INV[bgs[0]]["keys"]


def arm(z, key):
    return next(a for a in z["arms"] if a["key"] == key)


out = {"G1": {}, "G2": {}, "G3": {}}
print("\n-- G1 계통 성장의 출처 (X 팔 표지1 − ref 팔 표지1 · 공통 칸 · 3,000틱 적분)")
for c in CODES:
    db, dd, parts = [], [], {"death_age": [], "death_prestep": [], "death_other": []}
    rest = []
    for s in bgs:
        z = INV[s]; x, r = arm(z, c), arm(z, "racld")
        lim = min(len(x["bins"]), len(r["bins"]))
        per_tick = lambda a, k: sum(b[k][1] / max(b["stepped"][1] + b["not_stepped"][1], 1) * z["bin"] for b in a["bins"][:lim])
        db.append(per_tick(x, "births") - per_tick(r, "births"))
        dd.append(-(per_tick(x, "deaths") - per_tick(r, "deaths")))
        for k in parts:
            parts[k].append(-(per_tick(x, k) - per_tick(r, k)))
        pt0 = lambda a, k: sum(b[k][0] / max(b["stepped"][0] + b["not_stepped"][0], 1) * z["bin"] for b in a["bins"][:lim])
        rest.append((pt0(x, "births") - pt0(x, "deaths")) - (pt0(r, "births") - pt0(r, "deaths")))
    B, Dd = boot(db), boot(dd)
    main = "주로 출생 쪽" if abs(B[0]) > abs(Dd[0]) else "주로 사망 쪽"
    out["G1"][c] = {"d_birth": B[:3], "d_death": Dd[:3], "main": main, "death_parts": {k: boot(v)[:3] for k, v in parts.items()}}
    print("   %-7s Δ출생 %+.3f [%+.3f, %+.3f] · Δ사망 %+.3f [%+.3f, %+.3f] → **%s**" % (c, *B[:3], *Dd[:3], main))
    Rr = boot(rest)
    out["G1"][c]["rest_net"] = Rr[:3]
    print("           표지0(나머지) 순성장 차 %+.3f [%+.3f, %+.3f] (서술 — Δ ≈ 표지1 쪽 − 표지0 쪽)" % Rr[:3])
    print("           사망 쪽 원인: " + " · ".join("%s %+.3f [%+.3f, %+.3f]" % (k.replace("death_", ""), *out["G1"][c]["death_parts"][k]) for k in parts))


def tpl(bins, L, lo=None):
    bb = [b for b in bins if lo is None or b["t"] >= lo]
    ok = sum(b["copy_ok"][L] for b in bb); ph = sum(b["copy_phase"][L] for b in bb)
    return ph / ok if ok else float("nan")


print("\n-- G2 가로채기 — 글자당 시간 비 (섞인 접시 대 단일 접시)")
mono_ln = {}
for c in ("rascld", "rsacld"):
    v = [math.log(tpl(MO[(c, s)]["bins"], 0, 15000) / tpl(MO[("racld", s)]["bins"], 0, 15000))
         for s in MSEEDS if (c, s) not in mext and ("racld", s) not in mext]
    mono_ln[c] = boot(v)
    mix = [math.log(tpl(arm(INV[s], c)["bins"], 1) / tpl(arm(INV[s], "racld")["bins"], 1)) for s in bgs]
    mx = boot(mix)
    diff_bs = [a - b for a, b in zip(mx[3], mono_ln[c][3])]
    dm = mx[0] - mono_ln[c][0]
    lo, hi = pct(diff_bs, 0.025), pct(diff_bs, 0.975)
    v2 = ("섞이면 글자를 더 빨리 얻는다 — 가로채기" if hi < 0 else ("섞이면 글자를 더 늦게 얻는다" if lo > 0 else "섞여도 글자 얻는 속도 비는 같다"))
    out["G2"][c] = {"mix_ratio": [math.exp(x) for x in mx[:3]], "mono_ratio": [math.exp(x) for x in mono_ln[c][:3]], "ln_diff": [dm, lo, hi], "verdict": v2}
    tag = "**" if c == "rascld" else ""
    print("   %-7s 섞인 비 %.3f [%.3f, %.3f] · 단일 비 %.3f [%.3f, %.3f] · ln 차 %+.3f [%+.3f, %+.3f] → %s%s%s%s"
          % (c, *[math.exp(x) for x in mx[:3]], *[math.exp(x) for x in mono_ln[c][:3]], dm, lo, hi, tag, v2, tag,
             "" if c == "rascld" else "  (서술)"))

print("\n-- G3 서술 (배경 평균 · 표지1 / 단일은 15,000틱 이후)")
print("   %-12s %9s %9s %8s %9s %8s %9s %9s" % ("", "출생/틱", "사망/틱", "평균나이", "헛손질/출생", "막힘/출생", "글자당틱", "틱/출생"))
def row(name, bins_list, L, lo=None):
    agg = {k: 0 for k in KEYS}
    for bins in bins_list:
        for b in bins:
            if lo is None or b["t"] >= lo:
                for k in KEYS:
                    agg[k] += b[k][L]
    n = agg["stepped"] + agg["not_stepped"]
    r = {"birth_rate": agg["births"] / n, "death_rate": agg["deaths"] / n, "mean_age": agg["age_sum"] / n,
         "stall_pb": agg["copy_stall"] / agg["births"], "block_pb": agg["alloc_blocked"] / agg["births"],
         "tpl": agg["copy_phase"] / agg["copy_ok"], "tpb": n / agg["births"]}
    out["G3"][name] = r
    print("   %-12s %9.5f %9.5f %8.0f %9.1f %8.1f %9.1f %9.1f" % (name, r["birth_rate"], r["death_rate"], r["mean_age"], r["stall_pb"], r["block_pb"], r["tpl"], r["tpb"]))
for c in ["racld"] + CODES:
    row("침입 " + c, [arm(INV[s], c)["bins"] for s in bgs], 1)
for c in MONO:
    row("단일 " + c, [MO[(c, s)]["bins"] for s in MSEEDS if (c, s) not in mext], 0, 15000)

json.dump(out, open(os.path.join(BASE, os.environ.get("PETRI_OUT", "_RESULT_invbud.json")), "w"), ensure_ascii=False, indent=1)
print("\n저장: " + os.environ.get("PETRI_OUT", "_RESULT_invbud.json"))
