# -*- coding: utf-8 -*-
"""해부 7 ① 판정 — 시간 예산 부활. 사전등록 해부7-사전등록-2026-09-16.md §1 (판정선 B1~B3 · 띠는 해부6 §3 그대로).

  경위: 원판 budget_analyze.py(시드 601~630)는 '멸종 0' 에서, 개정판 budget_analyze2.py(631~660)는
  '낳은 뒤 멸종 0' 에서 멈췄다 — 둘 다 밀도 4 · 8 의 창업 실패(개체가 1 을 넘은 기록이 없음)였다.
  사용자 결정(09-16)으로 규칙을 한 번 더 고친다. **새 시드 661~690 결과를 보기 전에** 쓰고 해시를 고정했다.

  바뀐 것 하나: **기록된 개체 수(500틱 칸 시작 · 끝)가 한 번도 10 을 넘지 않은 채 멸종한 접시는 창업 실패**로 보고,
  그 칸에서 그 시드를 통째로 빼고 수를 적는다. 기록상 10 을 넘은 뒤 멸종한 접시가 하나라도 있으면 중단한다.
  (칸 사이에 잠깐 10 을 넘었다 사라진 접시는 기록에 안 남으므로 창업 실패로 분류된다 — 미리 적어 둔 한계)
"""
import glob
import json
import math
import os
import random
import sys

BASE = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1 else "~/petri")
ROOT = os.environ.get("PETRI_ROOT", "bud6c")
SEEDS = [int(x) for x in os.environ["PETRI_SEEDS"].split(",")] if os.environ.get("PETRI_SEEDS") else list(range(661, 691))
RUNGS = [("mat_d04", "밀도 4"), ("mat_d08", "밀도 8"), ("mat_d14", "밀도 14"), ("mat_d32", "밀도 32"), ("space", "space")]
CODES = ["racld", "rascld", "rsacld"]
BANDS = {("rsacld", "rascld"): ((0.55, 0.80), (0.90, 1.10)), ("racld", "rascld"): ((0.35, 0.55), (0.75, 0.95))}
N_BOOT = 2000
RNG = random.Random(20261005)
NUM = ["stepped", "not_stepped", "unknown", "copy_ok", "copy_stall", "copy_del", "copy_idle",
       "alloc_ok", "alloc_blocked", "alloc_idle", "alloc_dead", "div_with_child", "div_idle", "births"]


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
    bs = [sum(vals[RNG.randrange(len(vals))] for _ in vals) / len(vals) for _ in range(N_BOOT)]
    return m, pct(bs, 0.025), pct(bs, 0.975)


def fname(rung, code, s):
    return os.path.join(BASE, ROOT, "bud_%s_%s_s%05d.json" % (rung, code, s))


D, missing = {}, []
for rung, _ in RUNGS:
    for code in CODES:
        for s in SEEDS:
            f = fname(rung, code, s)
            if os.path.exists(f):
                D[(rung, code, s)] = json.load(open(f))
            else:
                missing.append((rung, code, s))


def window(z, which):
    t = {k: 0 for k in NUM}
    for b in z["bins"]:
        if (which == "grow" and b["pop0"] < 0.5 * z["pop_end"]) or (which == "steady" and b["t0"] >= 15000):
            for k in NUM:
                t[k] += b[k]
    return t


W = {k: {"grow": window(z, "grow"), "steady": window(z, "steady")} for k, z in D.items()}
space_stall = sum(b["copy_stall"] for k, z in D.items() if k[0] == "space" for b in z["bins"])
dels = sum(b["copy_del"] for z in D.values() for b in z["bins"])
unk = max(((sum(b["unknown"] + b["alloc_dead"] for b in z["bins"])) / max(sum(b["stepped"] for b in z["bins"]), 1)) for z in D.values()) if D else 1
def births_total(z):
    return sum(b["births"] for b in z["bins"])


def max_pop(z):
    return max([b["pop0"] for b in z["bins"]] + [z["pop_end"]])


founder = [k for k, z in D.items() if z["extinct_at"] >= 0 and max_pop(z) <= 10]
ext = [k for k, z in D.items() if z["extinct_at"] >= 0 and max_pop(z) > 10]
RS = {rung: [s for s in SEEDS if not any((rung, c, s) in founder for c in CODES)] for rung, _ in RUNGS}
chk_ok, chk_msg = False, "mono 대조 파일 없음"
mf = sorted(glob.glob(os.path.join(BASE, ROOT + "_mono", "mono_mat_racld_mu0_s%05d.json" % SEEDS[0])))
if mf and ("mat_d08", "racld", SEEDS[0]) in D:
    a, b = json.load(open(mf[0]))["checksum"], D[("mat_d08", "racld", SEEDS[0])]["checksum"]
    chk_ok, chk_msg = a == b, "%s vs %s" % (a, b)
print("== 해부 7 ① 시간 예산 — 판정 (시드 661~690)")
print("  창업 실패(기록상 개체 ≤ 10 인 채 멸종) %d접시 → 뺀 시드: %s" % (len(founder), {n: len(SEEDS) - len(RS[r]) for r, n in RUNGS}))
print("  파일 %d/%d · 10 을 넘은 뒤 멸종 %d · space 헛손질 %d · μ0 빠뜨림 %d · 모름 몫 최대 %.3f%% · budget/mono 체크섬 %s (%s)"
      % (len(D), len(RUNGS) * len(CODES) * len(SEEDS), len(ext), space_stall, dels, 100 * unk, chk_ok, chk_msg))
if missing or space_stall or dels or unk >= 0.01 or not chk_ok or ext:
    print("🚨 점검 실패 — 판정하지 않는다")
    for m in missing[:10]:
        print("    없음", m)
    sys.exit(1)

out = {"table": {}, "B1": {}, "B2": {}, "B3": {}}
print("\n-- 출생당 개체-틱 (시드 합) · 다 찬 뒤 / 채우는 동안")
print("   %-8s %-7s | %8s %8s | %8s %8s %8s %8s | %7s %7s" % ("칸", "코드", "틱/출생", "(채움)", "헛손질", "막힘", "복사성공", "그 밖", "헛손질률", "막힘률"))
for rung, name in RUNGS:
    for code in CODES:
        st = {k: sum(W[(rung, code, s)]["steady"][k] for s in RS[rung]) for k in NUM}
        gr = {k: sum(W[(rung, code, s)]["grow"][k] for s in RS[rung]) for k in NUM}
        B, G = st["births"], max(gr["births"], 1)
        other = st["stepped"] - st["copy_stall"] - st["alloc_blocked"] - st["copy_ok"]
        row = {"tpb": st["stepped"] / B, "tpb_grow": gr["stepped"] / G, "stall_pb": st["copy_stall"] / B,
               "block_pb": st["alloc_blocked"] / B, "ok_pb": st["copy_ok"] / B, "other_pb": other / B,
               "stall_rate": st["copy_stall"] / max(st["copy_stall"] + st["copy_ok"], 1),
               "block_rate": st["alloc_blocked"] / max(st["alloc_blocked"] + st["alloc_ok"], 1)}
        out["table"]["%s|%s" % (rung, code)] = row
        print("   %-8s %-7s | %8.1f %8.1f | %8.1f %8.1f %8.2f %8.1f | %7.3f %7.3f"
              % (name, code, row["tpb"], row["tpb_grow"], row["stall_pb"], row["block_pb"], row["ok_pb"], row["other_pb"],
                 row["stall_rate"], row["block_rate"]))

print("\n-- B1 (②) 다 찬 뒤 무엇이 낭비인가 — 헛손질 몫 − 막힘 몫 (racld · rascld 평균)")
for rung, name in RUNGS:
    v = []
    for s in RS[rung]:
        d = 0
        for code in ("racld", "rascld"):
            t = W[(rung, code, s)]["steady"]
            d += (t["copy_stall"] - t["alloc_blocked"]) / t["stepped"] / 2
        v.append(d)
    m, lo, hi = boot(v)
    verdict = "글자 찾기가 더 큰 낭비" if lo > 0 else ("자리 막힘이 더 큰 낭비" if hi < 0 else "비슷하다")
    out["B1"][name] = {"diff": [m, lo, hi], "verdict": verdict}
    print("   %-8s %+.3f [%+.3f, %+.3f] → **%s**" % (name, m, lo, hi, verdict))


def ratio(rung, a, b, key, which):
    v = []
    for s in RS[rung]:
        x, y = W[(rung, a, s)][which], W[(rung, b, s)][which]
        if key == "stall_pb":
            xa, yb = x["copy_stall"] / max(x["births"], 1), y["copy_stall"] / max(y["births"], 1)
        else:
            xa, yb = x["births"] / max(x["stepped"], 1), y["births"] / max(y["stepped"], 1)
        if xa > 0 and yb > 0:
            v.append(math.log(xa / yb))
    return boot(v)


print("\n-- B2 다시 시도 방식 — 출생당 헛손질 수의 비 (다 찬 뒤)")
for rung, name in (("mat_d04", "밀도 4"), ("mat_d08", "밀도 8")):
    sides = []
    for (a, b), (tb, ab) in BANDS.items():
        m, lo, hi = (math.exp(x) for x in ratio(rung, a, b, "stall_pb", "steady"))
        side = "시간" if tb[0] <= lo and hi <= tb[1] else ("시도" if ab[0] <= lo and hi <= ab[1] else "없음")
        sides.append(side)
        out["B2"]["%s|%s/%s" % (name, a, b)] = {"ratio": [m, lo, hi], "time_band": tb, "try_band": ab, "side": side}
        print("   %-6s %s/%s %.3f [%.3f, %.3f] · 시간 띠 %s · 시도 띠 %s → %s" % (name, a, b, m, lo, hi, tb, ab, side))
    v = "**시간이 묶는다 — 빨리 다시 시도해도 글자는 안 온다**" if sides == ["시간", "시간"] else (
        "**시도가 묶는다 — 한 번 시도할 때마다 같은 확률**" if sides == ["시도", "시도"] else "둘 다 아님")
    out["B2"][name] = v
    print("   %-6s → %s" % (name, v))

print("\n-- B3 채우는 속도의 기제 — 개체-틱당 출생 racld/rascld (채우는 동안)")
for rung, name in RUNGS:
    m, lo, hi = ratio(rung, "racld", "rascld", "birth_rate", "grow")
    v = "racld 가 더 자주 낳는다" if lo > 0 else ("씨앗이 더 자주 낳는다" if hi < 0 else "구별 안 됨")
    out["B3"][name] = {"ln_ratio": [m, lo, hi], "verdict": v}
    star = "  ← 판정 칸" if rung == "mat_d08" else ""
    print("   %-8s %+.4f [%+.4f, %+.4f] → %s%s" % (name, m, lo, hi, v, star))

json.dump(out, open(os.path.join(BASE, os.environ.get("PETRI_OUT", "_RESULT_budget3.json")), "w"), ensure_ascii=False, indent=1)
print("\n저장: " + os.environ.get("PETRI_OUT", "_RESULT_budget3.json"))
