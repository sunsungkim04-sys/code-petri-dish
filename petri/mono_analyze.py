"""해부 4 A · B 판정 — 단일 코드 접시. 사전등록 해부4-사전등록-2026-09-16.md §2 · §3 · §5. 결과 전에 썼다.

  A 성장 사다리: mat · racld 대 rascld · μ {0, 0.0025, 0.005, 0.01, 0.02, 0.04} · 시드 1~30 · 20,000틱
     D(μ) = ln P(10,000) 의 racld − rascld (같은 시드끼리 짝) · 시드 재추출 부트스트랩 2,000회 95%
     하한 > 0 → racld 가 더 잘 자란다 · 상한 < 0 → 씨앗 · 0 품음 → 구별 안 됨
     기울기: D 를 μ 에 최소제곱 회귀한 기울기의 구간 — 하한 > 0 이면 '오류가 많아질수록 racld 가 유리'
  B 자식 부하: 두 WT 의 한 글자 이웃 · μ 0 · 시드 1~3 · 4,000틱
     복제함 = 한 시드라도 끝 개체 ≥ 10 · U_dead = 흐름 가중(sub 1/11 · del 1/3 · ins 1/33) '못 함' 비율
  점검 ① 같은 시드 두 번 체크섬(A 자료 안에서 코드·μ·시드가 같은 파일이 둘이면 비교) ② 재료 보존 위반 0
       ④ μ=0 인데 출발 코드 몫 < 100% 면 중단 ⑤ 두 WT 가 '복제함' 이어야 한다
실행: python3 mono_analyze.py [petri 폴더, 기본 ~/petri]  →  _RESULT_mono.json
"""
import glob
import json
import math
import os
import random
import sys

BASE = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1 else "~/petri")
CODES = ["racld", "rascld"]
MUS = [0.0, 0.0025, 0.005, 0.01, 0.02, 0.04]
SEEDS = list(range(1, 31))
N_BOOT = 2000
RNG = random.Random(20260924)
LET = "nsracldjehx"
W = {"sub": 1 / 11, "del": 1 / 3, "ins": 1 / 33}
ALIVE_MIN = 10


def pct(xs, p):
    xs = sorted(xs)
    k = (len(xs) - 1) * p
    f = math.floor(k)
    c = min(f + 1, len(xs) - 1)
    return xs[f] + (xs[c] - xs[f]) * (k - f)


def load(folder):
    out = {}
    for f in sorted(glob.glob(os.path.join(BASE, folder, "mono_mat_*.json"))):
        z = json.load(open(f))
        out[(z["code"], z["mu"], z["seed"])] = z
    return out


def pop_at(z, t):
    for r in z["samples"]:
        if r[0] == t:
            return r[1], r[3]
    last = z["samples"][-1]
    return (0, 0) if z["extinct_at"] >= 0 else (last[1], last[3])


print("== A 성장 사다리")
A = load("mono_curve")
problems = []
for code in CODES:
    for mu in MUS:
        for s in SEEDS:
            if (code, mu, s) not in A:
                problems.append(f"없음 {code} μ{mu} s{s}")
viol = sum(z["conservation_violations"] for z in A.values())
mu0_bad = [k for k, z in A.items() if k[1] == 0.0 and z["samples"][-1][1] and z["samples"][-1][3] != z["samples"][-1][1]]
print(f"  파일 {len(A)}/{len(CODES) * len(MUS) * len(SEEDS)} · 재료 보존 위반 {viol} · 점검 ④ μ0 에서 출발 코드 아닌 개체가 있는 접시 {len(mu0_bad)}")
if problems or viol or mu0_bad:
    print("🚨 완결성/계측기 점검 실패 — A 를 판정하지 않는다")
    for p in problems[:10]:
        print("   ", p)
    A = None

out = {"A": {}, "B": {}}
if A:
    idxs = [[RNG.randrange(len(SEEDS)) for _ in range(len(SEEDS))] for _ in range(N_BOOT)]
    ext = {c: {mu: 0 for mu in MUS} for c in CODES}
    Dm = {}
    print(f"  {'μ':>7} | {'racld ln P(1만)':>16} | {'rascld':>16} | {'D = racld − rascld':>26} | 판정")
    for mu in MUS:
        d, lp = [], {c: [] for c in CODES}
        for s in SEEDS:
            vals = {}
            for c in CODES:
                z = A[(c, mu, s)]
                p10, _ = pop_at(z, 10000)
                if z["extinct_at"] >= 0:
                    ext[c][mu] += 1
                vals[c] = math.log(p10 + 0.5)
                lp[c].append(vals[c])
            d.append(vals["racld"] - vals["rascld"])
        m = sum(d) / len(d)
        bm = sorted(sum(d[i] for i in ix) / len(d) for ix in idxs)
        lo, hi = pct(bm, .025), pct(bm, .975)
        v = "racld 가 더 잘 자란다" if lo > 0 else ("씨앗이 더 잘 자란다" if hi < 0 else "구별 안 됨")
        Dm[mu] = (m, lo, hi)
        out["A"][str(mu)] = {"D": [m, lo, hi], "verdict": v, "mean_lnP": {c: sum(lp[c]) / len(lp[c]) for c in CODES},
                             "extinct": {c: ext[c][mu] for c in CODES}}
        print(f"  {mu:>7} | {sum(lp['racld']) / len(SEEDS):>16.2f} | {sum(lp['rascld']) / len(SEEDS):>16.2f} | "
              f"{m:>+10.3f} [{lo:+.3f}, {hi:+.3f}] | {v} · 멸종 {ext['racld'][mu]}/{ext['rascld'][mu]}")
    # μ 에 대한 D 의 기울기 (시드 재추출)
    xs = MUS
    def slope(dd):
        mx = sum(xs) / len(xs)
        my = sum(dd) / len(dd)
        num = sum((x - mx) * (y - my) for x, y in zip(xs, dd))
        den = sum((x - mx) ** 2 for x in xs)
        return num / den
    per_seed = {mu: [math.log(pop_at(A[("racld", mu, s)], 10000)[0] + 0.5) - math.log(pop_at(A[("rascld", mu, s)], 10000)[0] + 0.5) for s in SEEDS] for mu in MUS}
    obs = slope([sum(per_seed[mu]) / len(SEEDS) for mu in MUS])
    bs = sorted(slope([sum(per_seed[mu][i] for i in ix) / len(SEEDS) for mu in MUS]) for ix in idxs)
    slo, shi = pct(bs, .025), pct(bs, .975)
    sv = "오류가 많아질수록 racld 가 유리해진다" if slo > 0 else ("반대로 간다" if shi < 0 else "기울기 검출 안 됨")
    out["A"]["slope"] = [obs, slo, shi, sv]
    print(f"  D 의 μ 기울기: {obs:+.2f} [{slo:+.2f}, {shi:+.2f}] → {sv}")
    signs = [Dm[mu][0] > 0 for mu in MUS]
    cross = next((MUS[i] for i in range(1, len(MUS)) if signs[i] != signs[i - 1]), None)
    out["A"]["cross"] = cross
    print(f"  부호가 바뀌는 μ: {cross if cross is not None else '없음(사다리 전 구간 같은 방향)'}")

print("\n== B 자식 부하")
B = load("mono_nb")
codes = sorted({k[0] for k in B})
cls, split = {}, []
for c in codes:
    pops = [B[(c, 0.0, s)]["samples"][-1][1] if (c, 0.0, s) in B else -1 for s in (1, 2, 3)]
    if min(pops) < 0:
        continue
    ok = [p >= ALIVE_MIN for p in pops]
    cls[c] = any(ok)
    if any(ok) and not all(ok):
        split.append((c, pops))
print(f"  이웃 코드 {len(cls)} · 시드에 따라 갈린 코드 {len(split)}" + (" · " + " · ".join(f"{c}{p}" for c, p in split[:8]) if split else ""))
print(f"  점검 ⑤ 두 WT: " + " · ".join(f"{c} {'복제함' if cls.get(c) else '🚨 못 함'}" for c in CODES))


def neighbours(wt):
    out = {}
    def add(key, label):
        if key == wt or not (2 <= len(key) <= 48):
            return
        out.setdefault(key, []).append(label)
    for p in range(len(wt)):
        for ch in LET:
            if ch != wt[p]:
                add(wt[:p] + ch + wt[p + 1:], "sub")
    for p in range(len(wt)):
        add(wt[:p] + wt[p + 1:], "del")
    for p in range(len(wt) + 1):
        for ch in LET:
            add(wt[:p] + ch + wt[p:], "ins")
    return out


for wt in CODES:
    nb = neighbours(wt)
    tot = sum(sum(W[l] for l in labs) for labs in nb.values())
    dead = sum(sum(W[l] for l in labs) for k, labs in nb.items() if not cls.get(k, True))
    n_dead = sum(1 for k in nb if not cls.get(k, True))
    u = dead / tot
    out["B"][wt] = {"n_neighbours": len(nb), "n_dead": n_dead, "U_dead": u, "flux_total": tot}
    print(f"  {wt:>7}: 이웃 {len(nb)} · 혼자서 번식 못 하는 이웃 {n_dead} ({n_dead / len(nb):.0%}) · **U_dead {u:.3f}**")
if all(w in out["B"] for w in CODES):
    a, b = out["B"]["racld"]["U_dead"], out["B"]["rascld"]["U_dead"]
    v = "racld 쪽 자식이 덜 죽는다" if a < b else ("씨앗 쪽 자식이 덜 죽는다" if b < a else "같다")
    if len(split) > 5:
        v += " · 🚨 시드에 따라 갈린 이웃이 5개를 넘어 보류"
    out["B"]["verdict"] = v
    print(f"  → U_dead racld {a:.3f} 대 씨앗 {b:.3f} · 차이 {a - b:+.3f} → {v}")
json.dump(out, open(os.path.join(BASE, "_RESULT_mono.json"), "w"), ensure_ascii=False, indent=1)
print("\n저장: _RESULT_mono.json")
