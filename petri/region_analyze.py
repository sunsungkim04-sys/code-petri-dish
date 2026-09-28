"""해부 4C 판정 — 에너지는 왜 한 코드로 고정되지 않나. 사전등록 해부4-사전등록-2026-09-16.md §4 · §5. 결과 전에 썼다.

  C1 연결된 구름인가 (기존 본실험 기록만): 에너지 생존 접시의 끝 코드 전수에서 개체 0.2% 이상인 코드만 남기고,
     한 글자 편집(바뀜 · 빠뜨림 · 끼어듦)으로 이어지는 코드끼리 선을 그어 성분을 센다.
     최대 성분이 그 코드들의 개체 합에서 70% 이상이면 '하나의 연결된 구름'.
  C2 빛 얼룩이 다른 코드를 기르나 (region.js): 구역(L1 · L2 · 어두운 곳) 사이 Bray–Curtis 중앙값이
     구역 안(무작위 절반 두 벌)의 95 백분위보다 크면 '빛 얼룩이 서로 다른 코드를 기른다'.
     비교는 모두 **각 조성의 상위 40 코드**로 맞춘다(절반 기록이 상위 40까지라 잘림을 같게 두려는 것).
  점검: C2 는 다시 키운 50,000틱 체크섬이 본실험과 같아야 한다 — 하나라도 다르면 C2 중단.
실행: python3 region_analyze.py [petri 폴더, 기본 ~/petri]  →  _RESULT_region.json
"""
import glob
import json
import math
import os
import sys

BASE = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1 else "~/petri")
SHARE_MIN, COMP_MIN, TOPN = 0.002, 0.70, 40


def pct(xs, p):
    xs = sorted(x for x in xs if x == x)
    if not xs:
        return float("nan")
    k = (len(xs) - 1) * p
    f = math.floor(k)
    c = min(f + 1, len(xs) - 1)
    return xs[f] + (xs[c] - xs[f]) * (k - f)


def one_edit(a, b):
    if a == b:
        return False
    if abs(len(a) - len(b)) > 1:
        return False
    if len(a) == len(b):
        return sum(x != y for x, y in zip(a, b)) == 1
    short, long_ = (a, b) if len(a) < len(b) else (b, a)
    return any(short == long_[:i] + long_[i + 1:] for i in range(len(long_)))


def bray(a, b):
    keys = set(a) | set(b)
    den = sum(a.get(k, 0) + b.get(k, 0) for k in keys)
    return sum(abs(a.get(k, 0) - b.get(k, 0)) for k in keys) / den if den else float("nan")


def norm(pairs):
    pairs = pairs[:TOPN]
    tot = sum(c for _, c in pairs)
    return {k: c / tot for k, c in pairs} if tot else {}


print("== C1 연결된 구름인가 (본실험 에너지 생존 접시)")
shares, ncodes, comps = [], [], []
for f in sorted(glob.glob(os.path.join(BASE, "main", "energy_d8_mu0p01_s*.json"))):
    z = json.load(open(f))
    if z["extinct_at"] >= 0:
        continue
    tot = sum(c for _, c in z["final_counts"])
    sel = [(k, c) for k, c in z["final_counts"] if c / tot >= SHARE_MIN]
    if not sel:
        continue
    keys = [k for k, _ in sel]
    cnt = dict(sel)
    parent = {k: k for k in keys}
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for i in range(len(keys)):
        for j in range(i + 1, len(keys)):
            if one_edit(keys[i], keys[j]):
                a, b = find(keys[i]), find(keys[j])
                if a != b:
                    parent[a] = b
    groups = {}
    for k in keys:
        groups.setdefault(find(k), []).append(k)
    sizes = sorted((sum(cnt[k] for k in g) for g in groups.values()), reverse=True)
    sel_tot = sum(cnt.values())
    shares.append(sizes[0] / sel_tot)
    ncodes.append(len(keys))
    comps.append(len(groups))
v1 = "하나의 연결된 구름" if pct(shares, .5) >= COMP_MIN else "떨어진 섬이 여럿"
print(f"  접시 {len(shares)} · 0.2% 이상 코드 수 중앙 {pct(ncodes, .5):.0f} · 성분 수 중앙 {pct(comps, .5):.0f}")
print(f"  최대 성분이 차지하는 개체 비율 중앙 {pct(shares, .5):.1%} [5% {pct(shares, .05):.1%} · 95% {pct(shares, .95):.1%}] → {v1}")

print("\n== C2 빛 얼룩이 다른 코드를 기르나")
files = sorted(glob.glob(os.path.join(BASE, "region", "region_energy_s*.json")))
bad = [json.load(open(f))["seed"] for f in files if not json.load(open(f))["checksum_match"]]
print(f"  접시 {len(files)} · 체크섬 불일치 {len(bad)}")
out = {"C1": {"share_median": pct(shares, .5), "n_dishes": len(shares), "verdict": v1}}
if bad:
    print("🚨 점검 ③ 실패 — C2 중단")
    out["C2"] = {"verdict": "중단 — 체크섬 불일치"}
else:
    between, within, hsh = [], [], {"L1": [], "L2": [], "dark": []}
    for f in files:
        z = json.load(open(f))
        R = z["regions"]
        comp = {k: norm(R[k]["full"]) for k in R}
        for k in R:
            hsh[k].append(R[k]["h_share"])
            within.append(bray(norm(R[k]["half_a"]), norm(R[k]["half_b"])))
        for a, b in (("L1", "L2"), ("L1", "dark"), ("L2", "dark")):
            between.append(bray(comp[a], comp[b]))
    w95 = pct(within, .95)
    bmed = pct(between, .5)
    v2 = "빛 얼룩이 서로 다른 코드를 기른다" if bmed > w95 else "공간 분화 없음"
    print(f"  구역 안 Bray–Curtis 중앙 {pct(within, .5):.3f} · 95% {w95:.3f}")
    print(f"  구역 사이 중앙 {bmed:.3f} [5% {pct(between, .05):.3f} · 95% {pct(between, .95):.3f}] → {v2}")
    print("  구역별 h 든 비율 중앙: " + " · ".join(f"{k} {pct(v, .5):.0%}" for k, v in hsh.items()))
    print("  구역별 개체 수 중앙: " + " · ".join(f"{k} {pct([json.load(open(f))['regions'][k]['n'] for f in files], .5):,.0f}" for k in hsh))
    out["C2"] = {"within_median": pct(within, .5), "within_95": w95, "between_median": bmed,
                 "h_share": {k: pct(v, .5) for k, v in hsh.items()}, "verdict": v2}
json.dump(out, open(os.path.join(BASE, "_RESULT_region.json"), "w"), ensure_ascii=False, indent=1)
print("\n저장: _RESULT_region.json")
