"""이긴 코드 해부(가상 DMS) 판정 — 설계 노트 해부-DMS-설계-2026-09-15.md §3~§6 · §7 개정(09-15 파일럿 후, 돌연변이 팔 결과 0건 시점).

  팔의 퍼짐 s = ln((m_T + 0.5) / m_0) − ln((r_T + 0.5) / r_0)   m = 표지 계통 · r = 나머지 · T = 마지막 기록
  기준 = 같은 배경의 WT 팔 11개(ref + neutral 10) s 평균 · 가짜 돌연변이 neutral_k 는 자기를 뺀 10개 평균
  Δ = s − 기준 · 코드마다 배경 평균 · 배경 재추출 부트스트랩 2,000회 95% 구간
  해로움 = 상한 < 0 · 이로움 = 하한 > 0 · 구별 안 됨 = 0 을 품음 · 치명 = 쓴 배경 전부에서 표지 계통 소멸
  점검 ② 체크섬 불일치가 하나라도 있으면 전부 중단 · ① 복제 충실도 · ③ 재료 보존 · 배경 규칙(끝 우세 = WT) 실패 배경은 빼고 보고
  점검 ④ 가짜 돌연변이 10개 중 3개 이상이 해로움 · 이로움 → 분류 보류 · ⑤ d 빠뜨림이 해로움 아님 → 분류 보류
실행: python3 dms_analyze.py [petri 폴더, 기본 ~/petri] [main | pilot | mu1]
  mu1 = 해부 2A(팔 μ 1%, 폴더 dms_mu1) — main 과 같은 판정선 · 출력 _RESULT_dms_mu1.json (09-15 추가, main 출력은 그대로)
"""
import glob
import json
import math
import os
import random
import sys
from collections import Counter

BASE = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1 else "~/petri")
MODE = sys.argv[2] if len(sys.argv) > 2 else "main"
LET = "nsracldjehx"
WTS = [("space", "reasccld"), ("mat", "racld"), ("energy", "reascld")]
FULL = MODE != "pilot"
MU_EXPECT = 0.01 if MODE == "mu1" else 0.0
PLANNED = {"space": 20, "mat": 20, "energy": 60} if FULL else {"space": 1, "mat": 1, "energy": 1}
N_BOOT = 2000
RNG = random.Random(20260915)


def pct(xs, p):
    xs = sorted(xs)
    k = (len(xs) - 1) * p
    f = math.floor(k)
    c = min(f + 1, len(xs) - 1)
    return xs[f] + (xs[c] - xs[f]) * (k - f)


def s_of(arm):
    r = arm["series"]
    n0, m0 = r[0][1], r[0][2]
    nT, mT = r[-1][1], r[-1][2]
    return math.log((mT + 0.5) / m0) - math.log((nT - mT + 0.5) / (n0 - m0))


def sd(xs):
    mu = sum(xs) / len(xs)
    return math.sqrt(sum((x - mu) ** 2 for x in xs) / (len(xs) - 1))


out = {"mode": MODE, "wts": {}}
print(f"모드 {MODE} · 폴더 dms_{MODE}")
halt = False
for cond, wt in WTS:
    files = sorted(glob.glob(os.path.join(BASE, f"dms_{MODE}", f"dms_{cond}_{wt}_s*_{'all' if FULL else 'pilot'}.json")))
    print(f"\n==================== {cond} · WT {wt} · 배경 파일 {len(files)}/{PLANNED[cond]}")
    if len(files) != PLANNED[cond]:
        print("  🚨 계획한 배경이 다 모이지 않았다 — 판정하지 않는다")
        halt = True
        continue
    bgs, dropped = [], []
    for f in files:
        z = json.load(open(f))
        n_expect = 11 + (z["n_mutants_total"] if FULL else 1)
        why = []
        if abs(z.get("mu_assay", 0.0) - MU_EXPECT) > 1e-12:
            print(f"  🚨 팔 μ {z.get('mu_assay', 0.0)} ≠ 모드 {MODE} 의 {MU_EXPECT} — 시드 {z['seed']} · 다른 실험 파일이 섞였다 → 전부 중단")
            halt = True
        if not z["checksum_match"]:
            print(f"  🚨 점검 ② 체크섬 불일치 — 시드 {z['seed']} ({z['grow_checksum']} ≠ {z['main_checksum']}) · v0.3.1 이 궤적을 바꿨다 → 전부 중단")
            halt = True
        if not z["fidelity"]["same"]:
            why.append("① 복제 충실도")
        if not z["grow_top"] or z["grow_top"][0][0] != wt:
            why.append("배경 규칙")
        if len(z["arms"]) != n_expect:
            why.append(f"팔 수 {len(z['arms'])}≠{n_expect}")
        bad_cons = [a["key"] for a in z["arms"] if not a["cons_ok"]]
        if bad_cons:
            why.append(f"③ 재료 보존 {len(bad_cons)}팔")
        if why:
            dropped.append((z["seed"], why))
        else:
            bgs.append(z)
    print(f"  쓴 배경 {len(bgs)} · 뺀 배경 {len(dropped)} " + (" · ".join(f"시드 {s}: {'/'.join(w)}" for s, w in dropped) if dropped else ""))
    if halt or not bgs:
        continue

    n_mut_total = bgs[0]["n_mutants_total"]
    wt_sd = [sd([s_of(a) for a in z["arms"] if a["kind"] in ("ref", "neutral")]) for z in bgs]
    print(f"  WT 팔 11개의 s 표준편차(배경별): 중앙 {pct(wt_sd, .5):.3f} · 범위 {min(wt_sd):.3f}~{max(wt_sd):.3f}")
    grow_s = [z["grow_ms"] / 1000 for z in bgs]
    arm_s = [a["elapsed_ms"] / 1000 for z in bgs for a in z["arms"]]
    print(f"  시간: 키우기 중앙 {pct(grow_s, .5):.0f}s · 팔 하나 중앙 {pct(arm_s, .5):.1f}s · 끼운 개체 중앙 {pct([z['arms'][0]['inj']['done'] for z in bgs], .5):.0f}")

    deltas = {}
    labels_of = {}
    extinct = Counter()
    for z in bgs:
        wt_arms = [a for a in z["arms"] if a["kind"] in ("ref", "neutral")]
        s_wt = {a["labels"][0]: s_of(a) for a in wt_arms}
        base_all = sum(s_wt.values()) / len(s_wt)
        for a in z["arms"]:
            if a["kind"] == "ref":
                continue
            if a["kind"] == "neutral":
                lab = a["labels"][0]
                others = [v for k, v in s_wt.items() if k != lab]
                d = s_wt[lab] - sum(others) / len(others)
                code = lab
            else:
                d = s_of(a) - base_all
                code = a["key"]
                labels_of[code] = a["labels"]
            deltas.setdefault(code, []).append(d)
            if a["series"][-1][2] == 0:
                extinct[code] += 1

    nb = len(bgs)
    boot_idx = [[RNG.randrange(nb) for _ in range(nb)] for _ in range(N_BOOT)]

    def summarize(code):
        ds = deltas[code]
        if len(ds) != nb:
            return None
        mean = sum(ds) / nb
        if nb >= 2:
            bm = [sum(ds[i] for i in idx) / nb for idx in boot_idx]
            lo, hi = pct(bm, .025), pct(bm, .975)
        else:
            lo = hi = float("nan")
        if lo != lo:
            cls = "—"
        elif hi < 0:
            cls = "해로움"
        elif lo > 0:
            cls = "이로움"
        else:
            cls = "구별 안 됨"
        return {"code": code, "labels": labels_of.get(code, [code]), "n": nb, "mean": mean, "lo": lo, "hi": hi, "cls": cls,
                "extinct": extinct[code], "lethal": extinct[code] == nb}

    res = {code: summarize(code) for code in deltas}
    res = {k: v for k, v in res.items() if v is not None}

    neut = [res[f"neutral:{k}"] for k in range(1, 11)]
    print("  점검 ④ 가짜 돌연변이 10개: " + " · ".join(f"{r['mean']:+.2f}[{r['lo']:+.2f},{r['hi']:+.2f}]{'' if r['cls'] in ('구별 안 됨', '—') else ' ' + r['cls']}" for r in neut))
    n_false = sum(r["cls"] in ("해로움", "이로움") for r in neut)
    pos_label = f"del:{wt.rfind('d')}"
    pos = next(r for r in res.values() if pos_label in r["labels"])
    print(f"  점검 ⑤ 양성 대조 {pos_label}({pos['code']}): Δ {pos['mean']:+.2f} [{pos['lo']:+.2f}, {pos['hi']:+.2f}] · {pos['cls']} · 소멸 {pos['extinct']}/{nb}")

    if not FULL:
        ds = [d for k in (f"neutral:{i}" for i in range(1, 11)) for d in deltas[k]]
        print(f"  (파일럿) 가짜 돌연변이 Δ 범위 {min(ds):+.3f} ~ {max(ds):+.3f} · 분류는 배경 1개라 하지 않는다")
        out["wts"][wt] = {"cond": cond, "n_bg": nb, "wt_sd": wt_sd, "neutral": neut, "positive": pos}
        continue

    hold = []
    if n_false >= 3:
        hold.append(f"④ 가짜 돌연변이 {n_false}/10 이 효과 있음으로 찍힘")
    print(f"  → 점검 ④ {'🚨 발동 — ' + hold[-1] if n_false >= 3 else f'통과 ({n_false}/10)'}")
    muts = [r for k, r in res.items() if not k.startswith("neutral:")]
    if len(muts) != n_mut_total:
        print(f"  🚨 돌연변이 코드 수 {len(muts)} ≠ {n_mut_total}")
    cnt = Counter(r["cls"] for r in muts)
    lethal = [r for r in muts if r["lethal"]]
    print(f"  돌연변이 {len(muts)}개: 해로움 {cnt['해로움']} (그중 치명 {len(lethal)}) · 구별 안 됨 {cnt['구별 안 됨']} · 이로움 {cnt['이로움']}")
    out["wts"][wt] = {"cond": cond, "n_bg": nb, "dropped": dropped, "wt_sd": wt_sd, "neutral": neut, "n_false": n_false,
                      "positive": pos, "mutants": sorted(muts, key=lambda r: r["mean"]), "hold": hold}

    print("  자리별 — 무늬 바뀜(기능만 지움) · 빠뜨림(기능+길이) · 바뀜 10개 중 해로움 · 끼어듦 11개 중 해로움/이로움")
    by_label = {lab: r for r in muts for lab in r["labels"]}
    for p, ch in enumerate(wt):
        xs = by_label[f"sub:{p}:x"]
        dl = by_label[f"del:{p}"]
        subs = [by_label[f"sub:{p}:{c}"] for c in LET if c != ch]
        role = "핵심" if xs["cls"] == "해로움" else "군더더기 후보"
        print(f"    {p} {ch}  무늬 {xs['mean']:+.2f} [{xs['lo']:+.2f},{xs['hi']:+.2f}] {xs['cls']:<5} · 빠뜨림 {dl['mean']:+.2f} {dl['cls']:<5}{' 치명' if dl['lethal'] else ''} · "
              f"바뀜 해로움 {sum(r['cls'] == '해로움' for r in subs)}/10 · → {role}")
    for p in range(len(wt) + 1):
        ins = [by_label[f"ins:{p}:{c}"] for c in LET]
        print(f"    끼어듦 자리 {p}: 해로움 {sum(r['cls'] == '해로움' for r in ins)}/11 · 이로움 {sum(r['cls'] == '이로움' for r in ins)}/11")
    good = sorted((r for r in muts if r["cls"] == "이로움"), key=lambda r: -r["mean"])
    print("  이로움 전부: " + (" · ".join(f"{r['code']}({'/'.join(r['labels'])}) {r['mean']:+.2f}[{r['lo']:+.2f},{r['hi']:+.2f}]" for r in good) if good else "없음"))
    print(f"  🟡 여러 번 검정: 진짜 중립이면 약 5% 가 우연히 찍힌다 — 돌연변이 {len(muts)}개면 약 {0.05 * len(muts):.0f}개")

if FULL and not halt:
    pos_ok = all(out["wts"][wt]["positive"]["cls"] == "해로움" for _, wt in WTS if wt in out["wts"])
    print(f"\n점검 ⑤ 양성 대조 세 WT 모두 해로움: {'통과' if pos_ok else '🚨 발동 — 분류 보류'}")
    out["positive_all_harmful"] = pos_ok
    json.dump(out, open(os.path.join(BASE, f"_RESULT_dms_{MODE}.json"), "w"), ensure_ascii=False, indent=1)
if halt:
    print("\n🚨 판정 중단 사유가 있다 — 위를 볼 것")
    sys.exit(1)
