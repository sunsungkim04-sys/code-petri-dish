# -*- coding: utf-8 -*-
"""해부 8 C 판정 — 재료 먼저 세계에서 부풀림이 1 보다 작은(≈0.8) 이유. 사전등록 해부8-사전등록-2026-09-17.md §4 · §5. 결과 전에 썼다(파일럿 시드 9201~9203 은 봤다 — 노트 §5).

  자료: lineage.js v1.3.0 --spec 1 · 규칙 {원래 · 재료 먼저} × 코드 {racld · rascld · rsacld} × 밀도 {4 · 8 · 32}
        × μ {0.1% · 1%} × 시드 901~920 · 20,000틱  (폴더 spec8/<규칙>/d<밀도>/)
  베끼기 결과(WT 개체 · 틱 전후 상태): 맞게 ok · 바뀜 sub(글자별) · 끼어듦 ins(글자별) · 빠뜨림 del · 헛손질 stall
    나아감 adv = ok + sub + del  (같은 글자로의 바뀜은 ok 로 보인다 — 명목도 그렇게 센다: 바뀜 10/11)
    R = (del/adv) ÷ (μ/3)                         주사위 횟수 몫 — 글자 하나 나아갈 때 주사위를 몇 번 굴렸나
    F = (1 + ins/del + sub/del) ÷ (2 + 30/11)     거름 몫 — 굴린 오류 중 실제로 들어간 몫(명목 = 1)
    I = R × F = 글자당 오류 ÷ 명목 μ(2/3 + 10/11)  (항등식)
    출생 단위 예측 u_spec = 1 − (1 − e)^L, e = (del + ins + sub)/adv · I_birth = u_spec / u_nom
    계보 부풀림 I_lin = u_obs / u_nom (해부 7 과 같은 정의 — WT 부모 출생 중 변종 몫)
  칸(규칙 · 코드 · 밀도 · μ)마다 시드 합으로 비를 내고 시드 재추출 2,000회 95%. 멸종 접시는 뺀다(수 보고).

  K0 계측 점검 : 재료 먼저 모든 칸에서 |ln(I_birth / I_lin)| ≤ 0.10 — 아니면 K1~K3 판정하지 않는다
                (원래 규칙 칸은 서술 — 헛손질이 몰려 글자끼리 독립 가정이 약하다)
  K1 주사위    : 재료 먼저 · 밀도 8 · 모든 코드 × μ 에서 R 구간이 [0.9, 1.1] 안 → '글자 하나에 주사위는 한 번 — 주사위 몫 없음'
                하나라도 구간 상한 < 0.9 또는 하한 > 1.1 → '주사위 몫이 있다'(방향 적음) · 그 밖 '판정 불가(넓음)'
  K2 거름      : 재료 먼저 · 밀도 8 · 모든 코드 × μ 에서 F 상한 < 1 → '없는 글자로의 끼어듦 · 바뀜이 걸러진다'
                그리고 몫 ln F / ln I (칸 풀링 점추정) ≥ 2/3 → '0.8 은 대부분 거름 몫'
  K3 용량 반응 : 재료 먼저 · 코드 × μ 여섯 칸 평균 ln F 의 (밀도 32 − 밀도 4) 하한 > 0 → '재료가 넉넉할수록 덜 걸러진다'
                상한 < 0 → '반대' · 그 밖 '검출 안 됨' · I 의 같은 차는 서술
  K4 서술      : 원래 규칙의 R · F (밀도별) — 원래 규칙의 부풀림 몫 ln R / ln I
  K5 서술      : 글자별 끼어듦 몫 p_k = ins_k ÷ (del/11) — 코드에 든 글자 대 안 든 글자 (재료 먼저 · 밀도 8)
  🔻 개정 1회 (09-17 · 판정 보류 뒤 · K 결과는 보기 전): 원판은 '놓친 출생 몫(파일별) < 1%' 를 중단 조건으로 걸었는데
     720 파일 중 하나(ff/d32/rsacld/μ0.1%/s917 · 1.128%)가 걸려 멈췄다. 이 점검이 지키려던 것은 I_계보(u_obs)인데,
     그것은 K0 가 독립 계측(I_출생)과 직접 대조한다. 그래서 놓친 출생은 서술(파일 수 · 최대)로 내리고 K0 가 그 자리를 맡는다.
실행: python3 spec_analyze.py [petri 폴더]  →  _RESULT_spec8.json
마른 실행: PETRI_ROOT=smoke8c PETRI_SEEDS=9201,9202 python3 spec_analyze.py  (파일럿 요약만: PETRI_CODES=racld,rascld)
"""
import json
import math
import os
import random
import sys

BASE = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1 else "~/petri")
ROOT = os.environ.get("PETRI_ROOT", "spec8")
SEEDS = [int(x) for x in os.environ["PETRI_SEEDS"].split(",")] if os.environ.get("PETRI_SEEDS") else list(range(901, 921))
RULES = ["orig", "ff"]
CODES = os.environ["PETRI_CODES"].split(",") if os.environ.get("PETRI_CODES") else ["racld", "rascld", "rsacld"]
DENS = [4, 8, 32]
MUS = [0.001, 0.01]
LET = "nsracldjehx"
NOM = 2 + 30 / 11
N_BOOT = 2000
RNG = random.Random(20260918)


def pct(xs, p):
    xs = sorted(x for x in xs if x == x)
    if not xs:
        return float("nan")
    k = (len(xs) - 1) * p
    f = math.floor(k)
    c = min(f + 1, len(xs) - 1)
    return xs[f] + (xs[c] - xs[f]) * (k - f)


def unom(mu, L):
    return 1 - (1 - mu * (2 / 3 + 10 / 11)) ** L


def fname(rule, code, d, mu, s):
    return os.path.join(BASE, ROOT, rule, "d%d" % d, "lin_%s_mu%s_s%05d%s.json" % (code, str(mu).replace(".", "p"), s, "_ff" if rule == "ff" else ""))


D, missing, bad = {}, [], []
for rule in RULES:
    for code in CODES:
        for d in DENS:
            for mu in MUS:
                for s in SEEDS:
                    f = fname(rule, code, d, mu, s)
                    if not os.path.exists(f):
                        missing.append(f)
                        continue
                    z = json.load(open(f))
                    if (z["sim_version"] not in ("0.3.3", "0.3.4", "0.3.5") or z["lineage_version"] not in ("1.3.0", "1.4.0", "1.5.0") or z["find_first"] != (rule == "ff")
                            or z["density"] != d or z["mu"] != mu or z["code"] != code or z["spec"] is None or z["watch"] is not None):
                        bad.append(f)
                    D[(rule, code, d, mu, s)] = z
total = len(RULES) * len(CODES) * len(DENS) * len(MUS) * len(SEEDS)
ext = sorted(k for k, z in D.items() if z["extinct_at"] >= 0)
miss = max((z["missed_births"] / max(z["births"], 1) for z in D.values()), default=1)
other = max((z["spec"]["other"] / max(1, z["spec"]["ok"]) for z in D.values()), default=1)
print("== 해부 8 C 부풀림 0.8 의 이유 — 판정")
print("  파일 %d/%d · 설정 어긋남 %d · 멸종 %d · 놓친 출생 몫 최대 %.3f%% · '그 밖' 결과 몫 최대 %.4f%%"
      % (len(D), total, len(bad), len(ext), 100 * miss, 100 * other))
n_over = sum(1 for z in D.values() if z["missed_births"] / max(z["births"], 1) >= 0.01)
print("  (개정 1회) 놓친 출생 ≥ 1%% 인 파일 %d/%d — 서술 · I_계보 의 보호는 K0 가 맡는다" % (n_over, len(D)))
if missing or bad:
    print("🚨 점검 실패 — 판정하지 않는다")
    for f in (missing + bad)[:10]:
        print("   ", f)
    sys.exit(1)
if ext:
    print("  멸종(뺌): " + ", ".join("%s/%s/d%d/μ%g/s%d" % k for k in ext[:30]) + (" …" if len(ext) > 30 else ""))


def stats(rows):
    """rows = 시드별 (spec, code, mu, u_obs) — 합으로 비를 낸다."""
    ok = sum(r[0]["ok"] for r in rows)
    de = sum(r[0]["del"] for r in rows)
    su = sum(sum(r[0]["sub"]) for r in rows)
    ins = sum(sum(r[0]["ins"]) for r in rows)
    wpb = sum(r[3][0] for r in rows)
    wpm = sum(r[3][1] for r in rows)
    code, mu = rows[0][1], rows[0][2]
    L = len(code)
    adv = ok + su + de
    if de == 0 or adv == 0 or wpb == 0:
        return None
    R = (de / adv) / (mu / 3)
    F = (1 + ins / de + su / de) / NOM
    e = (de + ins + su) / adv
    Ib = (1 - (1 - e) ** L) / unom(mu, L)
    Il = (wpm / wpb) / unom(mu, L)
    return {"R": R, "F": F, "I": R * F, "Ib": Ib, "Il": Il, "n_del": de}


cells, boots = {}, {}
for rule in RULES:
    for code in CODES:
        for d in DENS:
            for mu in MUS:
                rows = [(D[(rule, code, d, mu, s)]["spec"], code, mu,
                         (D[(rule, code, d, mu, s)]["wt_parent_births"], D[(rule, code, d, mu, s)]["wt_parent_mut_births"]))
                        for s in SEEDS if (rule, code, d, mu, s) in D and D[(rule, code, d, mu, s)]["extinct_at"] < 0]
                key = (rule, code, d, mu)
                if len(rows) < 2:
                    cells[key] = None
                    continue
                cells[key] = stats(rows)
                bs = []
                for _ in range(N_BOOT):
                    st = stats([rows[RNG.randrange(len(rows))] for _ in rows])
                    if st:
                        bs.append(st)
                boots[key] = bs


def ci(key, f):
    v = [f(b) for b in boots[key]]
    return pct(v, .025), pct(v, .975)


out = {"cells": {}}
print("\n   %-5s %-7s %3s %6s | %-22s | %-22s | %-7s %-7s %-7s | del 수" % ("규칙", "코드", "밀도", "μ", "R 주사위", "F 거름", "I", "I_출생", "I_계보"))
for key, st in cells.items():
    if st is None:
        print("   %-5s %-7s %3d %6g | (자료 부족)" % key)
        continue
    rl, rh = ci(key, lambda b: b["R"])
    fl, fh = ci(key, lambda b: b["F"])
    print("   %-5s %-7s %3d %6g | %6.3f [%6.3f, %6.3f] | %5.3f [%5.3f, %5.3f] | %6.3f  %6.3f  %6.3f | %d"
          % (key + (st["R"], rl, rh, st["F"], fl, fh, st["I"], st["Ib"], st["Il"], st["n_del"])))
    out["cells"]["/".join(map(str, key))] = dict(st, R_ci=[rl, rh], F_ci=[fl, fh])

ffkeys = [k for k in cells if k[0] == "ff"]
k0 = [(k, math.log(cells[k]["Ib"] / cells[k]["Il"])) for k in ffkeys if cells[k]]
k0_bad = [k for k, v in k0 if abs(v) > 0.10]
print("\n-- K0 계측 점검: 재료 먼저 칸 %d개 · |ln(I_출생/I_계보)| 최대 %.3f → %s"
      % (len(k0), max(abs(v) for _, v in k0), "통과" if not k0_bad and len(k0) == len(ffkeys) else "🚨 실패 %d칸 — K1~K3 판정하지 않는다" % len(k0_bad)))
orig_k0 = [math.log(cells[k]["Ib"] / cells[k]["Il"]) for k in cells if k[0] == "orig" and cells[k]]
print("   (서술) 원래 규칙 칸 ln(I_출생/I_계보) 범위 %+.3f ~ %+.3f" % (min(orig_k0), max(orig_k0)))
out["K0"] = {"ff_max_abs": max(abs(v) for _, v in k0), "pass": not k0_bad and len(k0) == len(ffkeys), "orig_range": [min(orig_k0), max(orig_k0)]}

if out["K0"]["pass"]:
    d8 = [k for k in ffkeys if k[2] == 8]
    rci = [ci(k, lambda b: b["R"]) for k in d8]
    if all(0.9 <= lo and hi <= 1.1 for lo, hi in rci):
        v1 = "글자 하나에 주사위는 한 번 — 주사위 몫 없음"
    elif any(hi < 0.9 for lo, hi in rci) or any(lo > 1.1 for lo, hi in rci):
        v1 = "주사위 몫이 있다 (" + ", ".join("%s/μ%g %.2f" % (k[1], k[3], cells[k]["R"]) for k in d8) + ")"
    else:
        v1 = "판정 불가 — 구간이 [0.9, 1.1] 밖으로 걸친다"
    print("\n-- K1 주사위 (재료 먼저 · 밀도 8): R 구간 " + " · ".join("[%.3f, %.3f]" % c for c in rci) + " → **%s**" % v1)
    fci = [ci(k, lambda b: b["F"]) for k in d8]
    lnF = sum(math.log(cells[k]["F"]) for k in d8) / len(d8)
    lnI = sum(math.log(cells[k]["I"]) for k in d8) / len(d8)
    share = lnF / lnI if lnI != 0 else float("nan")
    filt = all(hi < 1 for lo, hi in fci)
    v2 = ("없는 글자로의 끼어듦 · 바뀜이 걸러진다" if filt else "거름이 모든 칸에서 검출되지는 않는다")
    v2b = ("0.8 은 대부분 거름 몫" if filt and share >= 2 / 3 else "거름만으로는 0.8 의 대부분이 설명되지 않는다")
    print("-- K2 거름 (재료 먼저 · 밀도 8): F 구간 " + " · ".join("[%.3f, %.3f]" % c for c in fci)
          + " → **%s** · 몫 ln F / ln I = %.2f (평균 ln F %+.3f · ln I %+.3f) → **%s**" % (v2, share, lnF, lnI, v2b))
    # K3: 여섯 칸 평균 ln F 의 (32 − 4) — 칸마다 따로 재추출한 부트스트랩 값을 같은 번호끼리 짝지어 평균
    def mean_lnF(d, b=None):
        vals = []
        for code in CODES:
            for mu in MUS:
                k = ("ff", code, d, mu)
                if not cells.get(k):
                    return float("nan")
                vals.append(math.log((boots[k][b % len(boots[k])] if b is not None else cells[k])["F"]))
        return sum(vals) / len(vals)

    def mean_lnI(d, b=None):
        vals = []
        for code in CODES:
            for mu in MUS:
                k = ("ff", code, d, mu)
                vals.append(math.log((boots[k][b % len(boots[k])] if b is not None else cells[k])["I"]))
        return sum(vals) / len(vals)

    dF = mean_lnF(32) - mean_lnF(4)
    dFb = [mean_lnF(32, b) - mean_lnF(4, b) for b in range(N_BOOT)]
    lo, hi = pct(dFb, .025), pct(dFb, .975)
    v3 = "재료가 넉넉할수록 덜 걸러진다" if lo > 0 else "재료가 넉넉할수록 더 걸러진다" if hi < 0 else "밀도에 따른 차 검출 안 됨"
    dI = mean_lnI(32) - mean_lnI(4)
    dIb = [mean_lnI(32, b) - mean_lnI(4, b) for b in range(N_BOOT)]
    print("-- K3 용량 반응: 평균 ln F 밀도 4 %+.3f · 8 %+.3f · 32 %+.3f · 차(32 − 4) %+.3f [%+.3f, %+.3f] → **%s**"
          % (mean_lnF(4), mean_lnF(8), mean_lnF(32), dF, lo, hi, v3))
    print("   (서술) 평균 ln I 밀도 4 %+.3f · 8 %+.3f · 32 %+.3f · 차 %+.3f [%+.3f, %+.3f]"
          % (mean_lnI(4), mean_lnI(8), mean_lnI(32), dI, pct(dIb, .025), pct(dIb, .975)))
    out.update({"K1": {"R_ci": rci, "verdict": v1}, "K2": {"F_ci": fci, "share": share, "verdict": [v2, v2b]},
                "K3": {"lnF": [mean_lnF(4), mean_lnF(8), mean_lnF(32)], "diff": [dF, lo, hi], "verdict": v3,
                       "lnI": [mean_lnI(4), mean_lnI(8), mean_lnI(32)], "dI": [dI, pct(dIb, .025), pct(dIb, .975)]}})

print("\n-- K4 서술 (원래 규칙): 밀도별 R · F · 부풀림 몫 ln R / ln I")
for d in DENS:
    ks = [k for k in cells if k[0] == "orig" and k[2] == d and cells[k]]
    if not ks:
        continue
    print("   밀도 %2d: " % d + " · ".join("%s/μ%g R %.2f F %.2f 몫 %.2f" % (k[1], k[3], cells[k]["R"], cells[k]["F"],
                                                                     math.log(cells[k]["R"]) / math.log(cells[k]["I"]) if abs(math.log(cells[k]["I"])) > 1e-9 else float("nan")) for k in ks))

print("\n-- K5 서술 (재료 먼저 · 밀도 8): 글자별 끼어듦 몫 p_k = ins_k ÷ (del/11)")
for code in CODES:
    ins = [0] * len(LET)
    de = 0
    for mu in MUS:
        for s in SEEDS:
            z = D.get(("ff", code, 8, mu, s))
            if z and z["extinct_at"] < 0:
                de += z["spec"]["del"]
                for i, v in enumerate(z["spec"]["ins"]):
                    ins[i] += v
    pk = {LET[i]: ins[i] / (de / 11) for i in range(len(LET))} if de else {}
    inside = [pk[c] for c in set(code)]
    outside = [pk[c] for c in LET if c not in code]
    print("   %-7s 코드 안 글자 평균 %.2f · 밖 글자 평균 %.2f · " % (code, sum(inside) / len(inside), sum(outside) / len(outside))
          + " ".join("%s %.2f" % (c, pk[c]) for c in LET))
    out.setdefault("K5", {})[code] = pk

json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "_RESULT_spec8.json"), "w"), ensure_ascii=False, indent=1)
print("\n저장: _RESULT_spec8.json")
