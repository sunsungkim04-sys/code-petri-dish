#!/usr/bin/env python3
"""p4win_analyze.py — 원고 Supp Table S10 · S11 의 P4(복사당 추첨 수 R 의 고리 간 비교)를 같은 시간 창으로 다시 센다
(2026-10-07 · 사후 · 판정 아님 — 사전등록 판정 _RESULT_mini26.txt · _RESULT_mini27.txt 는 그대로다)

왜
  동결 판정기는 mut 팔 안의 주입 L2(`1_2`)와 상주 L4(`0_4`)의 R 을 assay 3,000틱 합계(draws ÷ positions)로 비교한다
  (mini26_analyze.py:442–456 · mini27_analyze.py:420–426 · 607–636). 주입 L2 는 퍼지면서 복사가 assay 뒤로 몰리고
  상주 L4 는 앞으로 몰린다(해부26 census: L2 가임 개체-틱 앞·중·뒤 3분의 1 = 0.165 · 0.307 · 0.529 · L4 0.396 · 0.342 · 0.262).
  R 이 assay 동안 움직이면 합계 비교는 두 계통에 서로 다른 시간대를 싣는다(갈래 A 판정 노트 §5.2 — 다른 배경에서 본 것).
  여기서는 원고 배경(해부26 ρ8 26001–26020 · 해부27 ρ8 27001–27020)을 창별 계수기 래퍼(galA_window.py v1.0.0 · 무수정)로
  다시 돌린 출력에서 같은 창 값을 낸다. 실행은 p4win_launch.sh · p4win_job.sh.

쓰는 법 (petri/ 안에서)
  python3 p4win_analyze.py selftest GALA_STD_DIR          비대상 시험: 갈래 A 출력(ρ8 · 배경 41001–41020)으로 §5.2 수치를 다시 낸다
  python3 p4win_analyze.py ident RUN_DIR REF_DIR STUDY     가드 G1 · G2 만 — 파일럿 시드 시운전용 · STUDY = 26 | 27
  python3 p4win_analyze.py judge RUN26 RUN27 [--json OUT]   본 계산: 가드 G1–G4 → 같은 창 P4 (해부26 · 해부27)

정의 (결과 보기 전에 고정 — 이 파일의 sha256 을 p4win_frozen.sha256 에 적은 뒤 본 실행)
  창       250틱 = galA_window 의 50틱 칸 5 개(홀짝 합침) → 12 창. 칸 b = (t − 1) // 50 · 창 w = b // 5.
  R_w      계통별 draws_w / positions_w.
  창 자격   두 계통 모두 positions_w ≥ 20. 빠진 창은 세어 인쇄하고, 한 양에서 빠진 창이 전체의 5% 를 넘으면 그 양은 판정 불가.
  1차      같은 창 · 창마다 같은 무게: 배경 값 = 자격 창에 걸친 [ln R2_w − ln R4_w] 의 단순 평균.
  P4b 창    배경 값 = 자격 창(두 R_w > 1)에 걸친 [ln((R2_w − 1)·2) − ln((R4_w − 1)·4)] 의 단순 평균.
  2차(서술) L2 의 시간 무게로 맞춤: ln R2 − ln R4⁽²⁾, R4⁽²⁾ = Σ_w pos2_w · R4_w / Σ_w pos2_w (자격 창 · R2 도 자격 창 합).
  번인(서술) 앞 4 창(t ≤ 1000)을 뺀 1차 값 · 점추정만(구간으로 다시 시험하지 않는다).
  추세(서술) 계통별 R 을 앞·중·뒤 3분의 1(창 0–3 · 4–7 · 8–11)로 — 창 R 의 같은 무게 평균을 배경에 걸쳐 평균 ·
            positions 몫은 같은 셋의 합 ÷ 전체 합.
  갈래 A 꼴(서술) 판정 노트 §5.2 와 같은 집계 — 합계 = 배경을 합친 ln(ΣD2/ΣP2) − ln(ΣD4/ΣP4) · 같은 창 = 모든 (배경, 창)
            R_w 의 같은 무게 평균을 계통별로 낸 뒤 ln 차. 비대상 시험(selftest)이 §5.2 의 glob μ0 네 수치를 이 꼴로 재현한다.
            1차 정의를 배경별 값으로 둔 까닭: 원고 P4a 는 배경별 ln 차의 평균 · 배경 재추출 구간이다 — 그 꼴 그대로 창만 맞춘다.
  배경 자격 · 재추출 · 판정 꼴은 동결 판정기와 같다:
     자격 = 합계 positions ≥ 100 · R 유한 · R > 1 (두 계통) · 자격 배경 ≥ 15.
     해부26 재추출 = rng('haebu26', 'P4') · 해부27 = rng('haebu27', 'P4|loc') · rng('haebu27', 'P4|glob') 독립.
     → 합계 값과 같은 창 값이 같은 재추출 행 위에서 나온다(차이의 구간도 같은 행으로).
     P4a 꼴 = 하한 > 0 · P4b 꼴 = 비 구간이 [0.8, 1.25] 안 ✅ / 밖 ❌ / 걸침 🟡. 같은 창 값에 붙인 표시는 서술이다.

가드 (하나라도 실패하면 같은 창 숫자를 내지 않는다)
  G1 동일성  다시 돌린 mut 팔(roll · μ 0)과 성장이 원 출력(mini26_out · mini27_out)과 정규형 JSON 으로 같다.
  G2 창 합   모든 계통 · 열 계수 전부에서 창 합 = 팔 합계 · 칸 0–59 · 홀짝 0/1 · 계통 키 집합 같음.
  G3 입력    시드 = 원 판정 배경 20 개씩 · 곁 파일의 mini27 sha256 = 동결값 · 래퍼 판 · 설정(--mus 0 · --rules roll · --nref 0)
             · 내려받기 목록(RUN.manifest.sha256)과 파일 sha256 이 같다.
  G4 계기    재계산한 합계 값(점추정 · 구간)이 동결 판정기 JSON(_RESULT_mini26.json · _RESULT_mini27.json)의 P4 값과 같다(|차| ≤ 1e-9).
"""
import glob
import hashlib
import json
import math
import os
import sys

import numpy as np

VERSION = "p4win_analyze v1.0.0"
HERE = os.path.dirname(os.path.abspath(__file__))
N_BOOT = 2000
WBIN = 50
WIN_BINS = 5                    # 250틱 창 = 50틱 칸 5 개
N_WIN = 12                      # 3,000 / 250
N_BINS = 60                     # 3,000 / 50
MIN_POS_WIN = 20
MAX_DROP_FRAC = 0.05
MIN_POS = 100                   # 동결 판정기와 같다
MIN_BG_P4 = 15
S_BAND = (math.log(0.8), math.log(1.25))
BURN_WIN = 4                    # 창 0–3 = t ≤ 1000
THIRDS = ((0, 4), (4, 8), (8, 12))
CNAMES = ["attempts", "draws", "positions", "copy_ticks", "stalls", "filtered", "err_realized",
          "kids_ok", "kids_sterile", "aborted"]
I_DRAW, I_POS = 1, 2
MINI27_SHA = "882891f63d68c516b196fc0bbfb00528f89a9e96cda70d876d4c5acbf53e1111"
WRAPPER_VERSION = "galA_window v1.0.0"
CONF_KEYS = ("K", "death", "grow", "assay", "every", "n0", "frac", "min_pop")
TOL = 1e-9
STUDY = {
    "26": {"seeds": list(range(26001, 26021)), "rho": 8.0, "geoms": ["glob"], "boot": "haebu26",
           "ref": "mini26_out", "ref_fn": "mini26_r8_s%d.json", "frozen": "_RESULT_mini26.json"},
    "27": {"seeds": list(range(27001, 27021)), "rho": 8.0, "geoms": ["glob", "loc"], "boot": "haebu27",
           "ref": "mini27_out", "ref_fn": "mini27_r8_s%d.json", "frozen": "_RESULT_mini27.json"},
}
# 갈래 A 판정 노트 §5.2 의 수치(비대상 시험의 맞출 값 · 같은 창 정의가 갈래 A 와 같은지 본다)
GALA_52 = {("glob", 0.0): (+0.0519, -0.0055, 8.9081, 10.7006)}


def rng_of(*parts):
    key = "|".join(str(p) for p in parts).encode()
    h = int.from_bytes(hashlib.sha256(key).digest()[:16], "little")
    return np.random.Generator(np.random.PCG64(np.random.SeedSequence(h)))


def boot_idx(prefix, cell, n):
    return rng_of(prefix, cell).integers(0, n, size=(N_BOOT, n))


def ci(vals, lo=2.5, hi=97.5):
    v = np.asarray(vals, float)
    v = v[np.isfinite(v)]
    if len(v) == 0:
        return float("nan"), float("nan")
    return float(np.percentile(v, lo)), float(np.percentile(v, hi))


def sha_file(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()


def canon(x):
    return json.dumps(x, sort_keys=True)


def same_mu(a, b):
    return abs(float(a) - float(b)) < 1e-12


# ───────────────────────── 읽기 ─────────────────────────
def to_windows(rows):
    """곁 파일 행 [칸, 홀짝, 계수 10 개] → (12, 10) 창 배열(홀짝 합침). 범위 밖 칸 · 홀짝 · 길이는 예외."""
    W = np.zeros((N_WIN, len(CNAMES)), dtype=np.int64)
    for r in rows:
        b, p = r[0], r[1]
        if len(r) != 2 + len(CNAMES) or not (0 <= b < N_BINS) or p not in (0, 1):
            raise ValueError("창 행 이상: 칸 %r · 홀짝 %r · 길이 %d" % (b, p, len(r)))
        W[b // WIN_BINS] += np.asarray(r[2:], dtype=np.int64)
    return W


def side_arm(side, geom, rule, mu, kind):
    hits = [a for a in side["arms"] if a["geom"] == geom and a["rule"] == rule and same_mu(a["mu"], mu)
            and a["kind"] == kind]
    if len(hits) != 1:
        raise KeyError("곁 파일 팔 %s/%s/%g/%s 이 %d 개" % (geom, rule, mu, kind, len(hits)))
    return hits[0]


def main_arm(z, geom, rule, mu, kind):
    arms = z["by_geom"][geom]["arms"] if "by_geom" in z else z["arms"]
    hits = [a for a in arms if a["rule"] == rule and same_mu(a["mu"], mu) and a["kind"] == kind]
    if len(hits) != 1:
        raise KeyError("본 출력 팔 %s/%s/%g/%s 이 %d 개" % (geom, rule, mu, kind, len(hits)))
    return hits[0]


def main_grow(z, geom):
    return z["by_geom"][geom]["grow"] if "by_geom" in z else z["grow"]


# ───────────────────────── 가드 ─────────────────────────
def g2_window_sums(arm, sarm, tag, bad):
    """G2 — 창 합 = 팔 합계(계통 · 열 계수 전부). 창 배열 dict 를 돌려준다."""
    W = {}
    ck, wk = set(arm["counters"]), set(sarm["win"])
    if ck != wk:
        bad.append("G2 %s 계통 키 다름: 계수기 %s · 창 %s" % (tag, sorted(ck), sorted(wk)))
        return W
    for key in sorted(ck):
        try:
            W[key] = to_windows(sarm["win"][key])
        except ValueError as e:
            bad.append("G2 %s %s %s" % (tag, key, e))
            continue
        tot = W[key].sum(axis=0)
        for i, nm in enumerate(CNAMES):
            if int(tot[i]) != int(arm["counters"][key][nm]):
                bad.append("G2 %s %s %s 창 합 %d ≠ 합계 %d" % (tag, key, nm, int(tot[i]), arm["counters"][key][nm]))
    return W


def check_run(run_dir, ref_dir, study, seeds, bad, strict):
    """G1 · G2 · G3(설정 · 판) — 배경마다 {geom: (본 팔, 창 dict)} 를 돌려준다."""
    S = STUDY[study]
    out = {}
    n_id = 0
    for s in seeds:
        tag = "h%s s%d" % (study, s)
        pz = os.path.join(run_dir, "mini27_r%g_s%d.json" % (S["rho"], s))
        ps = os.path.join(run_dir, "galA_r%g_s%d.win.json" % (S["rho"], s))
        pr = os.path.join(ref_dir, S["ref_fn"] % s)
        if not (os.path.exists(pz) and os.path.exists(ps) and os.path.exists(pr)):
            bad.append("G3 %s 파일 없음(본 %s · 곁 %s · 원 %s)" % (tag, os.path.exists(pz), os.path.exists(ps), os.path.exists(pr)))
            continue
        z, side, ref = json.load(open(pz)), json.load(open(ps)), json.load(open(pr))
        # G3 — 설정 · 판
        if z["geoms"] != S["geoms"] or [float(m) for m in z["mus"]] != [0.0] or z["rules"] != ["roll"] \
                or int(z["config"]["nref"]) != 0 or float(z["config"]["rho"]) != S["rho"] or int(z["config"]["seed"]) != s:
            bad.append("G3 %s 설정 다름: geoms %s · mus %s · rules %s · nref %s" % (tag, z["geoms"], z["mus"], z["rules"],
                                                                              z["config"]["nref"]))
        for k in CONF_KEYS:
            if abs(float(z["config"][k]) - float(ref["config"][k])) > 1e-12:
                bad.append("G3 %s 설정 %s %s ≠ 원 %s" % (tag, k, z["config"][k], ref["config"][k]))
        if side.get("version") != WRAPPER_VERSION or side.get("mini27_sha256") != MINI27_SHA or side.get("design") != "std" \
                or int(side.get("wbin", -1)) != WBIN or float(side.get("wave", -1)) != 0.0 \
                or side.get("L_res") != 4 or side.get("L_inj_mut") != 2:
            bad.append("G3 %s 곁 파일 판 · 설정 다름: %s · %s · wbin %s · wave %s" % (
                tag, side.get("version"), str(side.get("mini27_sha256"))[:12], side.get("wbin"), side.get("wave")))
        if strict and (side.get("python"), side.get("numpy")) != ("3.13.12", "2.4.3"):
            bad.append("G3 %s 실행 환경 %s · numpy %s (기대 3.13.12 · 2.4.3)" % (tag, side.get("python"), side.get("numpy")))
        out[s] = {}
        for geom in S["geoms"]:
            gtag = "%s %s" % (tag, geom)
            try:
                arm, sarm = main_arm(z, geom, "roll", 0.0, "mut"), side_arm(side, geom, "roll", 0.0, "mut")
                rarm = main_arm(ref, geom, "roll", 0.0, "mut")
            except KeyError as e:
                bad.append("G1 %s %s" % (gtag, e))
                continue
            # G1 — 동일성(팔 · 성장 정규형 JSON)
            if canon(arm) != canon(rarm):
                diff = [k for k in sorted(set(arm) | set(rarm)) if canon(arm.get(k)) != canon(rarm.get(k))]
                bad.append("G1 %s mut 팔 다름: %s" % (gtag, diff))
            elif canon(main_grow(z, geom)) != canon(main_grow(ref, geom)):
                bad.append("G1 %s 성장 다름" % gtag)
            else:
                n_id += 1
            W = g2_window_sums(arm, sarm, gtag, bad)
            out[s][geom] = (arm, W)
    return out, n_id


def check_manifest(run_dir, bad):
    """G3 — 내려받기 목록(RUN_DIR.manifest.sha256: 'sha256  상대경로' 줄)과 디스크 파일이 같고, 목록 밖 json 이 없다."""
    mf = run_dir.rstrip("/") + ".manifest.sha256"
    if not os.path.exists(mf):
        bad.append("G3 목록 없음 %s" % os.path.basename(mf))
        return 0
    listed = {}
    for line in open(mf, encoding="utf-8"):
        if line.strip():
            h, p = line.split(None, 1)
            listed[os.path.basename(p.strip())] = h
    disk = {os.path.basename(p) for p in glob.glob(os.path.join(run_dir, "*.json"))}
    if set(listed) != disk:
        bad.append("G3 목록 · 디스크 파일 집합 다름: 목록만 %s · 디스크만 %s" % (sorted(set(listed) - disk)[:5],
                                                                     sorted(disk - set(listed))[:5]))
    for fn, h in listed.items():
        p = os.path.join(run_dir, fn)
        if os.path.exists(p) and sha_file(p) != h:
            bad.append("G3 sha256 다름 %s" % fn)
    return len(listed)


# ───────────────────────── 배경 하나의 양 ─────────────────────────
def bg_stats(arm, W):
    c2, c4 = arm["counters"].get("1_2"), arm["counters"].get("0_4")
    nan = float("nan")
    if not c2 or not c4 or c2["positions"] == 0 or c4["positions"] == 0 or "1_2" not in W or "0_4" not in W:
        return {"qual": False, "P2": c2["positions"] if c2 else 0, "P4": c4["positions"] if c4 else 0}
    R2, R4 = c2["draws"] / c2["positions"], c4["draws"] / c4["positions"]
    qual = c2["positions"] >= MIN_POS and c4["positions"] >= MIN_POS and math.isfinite(R2) and math.isfinite(R4) \
        and R2 > 1 and R4 > 1
    d2, p2 = W["1_2"][:, I_DRAW].astype(float), W["1_2"][:, I_POS].astype(float)
    d4, p4 = W["0_4"][:, I_DRAW].astype(float), W["0_4"][:, I_POS].astype(float)
    ok = (p2 >= MIN_POS_WIN) & (p4 >= MIN_POS_WIN)
    with np.errstate(divide="ignore", invalid="ignore"):
        r2 = np.where(p2 > 0, d2 / np.where(p2 > 0, p2, 1), nan)
        r4 = np.where(p4 > 0, d4 / np.where(p4 > 0, p4, 1), nan)
    okb = ok & (r2 > 1) & (r4 > 1)
    burn = ok & (np.arange(N_WIN) >= BURN_WIN)
    st = {"qual": bool(qual), "P2": int(c2["positions"]), "P4": int(c4["positions"]), "R2": R2, "R4": R4,
          "lnd_sum": math.log(R2) - math.log(R4),
          "lns_sum": (math.log((R2 - 1) * 2) - math.log((R4 - 1) * 4)) if (R2 > 1 and R4 > 1) else nan,
          "n_drop": int((~ok).sum()), "n_drop_b": int((~okb).sum())}
    if ok.any():
        st["lnd_win"] = float(np.mean(np.log(r2[ok]) - np.log(r4[ok])))
        st["R2_eq"] = float(np.exp(np.mean(np.log(r2[ok]))))
        st["R4_eq"] = float(np.exp(np.mean(np.log(r4[ok]))))
        R2v = d2[ok].sum() / p2[ok].sum()
        R4rw = float((p2[ok] * r4[ok]).sum() / p2[ok].sum())
        st["lnd_rw"] = math.log(R2v) - math.log(R4rw)
    else:
        st["lnd_win"] = st["R2_eq"] = st["R4_eq"] = st["lnd_rw"] = nan
    st["lns_win"] = float(np.mean(np.log((r2[okb] - 1) * 2) - np.log((r4[okb] - 1) * 4))) if okb.any() else nan
    st["lnd_burn"] = float(np.mean(np.log(r2[burn]) - np.log(r4[burn]))) if burn.any() else nan
    with np.errstate(invalid="ignore"):
        st["R2_th"] = [float(np.nanmean(r2[a:b])) if np.isfinite(r2[a:b]).any() else nan for a, b in THIRDS]
        st["R4_th"] = [float(np.nanmean(r4[a:b])) if np.isfinite(r4[a:b]).any() else nan for a, b in THIRDS]
    st["sh2_th"] = [float(p2[a:b].sum() / p2.sum()) for a, b in THIRDS]
    st["sh4_th"] = [float(p4[a:b].sum() / p4.sum()) for a, b in THIRDS]
    st["D2"], st["D4"] = int(c2["draws"]), int(c4["draws"])
    st["r2w"], st["r4w"] = [float(x) for x in r2], [float(x) for x in r4]
    return st


def gala_style(rows):
    """갈래 A §5.2 꼴(서술): (배경을 합친 합계 ln 차, 모든 (배경, 창) R_w 같은 무게 평균의 ln 차)"""
    D2, P2 = sum(r["D2"] for r in rows), sum(r["P2"] for r in rows)
    D4, P4 = sum(r["D4"] for r in rows), sum(r["P4"] for r in rows)
    with np.errstate(invalid="ignore"):
        g2 = np.nanmean(np.array([r["r2w"] for r in rows], float))
        g4 = np.nanmean(np.array([r["r4w"] for r in rows], float))
    return math.log(D2 / P2) - math.log(D4 / P4), math.log(g2) - math.log(g4)


def col(rows, k):
    return np.array([r[k] for r in rows], float)


def summarize(rows, idx):
    """자격 배경(시드 순) 목록 → 평균 · 같은 재추출 행 위 구간 · 서술 양"""
    out = {"k": len(rows)}
    for k in ("lnd_sum", "lnd_win", "lns_sum", "lns_win", "lnd_rw"):
        v = col(rows, k)
        out[k] = [float(v.mean()), *ci(v[idx].mean(axis=1))]
    dv = col(rows, "lnd_win") - col(rows, "lnd_sum")
    out["win_minus_sum"] = [float(dv.mean()), *ci(dv[idx].mean(axis=1))]
    out["lnd_burn"] = float(col(rows, "lnd_burn").mean())
    for k in ("R2", "R4", "R2_eq", "R4_eq"):
        out["med_" + k] = float(np.median(col(rows, k)))
    for k in ("R2_th", "R4_th", "sh2_th", "sh4_th"):
        out[k] = [float(x) for x in np.mean(np.array([r[k] for r in rows], float), axis=0)]
    out["gala_sum"], out["gala_win"] = gala_style(rows)
    nw = N_WIN * len(rows)
    out["drop"] = [int(sum(r["n_drop"] for r in rows)), nw]
    out["drop_b"] = [int(sum(r["n_drop_b"] for r in rows)), nw]
    out["und_win"] = out["drop"][0] > MAX_DROP_FRAC * nw
    out["und_win_b"] = out["drop_b"][0] > MAX_DROP_FRAC * nw
    return out


def band(lo, hi):
    if S_BAND[0] <= lo and hi <= S_BAND[1]:
        return "✅"
    if hi < S_BAND[0] or lo > S_BAND[1]:
        return "❌"
    return "🟡"


def fmt3(v):
    return "%+.4f [%+.4f, %+.4f]" % tuple(v)


def fmt_ratio(v):
    return "%.3f [%.3f, %.3f]" % tuple(math.exp(x) for x in v)


# ───────────────────────── 비대상 시험 ─────────────────────────
def selftest(d):
    print("== %s · 비대상 시험 — 갈래 A 출력 %s (ρ8 · roll · mut)" % (VERSION, d))
    bad = []
    seeds = sorted(int(os.path.basename(p).split("_s")[1].split(".")[0])
                   for p in glob.glob(os.path.join(d, "galA_r8_s*.win.json")))
    print("   배경 %d: %s … %s" % (len(seeds), seeds[0] if seeds else "-", seeds[-1] if seeds else "-"))
    cells = {}
    for s in seeds:
        z = json.load(open(os.path.join(d, "mini27_r8_s%d.json" % s)))
        side = json.load(open(os.path.join(d, "galA_r8_s%d.win.json" % s)))
        for sa in side["arms"]:
            if sa["rule"] != "roll" or sa["kind"] != "mut":
                continue
            arm = main_arm(z, sa["geom"], "roll", sa["mu"], "mut")
            W = g2_window_sums(arm, sa, "selftest s%d %s %g" % (s, sa["geom"], sa["mu"]), bad)
            st = bg_stats(arm, W)
            if "lnd_win" in st:
                cells.setdefault((sa["geom"], float(sa["mu"])), []).append(st)
    print("-- G2(창 합 = 합계) 문제 %d" % len(bad))
    for b in bad[:10]:
        print("   ", b)
    print("-- 기하 · μ | 배경 | 1차 정의(배경별 · 평균): 합계 · 같은 창 · L2 무게 | 갈래 A 꼴: 합계 · 같은 창"
          " | L4 R 앞→뒤 3분의 1 | L2 positions 몫 앞·중·뒤 | 빠진 창")
    miss = 0
    for (geom, mu) in sorted(cells):
        rows = cells[(geom, mu)]
        r4 = np.mean([r["R4_th"] for r in rows], axis=0)
        sh = np.mean([r["sh2_th"] for r in rows], axis=0)
        gs, gw = gala_style(rows)
        line = "   %-5s %7.4f%% | %2d | %+.4f · %+.4f · %+.4f | %+.4f · %+.4f | %.4f → %.4f | %.2f · %.2f · %.2f | %d" % (
            geom, 100 * mu, len(rows), col(rows, "lnd_sum").mean(), np.nanmean(col(rows, "lnd_win")),
            np.nanmean(col(rows, "lnd_rw")), gs, gw, r4[0], r4[2], sh[0], sh[1], sh[2], sum(r["n_drop"] for r in rows))
        if (geom, mu) in GALA_52:
            e = GALA_52[(geom, mu)]
            got = (round(gs, 4), round(gw, 4), round(float(r4[0]), 4), round(float(r4[2]), 4))
            same = all(abs(a - b) < 5e-9 for a, b in zip(got, e))
            miss += 0 if same else 1
            line += "\n      ← §5.2 갈래 A 꼴: 합계 %+.4f · 같은 창 %+.4f · L4 R %.4f → %.4f · 재현 %s" % (*e, "✅" if same else "❌")
        print(line)
    print("-- 결과: G2 문제 %d · §5.2 재현 실패 %d" % (len(bad), miss))
    return 1 if (bad or miss) else 0


# ───────────────────────── 가드만 ─────────────────────────
def ident(run_dir, ref_dir, study):
    S = STUDY[study]
    seeds = sorted(int(os.path.basename(p).split("_s")[1].split(".")[0])
                   for p in glob.glob(os.path.join(run_dir, "mini27_r%g_s*.json" % S["rho"])))
    bad = []
    out, n_id = check_run(run_dir, ref_dir, study, seeds, bad, strict=False)
    print("== %s · 가드만 h%s · 시드 %s · G1 같은 팔 %d · 문제 %d" % (VERSION, study, seeds, n_id, len(bad)))
    for b in bad[:20]:
        print("   ", b)
    for s in out:
        for geom, (arm, W) in out[s].items():
            st = bg_stats(arm, W)
            print("   s%d %-4s 자격 %s · 빠진 창 %d · 창 합 R2 %.4f R4 %.4f (계수기와 같음: G2)" % (
                s, geom, st["qual"], st.get("n_drop", -1), st.get("R2", float("nan")), st.get("R4", float("nan"))))
    return 1 if bad else 0


# ───────────────────────── 본 계산 ─────────────────────────
def judge(run26, run27, json_out):
    print("== %s — 원고 Table S10 · S11 P4 의 같은 창 재계수 (사후 · 판정 아님)" % VERSION)
    print("   분석 환경 python %s · numpy %s · 이 파일 sha256 %s" % (sys.version.split()[0], np.__version__,
                                                                sha_file(os.path.abspath(__file__))))
    bad = []
    data = {}
    n_files = {}
    for study, run in (("26", run26), ("27", run27)):
        S = STUDY[study]
        n_files[study] = check_manifest(run, bad)
        data[study], n_id = check_run(run, os.path.join(HERE, S["ref"]), study, S["seeds"], bad, strict=True)
        print("-- h%s: 출력 %s · 목록 파일 %d · G1 같은 팔 %d / %d" % (study, run, n_files[study], n_id,
                                                                len(S["seeds"]) * len(S["geoms"])))
    # 배경 순서 = 시드 순 · 자격 = 동결 판정기와 같다
    stats = {}
    for study in ("26", "27"):
        for geom in STUDY[study]["geoms"]:
            rows = []
            for s in STUDY[study]["seeds"]:
                if s in data[study] and geom in data[study][s]:
                    arm, W = data[study][s][geom]
                    st = bg_stats(arm, W)
                    st["seed"] = s
                    rows.append(st)
            stats[(study, geom)] = rows
    # G4 — 합계 값이 동결 판정기와 같나
    q26 = [r for r in stats[("26", "glob")] if r["qual"]]
    q27l = [r for r in stats[("27", "loc")] if r["qual"]]
    q27g = [r for r in stats[("27", "glob")] if r["qual"]]
    if min(len(q26), len(q27l), len(q27g)) < MIN_BG_P4:
        bad.append("자격 배경 부족: h26 %d · h27 loc %d · glob %d (최소 %d)" % (len(q26), len(q27l), len(q27g), MIN_BG_P4))
    if not bad:
        f26 = json.load(open(os.path.join(HERE, STUDY["26"]["frozen"])))["P4"]
        f27 = json.load(open(os.path.join(HERE, STUDY["27"]["frozen"])))["P4"]
        i26 = boot_idx("haebu26", "P4", len(q26))
        ia, ib = boot_idx("haebu27", "P4|loc", len(q27l)), boot_idx("haebu27", "P4|glob", len(q27g))
        mine = {
            "h26 lnR": [col(q26, "lnd_sum").mean(), *ci(col(q26, "lnd_sum")[i26].mean(axis=1))],
            "h26 lnS": [col(q26, "lns_sum").mean(), *ci(col(q26, "lns_sum")[i26].mean(axis=1))],
            "h27 dlnR": [col(q27l, "lnd_sum").mean() - col(q27g, "lnd_sum").mean(),
                         *ci(col(q27l, "lnd_sum")[ia].mean(axis=1) - col(q27g, "lnd_sum")[ib].mean(axis=1))],
            "h27 lnS_loc": [col(q27l, "lns_sum").mean(), *ci(col(q27l, "lns_sum")[ia].mean(axis=1))]}
        froz = {"h26 lnR": f26["lnR"], "h26 lnS": f26["lnS"], "h27 dlnR": f27["dlnR"], "h27 lnS_loc": f27["lnS_loc"]}
        worst = 0.0
        for k in mine:
            dmax = max(abs(float(a) - float(b)) for a, b in zip(mine[k], froz[k]))
            worst = max(worst, dmax)
            if dmax > TOL:
                bad.append("G4 %s 재계산 %s ≠ 동결 %s (최대 차 %.2e)" % (k, [round(float(x), 6) for x in mine[k]],
                                                                 [round(float(x), 6) for x in froz[k]], dmax))
        print("-- G4 합계 값 = 동결 판정기 JSON: 최대 차 %.2e (허용 %.0e) · 자격 배경 h26 %d · h27 loc %d · glob %d"
              % (worst, TOL, len(q26), len(q27l), len(q27g)))
    print("-- 가드 문제 %d" % len(bad))
    for b in bad[:40]:
        print("   ", b)
    if bad:
        print("** 가드 실패 — 같은 창 숫자를 내지 않는다")
        return 2

    res = {"version": VERSION, "analyzer_sha256": sha_file(os.path.abspath(__file__)),
           "python": sys.version.split()[0], "numpy": np.__version__,
           "definitions": {"window_ticks": WBIN * WIN_BINS, "n_win": N_WIN, "min_pos_win": MIN_POS_WIN,
                           "max_drop_frac": MAX_DROP_FRAC, "burn_win": BURN_WIN},
           "inputs": {"run26": run26, "run27": run27, "n_files": n_files}}

    # ── 해부26 (Table S10 · 최소 모형 · 전역 풀) ──
    s26 = summarize(q26, i26)
    res["h26"] = s26
    print("\n-- 해부26 (Table S10) · draw-first · ρ8 · μ 0 · mut 팔 · 자격 배경 %d / %d" % (s26["k"], len(STUDY["26"]["seeds"])))
    print("   R 중앙(합계)      L2 %.2f · L4 %.2f   ← 원고 9.86 · 9.30" % (s26["med_R2"], s26["med_R4"]))
    print("   R 중앙(같은 창 기하평균) L2 %.2f · L4 %.2f" % (s26["med_R2_eq"], s26["med_R4_eq"]))
    print("   P4a ln R(2) − ln R(4)  합계 %s ← 동결 판정기 값(원고 +0.055 [+0.041, +0.071])" % fmt3(s26["lnd_sum"]))
    print("                         같은 창 %s → 꼴(하한 > 0): %s%s" % (
        fmt3(s26["lnd_win"]), "✅" if s26["lnd_win"][1] > 0 else "❌", " · 판정 불가(빠진 창 > 5%)" if s26["und_win"] else ""))
    print("                         같은 창 − 합계 %s (같은 재추출 행)" % fmt3(s26["win_minus_sum"]))
    print("                         2차 · L2 무게로 맞춤 %s (서술)" % fmt3(s26["lnd_rw"]))
    print("                         번인(앞 4 창 뺌) 같은 창 점추정 %+.4f (서술)" % s26["lnd_burn"])
    print("                         갈래 A §5.2 꼴(서술): 배경 합친 합계 %+.4f · 창 R 전체 평균 %+.4f" % (s26["gala_sum"], s26["gala_win"]))
    print("   P4b (R−1)·L 비      합계 %s ← 원고 0.532 [0.523, 0.541] · 꼴 %s" % (fmt_ratio(s26["lns_sum"]), band(*s26["lns_sum"][1:])))
    print("                         같은 창 %s · 꼴 %s%s" % (fmt_ratio(s26["lns_win"]), band(*s26["lns_win"][1:]),
                                                    " · 판정 불가(빠진 창 > 5%)" if s26["und_win_b"] else ""))
    print("   추세(배경 평균): L4 R 앞·중·뒤 %.3f · %.3f · %.3f | L2 R %.3f · %.3f · %.3f" % (*s26["R4_th"], *s26["R2_th"]))
    print("                   L2 positions 몫 %.3f · %.3f · %.3f | L4 %.3f · %.3f · %.3f" % (*s26["sh2_th"], *s26["sh4_th"]))
    print("   빠진 창 %d / %d (P4b %d / %d)" % (*s26["drop"], *s26["drop_b"]))

    # ── 해부27 (Table S11 · 국소 결핍) ──
    sl, sg = summarize(q27l, ia), summarize(q27g, ib)
    dd_sum = col(q27l, "lnd_sum")[ia].mean(axis=1) - col(q27g, "lnd_sum")[ib].mean(axis=1)
    dd_win = col(q27l, "lnd_win")[ia].mean(axis=1) - col(q27g, "lnd_win")[ib].mean(axis=1)
    p4a_sum = [float(col(q27l, "lnd_sum").mean() - col(q27g, "lnd_sum").mean()), *ci(dd_sum)]
    p4a_win = [float(col(q27l, "lnd_win").mean() - col(q27g, "lnd_win").mean()), *ci(dd_win)]
    res["h27"] = {"loc": sl, "glob": sg, "P4a_sum": p4a_sum, "P4a_win": p4a_win}
    print("\n-- 해부27 (Table S11) · draw-first · ρ8 · μ 0 · mut 팔 · 자격 배경 loc %d · glob %d" % (sl["k"], sg["k"]))
    for nm, st in (("loc", sl), ("glob", sg)):
        print("   %-4s R 중앙(합계) L2 %.2f · L4 %.2f | 같은 창 기하평균 L2 %.2f · L4 %.2f" % (
            nm, st["med_R2"], st["med_R4"], st["med_R2_eq"], st["med_R4_eq"]))
        print("        ln R2 − ln R4 합계 %s · 같은 창 %s · 차 %s" % (fmt3(st["lnd_sum"]), fmt3(st["lnd_win"]),
                                                             fmt3(st["win_minus_sum"])))
        print("        2차 · L2 무게 %s · 번인 같은 창 %+.4f · 갈래 A 꼴 합계 %+.4f · 같은 창 %+.4f (서술)" % (
            fmt3(st["lnd_rw"]), st["lnd_burn"], st["gala_sum"], st["gala_win"]))
        print("        추세 L4 R 앞·중·뒤 %.3f · %.3f · %.3f | L2 R %.3f · %.3f · %.3f | L2 몫 %.3f · %.3f · %.3f · 빠진 창 %d / %d"
              % (*st["R4_th"], *st["R2_th"], *st["sh2_th"], *st["drop"]))
    print("   P4a [loc] − [glob]  합계 %s ← 동결 판정기 값(원고 +0.377 [+0.354, +0.400])" % fmt3(p4a_sum))
    print("                       같은 창 %s → 꼴(하한 > 0): %s%s" % (
        fmt3(p4a_win), "✅" if p4a_win[1] > 0 else "❌",
        " · 판정 불가(빠진 창 > 5%)" if (sl["und_win"] or sg["und_win"]) else ""))
    print("   P4b loc (R−1)·L 비  합계 %s ← 원고 0.783 [0.766, 0.801] · 꼴 %s" % (fmt_ratio(sl["lns_sum"]), band(*sl["lns_sum"][1:])))
    print("                       같은 창 %s · 꼴 %s · 상한 − ln 0.8 = %+.4f (합계 %+.4f)%s" % (
        fmt_ratio(sl["lns_win"]), band(*sl["lns_win"][1:]), sl["lns_win"][2] - S_BAND[0], sl["lns_sum"][2] - S_BAND[0],
        " · 판정 불가(빠진 창 > 5%)" if sl["und_win_b"] else ""))

    # ── 배경별 원장 ──
    print("\n-- 배경별(시드 · 자격 · R2 · R4 합계 · ln 차 합계 · 같은 창 · 빠진 창)")
    for (study, geom), rows in stats.items():
        for r in rows:
            print("   h%s %-4s s%d %s | %.3f %.3f | %+.4f %+.4f | %d" % (
                study, geom, r["seed"], "자격" if r["qual"] else "제외", r.get("R2", float("nan")), r.get("R4", float("nan")),
                r.get("lnd_sum", float("nan")), r.get("lnd_win", float("nan")), r.get("n_drop", -1)))
    res["per_bg"] = {"h%s_%s" % k: v for k, v in stats.items()}
    out = json_out or os.path.join(HERE, "_RESULT_p4win.json")
    with open(out, "w") as f:
        json.dump(res, f, indent=1, default=lambda x: None if (isinstance(x, float) and not math.isfinite(x)) else str(x))
    print("\n결과 JSON: %s" % os.path.basename(out))
    return 0


def main():
    a = sys.argv[1:]
    if len(a) >= 2 and a[0] == "selftest":
        return selftest(a[1])
    if len(a) >= 4 and a[0] == "ident":
        return ident(a[1], a[2], a[3])
    if len(a) >= 3 and a[0] == "judge":
        jo = a[a.index("--json") + 1] if "--json" in a else None
        return judge(a[1], a[2], jo)
    print(__doc__)
    return 1


if __name__ == "__main__":
    sys.exit(main())
