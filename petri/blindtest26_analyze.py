"""해부 26 판정기 개정 눈가림 시험 — 완전히 지어낸 가짜 JSON(시뮬레이션 아님 · 시드 90101~ 버리는 번호).
출력은 계기 줄만(rc · G0 문제 수 · '판정 불가' 표시 · 결과 문장의 상태 낱말)."""
import io, json, math, os, random, shutil, sys, contextlib, importlib.util
import numpy as np
P = os.path.join(os.path.dirname(os.path.abspath(__file__)), "mini26_analyze.py")
spec = importlib.util.spec_from_file_location("an", P); an = importlib.util.module_from_spec(spec); spec.loader.exec_module(an)
an.JUDGE_SEEDS = list(range(90101, 90131))
import tempfile
W = tempfile.mkdtemp(prefix="blind26_") + "/fake"   # 볼트 밖 임시 폴더 · judge 의 --json 도 그 안

def arm(rule, mu, kind, rnd, ext=False, pos=500):
    L = 2 if kind == "mut" else 4
    ser = []
    for i in range(13):
        n = 600 + rnd.randint(-20, 20); m = 60 + rnd.randint(-30, 30)
        ser.append([250 * i, n, m, m, n - 5])
    ser[0] = [0, 600, 60, 60, 600]
    if ext:
        ser[-1] = [3000, 0, 0, 0, 0]
    c = lambda P_: {"attempts": 3 * P_, "draws": 3 * P_ if L == 2 else 2 * P_, "positions": P_, "copy_ticks": 3 * P_ * L,
                    "stalls": P_, "filtered": 0, "err_realized": int(P_ * mu), "kids_ok": 50, "kids_sterile": 1 + rnd.randint(0, 9),
                    "aborted": 0}
    return {"rule": rule, "mu": mu, "kind": kind, "L_inj": L, "n_inject": 60, "skipped": 0, "series": ser,
            "extinct_tick": -1, "counters": {"0_4": c(pos), "1_%d" % L: c(pos)}, "cons_ok": True, "cons_checks": 14,
            "cons_bad": [], "end_free": [0] * 11}

def make(d, seeds, mod=None):
    os.makedirs(d, exist_ok=True)
    for rho in (8.0, 32.0):
        for s in seeds:
            rnd = random.Random("%s|%s" % (rho, s))
            z = {"version": an.VERSION, "config": dict(an.CONFIG, rho=rho, seed=s), "mus": an.MUS, "rules": an.RULES, "M": 8000,
                 "grow": {"series": [], "extinct_tick": -1, "n_end": 600, "qualified": True, "checksum": "x",
                          "counters": {"0_4": {"draws": 30, "positions": 10, "attempts": 30, "stalls": 20}},
                          "cons_ok": True, "cons_checks": 26, "cons_bad": [], "end_free": [0] * 11, "occ_end": 990},
                 "arms": [arm(r, mu, k, rnd) for r in an.RULES for mu in an.MUS for k in ["mut"] + an.REF_KINDS], "seconds": 1.0}
            am = {(x["rule"], x["mu"], x["kind"]): x for x in z["arms"]}
            for k in ["mut"] + an.REF_KINDS:          # μ0 에서 두 규칙 같게(실제 모형의 성질)
                am[("find", 0.0, k)]["series"] = json.loads(json.dumps(am[("roll", 0.0, k)]["series"]))
                am[("find", 0.0, k)]["counters"] = json.loads(json.dumps(am[("roll", 0.0, k)]["counters"]))
            if mod:
                mod(rho, s, z)
            json.dump(z, open(os.path.join(d, "mini26_r%g_s%d.json" % (rho, s)), "w"))

def run(name, mod=None, repro=True):
    shutil.rmtree(W, ignore_errors=True)
    seeds = list(range(90101, 90121))
    make(W + "/out", seeds, mod)
    if repro:
        os.makedirs(W + "/rep", exist_ok=True)
        for rho in (8, 32):
            f = "mini26_r%d_s90101.json" % rho
            shutil.copy(W + "/out/" + f, W + "/rep/" + f)
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        try:
            rc = an.judge(W + "/out", W + "/rep" if repro else None, W + "/res.json")
        except Exception as e:
            rc = "EXC %s: %s" % (type(e).__name__, e)
    out = buf.getvalue().splitlines()
    keep = [l for l in out if ("G0 점검 문제" in l or "판정 불가" in l or "G0 실패" in l or "재현" in l or "다름" in l
                               or "키 집합" in l or l.startswith("결과 문장"))]
    # 결과 문장은 상태 낱말만
    keep = [("결과 문장(상태): " + " ".join(w for w in l.split() if w in ("✅", "❌", "🟡") or "불가" in w or "보류" in w))
            if l.startswith("결과 문장") else l.strip()[:150] for l in keep]
    print("[%s] rc=%s" % (name, rc))
    for l in keep:
        print("    ", l)

KIND = "mut"
def ext_at(rule, rho_t, mu_t, nbg):
    def m(rho, s, z):
        if rho == rho_t and s < 90101 + nbg:
            for a in z["arms"]:
                if a["rule"] == rule and a["mu"] == mu_t and a["kind"] == KIND:
                    a["series"][-1] = [3000, 0, 0, 0, 0]
    return m

def lowpos(nbg):
    def m(rho, s, z):
        if rho == 8.0 and s < 90101 + nbg:
            for a in z["arms"]:
                if a["rule"] == "roll" and a["mu"] == 0.0 and a["kind"] == "mut":
                    a["counters"]["1_2"]["positions"] = 50
                if a["rule"] == "find" and a["mu"] == 0.0 and a["kind"] == "mut":
                    a["counters"]["1_2"]["positions"] = 50
    return m

def badinj(rho, s, z):
    if rho == 8.0 and s == 90105:
        z["arms"][7]["n_inject"] = 59

def badkeys(rho, s, z):
    if rho == 32.0 and s == 90103:
        for a in z["arms"]:
            if a["rule"] == "find" and a["mu"] == 0.0 and a["kind"] == "mut":
                a["counters"]["9_9"] = dict(a["counters"]["0_4"])

run("T0 온전")
run("T1 --repro 없음", repro=False)
run("T2 roll·8·μ1% 칸 뺀 배경 3", ext_at("roll", 8.0, 0.01, 3))
run("T2b roll·8·μ1% 칸 뺀 배경 2(경계 · 판정 가능해야)", ext_at("roll", 8.0, 0.01, 2))
run("T3 find·8·μ0.5% 뺀 배경 3 (P2 만 불가)", ext_at("find", 8.0, 0.005, 3))
run("T3b roll·32·μ0.2% 뺀 배경 3 (P3 만 불가)", ext_at("roll", 32.0, 0.002, 3))
run("T4 P4 자격 14/20", lowpos(6))
run("T4b P4 자격 0/20", lowpos(20))
KIND = "ref2"
run("T2c roll·8·μ1% 기준 팔 전멸 배경 3 → G1 멈춤", ext_at("roll", 8.0, 0.01, 3))
KIND = "mut"
run("T6 팔 하나 n_inject 다름", badinj)
run("T7 μ0 계수기 키 집합 다름", badkeys)
# T5 단위: NaN 곡선
try:
    an.cens_flip([0.1, float("nan")] + [0.0] * 6); print("[T5a] NaN 곡선 → 예외 없음 ✗")
except an.Undecidable:
    print("[T5a] NaN 곡선 → Undecidable ✓ (예전엔 조용히 +inf)")
D = {mu: np.array([0.1] * 20) for mu in an.MUS}; D[0.002] = np.array([float("nan")] * 20)
cb = an.curve_boot(D, an.boot_idx("x", 20))
print("[T5b] 전부 NaN 인 μ 가 있는 곡선의 재추출 → NaN 행 %d/%d · inf 행 %d" % (np.isnan(cb).sum(), len(cb), (cb == math.inf).sum()))
