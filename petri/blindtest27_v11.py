# 기록용 사본 — 10-02 실제 실행 위치는 볼트 밖 임시 폴더(가짜 자료 out · rep · idbad · repbad 가 그 옆에 있어야 돈다). 인자: petri 디렉터리 경로.
"""해부 27 판정기 v1.1 눈가림 시험 — 버리는 시드 90201~90203 · 짧은 설정(성장 400 · 침입 200 · 기록 50) · 계기 줄만.
S4: 결과 · 해석 · 갈래 · 열쇠 · 판정 기호(✅❌🟡) 줄은 인쇄하지 않는다. 남는 것 = rc · G0/G1 머리 줄(숫자 지움) · '판정 불가' 표지의 유무.
interp() 단위 시험은 지어낸 입력(자료 아님)으로 갈래 머리만 본다."""
import sys, io, re, os, json, shutil, contextlib, hashlib
P = sys.argv[1]
sys.path.insert(0, P)
import mini27_analyze as A
S = os.path.dirname(os.path.abspath(__file__))
A.JUDGE_SEEDS = [90201, 90202, 90203]; A.N_BG = 3; A.MIN_BG_P4 = 2; A.N_BOOT = 200
A.CONFIG = dict(A.CONFIG, grow=400, assay=200, every=50)
ID, M26 = os.path.join(P, "ident27"), os.path.join(P, "mini26_out")
SHA = hashlib.sha256(open(os.path.join(P, "mini27.py"), "rb").read()).hexdigest()
BANNED = ("결과", "해석", "갈래", "열쇠", "✅", "❌", "🟡", "단서", "교차", "Δ", "R 중앙", "P(", "반전")


def write_log(d, sha=SHA, ver="python 3.9.6 numpy 2.0.2", err=0):
    with open(d.rstrip("/") + ".log", "w", encoding="utf-8") as f:
        f.write("===== 시작 x · 작업 0 · 병렬 0 · %s\n%s  mini27.py\n===== 끝 x · DONE 0 · ERROR %d\n" % (ver, sha, err))


def run(d, tag, rep=None, ident=ID):
    rep = rep or S + "/rep"
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        try:
            rc = A.judge(d, rep, ident, M26, S + "/res_v11_%s.json" % tag)
        except Exception as e:
            rc = "EXC %s: %s" % (type(e).__name__, str(e)[:80])
    out = buf.getvalue().splitlines()
    heads = []
    for l in out:
        l = l.strip()
        if not (l.startswith("-- G0") or l.startswith("**") or l.startswith("M3") or "M3" in l[:40]):
            continue
        if any(b in l for b in BANNED):
            continue
        heads.append(re.sub(r"[-+]?\d[\d.]*%?", "#", l)[:90])
    und = {k: any(l.startswith(k) and "판정 불가" in l for l in out) for k in ("-- P1", "-- P3", "해석")}
    print("[%s] rc=%s · 줄 %d · 판정 불가 표지 P1 %s · P3 %s · 해석 %s" % (tag, rc, len(out), und["-- P1"], und["-- P3"], und["해석"]))
    for h in heads[:6]:
        print("     ", h)
    return rc


write_log(S + "/out"); write_log(S + "/rep")
rcs = {}
rcs["온전"] = run(S + "/out", "온전")


def corrupt(tag, fn, base="out"):
    d = S + "/c_" + tag
    shutil.rmtree(d, ignore_errors=True); shutil.copytree(S + "/" + base, d)
    write_log(d)
    p = d + "/mini27_r8_s90202.json"; z = json.load(open(p)); fn(z); json.dump(z, open(p, "w"))
    return d


def c1(z):
    a = [x for x in z["by_geom"]["loc"]["arms"] if x["rule"] == "find" and x["mu"] == 0.0 and x["kind"] == "ref2"][0]
    a["series"][3][1] += 1
def c2(z): z["by_geom"]["place"]["arms"][5]["cons_ok"] = False
def c3(z): z["by_geom"]["glob"]["arms"].pop(17)
def c4(z): z["by_geom"]["loc"]["arms"][9]["n_inject"] += 1
def c5(z): z["by_geom"]["loc"]["grow"]["qualified"] = False
for tag, fn in (("μ0궤적", c1), ("보존깃발", c2), ("팔삭제", c3), ("끼우기수", c4), ("자격없음", c5)):
    rcs[tag] = run(corrupt(tag, fn), tag)

# 재현 · G2a 대조 망가뜨림
write_log(S + "/repbad"); rcs["재현불일치"] = run(S + "/out", "재현불일치", rep=S + "/repbad")
write_log(S + "/idbad"); rcs["G2a불일치"] = run(S + "/out", "G2a불일치", ident=S + "/idbad")

# M3 — 로그 묶기
d = S + "/c_m3hash"; shutil.rmtree(d, ignore_errors=True); shutil.copytree(S + "/out", d); write_log(d, sha="0" * 64)
rcs["M3해시"] = run(d, "M3해시")
d = S + "/c_m3np"; shutil.rmtree(d, ignore_errors=True); shutil.copytree(S + "/out", d); write_log(d, ver="python 3.13.12 numpy 2.4.3")
rcs["M3판다름"] = run(d, "M3판다름")
d = S + "/c_m3err"; shutil.rmtree(d, ignore_errors=True); shutil.copytree(S + "/out", d); write_log(d, err=1)
rcs["M3에러"] = run(d, "M3에러")
d = S + "/c_m3nolog"; shutil.rmtree(d, ignore_errors=True); shutil.copytree(S + "/out", d)
if os.path.exists(d + ".log"):
    os.remove(d + ".log")
rcs["M3로그없음"] = run(d, "M3로그없음")


# 판정 불가 전파 — loc · roll · ρ8 · 1.625% 끼운 팔을 세 배경 모두 전멸로
d = S + "/c_und"; shutil.rmtree(d, ignore_errors=True); shutil.copytree(S + "/out", d); write_log(d)
for s in (90201, 90202, 90203):
    p = d + "/mini27_r8_s%d.json" % s; z = json.load(open(p))
    for a in z["by_geom"]["loc"]["arms"]:
        if a["rule"] == "roll" and abs(a["mu"] - A.MU_TOP) < 1e-12 and a["kind"] == "mut":
            a["extinct_tick"] = 100; a["series"][-1][1] = 0; a["series"][-1][2] = 0
    json.dump(z, open(p, "w"))
r = S + "/rep_und"; shutil.rmtree(r, ignore_errors=True); os.makedirs(r)
for f in ("mini27_r8_s90201.json", "mini27_r32_s90201.json"):
    shutil.copy(d + "/" + f, r + "/" + f)      # 재현 대조는 같은 망가뜨린 자료끼리(이 시험은 전파만 본다)
write_log(r)
rcs["판정불가전파"] = run(d, "판정불가전파", rep=r)

print("\n== interp() 단위 시험(지어낸 입력) — 갈래 머리만")
cases = [(True, "none", "none", 1.0, "갈래 A"), (True, "reversal", "none", 1.0, "갈래 B"), (True, "none", "reversal", 1.0, "갈래 C"),
         (True, "indeterminate", "none", 1.0, "P1 ✅ · 열쇠 불확정"), (True, "none", "indeterminate", 1.0, "P1 ✅ · 열쇠 불확정"),
         (True, "undecidable", "none", 1.0, "P1 ✅ · 열쇠 판정 불가"), (False, "reversal", "none", 1.0, "갈래 E"),
         (False, "undecidable", "none", 1.0, "P1 ❌ · K1 판정 불가"), (False, "indeterminate", "none", -1.0, "갈래 D(가)"),
         (False, "none", "none", 1.0, "갈래 D(나)"), ("undecidable", "none", "none", 1.0, "판정 불가")]
okn = 0
for p1, k1, k2, r0, want in cases:
    got = A.interp(p1, k1, k2, r0)
    ok = got.startswith(want)
    okn += ok
    print("   P1 %-11s K1 %-13s K2 %-13s → %s %s" % (p1, k1, k2, want, "ok" if ok else "틀림: " + got[:30]))
print("   interp 맞음 %d/%d" % (okn, len(cases)))
print("\n== 요약 rc:", rcs)
