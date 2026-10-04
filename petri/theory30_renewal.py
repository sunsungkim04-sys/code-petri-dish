#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""[10-03 정본 사본 — scratchpad/theory/renewal_mustar.py 에서 옮김 · 독립 감사 scratchpad/audit_theory 반영 · 끝에 수렴(collapse) 블록 추가 · 출력 _RESULT_theory30.json]
Renewal-process account of the crossing error rate mu* (journal items M4 + M7).  POST HOC.

Reads ONLY existing result files in petri/ (read-only).  Every input is loaded from a file and
printed with its source key/line.  Nothing in the vault is written.

Theory (see report):
  draws per letter      R(q) = 1 + q*h                    (Wald identity; h = stalls per letter written)
  ins/sub survival      Pe(q) = b / (b + (1-b) q)         (geometric availability b of the random letter)
  filtering             F(q) = F0 * (11 + 41 Pe(q)) / 52  (letter-change weights: skip 1/3, ins 1/3, sub 10/11 of 52/33)
  inflation             I(q) = R(q) F(q)
  genomic rate          U_c = n_c (52/33) mu I_c          (letter changes per offspring)
  lineage penalty       Delta(mu) = Delta0 + G [ln(1 - lam_f v_f) - ln(1 - lam_r v_r)],  v = 1 - exp(-U)
  crossing              Delta(mu*) = 0;  linearised  mu* = s0 / ((52/33)(lam_f n_f I_f - lam_r n_r I_r)),  s0 = Delta0/G
"""
import json, math, os, sys, re
import numpy as np

P = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1 else
    "/Users/minseokim/Documents/Obsidian Vault/03_Projects/Computationial biology/Code-Petri-Dish/petri")
J = lambda f: json.load(open(os.path.join(P, f), encoding="utf-8"))
T = lambda f: open(os.path.join(P, f), encoding="utf-8").read().splitlines()
SRC = []
def src(name, val, where):
    SRC.append((name, val, where)); return val

# ---------------------------------------------------------------- inputs
inv28 = J("_RESULT_inv28.json")
QS = [0.0, 0.25, 0.5, 0.75, 1.0]
qk = lambda q: {0.0: "0", 0.25: "0.25", 0.5: "0.5", 0.75: "0.75", 1.0: "1"}[q]
MUS = [0.0, 0.001, 0.003, 0.005, 0.01]
DEL = {q: {mu: src("Delta q=%s mu=%s" % (qk(q), mu), inv28["delta"]["%s|rascld|%s" % (qk(q), {0.0: "0", 0.001: "0.001", 0.003: "0.003", 0.005: "0.005", 0.01: "0.01"}[mu])],
              "_RESULT_inv28.json delta['%s|rascld|%s']" % (qk(q), mu)) for mu in MUS} for q in QS}
SLOPE = {q: src("slope q=%s" % qk(q), inv28["slopes"]["%s|rascld" % qk(q)], "_RESULT_inv28.json slopes['%s|rascld']" % qk(q)) for q in QS}
CROSS = {q: src("crossing q=%s" % qk(q), inv28["D4"]["crossing"][qk(q)], "_RESULT_inv28.json D4.crossing['%s']" % qk(q)) for q in QS}
inv21 = J("_RESULT_inv21.json")
X_OLD = src("crossing SetA 8-pt", (inv21["R4"]["flip_old"], *inv21["R4"]["flip_old_ci"]), "_RESULT_inv21.json R4.flip_old(+_ci)")
X_NEW = src("crossing new 40 bg", (inv21["R4"]["flip_new"], *inv21["R4"]["flip_new_ci"]), "_RESULT_inv21.json R4.flip_new(+_ci)")
D0_NEW = src("Delta0 new 40 bg", inv21["delta_new"]["rascld"]["0.0"], "_RESULT_inv21.json delta_new.rascld['0.0']")
G = src("G generations/assay", J("_RESULT_model12.json")["G"], "_RESULT_model12.json G (=3000/420; Supp S8 M1)")
NF, NR = 6, 5   # len('rascld'), len('racld')  (Table 3)
LF, LR = 2, 4   # loop ticks (Table 3)

# stalls per letter in the invasion arms at mu = 0 (q-independent; D5 of inv28 shows q changes nothing at mu 0)
ib = J("_RESULT_invbud.json")["G3"]
TPL_F = src("T founder, invasion arm mu0", ib["침입 rascld"]["tpl"], "_RESULT_invbud.json G3['침입 rascld'].tpl (txt line 22)")
TPL_R = src("T resident, invasion arm mu0", ib["침입 racld"]["tpl"], "_RESULT_invbud.json G3['침입 racld'].tpl (txt line 21)")
p14 = J("_RESULT_pairs14.json")
def tr_mean(L):
    v = [x["T"] / x["R"] for k, x in p14.items() if k.startswith("code|") and x["L"] == L]
    return sum(v) / len(v), len(v)
TR2, n2 = tr_mean(2); TR4, n4 = tr_mean(4)
src("T/R loop-2 class mean (%d codes)" % n2, TR2, "_RESULT_pairs14.json code|*|* T/R where L=2")
src("T/R loop-4 class mean (%d codes)" % n4, TR4, "_RESULT_pairs14.json code|*|* T/R where L=4")
R_ARM = {"f": TPL_F / TR2, "r": TPL_R / TR4}
# alternative: R of same-loop family members directly counted (pairs14), the founder/resident themselves not counted there
R_FAM = {"f": np.mean([p14[k]["R"] for k in p14 if k.startswith("code|rascld|") and p14[k]["L"] == 2]),
         "r": np.mean([p14[k]["R"] for k in p14 if k.startswith("code|racld|") and p14[k]["L"] == 4])}
src("R loop-2 rascld-family mean (7-letter codes)", R_FAM["f"], "_RESULT_pairs14.json code|rascld|* R where L=2")
src("R loop-4 racld-family mean (6-letter codes)", R_FAM["r"], "_RESULT_pairs14.json code|racld|* R where L=4")

sp = J("_RESULT_spec8.json")["cells"]
F1 = {"f": src("F1 founder d8 mu0.1%", (sp["orig/rascld/8/0.001"]["F"], *sp["orig/rascld/8/0.001"]["F_ci"]), "_RESULT_spec8.json cells['orig/rascld/8/0.001'].F"),
      "r": src("F1 resident d8 mu0.1%", (sp["orig/racld/8/0.001"]["F"], *sp["orig/racld/8/0.001"]["F_ci"]), "_RESULT_spec8.json cells['orig/racld/8/0.001'].F")}
RS = {c: {mu: src("R single-code d8 %s mu%s" % (code, mu), sp["orig/%s/8/%s" % (code, mu)]["R"], "_RESULT_spec8.json cells['orig/%s/8/%s'].R" % (code, mu))
          for mu in ("0.001", "0.01")} for c, code in (("f", "rascld"), ("r", "racld"))}
IS = {c: {mu: sp["orig/%s/8/%s" % (code, mu)]["Il"] for mu in ("0.001", "0.01")} for c, code in (("f", "rascld"), ("r", "racld"))}
# draw-once census (one seed, mu 1%, d8) — only in the regression-check text
reg9 = T("_RESULT_reg9_check.txt")
def grab(code):
    for i, l in enumerate(reg9):
        if "lin_%s_mu0p01" % code in l and "remember_die True" in l:
            R = float(re.search(r"R ([\d.]+)", l).group(1)); F = float(re.search(r"F ([\d.]+)", l).group(1))
            h = float(re.search(r"헛손질/맞게 ([\d.]+)", l).group(1))
            return R, F, h, i + 1
F0, R0 = {}, {}
for c, code in (("f", "rascld"), ("r", "racld")):
    R, F, h, ln = grab(code); F0[c] = src("F0 draw-once %s" % code, F, "_RESULT_reg9_check.txt line %d" % ln); R0[c] = src("R0 draw-once %s" % code, R, "_RESULT_reg9_check.txt line %d" % ln)
m12 = J("_RESULT_model12.json")
SIG = {"f": src("sterile share of variants, founder", m12["flow|rascld"]["cls"]["불임"], "_RESULT_model12.json flow|rascld.cls['불임'] (txt line 7)"),
       "r": src("sterile share of variants, resident", m12["flow|racld"]["cls"]["불임"], "_RESULT_model12.json flow|racld.cls['불임'] (txt line 6)")}
VFLOW = {"f": src("variant fraction/birth founder, arms mu0.3%", m12["flow|rascld"]["U"], "_RESULT_model12.json flow|rascld.U"),
         "r": src("variant fraction/birth resident, arms mu0.3%", m12["flow|racld"]["U"], "_RESULT_model12.json flow|racld.U")}
lin = J("_RESULT_lineage.json")["table"]
ULIN = {c: {mu: src("u_obs lineage census d8 %s mu%s" % (code, mu), lin["d08|%s|%s" % (code, mu)]["u"][0], "_RESULT_lineage.json table['d08|%s|%s'].u" % (code, mu))
            for mu in ("0.001", "0.005", "0.01")} for c, code in (("f", "rascld"), ("r", "racld"))}

K = 52 / 33  # letter-change probability per draw / mu (Table 1)

# ---------------------------------------------------------------- theory
def b_from_F1(F1v, F0v):
    # F1 = F0 (11 + 41 b)/52  ->  b
    return min(1.0, max(1e-6, (52 * F1v / F0v - 11) / 41))

def Pe(q, b): return b / (b + (1 - b) * q) if q > 0 else 1.0

def infl(c, q, mu, mode, F1v=None):
    """I_c(q, mu).  mode: 'A' (R at mu0 in arms) · 'B' (R scaled by nominal mu, single-code census) · 'C' (scaled by realised rate)."""
    F1v = F1[c][0] if F1v is None else F1v
    b = b_from_F1(F1v, F0[c])
    F = F0[c] * (11 + 41 * Pe(q, b)) / 52
    h0 = R_ARM[c] - 1
    if mode == "FAM": h0 = R_FAM[c] - 1
    if mode in ("A", "FAM"):
        return (1 + q * h0) * F
    r1, r2 = RS[c]["0.001"], RS[c]["0.01"]
    if mode == "B":
        x = min(max(mu, 0.001), 0.01)
        g = math.exp(math.log(r2 / r1) * (math.log(x) - math.log(0.001)) / (math.log(0.01) - math.log(0.001)))
        return (1 + q * (h0 + 1) * g - q) * F if q > 0 else F
    if mode == "C":   # feedback through the realised per-letter rate x = mu * I  (single-code calibration of ln R on ln x)
        x1, x2 = 0.001 * IS[c]["0.001"], 0.01 * IS[c]["0.01"]
        I = (1 + q * h0) * F
        for _ in range(60):
            x = min(max(mu * I, x1), x2) if mu > 0 else x1
            g = math.exp(math.log(r2 / r1) * (math.log(x) - math.log(x1)) / (math.log(x2) - math.log(x1)))
            I = (1 + q * ((h0 + 1) * g - 1)) * F
        return I
    raise ValueError(mode)

def delta_model(q, mu, d0, lam, mode, F1f=None, F1r=None):
    Uf = NF * K * mu * infl("f", q, mu, mode, F1f); Ur = NR * K * mu * infl("r", q, mu, mode, F1r)
    vf, vr = 1 - math.exp(-Uf), 1 - math.exp(-Ur)
    return d0 + G * (math.log(1 - lam[0] * vf) - math.log(1 - lam[1] * vr))

def cross_exact(q, d0, lam, mode, **kw):
    # first sign change scanning upward (the lambda<1 curve can turn back up at large mu), searched up to 5%
    f = lambda m: delta_model(q, m, d0, lam, mode, **kw)
    grid = np.geomspace(1e-5, 0.05, 400); lo = hi = None
    for a_, b_ in zip(grid[:-1], grid[1:]):
        if f(a_) > 0 >= f(b_): lo, hi = a_, b_; break
    if lo is None: return None
    for _ in range(80):
        m = (lo + hi) / 2
        lo, hi = (m, hi) if f(m) > 0 else (lo, m)
    return (lo + hi) / 2

def cross_grid(vals, mus=MUS):
    for i in range(1, len(mus)):
        if vals[i - 1] > 0 >= vals[i]:
            return mus[i - 1] + (mus[i] - mus[i - 1]) * vals[i - 1] / (vals[i - 1] - vals[i])
    return None

def slope_ls(vals, mus=MUS):
    x = np.array(mus); y = np.array(vals); return float(np.polyfit(x, y, 1)[0])

def Kq(q, mode, mu=0.002, lam=(1, 1)):
    return lam[0] * NF * infl("f", q, mu, mode) - lam[1] * NR * infl("r", q, mu, mode)

D0 = {q: DEL[q][0.0][0] for q in QS}
D0_POOL = float(np.mean([D0[q] for q in QS]))

def fit_lambda(q, target, mode, d0=None, grid=False):
    d0 = D0[q] if d0 is None else d0
    lo, hi = 0.02, 0.9999   # common lambda, 0 < lambda < 1
    def xing(l):
        if grid:
            return cross_grid([delta_model(q, m, d0, (l, l), mode) for m in MUS]) or 1.0
        return cross_exact(q, d0, (l, l), mode)
    for _ in range(70):
        m = (lo + hi) / 2
        x = xing(m)
        lo, hi = (m, hi) if (x is None or x > target) else (lo, m)
    return (lo + hi) / 2

def cross_lin(q, d0, l, mode):
    # linearised crossing, used only for lambda > 1 (penalty G*lam*dU)
    lo, hi = 1e-7, 0.05
    f = lambda m: d0 - G * l * (NF * K * m * infl("f", q, m, mode) - NR * K * m * infl("r", q, m, mode))
    if f(hi) > 0: return None
    for _ in range(80):
        m = (lo + hi) / 2; lo, hi = (m, hi) if f(m) > 0 else (lo, m)
    return (lo + hi) / 2

# ---------------------------------------------------------------- measured crossing uncertainty (approximate)
rng = np.random.default_rng(20261003)
def meas_cross_ci(q, n=4000):
    pts = DEL[q]; xs = []
    for _ in range(n):
        v = [rng.normal(pts[m][0], (pts[m][2] - pts[m][1]) / 3.92) for m in MUS]
        xs.append(cross_grid(v))
    ok = [x for x in xs if x is not None]
    return (np.percentile(ok, 2.5), np.percentile(ok, 97.5), 1 - len(ok) / n) if ok else (None, None, 1.0)

# ================================================================= output
out = {}
print("=" * 100)
print("RENEWAL THEORY OF mu* — post hoc · inputs read from petri/ (read-only)")
print("=" * 100)
print("\n-- derived inputs")
for c in ("f", "r"):
    b = b_from_F1(F1[c][0], F0[c])
    print("  %s: R_arm(mu0) = T/(T/R)_L = %.3f  -> h = %.3f · W = hL = %.2f ticks · R_family(pairs14) %.3f · b (random-letter availability from F1/F0) = %.3f"
          % ("founder rascld" if c == "f" else "resident racld ", R_ARM[c], R_ARM[c] - 1, (R_ARM[c] - 1) * (LF if c == "f" else LR), R_FAM[c], b))
print("  s0 = Delta0/G : q=1 %.4f · pooled-q Delta0 %.4f -> s0 %.4f per generation" % (D0[1.0] / G, D0_POOL, D0_POOL / G))

print("\n-- (1) draws per letter and inflation vs q  (mode A: R from invasion arms at mu 0)")
print("   q     R_f     R_r     F_f    F_r     I_f     I_r    K(q)=n_f I_f - n_r I_r")
for q in QS:
    Rf, Rr = 1 + q * (R_ARM["f"] - 1), 1 + q * (R_ARM["r"] - 1)
    If, Ir = infl("f", q, 0, "A"), infl("r", q, 0, "A")
    print("  %4.2f  %6.2f  %6.2f   %.3f  %.3f  %6.2f  %6.2f   %7.2f" % (q, Rf, Rr, If / Rf, Ir / Rr, If, Ir, NF * If - NR * Ir))
print("   check q=0 R: theory 1.00 vs draw-once census founder %.2f · resident %.2f (reg9, one seed, mu 1%%)" % (R0["f"], R0["r"]))

# ---- (2) lambda-free prediction: ratio mu*(q)/mu*(1)
print("\n-- (2) lambda-, G-, Delta0-free prediction: mu*(q)/mu*(1) = K(1)/K(q)  (linearised, equal lambda)")
print("   q    meas mu*   meas ratio | ratio A (arms mu0) · FAM (pairs14 family R) · A w/ F const · B (nominal-mu R) · C (realised-rate R)")
meas_ratio = {}
for q in QS[1:]:
    mr = CROSS[q] / CROSS[1.0]; meas_ratio[q] = mr
    rA = Kq(1, "A") / Kq(q, "A"); rF = Kq(1, "FAM") / Kq(q, "FAM")
    # F constant (= F1) variant
    Kc = lambda qq: NF * (1 + qq * (R_ARM["f"] - 1)) * F1["f"][0] - NR * (1 + qq * (R_ARM["r"] - 1)) * F1["r"][0]
    rFc = Kc(1) / Kc(q)
    rB = Kq(1, "B", mu=CROSS[1.0]) / Kq(q, "B", mu=CROSS[q])
    rC = Kq(1, "C", mu=CROSS[1.0]) / Kq(q, "C", mu=CROSS[q])
    print("  %4.2f  %.3f%%   %5.2f     | %5.2f · %5.2f · %5.2f · %5.2f · %5.2f" % (q, 100 * CROSS[q], mr, rA, rF, rFc, rB, rC))
    out["ratio|%s" % q] = dict(meas=mr, A=rA, FAM=rF, Fconst=rFc, B=rB, C=rC)
print("   q=0: K(0) = %.2f vs K(1) = %.2f -> predicted mu*(0)/mu*(1) = %.1f (i.e. ~%.2f%%, beyond the 1%% grid); measured: no crossing within 1%%"
      % (Kq(0, "A"), Kq(1, "A"), Kq(1, "A") / Kq(0, "A"), 100 * CROSS[1.0] * Kq(1, "A") / Kq(0, "A")))

# ---- (3) absolute level, no fitted parameter: lambda bounds
print("\n-- (3) absolute mu*, NO fitted parameter: lambda = 1 (every variant lost) and lambda = sterile share (only sterile variants lost)")
print("   exact = root of the model curve; grid = model evaluated at the 5 measured mu and interpolated like the measurement")
print("   q    measured [approx 95%]            | lam=1  exact/grid (A)  | lam=sterile exact/grid (A) | lam=1 grid (C) | lam=sterile grid (C)")
for q in QS:
    lo, hi, pnone = meas_cross_ci(q)
    row = []
    for lam, mode in (((1, 1), "A"), ((SIG["f"], SIG["r"]), "A"), ((1, 1), "C"), ((SIG["f"], SIG["r"]), "C")):
        lam2 = (min(lam[0], 0.9999), min(lam[1], 0.9999))
        ex = cross_exact(q, D0[q], lam2, mode)
        gr = cross_grid([delta_model(q, m, D0[q], lam2, mode) for m in MUS])
        row.append((ex, gr))
    f = lambda x: ("%.3f%%" % (100 * x)) if x else "none<1%"
    ms = f(CROSS[q]) + (" [%.3f, %.3f]" % (100 * lo, 100 * hi) if lo else "") + (" P(none)=%.2f" % pnone if pnone > 0.001 else "")
    print("  %4.2f  %-30s | %s / %s | %s / %s | %s | %s" % (q, ms, f(row[0][0]), f(row[0][1]), f(row[1][0]), f(row[1][1]), f(row[2][1]), f(row[3][1])))
    out["abs|%s" % q] = dict(meas=CROSS[q], meas_ci=(lo, hi), lam1_A=row[0], lamS_A=row[1], lam1_C=row[2], lamS_C=row[3])

# ---- (4) flow route (draw-first only): variant fractions measured inside the invasion arms at mu 0.3%
print("\n-- (4) draw-first, independent flow instrument (variant fraction per birth inside the invasion arms, mu 0.3%)")
Uf3, Ur3 = -math.log(1 - VFLOW["f"]), -math.log(1 - VFLOW["r"])
dU3 = Uf3 - Ur3
print("   U_g(0.3%%) founder %.4f · resident %.4f · dU %.4f  -> per unit mu %.1f" % (Uf3, Ur3, dU3, dU3 / 0.003))
print("   implied inflation at 0.3%%: founder %.2f · resident %.2f (vs census-A at mu0: %.2f · %.2f)"
      % (Uf3 / (NF * K * 0.003), Ur3 / (NR * K * 0.003), infl("f", 1, 0, "A"), infl("r", 1, 0, "A")))
for nm, lam in (("lam=1", (1, 1)), ("lam=sterile", (SIG["f"], SIG["r"]))):
    # linear-in-mu scaling of U from 0.3%
    lo, hi = 1e-6, 0.05
    f = lambda m: D0[1.0] + G * (math.log(1 - lam[0] * (1 - math.exp(-Uf3 * m / 0.003))) - math.log(1 - lam[1] * (1 - math.exp(-Ur3 * m / 0.003))))
    for _ in range(80):
        m = (lo + hi) / 2; lo, hi = (m, hi) if f(m) > 0 else (lo, m)
    # shape from lineage census ln I(mu) (as in f4prime): U(mu) = U(0.3%) * (mu/0.003) * I(mu)/I(0.003)
    def Ilin(c, mu):
        mus = [0.001, 0.005, 0.01]; I = [-math.log(1 - ULIN[c]["%s" % ({0.001: "0.001", 0.005: "0.005", 0.01: "0.01"}[x])]) / ((NF if c == "f" else NR) * K * x) for x in mus]
        x = min(max(mu, 0.001), 0.01)
        return float(np.exp(np.interp(x, mus, np.log(I))))
    g = lambda m: D0[1.0] + G * (math.log(1 - lam[0] * (1 - math.exp(-Uf3 * m / 0.003 * Ilin("f", m) / Ilin("f", 0.003))))
                               - math.log(1 - lam[1] * (1 - math.exp(-Ur3 * m / 0.003 * Ilin("r", m) / Ilin("r", 0.003)))))
    lo2, hi2 = 1e-6, 0.05
    for _ in range(80):
        m2 = (lo2 + hi2) / 2; lo2, hi2 = (m2, hi2) if g(m2) > 0 else (lo2, m2)
    print("   %-11s mu* (U linear in mu) %.3f%% · (U shaped by lineage-census I(mu)) %.3f%%  | measured q=1 %.3f%% (5-pt) · Set A 8-pt %.3f%% [%.3f, %.3f] · new %.3f%% [%.3f, %.3f]"
          % (nm, 100 * m, 100 * m2, 100 * CROSS[1.0], 100 * X_OLD[0], 100 * X_OLD[1], 100 * X_OLD[2], 100 * X_NEW[0], 100 * X_NEW[1], 100 * X_NEW[2]))
    out["flow|%s" % nm] = (m, m2)

# ---- (5) M7: U/s at the measured crossing
print("\n-- (5) M7 — excess realised genomic rate at the measured crossing, in units of s0")
s0 = D0[1.0] / G
for nm, dU in (("flow (arms, 0.3%, scaled linearly to mu*)", dU3 * CROSS[1.0] / 0.003),
               ("census A (arms R at mu0, F1)", K * CROSS[1.0] * Kq(1, "A")),
               ("census C (R falling with realised rate)", K * CROSS[1.0] * Kq(1, "C", mu=CROSS[1.0]))):
    print("   %-45s dU* = %.4f per offspring · s0 = %.4f per generation · dU*/s0 = %.2f  (theory: 1/lambda in [1, %.2f])" % (nm, dU, s0, dU / s0, 1 / min(SIG.values())))
    out["Us|%s" % nm] = dU / s0
Uf_star = NF * K * CROSS[1.0] * infl("f", 1, CROSS[1.0], "C"); Ur_star = NR * K * CROSS[1.0] * infl("r", 1, CROSS[1.0], "C")
print("   (absolute, census C at mu*: U_founder %.3f · U_resident %.3f letter changes per offspring; nominal n*52/33*mu = %.4f · %.4f)"
      % (Uf_star, Ur_star, NF * K * CROSS[1.0], NR * K * CROSS[1.0]))

# ---- (6) one fitted parameter: common lambda, fit on one condition, predict the rest (leave-one-out)
print("\n-- (6) ONE fitted parameter (common lambda), fitted on the q = 1 crossing (grid-interpolated, like the measurement); predict the rest")
for mode in ("A", "C"):
    lam_fit = fit_lambda(1.0, CROSS[1.0], mode, grid=True)
    preds = {}
    for q in QS:
        l2 = (min(lam_fit, 0.9999),) * 2
        preds[q] = cross_grid([delta_model(q, m, D0[q], l2, mode) for m in MUS])
    sl = {q: slope_ls([delta_model(q, m, D0[q], (min(lam_fit, .9999),) * 2, mode) for m in MUS]) for q in QS}
    print("   mode %s: lambda_fit = %.3f" % (mode, lam_fit))
    for q in QS:
        print("     q %.2f  mu* pred %-9s meas %-9s | slope pred %7.1f  meas %7.1f [%.1f, %.1f]"
              % (q, ("%.3f%%" % (100 * preds[q])) if preds[q] else "none", ("%.3f%%" % (100 * CROSS[q])) if CROSS[q] else "none", sl[q], *SLOPE[q]))
    # new backgrounds: same lambda, their own Delta0 (8-pt grid of inv21)
    MU8 = [0, 0.001, 0.002, 0.003, 0.004, 0.005, 0.0075, 0.01]
    pn = cross_grid([delta_model(1.0, m, D0_NEW[0], (min(lam_fit, .9999),) * 2, mode) for m in MU8], MU8)
    po = cross_grid([delta_model(1.0, m, D0[1.0], (min(lam_fit, .9999),) * 2, mode) for m in MU8], MU8)
    print("     draw-first, 8-pt grid: Set A pred %.3f%% meas %.3f%% [%.3f, %.3f] · new 40 bg pred %.3f%% meas %.3f%% [%.3f, %.3f]"
          % (100 * po, 100 * X_OLD[0], 100 * X_OLD[1], 100 * X_OLD[2], 100 * pn, 100 * X_NEW[0], 100 * X_NEW[1], 100 * X_NEW[2]))
    out["fit1|%s" % mode] = dict(lam=lam_fit, preds=preds, slopes=sl, new=pn, old=po)
    # leave-one-out over the four crossing q's
    errs = []
    for qfit in QS[1:]:
        lf = fit_lambda(qfit, CROSS[qfit], mode, grid=True)
        for q in QS[1:]:
            if q == qfit: continue
            p = cross_grid([delta_model(q, m, D0[q], (min(lf, .9999),) * 2, mode) for m in MUS])
            errs.append((qfit, q, lf, p, CROSS[q], math.log(p / CROSS[q]) if p else float("nan")))
    le = [abs(e[5]) for e in errs]
    print("     LOO (fit on one q, predict the other three): lambda range %.3f–%.3f · |ln(pred/meas)| median %.3f · max %.3f"
          % (min(e[2] for e in errs), max(e[2] for e in errs), float(np.median(le)), max(le)))
    out["loo|%s" % mode] = errs

# ---- (7) Monte Carlo on input uncertainty (Delta0, F1, measured q=1 crossing for the fit)
print("\n-- (7) Monte Carlo (n=1000): Delta0 ~ N(inv28 CI), F1 ~ N(spec8 CI), lambda refit on a resampled q=1 crossing; mode C")
N = 1000
sims = {q: [] for q in QS}; lams = []
for i in range(N):
    d0s = {q: rng.normal(D0[q], (DEL[q][0.0][2] - DEL[q][0.0][1]) / 3.92) for q in QS}
    f1f = rng.normal(F1["f"][0], (F1["f"][2] - F1["f"][1]) / 3.92); f1r = rng.normal(F1["r"][0], (F1["r"][2] - F1["r"][1]) / 3.92)
    v = [rng.normal(DEL[1.0][m][0], (DEL[1.0][m][2] - DEL[1.0][m][1]) / 3.92) for m in MUS]
    target = cross_grid(v)
    # fit lambda on this draw
    lo, hi = 0.05, 0.9999
    for _ in range(40):
        m = (lo + hi) / 2
        x = cross_grid([delta_model(1.0, mu, d0s[1.0], (m, m), "C", F1f=f1f, F1r=f1r) for mu in MUS]) or 1.0
        lo, hi = (m, hi) if x > target else (lo, m)
    lam = (lo + hi) / 2; lams.append(lam)
    for q in QS:
        sims[q].append(cross_grid([delta_model(q, mu, d0s[q], (lam, lam), "C", F1f=f1f, F1r=f1r) for mu in MUS]))
print("   lambda_fit %.3f [%.3f, %.3f]" % (np.median(lams), np.percentile(lams, 2.5), np.percentile(lams, 97.5)))
for q in QS:
    ok = [x for x in sims[q] if x]
    lo, hi, pn = meas_cross_ci(q)
    if ok:
        print("   q %.2f pred %.3f%% [%.3f, %.3f] (none in %.0f%%) · meas %s"
              % (q, 100 * np.median(ok), 100 * np.percentile(ok, 2.5), 100 * np.percentile(ok, 97.5), 100 * (1 - len(ok) / N),
                 ("%.3f%% [%.3f, %.3f]" % (100 * CROSS[q], 100 * lo, 100 * hi)) if CROSS[q] else "none"))
    else:
        print("   q %.2f pred: no crossing within 1%% in %d/%d draws · meas none" % (q, N - len(ok), N))
    out["mc|%s" % q] = [x for x in sims[q]]

# ---- (8) out-of-sample: lifespan (second-founder world, assay scaled with lifespan so G fixed)
print("\n-- (8) out-of-sample check — lifespan constant a = 300, 600, 1200 (an11 B; second-founder world)")
an = J("_RESULT_an11.json")
anB = T("_RESULT_an11.txt")
dl = {}
for i, l in enumerate(anB):
    mm = re.match(r"\s+(0|0\.001|0\.003|0\.005|0\.01) \|\s+([+-][\d.]+)\s+\|\s+([+-][\d.]+)\s+\|\s+([+-][\d.]+)", l)
    if mm: dl[float(mm.group(1))] = (float(mm.group(2)), float(mm.group(3)), float(mm.group(4)), i + 1)
src("an11 B Delta table (point)", dl, "_RESULT_an11.txt lines %d-%d" % (min(v[3] for v in dl.values()), max(v[3] for v in dl.values())))
sA = {(c, a): an["A_cells"]["%s|8|%d" % (code, a)]["s"] for c, code in (("f", "rascld"), ("r", "racld")) for a in (300, 600, 1200)}
for k, v in sA.items(): src("h single-code mu0 d8 %s a=%d" % k, v, "_RESULT_an11.json A_cells['%s|8|%d'].s" % ({"f": "rascld", "r": "racld"}[k[0]], k[1]))
LIFE = {a: src("mean age at death d8 a=%d (rascld)" % a, an["A_cells"]["rascld|8|%d" % a]["life"], "_RESULT_an11.json A_cells['rascld|8|%d'].life" % a) for a in (300, 600, 1200)}
ILc = {(c, a): src("I lineage census mu0.1%% d8 %s a=%d" % (code, a), an["L_cell|%s|8|%d" % (code, a)]["I"], "_RESULT_an11.json L_cell|%s|8|%d.I" % (code, a))
       for c, code in (("f", "rascld"), ("r", "racld")) for a in (300, 1200)}
def Ilife(c, a):  # ln-interpolate in ln a (a = 600 not measured)
    return math.exp(np.interp(math.log(a), [math.log(300), math.log(1200)], [math.log(ILc[(c, 300)]), math.log(ILc[(c, 1200)])]))
print("   G(a) = assay window / mean age at death, window = 3000 a/300 (an11 B: lengths scaled by the lifespan multiple)")
for mode in ("A", "C"):
    lamA = out["fit1|%s" % mode]["lam"]
    print("   mode %s, lambda fitted on draw-first Set A = %.3f" % (mode, lamA))
    print("     a   G(a)   meas mu*(grid)  slope meas | pred (i) h x single-code h(a)/h(300) · slope | pred (ii) I x lineage-census I(a)/I(300) · slope")
    for j, a in enumerate((300, 600, 1200)):
        vals = [dl[m][j] for m in MUS]; meas = cross_grid(vals); d0a = vals[0]
        G0 = G; G = 3000 * (a / 300) / LIFE[a]
        saveR, saveF0 = dict(R_ARM), dict(F0)
        for c in ("f", "r"): R_ARM[c] = 1 + (saveR[c] - 1) * sA[(c, a)] / sA[(c, 300)]
        p1 = cross_grid([delta_model(1.0, m, d0a, (lamA,) * 2, mode) for m in MUS]); s1 = slope_ls([delta_model(1.0, m, d0a, (lamA,) * 2, mode) for m in MUS])
        R_ARM.update(saveR)
        for c in ("f", "r"): F0[c] = saveF0[c] * Ilife(c, a) / Ilife(c, 300)   # scales I(q=1) by the measured ratio (F0 enters I multiplicatively)
        p2 = cross_grid([delta_model(1.0, m, d0a, (lamA,) * 2, mode) for m in MUS]); s2 = slope_ls([delta_model(1.0, m, d0a, (lamA,) * 2, mode) for m in MUS])
        F0.update(saveF0); Gs = G; G = G0
        print("  %5d  %.2f  %.3f%%  %s | %.3f%% · %.1f | %.3f%% · %.1f" % (a, Gs, 100 * meas, {300: "-111.9", 600: "-129.6", 1200: "-135.5"}[a], 100 * p1, s1, 100 * p2, s2))
        out["life|%s|%d" % (mode, a)] = (meas, p1, p2)

print("\n-- sources")
for n, v, w in SRC:
    if isinstance(v, dict): v = "(table)"
    print("   %-52s %-40s %s" % (n, (", ".join("%.5g" % x for x in v) if isinstance(v, (list, tuple)) else ("%.5g" % v if isinstance(v, float) else v)), w))
# ---------------------------------------------------------------- collapse (10-03 · 독립 감사 collapse.py 와 같은 축)
SIGQ = {q: NF * infl("f", q, 0.0, "A") - NR * infl("r", q, 0.0, "A") for q in QS}
pts = []
for q in QS:
    d0 = DEL[q][0.0]
    for mu in MUS[1:]:
        d = DEL[q][mu]
        pts.append(dict(q=q, mu=mu, x=K * SIGQ[q] * mu, dd=d[0] - d0[0], d=list(d), d0=list(d0)))
out["collapse"] = dict(sigma={str(q): SIGQ[q] for q in QS}, points=pts, note="x = (52/33) * Sigma(q) * mu, Sigma(q) = n_f I_f(q) - n_r I_r(q), mode A; dd = Delta(mu) - Delta(0)")
print("\n-- collapse (mode A)")
for p_ in sorted(pts, key=lambda z: z["x"]): print("   q=%-5s mu=%-6s x=%.4f dd=%+.3f" % (p_["q"], p_["mu"], p_["x"], p_["dd"]))
json.dump(out, open(os.path.join(P, "_RESULT_theory30.json"), "w"), indent=1, default=str)
