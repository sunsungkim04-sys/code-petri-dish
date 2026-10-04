#!/usr/bin/env python3
"""해부 30 — decision script for the Avida conserved-material replication (frozen with the prereg note).

  python3 a30_analyze.py e1 <e1_rows.json>
  python3 a30_analyze.py e2 <arms_table.json>

All thresholds below are fixed before any confirmatory run (사전등록/해부30 §5). Intervals: 95% percentile
bootstrap, 2,000 resamples of the independent units (E1: seeds; E2: backgrounds), the same resampled list
for every cell inside one resample (paired), RNG random.Random("haebu30-avida").
"""
import json, math, random, sys
from collections import defaultdict

B = 2000
RNG_KEY = 'haebu30-avida'

# ---- E1 thresholds ----
E1_MU_MID = 0.001          # error rate of the P2 comparison (chosen from the pilot, 해부 30 §7)
E1_LS = (3, 4, 5, 6)
W_STRONG, W_APPROX = 1.25, 1.5     # P3 band: max/min of W over L
R_PLENTY_MAX = 1.10                # P4
# ---- E2 thresholds ----
E2_SLOPE_MUS = (0.0, 0.0003, 0.001, 0.003, 0.01)   # slope fitted over this ladder (μ as a proportion)
SKIP_MAX = 0.05


def pct(xs, lo=0.025, hi=0.975):
    xs = sorted(x for x in xs if x == x)
    if not xs:
        return (float('nan'), float('nan'))
    def q(p):
        k = p * (len(xs) - 1)
        f = math.floor(k); c = math.ceil(k)
        return xs[f] if f == c else xs[f] + (xs[c] - xs[f]) * (k - f)
    return (q(lo), q(hi))


def mean(xs):
    xs = [x for x in xs if x == x]
    return sum(xs) / len(xs) if xs else float('nan')


def boot(units, stat):
    """stat(list_of_units) -> float or dict of floats. Returns point, and per-key (lo, hi)."""
    rng = random.Random(RNG_KEY)
    point = stat(units)
    draws = []
    for _ in range(B):
        draws.append(stat([units[rng.randrange(len(units))] for _ in units]))
    if isinstance(point, dict):
        return {k: (point[k], pct([d[k] for d in draws])) for k in point}
    return point, pct(draws)


def fmt(p, ci):
    return f'{p:+.3f} [{ci[0]:+.3f}, {ci[1]:+.3f}]'


# =========================== E1 ===========================
def e1(path):
    rows = json.load(open(path))
    errs = [r for r in rows if 'error' in r]
    ok = [r for r in rows if 'error' not in r and not r.get('extinct') and r.get('written', 0) > 0]
    print(f'# E1 runs {len(rows)} · errors {len(errs)} · extinct/empty-window {len(rows) - len(errs) - len(ok)}')
    viol = sum(1 for r in rows if r.get('audit_viol'))
    rc = sum(1 for r in rows if r.get('rc') not in (0, None))
    tr_bad = [r for r in ok if r['mu'] == 0 and abs(r['TR'] - r['L']) > 0.02]   # single code only at zero error
    print(f'G1 audit-violation runs = {viol} · nonzero rc = {rc} · T/R off by > 0.02 from L = {len(tr_bad)} '
          f'(μ = 0 runs; max |T/R - L| = {max((abs(r["TR"] - r["L"]) for r in ok if r["mu"] == 0), default=float("nan")):.4f})')
    cell = defaultdict(dict)       # (density, rule, mu, L) -> seed -> row
    for r in ok:
        cell[(r['density'], r['rule'], r['mu'], r['L'])][r['seed']] = r
    seeds = sorted({r['seed'] for r in rows})

    def lnR(units, key, field='R'):
        return mean([math.log(cell[key][s][field]) for s in units if s in cell[key] and cell[key][s][field] > 0])

    # P1: adjacent ln R differences under draw-first, scarce, mu = 0
    def p1(units):
        out = {}
        for L in (3, 4, 5):
            out[f'd{L}{L + 1}'] = lnR(units, ('scarce', 'df', 0.0, L)) - lnR(units, ('scarce', 'df', 0.0, L + 1))
        return out
    r1 = boot(seeds, p1)
    p1_pass = all(v[1][0] > 0 for v in r1.values())
    print('\n## P1 (primary) ln R(L) − ln R(L+1), draw-first, scarce, μ = 0')
    for k, (p, ci) in r1.items():
        print(f'  {k}: {fmt(p, ci)}')
    print(f'  P1 = {"PASS" if p1_pass else "FAIL"} (all three lower bounds > 0)')
    for rn in ('ff', 'do'):
        vals = [mean([cell[("scarce", rn, 0.0, L)][s]["R"] for s in seeds if s in cell[("scarce", rn, 0.0, L)]]) for L in E1_LS]
        print(f'  check {rn}: mean R by L = ' + ' '.join(f'{v:.4f}' for v in vals))

    # P2: realised substitution inflation Ic, draw-first vs controls, scarce, mu = mid
    def lnIc(units, key):
        return mean([math.log(cell[key][s]['Ic']) for s in units if s in cell[key] and cell[key][s]['Ic'] > 0])
    def p2(units):
        out = {}
        for L in E1_LS:
            for rn in ('ff', 'do'):
                out[f'g_{rn}_L{L}'] = lnIc(units, ('scarce', 'df', E1_MU_MID, L)) - lnIc(units, ('scarce', rn, E1_MU_MID, L))
        for rn in ('ff', 'do'):
            out[f'G_{rn}'] = out[f'g_{rn}_L3'] - out[f'g_{rn}_L6']
        return out
    r2 = boot(seeds, p2)
    p2a = r2['g_ff_L3'][1][0] > 0 and r2['g_do_L3'][1][0] > 0
    p2b = r2['G_ff'][1][0] > 0 and r2['G_do'][1][0] > 0
    print(f'\n## P2 (primary) ln Ic(draw-first) − ln Ic(control), scarce, μ = {E1_MU_MID}')
    for k, (p, ci) in r2.items():
        print(f'  {k}: {fmt(p, ci)}')
    p2all = all(r2[f'g_{rn}_L{L}'][1][0] > 0 for rn in ('ff', 'do') for L in E1_LS)
    print(f'  P2a (primary: L3, both controls, lower > 0) = {"PASS" if p2a else "FAIL"} · every L = {"yes" if p2all else "no"}'
          f' · P2b (secondary: gap larger at L3 than L6, both controls) = {"PASS" if p2b else "FAIL"}')
    p2_pass = p2a
    verdict = '✅' if (p1_pass and p2_pass) else ('🟡' if (p1_pass or p2_pass) else '❌')
    print(f'\n## E1 primary verdict: {verdict}  (✅ P1 ∧ P2a · 🟡 one of them · ❌ neither)')

    # P3: W = (A - 1) L independent of L within a band
    def p3(units):
        ws = [mean([cell[('scarce', 'df', 0.0, L)][s]['W'] for s in units if s in cell[('scarce', 'df', 0.0, L)]]) for L in E1_LS]
        return {'Wratio': max(ws) / min(ws)}
    r3 = boot(seeds, p3)
    wr, wci = r3['Wratio']
    grade = 'strong (≤ 1.25)' if wr <= W_STRONG else ('approximate (≤ 1.5)' if wr <= W_APPROX else 'not supported (> 1.5)')
    ws = [mean([cell[('scarce', 'df', 0.0, L)][s]['W'] for s in seeds if s in cell[('scarce', 'df', 0.0, L)]]) for L in E1_LS]
    print(f'\n## P3 (secondary) W by L (draw-first, scarce, μ = 0): ' + ' '.join(f'L{L} {w:.2f}' for L, w in zip(E1_LS, ws)))
    print(f'  max/min = {wr:.3f} [{wci[0]:.3f}, {wci[1]:.3f}] → {grade}')
    for mu in sorted({k[2] for k in cell if k[2] > 0}):
        wsm = [mean([cell[('scarce', 'df', mu, L)][s]['W'] for s in seeds if s in cell[('scarce', 'df', mu, L)]]) for L in E1_LS]
        rsm = [mean([cell[('scarce', 'df', mu, L)][s]['R'] for s in seeds if s in cell[('scarce', 'df', mu, L)]]) for L in E1_LS]
        print(f'  (descriptive) μ = {mu}: W by L ' + ' '.join(f'{w:.2f}' for w in wsm) + ' · R by L ' + ' '.join(f'{r:.3f}' for r in rsm))

    # P4: plentiful removes the redraw factor
    rp = [mean([cell[('plentiful', 'df', 0.0, L)][s]['R'] for s in seeds if s in cell[('plentiful', 'df', 0.0, L)]]) for L in E1_LS]
    print(f'\n## P4 (secondary) draw-first plentiful μ = 0 mean R by L: ' + ' '.join(f'{v:.4f}' for v in rp)
          + f' → {"PASS" if all(v <= R_PLENTY_MAX for v in rp) else "FAIL"} (all ≤ {R_PLENTY_MAX})')

    print('\n## descriptives: density rule μ L | n | A R T/R W F Ic supply alive (seed means)')
    for key in sorted(cell):
        rs = list(cell[key].values())
        print(' '.join(map(str, key)), '|', len(rs), '|',
              ' '.join(f'{mean([r[f] for r in rs]):.4g}' for f in ('A', 'R', 'TR', 'W', 'F', 'Ic', 'supply', 'alive_mean')))


# =========================== E2 ===========================
def ols_slope(xs, ys):
    pts = [(x, y) for x, y in zip(xs, ys) if y == y]
    if len(pts) < 2:
        return float('nan')
    mx = sum(p[0] for p in pts) / len(pts); my = sum(p[1] for p in pts) / len(pts)
    sxx = sum((p[0] - mx) ** 2 for p in pts)
    return sum((p[0] - mx) * (p[1] - my) for p in pts) / sxx if sxx else float('nan')


def e2(path):
    arms = json.load(open(path))
    errs = [a for a in arms if 'error' in a]
    good = [a for a in arms if 'error' not in a]
    tgt = sum(a['target'] for a in good); skp = sum(a['skipped'] for a in good)
    incomplete = sum(1 for a in good if not a['complete'])
    viol = sum(1 for a in good if a['audit_viol_lines'] or a['audit_after_inject'])
    lost = defaultdict(int)
    for a in good:
        if a['lost']:
            lost[(a['code'], a['rule'], a['mu'])] += 1
    print(f'# E2 arms {len(arms)} · errors {len(errs)} · incomplete {incomplete} · audit-violation arms {viol} · '
          f'skipped/target {skp}/{tgt} = {skp / tgt if tgt else float("nan"):.4f} (guard ≤ {SKIP_MAX})')
    print('  lost injected lineages (code, rule, μ): ' + (', '.join(f'{k}: {v}' for k, v in sorted(lost.items())) or 'none'))
    # per background Δ
    by = defaultdict(list)    # (bg, code, rule, mu) -> arm deltas
    for a in good:
        if a['complete']:
            by[(a['bg'], a['code'], a['rule'], a['mu'])].append(a['delta_arm'])
    bgs = sorted({a['bg'] for a in arms})
    mus = sorted({a['mu'] for a in arms})
    rules = ('df', 'ff', 'do')

    def D(bg, code, rn, mu):
        rn_eff = 'df' if mu == 0 else rn     # zero error: rules identical (E0c), the df arms stand for all
        x = by.get((bg, code, rn_eff, mu)); r = by.get((bg, 'ref', rn_eff, mu))
        return mean(x) - mean(r) if x and r else float('nan')

    def curve(units, code, rn):
        return [mean([D(b, code, rn, mu) for b in units]) for mu in mus]

    def stat(units):
        out = {}
        for rn in rules:
            c3 = curve(units, 'L3', rn); c5 = curve(units, 'L5', rn)
            sel = [i for i, mu in enumerate(mus) if mu in E2_SLOPE_MUS]
            xs = [mus[i] for i in sel]
            out[f'b3_{rn}'] = ols_slope(xs, [c3[i] for i in sel])
            out[f'b5_{rn}'] = ols_slope(xs, [c5[i] for i in sel])
            out[f'c_{rn}'] = ols_slope(xs, [c5[i] - c3[i] for i in sel])
        out['D3_mu0'] = curve(units, 'L3', 'df')[mus.index(0.0)]
        out['D5_mu0'] = curve(units, 'L5', 'df')[mus.index(0.0)]
        for rn in ('ff', 'do'):
            out[f'Q1_{rn}'] = out['b3_df'] - out[f'b3_{rn}']
            out[f'Q2_{rn}'] = out['c_df'] - out[f'c_{rn}']
        return out
    r = boot(bgs, stat)
    print(f'\n## Q0 (precondition) Δ(L3) at μ = 0: {fmt(*r["D3_mu0"])} · Δ(L5) at μ = 0: {fmt(*r["D5_mu0"])}')
    q0 = r['D3_mu0'][1][0] > 0
    print(f'  Q0 = {"PASS" if q0 else "FAIL"} (lower > 0: the shorter loop invades at zero error)')
    print(f'\n## Q1 (primary) slope of Δ(L3) on μ over μ ∈ {E2_SLOPE_MUS}, per rule (per unit μ)')
    for rn in rules:
        print(f'  b3_{rn}: {fmt(*r[f"b3_{rn}"])}')
    q1 = r['Q1_ff'][1][1] < 0 and r['Q1_do'][1][1] < 0
    print(f'  b3_df − b3_ff: {fmt(*r["Q1_ff"])} · b3_df − b3_do: {fmt(*r["Q1_do"])} → Q1 = {"PASS" if q1 else "FAIL"} (both upper < 0)')
    print(f'\n## Q2 (secondary) slope of Δ(L5) − Δ(L3) on μ (fidelity component)')
    for rn in rules:
        print(f'  c_{rn}: {fmt(*r[f"c_{rn}"])}')
    q2 = r['Q2_ff'][1][0] > 0 and r['Q2_do'][1][0] > 0
    print(f'  c_df − c_ff: {fmt(*r["Q2_ff"])} · c_df − c_do: {fmt(*r["Q2_do"])} → Q2 = {"PASS" if q2 else "FAIL"} (both lower > 0)')
    verdict = '✅' if (q0 and q1) else ('🟡' if q1 else '❌')
    print(f'\n## E2 primary verdict: {verdict}  (✅ Q0 ∧ Q1 · 🟡 Q1 without Q0 · ❌ Q1 fails)')
    print('\n## Q3 (open) mean Δ by μ and crossing (sign change, linear interpolation)')
    for code in ('L3', 'L5'):
        for rn in rules:
            c = curve(bgs, code, rn)
            cross = 'none in range'
            for i in range(len(mus) - 1):
                if c[i] == c[i] and c[i + 1] == c[i + 1] and (c[i] > 0) != (c[i + 1] > 0):
                    cross = f'{mus[i] + (mus[i + 1] - mus[i]) * c[i] / (c[i] - c[i + 1]):.5f}'
                    break
            print(f'  {code} {rn}: ' + ' '.join(f'μ{mu:g} {v:+.3f}' for mu, v in zip(mus, c)) + f' · crossing {cross}')
    # noise floor: two reference arms of one background, rule and μ
    nf = []
    for (bg, code, rn, mu), xs in by.items():
        if code == 'ref' and len(xs) >= 2:
            nf.append(abs(xs[0] - xs[1]))
    if nf:
        print(f'\n## noise floor: 95th percentile |Δ_arm(ref k0) − Δ_arm(ref k1)| = {pct(nf, 0.95, 0.95)[0]:.3f} (n = {len(nf)})')
    # draws counted directly in the injected lineage
    print('\n## injected-lineage draws per letter (mean over arms): code rule μ | R_injected R_resident supply_injected')
    agg = defaultdict(list)
    for a in good:
        agg[(a['code'], a['rule'], a['mu'])].append(a)
    for k in sorted(agg):
        xs = agg[k]
        print(' '.join(map(str, k)), '|', f"{mean([a['inj_R'] for a in xs if a['inj_R']]):.3f}",
              f"{mean([a['res_R'] for a in xs if a['res_R']]):.3f}",
              f"{mean([a['inj_supply'] for a in xs if a['inj_supply'] is not None]):.4f}")


if __name__ == '__main__':
    {'e1': e1, 'e2': e2}[sys.argv[1]](sys.argv[2])
