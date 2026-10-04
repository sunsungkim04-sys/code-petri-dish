#!/usr/bin/env python3
"""해부 30 — E1 single-code census in Avida with conserved material.

Grid: L in {3,4,5,6} (family S, matched) x density {scarce, plentiful} x rule {df, ff, do} x mu x seeds.
Each dish: one ancestor at the centre (letters taken anywhere), T updates, census every 100 updates;
measures over the window (T/2, T] for census bucket 0 = organisms still carrying the founder code
(MATERIAL_CENSUS_FOUNDER; their offspring are counted whether or not they changed).

  python3 a30_e1.py run          run the grid (skips finished runs)       env: A30_SEEDS, A30_PAR, A30_ROOT
  python3 a30_e1.py summarize    per-run measures -> <root>/e1_rows.json
"""
import json, os, sys
from multiprocessing import Pool
import a30lib as A

ROOT = os.environ.get('A30_ROOT', os.path.expanduser('~/projects/avida-port/e30/e1'))
SEEDS = [int(x) for x in os.environ.get('A30_SEEDS', '1 2 3 4 5 6 7 8').split()]
PAR = int(os.environ.get('A30_PAR', '16'))
T = int(os.environ.get('A30_T', '12000'))
DENS = {'scarce': 200, 'plentiful': 640}
MUS = [0.0, 0.0003, 0.001, 0.003]
LS = (3, 4, 5, 6)


def jobs():
    return [(L, dn, rn, mu, s) for s in SEEDS for L in LS for dn in DENS for rn in A.RULES for mu in MUS]


def prefix(L, dn, rn, mu, s):
    return f'{ROOT}/runs/L{L}_{dn}_{rn}_mu{mu:g}_s{s}'


def one(a):
    L, dn, rn, mu, s = a
    p = prefix(*a)
    if os.path.exists(p + '.run.json') and json.load(open(p + '.run.json'))['rc'] == 0:
        return p, 'skip'
    ev = [f'u begin InjectSequence {A.FAMILY_S[L]} {A.CENTER} {A.CENTER + 1} -1 0', f'u {T + 1} Exit']
    info = A.run(p, ev, {'MATERIAL_DENSITY': DENS[dn], 'RANDOM_SEED': s, 'COPY_DRAW_RULE': A.RULES[rn],
                         'COPY_MAT_SUB_PROB': mu, 'MATERIAL_CENSUS_FOUNDER': A.FAMILY_S[L]})
    return p, info['rc']


def summarize():
    rows = []
    for a in jobs():
        L, dn, rn, mu, s = a
        p = prefix(*a)
        rec = {'L': L, 'density': dn, 'rule': rn, 'mu': mu, 'seed': s}
        try:
            info = json.load(open(p + '.run.json'))
            rec['rc'] = info['rc']
            cen = A.read_census(p)
            last = max(r['update'] for r in cen)
            rec['last_update'] = last
            rec['extinct'] = last < T or all(r['alive'] == 0 for r in cen if r['update'] == last and r['label'] == 0)
            s0 = A.window_sums(cen, 0, T // 2, T)          # bucket 0 = organisms still carrying the founder code
            m = A.measures(s0, L, 18, mu)
            rec.update(m)
            s4 = A.window_sums(cen, 4, T // 2, T)          # bucket 4 = all others (descendants that changed)
            rec['alive_other'] = s4['alive_mean']
            rec['TR_other'] = s4['gap_sum'] / s4['gap_n'] if s4['gap_n'] else float('nan')
            rec['sub_changed_per_written'] = s0['sub_changed'] / s0['written'] if s0['written'] else float('nan')
            rec['Ic'] = (rec['sub_changed_per_written'] / (mu * (A.N_KINDS - 1) / A.N_KINDS)) if mu > 0 else float('nan')
            rec['audit_viol'] = A.audit_violations(p)
        except Exception as e:
            rec['error'] = repr(e)
        rows.append(rec)
    with open(f'{ROOT}/e1_rows.json', 'w') as f:
        json.dump(rows, f)
    print('runs', len(rows), 'errors', sum('error' in r for r in rows), 'extinct', sum(bool(r.get('extinct')) for r in rows))


if __name__ == '__main__':
    os.makedirs(f'{ROOT}/runs', exist_ok=True)
    if sys.argv[1] == 'run':
        with Pool(PAR) as P:
            for p, rc in P.imap_unordered(one, jobs()):
                print(rc, p, flush=True)
    elif sys.argv[1] == 'summarize':
        summarize()
