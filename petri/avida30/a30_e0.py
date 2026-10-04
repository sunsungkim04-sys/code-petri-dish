#!/usr/bin/env python3
"""해부 30 — E0 validation of the Avida conserved-material port (pilot seeds >= 9000 only).

  python3 a30_e0.py b      conservation audit over long runs (+ negative control: no return on death)
  python3 a30_e0.py c      zero error: the three rules give identical trajectories
  python3 a30_e0.py de     pilot grid: T/R = L (d) and R by rule x density x L (e)
Outputs under $A30_ROOT (default ~/projects/avida-port/e30/e0) and a text summary on stdout.
"""
import gzip, json, math, os, sys
from multiprocessing import Pool
import a30lib as A

ROOT = os.environ.get('A30_ROOT', os.path.expanduser('~/projects/avida-port/e30/e0'))
PAR = int(os.environ.get('A30_PAR', '8'))
D_SCARCE, D_PLENTY = 200, 640


def dish(L, seed, T, sets, extra_events=()):
    ev = [f'u begin InjectSequence {A.FAMILY_S[L]} {A.CENTER} {A.CENTER + 1} -1 0'] + list(extra_events) + [f'u {T} Exit']
    return ev


# ---------------- (b) conservation ----------------
def job_b(a):
    name, L, seed, T, sets = a
    p = f'{ROOT}/b/{name}'
    info = A.run(p, dish(L, seed, T, sets), dict(sets, RANDOM_SEED=seed))
    lines = open(p + '.audit.txt').read().splitlines() if os.path.exists(p + '.audit.txt') else []
    viol_lines = [l for l in lines if 'kind_or_neg_violations=' in l and int(l.split('kind_or_neg_violations=')[1].split()[0]) > 0]
    rows = A.read_census(p)
    last = max(r['update'] for r in rows)
    births = sum(r['placed'] for r in rows)
    return name, info['rc'], info['elapsed_s'], last, births, len(viol_lines), (viol_lines[0] if viol_lines else '')


def part_b():
    T = 30000
    base = {'MATERIAL_DENSITY': D_SCARCE, 'MATERIAL_AUDIT_INTERVAL': 10, 'MATERIAL_AUDIT_ABORT': 1}
    stress = {'COPY_MAT_SUB_PROB': 0.01, 'COPY_MAT_INS_PROB': 0.0033, 'COPY_MAT_DEL_PROB': 0.0033}
    jobs = []
    for rn, r in A.RULES.items():
        jobs.append((f'stress_{rn}', 3, 9201, T, dict(base, COPY_DRAW_RULE=r, **stress)))
        jobs.append((f'sub_{rn}', 4, 9202, T, dict(base, COPY_DRAW_RULE=r, COPY_MAT_SUB_PROB=0.003)))
    jobs.append(('stress_do_q05', 3, 9203, T, dict(base, COPY_DRAW_RULE=2, COPY_REDRAW_Q=0.5, **stress)))
    jobs.append(('stress_df_plenty', 6, 9204, T, dict(base, MATERIAL_DENSITY=D_PLENTY, COPY_DRAW_RULE=0, **stress)))
    # negative control: letters are not returned on death; the audit must catch it (abort off, count)
    jobs.append(('NEG_no_death_return', 3, 9205, 5000, dict(base, COPY_DRAW_RULE=0, COPY_MAT_SUB_PROB=0.003,
                                                             MATERIAL_DEBUG_NO_DEATH_RETURN=1, MATERIAL_AUDIT_ABORT=0)))
    print('# E0b conservation: name rc elapsed_s last_update births audit_violation_lines first_violation')
    with Pool(PAR) as P:
        for r in P.imap(job_b, jobs):
            print(*r, flush=True)


# ---------------- (c) identical trajectories at mu = 0 ----------------
def job_c(a):
    tag, L, d, rn, seed, T = a
    p = f'{ROOT}/c/{tag}_{rn}'
    ev = dish(L, seed, T, {}, extra_events=[f'u {T} SavePopulation'])
    info = A.run(p, ev, {'MATERIAL_DENSITY': d, 'RANDOM_SEED': seed, 'COPY_DRAW_RULE': A.RULES[rn],
                         'MATERIAL_AUDIT_INTERVAL': 50}, keep_data=True)
    spop = f'{p}.data/detail-{T}.spop'
    body = [l for l in open(spop) if not l.startswith('#')]
    os.makedirs(f'{ROOT}/c/spop', exist_ok=True)
    with gzip.open(f'{ROOT}/c/spop/{tag}_{rn}.spop.gz', 'wt') as f:
        f.writelines(body)
    for fn in os.listdir(p + '.data'):
        os.remove(os.path.join(p + '.data', fn))
    os.rmdir(p + '.data')
    rows = A.read_census(p)
    drop = {'nodraw_stalls', 'draws', 'reused'}
    traj = [tuple(r[c] for c in A.COLS if c not in drop) for r in rows]
    return tag, rn, info['rc'], body, traj, sum(r['draws'] for r in rows), sum(r['written'] for r in rows)


def part_c():
    jobs = []
    for tag, L, d, seed in (('L3_d200', 3, D_SCARCE, 9301), ('L6_d200', 6, D_SCARCE, 9302), ('L4_d640', 4, D_PLENTY, 9303)):
        for rn in A.RULES:
            jobs.append((tag, L, d, rn, seed, 8000))
    with Pool(PAR) as P:
        res = P.map(job_c, jobs)
    by = {}
    for tag, rn, rc, body, traj, draws, written in res:
        by.setdefault(tag, {})[rn] = (rc, body, traj, draws, written)
    print('# E0c zero error: tag | spop identical df=ff=do | census (excluding draw columns) identical | draws/written per rule | organisms in spop')
    for tag, d in by.items():
        sp = d['df'][1] == d['ff'][1] == d['do'][1]
        ce = d['df'][2] == d['ff'][2] == d['do'][2]
        rcs = [d[r][0] for r in A.RULES]
        rr = ' '.join(f'{r}={d[r][3] / d[r][4]:.4f}' for r in A.RULES)
        norg = sum(int(l.split()[4]) for l in d['df'][1] if l.strip()) if d['df'][1] else 0
        print(tag, 'rc', rcs, 'spop_identical', sp, 'census_identical', ce, rr, 'n_lines', len(d['df'][1]), flush=True)


# ---------------- (d, e) pilot grid ----------------
def job_de(a):
    L, d, rn, mu, seed, T = a
    p = f'{ROOT}/de/L{L}_d{d}_{rn}_mu{mu}_s{seed}'
    info = A.run(p, dish(L, seed, T, {}), {'MATERIAL_DENSITY': d, 'RANDOM_SEED': seed, 'COPY_DRAW_RULE': A.RULES[rn],
                                           'COPY_MAT_SUB_PROB': mu, 'MATERIAL_CENSUS_FOUNDER': A.FAMILY_S[L]})
    rows = A.read_census(p)
    s = A.window_sums(rows, 0, T // 2, T)          # bucket 0: organisms carrying the founder code
    m = A.measures(s, L, 18, mu)
    s4 = A.window_sums(rows, 4, T // 2, T)         # bucket 4: all other organisms
    return dict(L=L, d=d, rule=rn, mu=mu, seed=seed, rc=info['rc'], elapsed=info['elapsed_s'],
                viol=A.audit_violations(p), alive_other=s4['alive_mean'], **m)


def part_de():
    T = 10000
    jobs = [(L, d, rn, mu, seed, T) for L in (3, 4, 5, 6) for d in (D_SCARCE, D_PLENTY) for rn in A.RULES
            for mu in (0.0, 0.001, 0.003) for seed in (9001, 9002)]
    with Pool(PAR) as P:
        res = P.map(job_de, jobs)
    with open(f'{ROOT}/de_rows.json', 'w') as f:
        json.dump(res, f, indent=0)
    print('# E0d/e pilot (seeds 9001, 9002; window 5000-10000): L d rule mu | A R T/R W F supply I alive(founder) alive(other) | rc viol')
    for r in res:
        print(r['L'], r['d'], r['rule'], r['mu'], r['seed'], '|', *(f"{r[k]:.4g}" for k in ('A', 'R', 'TR', 'W', 'F', 'supply', 'I', 'alive_mean', 'alive_other')),
              '|', r['rc'], r['viol'])


if __name__ == '__main__':
    os.makedirs(ROOT, exist_ok=True)
    for part in sys.argv[1:]:
        {'b': part_b, 'c': part_c, 'de': part_de}[part]()
