#!/usr/bin/env python3
"""해부 30 — E2 invasion assays in Avida with conserved material.

One Avida process per background: the resident code (family S, loop 4) is grown alone at zero error
(the three rules give identical trajectories there, E0c), and at update G the process forks one child
per arm (MaterialForkArms). Each arm clones the grown dish, re-seeds every random stream from its own
arm seed, replaces a random 10% of the organisms (same set for every arm of a background) with the
injected code under lineage label 1 (material-conserving: the replaced organism's letters return to
its cell, the injected letters are taken from its cell or nearby, else anywhere; failures are skipped
and counted), sets the arm's error rate and copying rule, and runs A further updates.

  python3 a30_e2.py run  <bg_seed> [<bg_seed> ...]     grow + fork arms (env: see DESIGN)
  python3 a30_e2.py table <root>                        per-arm Delta table (json) from finished runs

Δ (manuscript §3.1): per arm, x = ln((n1 + 0.5) / (n0 + 0.5)) with n1 the injected lineage and n0 the
resident lineage; Δ_arm = x(end) − x(injection). For a code in one background, rule and μ:
Δ = mean over its arms − mean over the reference arms (resident code injected) of the same background,
rule and μ.
"""
import gzip, hashlib, json, math, os, sys
import a30lib as A

ROOT = os.environ.get('A30_ROOT', os.path.expanduser('~/projects/avida-port/e30/e2'))
RESIDENT = A.FAMILY_S[4]
DESIGN = {
    'density': float(os.environ.get('A30_DENSITY', '200')),
    'G': int(os.environ.get('A30_G', '8000')),                   # growth updates (zero error)
    'A': int(os.environ.get('A30_A', '10000')),                   # assay updates
    'mus': [float(x) for x in os.environ.get('A30_MUS', '0 0.0003 0.001 0.003 0.01 0.03').split()],
    'rules': os.environ.get('A30_RULES', 'df ff do').split(),
    'codes': os.environ.get('A30_CODES', 'ref:4 L3:3 L5:5').split(),   # name:loop (ref = resident code)
    'n_arms': int(os.environ.get('A30_NARMS', '3')),             # arms per code x rule x mu
    'frac': 0.1,
    'parallel': int(os.environ.get('A30_PAR', '6')),
    'q': float(os.environ.get('A30_Q', '-1')),
}


def seed_of(*parts):
    h = hashlib.sha256('|'.join(map(str, parts)).encode()).hexdigest()
    return int(h[:7], 16) + 1     # < 2^28, positive


def arms_for(bg):
    lines, meta = [], []
    inj_seed = seed_of('inject', bg)
    for mu in DESIGN['mus']:
        rules = ['df'] if mu == 0 else DESIGN['rules']     # at zero error the rules are identical (E0c)
        for rn in rules:
            for code in DESIGN['codes']:
                name, L = code.split(':')
                for k in range(DESIGN['n_arms']):
                    arm = f'{name}_mu{mu:g}_{rn}_k{k}'
                    prefix = f'{ROOT}/bg{bg}/arms/{arm}'
                    ms = seed_of('arm', bg, arm)
                    lines.append(f"{arm} {prefix} {ms} {A.FAMILY_S[int(L)]} {DESIGN['frac']} {inj_seed} 1 "
                                 f"{A.RULES[rn]} {DESIGN['q']} {mu} 0 0")
                    meta.append({'arm': arm, 'code': name, 'L': int(L), 'mu': mu, 'rule': rn, 'k': k, 'main_seed': ms})
    return lines, meta, inj_seed


def run_bg(bg):
    d = f'{ROOT}/bg{bg}'
    os.makedirs(f'{d}/arms', exist_ok=True)
    lines, meta, inj_seed = arms_for(bg)
    with open(f'{d}/arms.txt', 'w') as f:
        f.write('# name prefix main_seed sequence frac inj_seed label rule q psub pins pdel\n')
        f.write('\n'.join(lines) + '\n')
    with open(f'{d}/arms_meta.json', 'w') as f:
        json.dump({'design': DESIGN, 'bg': bg, 'inj_seed': inj_seed, 'arms': meta}, f, indent=0)
    G, Aa = DESIGN['G'], DESIGN['A']
    ev = [f'u begin InjectSequence {RESIDENT} {A.CENTER} {A.CENTER + 1} -1 0',
          f"u {G} MaterialForkArms {d}/arms.txt {DESIGN['parallel']}",
          f'u {G + Aa + 1} Exit']
    info = A.run(f'{d}/grow', ev, {'MATERIAL_DENSITY': DESIGN['density'], 'RANDOM_SEED': bg, 'COPY_DRAW_RULE': 0,
                                   'COPY_MAT_SUB_PROB': 0.0, 'MATERIAL_AUDIT_INTERVAL': 100})
    done = open(f'{d}/arms.txt.done').read().strip() if os.path.exists(f'{d}/arms.txt.done') else 'missing'
    print(bg, 'rc', info['rc'], 'elapsed_s', info['elapsed_s'], done, flush=True)
    return info


def arm_record(prefix, G, Aa):
    info = {}
    for line in open(prefix + '.info'):
        k, v = line.split(None, 1)
        info[k] = v.strip()
    rows = A.read_census(prefix)
    ups = sorted({r['update'] for r in rows})
    u0, u1 = ups[0], ups[-1]
    def n(label, u):
        for r in rows:
            if r['update'] == u and r['label'] == label:
                return r['alive']
        return 0
    n0a, n1a, n0b, n1b = n(0, u0), n(1, u0), n(0, u1), n(1, u1)
    x0 = math.log((n1a + 0.5) / (n0a + 0.5))
    x1 = math.log((n1b + 0.5) / (n0b + 0.5))
    s1 = A.window_sums(rows, 1, u0, u1)
    s0 = A.window_sums(rows, 0, u0, u1)
    return {'u0': u0, 'u1': u1, 'n0_start': n0a, 'n1_start': n1a, 'n0_end': n0b, 'n1_end': n1b,
            'delta_arm': x1 - x0, 'lost': n1b == 0, 'complete': u1 >= u0 + Aa - 100,
            'done': int(info['done']), 'target': int(info['target']), 'skipped': int(info['skipped']),
            'audit_after_inject': int(info['audit_after_inject']), 'audit_viol_lines': A.audit_violations(prefix),
            'inj_A': s1['attempts'] / s1['written'] if s1['written'] else None,
            'inj_R': s1['draws'] / s1['written'] if s1['written'] else None,
            'inj_supply': s1['placed_changed'] / s1['placed'] if s1['placed'] else None,
            'res_R': s0['draws'] / s0['written'] if s0['written'] else None,
            'res_supply': s0['placed_changed'] / s0['placed'] if s0['placed'] else None}


def table(root):
    out = []
    for dn in sorted(os.listdir(root)):
        if not dn.startswith('bg'):
            continue
        meta = json.load(open(f'{root}/{dn}/arms_meta.json'))
        G, Aa = meta['design']['G'], meta['design']['A']
        for m in meta['arms']:
            p = f"{root}/{dn}/arms/{m['arm']}"
            rec = dict(m, bg=meta['bg'])
            try:
                rec.update(arm_record(p, G, Aa))
            except Exception as e:   # missing or broken arm: recorded, never silently dropped
                rec['error'] = repr(e)
            out.append(rec)
    with open(f'{root}/arms_table.json', 'w') as f:
        json.dump(out, f)
    print('arms', len(out), 'errors', sum('error' in r for r in out))
    return out


if __name__ == '__main__':
    cmd = sys.argv[1]
    if cmd == 'run':
        for bg in sys.argv[2:]:
            run_bg(int(bg))
    elif cmd == 'table':
        table(sys.argv[2] if len(sys.argv) > 2 else ROOT)
