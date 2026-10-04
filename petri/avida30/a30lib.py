#!/usr/bin/env python3
"""Avida conserved-material port (해부 30) — shared runner and census reader.

Runs the patched Avida (branch material-draw) with cfg/avida.cfg plus -set overrides.
Every run writes <prefix>.census.gz (per lineage label, every MATERIAL_CENSUS_INTERVAL updates),
<prefix>.audit.txt (conservation audit lines), <prefix>.log and <prefix>.run.json.
"""
import gzip, json, os, subprocess, time

HOME = os.path.expanduser('~')
AVIDA = os.environ.get('A30_AVIDA', f'{HOME}/projects/avida-port/avida/cbuild/bin/avida')
CFG = os.environ.get('A30_CFG', f'{HOME}/projects/avida-port/e30/cfg')

# Matched ancestors (heads instruction set, h-copy -> h-copy-mat at the same letter 'v').
# Family S (18 letters): one movable nop-C moves from before h-search into the copy loop
# (between h-copy 'v' and if-label 'f'); loop L = 3 + number of nop-C inside the loop.
# Same length, same composition, same merit (constant), loop L 3..6.
FAMILY_S = {
    3: 'wzcagcccczvfcaxgab',
    4: 'wzcagccczvcfcaxgab',
    5: 'wzcagcczvccfcaxgab',
    6: 'wzcagczvcccfcaxgab',
}
N_KINDS = 26          # heads instruction set size (letters a..z)
CENTER = 30 * 60 + 30  # 60 x 60 world

# World constants fixed for 해부 30 (also in cfg/avida.cfg; repeated here so every run records them)
WORLD = {
    'WORLD_X': 60, 'WORLD_Y': 60, 'WORLD_GEOMETRY': 2,
    'SPECULATIVE': 0, 'SLICING_METHOD': 0, 'AVE_TIME_SLICE': 1,
    'BASE_MERIT_METHOD': 0, 'BIRTH_METHOD': 3, 'PREFER_EMPTY': 1, 'ALLOW_PARENT': 0,
    'DEATH_METHOD': 0, 'MATERIAL_LIFESPAN': 1200, 'MATERIAL_DEATH_RATE': 0.0003,
    'MATERIAL_MODE': 2, 'MATERIAL_REACH': 4, 'MATERIAL_DIFFUSE': 0.15,
    'MATERIAL_INJECT_TAKE': 0, 'MATERIAL_CENSUS_INTERVAL': 100,
}

RULES = {'df': 0, 'ff': 1, 'do': 2}


def run(prefix, events, sets, timeout=None, keep_data=False):
    """Run one Avida process. events: list of event lines. sets: dict of config overrides."""
    os.makedirs(os.path.dirname(prefix), exist_ok=True)
    ev = prefix + '.events'
    with open(ev, 'w') as f:
        f.write('\n'.join(events) + '\n')
    allsets = dict(WORLD)
    allsets.update(sets)
    allsets['EVENT_FILE'] = ev
    allsets['DATA_DIR'] = prefix + '.data'
    allsets['MATERIAL_CENSUS_FILE'] = prefix
    allsets.setdefault('VERBOSITY', 0)
    cmd = [AVIDA, '-c', 'avida.cfg']
    for k, v in allsets.items():
        cmd += ['-set', k, str(v)]
    t0 = time.time()
    with open(prefix + '.log', 'w') as log:
        rc = subprocess.run(cmd, cwd=CFG, stdout=log, stderr=subprocess.STDOUT, timeout=timeout).returncode
    el = time.time() - t0
    info = {'prefix': prefix, 'rc': rc, 'elapsed_s': round(el, 2), 'sets': allsets, 'events': events}
    with open(prefix + '.run.json', 'w') as f:
        json.dump(info, f, indent=1)
    if not keep_data:
        d = prefix + '.data'
        if os.path.isdir(d):
            for fn in os.listdir(d):
                os.remove(os.path.join(d, fn))
            os.rmdir(d)
    return info


COLS = ('update label alive attempts stalls nodraw_stalls draws reused err_draws err_stalls del ins_written '
        'sub_draws sub_written sub_changed written gap_sum gap_n divides div_fail_mat uncopied_mat extra_ret '
        'overwrite_ret placed dropped placed_changed placed_len free_total').split()


def read_census(prefix):
    rows = []
    with gzip.open(prefix + '.census.gz', 'rt') as f:
        head = f.readline().split()
        assert head == COLS, ('census header changed', head)
        for line in f:
            v = line.split()
            if len(v) != len(COLS):
                continue
            rows.append(dict(zip(COLS, map(int, v))))
    return rows


def audit_violations(prefix):
    """Number of audit lines reporting a violation (forced audit lines with 0 violations are not counted)."""
    n = 0
    p = prefix + '.audit.txt'
    if not os.path.exists(p):
        return None
    for line in open(p):
        if 'kind_or_neg_violations=' in line:
            v = int(line.split('kind_or_neg_violations=')[1].split()[0])
            if v > 0:
                n += 1
    return n


def window_sums(rows, label, u0, u1):
    """Sum census counters of one label over updates (u0, u1]."""
    s = {c: 0 for c in COLS[3:]}
    alive = []
    for r in rows:
        if r['label'] == label and u0 < r['update'] <= u1:
            for c in COLS[3:]:
                s[c] += r[c]
            alive.append(r['alive'])
    s['alive_mean'] = sum(alive) / len(alive) if alive else 0.0
    s['n_rows'] = len(alive)
    return s


def measures(s, L, n_letters, mu_sub):
    """Copy-outcome census measures of one window (see 해부 30 §5)."""
    w = s['written']
    out = {
        'written': w,
        'A': s['attempts'] / w if w else float('nan'),          # attempts per letter written
        'R': s['draws'] / w if w else float('nan'),             # draws per letter written
        'TR': s['gap_sum'] / s['gap_n'] if s['gap_n'] else float('nan'),  # cycles per attempt
        'h': s['stalls'] / w if w else float('nan'),            # stalls per letter written
        'F': s['sub_written'] / s['sub_draws'] if s['sub_draws'] else float('nan'),
        'placed': s['placed'],
        'supply': s['placed_changed'] / s['placed'] if s['placed'] else float('nan'),
        'alive_mean': s['alive_mean'],
    }
    out['W'] = (out['A'] - 1.0) * L if w else float('nan')     # stall cycles per letter
    if mu_sub > 0 and s['placed']:
        nominal = 1.0 - (1.0 - mu_sub * (N_KINDS - 1) / N_KINDS) ** n_letters
        out['nominal'] = nominal
        out['I'] = out['supply'] / nominal
    else:
        out['nominal'] = float('nan')
        out['I'] = float('nan')
    return out
