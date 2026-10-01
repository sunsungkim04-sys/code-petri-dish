#!/usr/bin/env python3
"""해부 25 독립 재계수 (R13 · §12d).

판정기(pairs25_analyze.py 등)와 _RESULT_* 판정 출력을 열지 않고, 사전등록 노트
(해부25 §1·§4·§5 · 해부20 §6 · 해부21 §4 · 해부14 R 정의)의 정의만으로 원자료 JSON 에서
판정 통계량을 새로 계산한다. 난수 흐름은 판정기와 다르다(numpy default_rng, 칸 이름 해시).

사용: python3 recount25_independent.py <원자료 루트>
  루트 아래: main25/ inv25P_<b>/mu0/ inv25L/mu*/ inv25ib_racld/mu0/
             inv6/mu*/ inv21/mu*/ inv15E_<b>/mu0/ inv20P_<b>/mu0/ inv15Eib_racld/mu0/
  (서버 ~/petri 에서 rsync/tar 로 받은 사본)
"""
import json, glob, os, sys, math, hashlib
import numpy as np

ROOT = sys.argv[1] if len(sys.argv) > 1 else '.'
NB = 2000
MAIN5 = ['racld', 'rascld', 'rsacld', 'racldx', 'rascled']
DESC4 = ['acld', 'rcald', 'crald', 'reasccld']
MUS = [('mu0', 0.0), ('mu0p001', 0.001), ('mu0p002', 0.002), ('mu0p003', 0.003),
       ('mu0p004', 0.004), ('mu0p005', 0.005), ('mu0p0075', 0.0075), ('mu0p01', 0.01)]


def rng(name):
    h = int(hashlib.sha256(('recount25|' + name).encode()).hexdigest()[:16], 16)
    return np.random.default_rng(h)


def ci(x):
    return float(np.percentile(x, 2.5)), float(np.percentile(x, 97.5))


def J(p):
    with open(p) as f:
        return json.load(f)


# ---------- 고리 틱 (해부20 11d 규칙: 표지 s 뒤 ~ 첫 l 포함 · x · s 제외) ----------
def loop_ticks(code):
    start = code.index('s') + 1 if 's' in code else 0
    seg = code[start:code.index('l', start) + 1]
    return sum(1 for ch in seg if ch not in 'xs')


# ---------- s (짝: 해부20 §6) ----------
def s_pair(arm):
    ser = arm['series']
    t0, tT = ser[0], ser[-1]
    # 열: t · pop · tag1 · tag1_exact_injected · rest_exact_wt — m = tag1(끼운 표지 계통), n = pop
    n0, m0 = t0[1], t0[2]
    nT, mT = tT[1], tT[2]
    return math.log((mT + 0.5) / m0) - math.log((nT - mT + 0.5) / (n0 - m0))


def load_pairs(dirname):
    """반환: {seed: {code: s}} (mut 팔만) · 메타"""
    out, meta = {}, {}
    for p in sorted(glob.glob(os.path.join(ROOT, dirname, 'mu0', '*.json'))):
        d = J(p)
        seed = d['seed']
        out[seed] = {a['key']: s_pair(a) for a in d['arms'] if a['kind'] == 'mut'}
        collapsed = sum(1 for a in d['arms'] if a['kind'] == 'mut'
                        and a['series'][-1][2] / a['series'][0][2] < 0.10)
        meta[seed] = dict(gc=d['grow_checksum'], match=d.get('checksum_match'),
                          mu=d['mu_assay'], dens=d.get('density', d['opts'].get('density')),
                          top=d.get('main_top'), collapsed=collapsed,
                          fid=d['fidelity']['same'], cons=all(a.get('cons_ok', True) for a in d['arms']))
    return out, meta


def pair_list(base, codes):
    g = {}
    for c in codes:
        g.setdefault(loop_ticks(c), []).append(c)
    pairs, negs = [], []
    ks = sorted(g)
    for L in ks:
        if L + 1 in g:
            for a in g[L]:
                for b in g[L + 1]:
                    pairs.append((a, b, 1))  # ΔΔ = s(b: L+1) - s(a: L)
        cs = g[L]
        for i in range(len(cs)):
            for j in range(i + 1, len(cs)):
                negs.append((cs[i], cs[j]))
    return pairs, negs, {k: len(v) for k, v in g.items()}


def pair_matrix(S, base, seeds):
    codes = sorted(S[seeds[0]].keys())
    pairs, negs, grades = pair_list(base, codes)
    M = np.array([[(S[sd][b] - S[sd][a]) / dl for (a, b, dl) in pairs] for sd in seeds])
    return pairs, negs, grades, M


def boot_mean(X, name):
    """X: (n_bg, ...) → 배경 재추출 평균"""
    r = rng(name)
    n = X.shape[0]
    idx = r.integers(0, n, size=(NB, n))
    return X[idx].mean(axis=1)


# ---------- 사다리 (해부21 §4: s = ln((m_T+.5)/m_0) - ln((r_T+.5)/r_0) · 기준 = WT 팔 11 평균) ----------
def s_lad(arm):
    # m = tag1(열 2) · r = pop − tag1 (나머지). 다른 열 조합(exact 열 3 · rest_exact_wt 열 4)은
    # 6C 정본(Δ(1%) −0.753 · 교차 0.1827%)을 재현하지 못했다 — 이 조합만 재현(자기 검산 줄 참조).
    a, b = arm['series'][0], arm['series'][-1]
    return math.log((b[2] + 0.5) / a[2]) - math.log((b[1] - b[2] + 0.5) / (a[1] - a[2]))


def load_ladder(dirname):
    """반환: {seed: {mu: {code: Δ}}}, meta {seed: set(grow_checksum), top, extinct}"""
    out, meta = {}, {}
    for md, mu in MUS:
        for p in glob.glob(os.path.join(ROOT, dirname, md, '*.json')):
            d = J(p)
            assert abs(d['mu_assay'] - mu) < 1e-12, p
            sd = d['seed']
            m = meta.setdefault(sd, dict(gc=set(), top=set(), ext=set(), match=set()))
            m['gc'].add(d['grow_checksum']); m['match'].add(d.get('checksum_match'))
            m['top'].add(d['grow_top'][0][0] if d.get('grow_top') else d.get('main_top'))
            m['ext'].add(d.get('grow_extinct'))
            if not d['arms']:  # 배경 멸종 → 팔 없음 (자격에서 빠진다)
                continue
            wt = [s_lad(a) for a in d['arms'] if a['kind'] in ('ref', 'neutral')]
            assert len(wt) == 11 and len(d['arms']) == 14, p
            ref = sum(wt) / len(wt)
            out.setdefault(sd, {})[mu] = {a['key']: s_lad(a) - ref for a in d['arms'] if a['kind'] == 'mut'}
    return out, meta


MUV = np.array([m for _, m in MUS])


def crossing(curve):
    """curve: Δ 8점 (μ 순) · 검열: Δ(0) ≤ 0 → 0 · 1% 까지 + → +inf"""
    if curve[0] <= 0:
        return 0.0
    for j in range(1, len(curve)):
        if curve[j] <= 0:
            a, b = curve[j - 1], curve[j]
            return MUV[j - 1] + a * (MUV[j] - MUV[j - 1]) / (a - b)
    return math.inf


def lad_array(L, seeds, code='rascld'):
    return np.array([[L[sd][mu][code] for _, mu in MUS] for sd in seeds])


def boot_cross(A, name):
    r = rng(name)
    n = A.shape[0]
    idx = r.integers(0, n, size=(NB, n))
    means = A[idx].mean(axis=1)
    return np.array([crossing(c) for c in means]), means


def slope(A):
    # 기울기: 배경 평균 곡선의 μ 에 대한 최소제곱 기울기 (서술용)
    y = A.mean(axis=0)
    return float(np.polyfit(MUV, y, 1)[0])


# ---------- R (해부14: (copy_ok+copy_stall+copy_del)/copy_ok · 끼운 계통 = 인덱스 1 · 칸 합) ----------
def load_R(dirname):
    out = {}
    for p in glob.glob(os.path.join(ROOT, dirname, 'mu0', '*.json')):
        d = J(p)
        rr = {}
        for a in d['arms']:
            if a['kind'] != 'mut':
                continue
            ok = sum(b['copy_ok'][1] for b in a['bins'])
            st = sum(b['copy_stall'][1] for b in a['bins'])
            de = sum(b['copy_del'][1] for b in a['bins'])
            rr[a['key']] = (ok + st + de) / ok
        out[d['seed']] = dict(R=rr, gc=d['grow_checksum'], match=d['checksum_match'])
    return out


def fmt(x, p=3):
    return f'{x:+.{p}f}'


def main():
    pr = print
    problems = []
    # ===== 성장 · 자격 =====
    grow = {}
    for p in glob.glob(os.path.join(ROOT, 'main25', '*.json')):
        d = J(p)
        grow[d['opts']['seed']] = d
    seeds_all = sorted(grow)
    cons_bad = [s for s in seeds_all if grow[s]['conservation_violations'] != 0]
    qual = [s for s in seeds_all if grow[s]['extinct_at'] == -1 and grow[s]['final_counts'][0][0] == 'racld'
            and grow[s]['ticks_run'] == 50000 and grow[s]['opts']['density'] == 16]
    BG = qual[:20]
    pr(f'[G] 성장 기록 {len(seeds_all)} ({seeds_all[0]}~{seeds_all[-1]}) · 재료 보존 위반 {len(cons_bad)} · 자격 {len(qual)}/{len(seeds_all)} · 앞 20 = {BG}')
    if cons_bad:
        problems.append('cons')
    # ===== 짝 P =====
    S16, M16 = {}, {}
    for b in MAIN5 + DESC4:
        S16[b], M16[b] = load_pairs('inv25P_' + b)
        if sorted(S16[b]) != BG:
            problems.append(f'P {b} 배경 목록 다름')
        for sd, m in M16[b].items():
            if not (m['match'] and m['mu'] == 0 and m['dens'] == 16 and m['fid'] and m['cons']):
                problems.append(f'P {b} {sd} 메타 {m}')
    Lad16, LM16 = load_ladder('inv25L')
    if sorted(Lad16) != BG or any(len(Lad16[s]) != 8 for s in BG):
        problems.append('L 배경/μ 목록')
    IB16 = load_R('inv25ib_racld')
    if sorted(IB16) != BG:
        problems.append('IB 배경 목록')
    nP = sum(len(S16[b]) for b in S16); nL = sum(len(Lad16[s]) for s in Lad16)
    pr(f'[G] 파일 P {nP} · L {nL} · IB {len(IB16)}')
    # 같은 시드 18 파일 grow_checksum
    gc_bad = 0
    for sd in BG:
        g = {M16[b][sd]['gc'] for b in S16} | LM16[sd]['gc'] | {IB16[sd]['gc']}
        mm = {M16[b][sd]['match'] for b in S16} | LM16[sd]['match'] | {IB16[sd]['match']}
        if len(g) != 1 or mm != {True}:
            gc_bad += 1
    pr(f'[G] 같은 시드 18 파일 grow_checksum 다름(또는 MATCH 아님) {gc_bad}/20')
    # 주 다섯 짝
    tot_pairs = 0; per_pair = []
    allM = []
    for b in MAIN5:
        pairs, negs, grades, M = pair_matrix(S16[b], b, BG)
        tot_pairs += len(pairs)
        allM.append(M)
        for k, (a, c, _) in enumerate(pairs):
            per_pair.append((b, a, c, M[:, k]))
    pr(f'[G] 주 다섯 짝 수 {tot_pairs} (예상 54) · 고리 등급 ' +
       ' · '.join(f"{b} {pair_list(b, sorted(S16[b][BG[0]]))[2]}" for b in MAIN5))
    if tot_pairs != 54:
        problems.append('짝 수')
    pr(f'[G] 점검 문제 {len(problems)} {problems}')

    # ===== P1 =====
    X16 = np.hstack(allM)  # 20 × 54
    w16 = X16.mean(axis=1)
    b16 = boot_mean(w16, 'P1')
    p1 = (w16.mean(), *ci(b16))
    pr(f'\nP1 주 다섯 합친 틱당 ΔΔ {fmt(p1[0])} [{fmt(p1[1])}, {fmt(p1[2])}] → {"✅" if p1[2] < 0 else "❌"} (선: 상한 < 0)')
    # ===== P2 / N1 =====
    n_ok = 0; n_neg_pt = 0; fails = []
    for (b, a, c, v) in per_pair:
        bb = boot_mean(v, f'pair|{b}|{a}|{c}')
        lo, hi = ci(bb)
        if hi < 0:
            n_ok += 1
        else:
            fails.append((b, a, c, v.mean(), lo, hi))
        if v.mean() < 0:
            n_neg_pt += 1
    pr(f'P2 짝 상한 < 0: {n_ok}/54 → {"✅" if n_ok >= 49 else "❌"} (선 ≥ 49) · 점추정 음수 {n_neg_pt}/54 · 실패 {fails}')
    # 짝 최대 상한(경계 기록)
    his = []
    for (b, a, c, v) in per_pair:
        bb = boot_mean(v, f'pair|{b}|{a}|{c}')
        his.append((ci(bb)[1], b, a, c))
    his.sort(reverse=True)
    pr(f'   가장 0 에 가까운 상한 셋: ' + ' · '.join(f'{b} {a}->{c} {h:+.3f}' for h, b, a, c in his[:3]))
    for b in MAIN5:
        pairs, negs, grades, M = pair_matrix(S16[b], b, BG)
        w = M.mean(axis=1); bb = boot_mean(w, f'W1|{b}'); lo, hi = ci(bb)
        coll = sum(M16[b][s]['collapsed'] for s in BG); collbg = sum(1 for s in BG if M16[b][s]['collapsed'])
        pr(f'   {b:8s} W1 {fmt(w.mean())} [{fmt(lo)}, {fmt(hi)}] · 짝 {len(pairs)} · 무너진 mut 팔 {coll} (배경 {collbg})')
    for b in DESC4:
        pairs, negs, grades, M = pair_matrix(S16[b], b, BG)
        w = M.mean(axis=1); bb = boot_mean(w, f'W1|{b}'); lo, hi = ci(bb)
        his_b = [ci(boot_mean(M[:, k], f'pair|{b}|{k}'))[1] for k in range(M.shape[1])]
        n1 = sum(1 for h in his_b if h < 0)
        coll = sum(M16[b][s]['collapsed'] for s in BG); collbg = sum(1 for s in BG if M16[b][s]['collapsed'])
        pr(f'   (서술) {b:8s} W1 {fmt(w.mean())} [{fmt(lo)}, {fmt(hi)}] · N1 {n1}/{len(pairs)} (최대 상한 {max(his_b):+.3f}) · 무너진 mut 팔 {coll} (배경 {collbg})')

    # ===== P4 : 밀도 8 · 40 배경 (inv15E 20 + inv20P 20) =====
    S8 = {b: {} for b in MAIN5}
    for b in MAIN5:
        for dn in ('inv15E_' + b, 'inv20P_' + b):
            s, _ = load_pairs(dn)
            S8[b].update(s)
    BG8 = sorted(S8['racld'])
    assert all(sorted(S8[b]) == BG8 for b in MAIN5) and len(BG8) == 40, len(BG8)
    X8 = np.hstack([pair_matrix(S8[b], b, BG8)[3] for b in MAIN5])
    assert X8.shape == (40, 54)
    w8 = X8.mean(axis=1)
    b8 = boot_mean(w8, 'P4|d8')
    b16b = boot_mean(w16, 'P4|d16')
    dP4 = b16b - b8
    p4 = (w16.mean() - w8.mean(), *ci(dP4))
    read4 = '줄었다(예측)' if p4[1] > 0 else ('커졌다' if p4[2] < 0 else '구별 안 됨')
    pr(f'\nP4 밀도 8(40) 합친 {fmt(w8.mean())} [{fmt(ci(b8)[0])}, {fmt(ci(b8)[1])}] · 차 D16 − D8 {fmt(p4[0])} [{fmt(p4[1])}, {fmt(p4[2])}] → {read4}')
    # 밀도 8 자기 검산: 해부20 주 다섯 W1 (기준 표 −0.581 · −0.602 · −0.561 · −0.633 · −0.892)
    S20 = {b: load_pairs('inv20P_' + b)[0] for b in MAIN5}
    bg20 = sorted(S20['racld'])
    pr('   [자기 검산] 해부20 주 다섯 W1 점추정: ' + ' · '.join(
        f"{b} {pair_matrix(S20[b], b, bg20)[3].mean():+.3f}" for b in MAIN5))

    # ===== 사다리 =====
    A16 = lad_array(Lad16, BG)
    L6, LM6 = load_ladder('inv6')
    bg6 = sorted(L6)
    L21, LM21 = load_ladder('inv21')
    q21 = [s for s in sorted(LM21) if LM21[s]['top'] == {'racld'} and LM21[s]['ext'] == {-1}]
    bg21 = q21[:40]
    pr(f'\n[사다리] 밀도 8: 해부6C 배경 {len(bg6)} ({bg6[0]}~{bg6[-1]}) · 해부21 자격 {len(q21)}/{len(L21)} → 앞 40')
    A6 = lad_array(L6, bg6); A21 = lad_array(L21, bg21)
    A8 = np.vstack([A6, A21])
    pr(f'   [자기 검산] 6C 교차 {crossing(A6.mean(0))*100:.4f}% (정본 0.1827%) · Δ(0) {A6.mean(0)[0]:+.3f} (정본 +0.388) · Δ(1%) {A6.mean(0)[-1]:+.3f} (정본 −0.753)')
    pr(f'   [자기 검산] 21 교차 {crossing(A21.mean(0))*100:.4f}% (정본 0.170) · Δ(0) {A21.mean(0)[0]:+.3f} · Δ(1%) {A21.mean(0)[-1]:+.3f} (정본 +0.383 · −0.793)')
    c8 = crossing(A8.mean(0))
    bc8, _ = boot_cross(A8, 'L2|d8')
    # 층화 재추출(6C 20 · 21 40 따로) — 민감도
    r = rng('L2|d8strat')
    i6 = r.integers(0, 20, size=(NB, 20)); i21 = r.integers(0, 40, size=(NB, 40))
    bc8s = np.array([crossing(np.vstack([A6[i6[k]], A21[i21[k]]]).mean(0)) for k in range(NB)])
    pr(f'   합친 60 교차 {c8*100:.4f}% [{ci(bc8)[0]*100:.3f}, {ci(bc8)[1]*100:.3f}] (정본 0.173 [0.163, 0.185]) · 층화 재추출 [{ci(bc8s)[0]*100:.3f}, {ci(bc8s)[1]*100:.3f}]')

    m16 = A16.mean(0)
    bc16, bm16 = boot_cross(A16, 'L2|d16')
    pr(f'\n밀도 16 씨앗 rascld Δ 곡선: ' + ' · '.join(f'{mu*100:g}% {v:+.3f}' for mu, v in zip(MUV, m16)))
    d0 = ci(bm16[:, 0]); d1 = ci(bm16[:, -1])
    L1 = d0[0] > 0 and d1[1] < 0
    pr(f'L1 Δ(0) {m16[0]:+.3f} [{d0[0]:+.3f}, {d0[1]:+.3f}] · Δ(1%) {m16[-1]:+.3f} [{d1[0]:+.3f}, {d1[1]:+.3f}] → {"✅" if L1 else "❌"}')
    for j in (1, 3):
        cc = ci(bm16[:, j]); pr(f'   (서술) Δ({MUV[j]*100:g}%) {m16[j]:+.3f} [{cc[0]:+.3f}, {cc[1]:+.3f}]')
    c16 = crossing(m16)
    cens_hi = float(np.mean(np.isinf(bc16))); cens_lo = float(np.mean(bc16 == 0))
    d = bc16 - bc8
    infinf = int(np.sum(np.isinf(bc16) & np.isinf(bc8)))
    with np.errstate(invalid='ignore'):
        n_up = int(np.sum(d > 0)); n_dn = int(np.sum(d < 0)); n_tie = int(np.sum(d == 0))
    P_up, P_dn = n_up / NB, n_dn / NB
    fin = d[np.isfinite(d)]
    dci = ci(fin)
    L2 = '✅ 위로' if P_up >= 0.975 else ('❌ 아래로' if P_dn >= 0.975 else '🟡 구별 안 됨')
    pr(f'L2 교차 D16 {c16*100:.4f}% [{ci(bc16)[0]*100:.3f}, {ci(bc16)[1]*100:.3f}] 대 D8(60) {c8*100:.4f}% · d {(c16-c8)*100:+.4f}%p [{dci[0]*100:+.4f}, {dci[1]*100:+.4f}] · '
       f'P(위로) {P_up:.4f} · P(아래로) {P_dn:.4f} · 분모 {NB} · ∞−∞ {infinf} · 동률 {n_tie} · cens_hi {cens_hi} · cens_lo {cens_lo} → {L2}')
    # 층화 기준으로도
    ds = bc16 - bc8s
    pr(f'   (민감도) 밀도 8 층화 재추출: P(위로) {np.sum(ds > 0)/NB:.4f} · d [{ci(ds)[0]*100:+.4f}, {ci(ds)[1]*100:+.4f}]')
    # 난수 흐름 민감도: 다른 시드 10 개
    pus = []
    for k in range(100):
        a, _ = boot_cross(A16, f'L2|d16|rep{k}'); b, _ = boot_cross(A8, f'L2|d8|rep{k}')
        pus.append(np.sum((a - b) > 0) / NB)
    pr(f'   (민감도) 난수 흐름 100 벌 P(위로): min {min(pus):.4f} · 중앙 {np.median(pus):.4f} · max {max(pus):.4f} · ≥0.975 인 벌 {sum(p >= 0.975 for p in pus)}/100')
    band = 0.0021 <= c16 <= 0.0065
    pr(f'   R2 서술: §3 띠 [0.21%, 0.65%] 안 {"예" if band else "아니오"} · d 상한 < +0.037%p(하단 배제) {"예" if dci[1] < 0.00037 else "아니오"}')
    pr(f'   (서술) 기울기(최소제곱, 배경 평균 곡선) D16 {slope(A16):.1f} · D8(60) {slope(A8):.1f} · 6C {slope(A6):.1f}')

    # ===== 조작 확인 R4 =====
    IB8 = load_R('inv15Eib_racld')
    six = ['nracld', 'rnacld', 'rancld', 'racnld', 'raclnd', 'racldn']
    def lnR_bg(IB, seeds):
        return np.array([np.mean([math.log(IB[s]['R'][c]) for c in six]) for s in seeds])
    bg8ib = sorted(IB8)
    l16 = lnR_bg(IB16, BG); l8 = lnR_bg(IB8, bg8ib)
    dR = boot_mean(l16, 'R4|d16') - boot_mean(l8, 'R4|d8')
    dRpt = l16.mean() - l8.mean(); dRci = ci(dR)
    rule = ('§3 전제만큼 줄었다' if dRci[1] < math.log(0.5) else
            ('줄었지만 §3 전제(비 ≤ 0.5)에 못 미친다' if dRci[1] < 0 else '줄지 않았다'))
    pr(f'\nR4 조작 확인: d_R {dRpt:+.3f} [{dRci[0]:+.3f}, {dRci[1]:+.3f}] · 비 e^d_R {math.exp(dRpt):.2f} → "{rule}"')
    def g45(IB, seeds):
        return np.array([np.mean([math.log(IB[s]['R'][c]) for c in ('raclnd', 'racldn')]) -
                         np.mean([math.log(IB[s]['R'][c]) for c in ('nracld', 'rnacld', 'rancld', 'racnld')]) for s in seeds])
    for nm, IB, ss in (('D16', IB16, BG), ('D8', IB8, bg8ib)):
        g = g45(IB, ss); bb = ci(boot_mean(g, 'g45|' + nm))
        Rmed = np.median([IB[s]['R'][c] for s in ss for c in six])
        pr(f'   {nm} ln R(4) − ln R(5) {g.mean():+.3f} [{bb[0]:+.3f}, {bb[1]:+.3f}] · R 중앙 {Rmed:.2f}' + ('  (정본 +0.211 [+0.205, +0.216])' if nm == 'D8' else ''))
    pr(f'\n판정 요약: G0 문제 {len(problems)} · P1 {"✅" if p1[2] < 0 else "❌"} · P2 {n_ok}/54 · L1 {"✅" if L1 else "❌"} · L2 {L2} · P4 {read4} · R4 "{rule}"')


if __name__ == '__main__':
    main()
