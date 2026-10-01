#!/usr/bin/env python3
"""해부 24 독립 재계수 (2026-10-01)

판정기(pairs24_analyze.py · pairs14/15/20_analyze.py)와 _RESULT_pairs24_* 를 열지 않고,
사전등록 노트(해부24-사전등록-2026-10-01.md §1 · §2 · §4)와 해부 14 노트 §1 · §4 · §4b 의 정의만으로
원자료 JSON 에서 B1 · B2 · B3 · B4 · B5 · 자기 검산 · 판정 가능 조건 · G0a 를 다시 센다.

실행(lab101 · 서버에 파일을 쓰지 않음):
    ssh lab101 'cd ~/petri && python3 - --set 2' < recount24_independent.py
    ssh lab101 'cd ~/petri && python3 - --set 3' < recount24_independent.py

정의(노트에서 옮김 · 해석이 필요한 곳은 ☆ 표시):
- 끼운 계통 = bins 의 각 칸 [배경, 표지1] 중 index 1 · 팔 3,000틱 bins 합 · 배경마다 하나.
- R = (copy_ok + copy_stall + copy_del)/copy_ok · T = copy_phase/copy_ok.
- L = 고리 한 바퀴 틱 수: sim.js 의 'l'(되돌기)은 l 바로 앞에서 뒤로 첫 's'(표지)를 찾아 그 다음으로,
  없으면 0 으로 간다. L = 그 구간[target .. l] 중 틱을 쓰는 글자 수('x' 는 틱 안 씀) ·
  k = 그 구간의 'c' 수. (sim.js 의 step/copyOne 을 읽어 직접 유도 — pairs13/20 의 함수를 보지 않음.)
- 코드 = inserts(b,'n') = 밑 코드의 모든 자리에 n 하나를 넣은 것(중복 제거).
- 같은-고리 쌍 = 같은 밑 코드 안에서 L 이 같은 두 코드.
- 멸종 = 끼운 계통 series 마지막 m == 0 (stopped != -1 과 교차 확인).
- ☆ 이웃 등급 효과 = 배경마다 ln(median R of 등급 L) − ln(median R of 등급 L+1)(짧은 − 긴),
  밑 코드 공통 배경(모든 코드 생존) 평균. 해부 14 노트 §4 A1 의 '등급 중앙값 R' 문구를 따른 기본값.
  `--grade-meanln` 이면 등급 대표값 = 등급 안 코드들의 ln R 평균(민감도 · 판정기 출력 머리줄 '등급 평균 ln R' 과 대조용).
- 같은-고리 쌍 ln R 비 = 배경마다 ln(R_i/R_j), 두 코드 공통 생존 배경 평균(B2).
- B4 Q = max 같은-고리 |평균 ln R 비| ÷ min 이웃 등급 효과 — 밑 코드 공통 배경에서,
  재추출마다 공통 배경을 한꺼번에 뽑아 둘 다 다시 낸다. 판정: 95% 백분위 상한 < 0.5.
- B5 = N3 실패 쌍 (a, c) 마다 배경별 ln(T/R)_a − ln(T/R)_c, 두 코드 공통 생존 배경 평균 · 95% 구간 ⊂ [−0.04, +0.04].
- B3 = |pooled T/R_i − pooled T/R_j| ≤ 0.10. ☆ pooled T/R = Σcopy_phase / Σ시도(생존 배경 합).
- 자기 검산 = pooled T/R − L/k ∈ [−0.5, +1.5].
- 재추출 2,000 · 백분위 95% · 난수 = random.Random("recount24|<통계량 이름>") (판정기와 다른 키 → 구간 끝은 조금 다를 수 있음).
"""
import json, math, os, random, statistics, sys, hashlib
from itertools import combinations

N_BOOT = 2000
TAU, TR_PAIR, Q_MAX, D_MAX = 0.10, 0.10, 0.5, 0.04
SELF_LO, SELF_HI = -0.5, 1.5
N_MIN, DEAD_MAX = 10, 4
EXPECT = {2: (40, 56, 6), 3: (56, 84, 12)}
BASES = {2: ['racld', 'rascld', 'rsacld', 'racldx', 'acld', 'rascled'],
         3: ['rcald', 'crald', 'reasccld', 'racld', 'rascld', 'rsacld', 'racldx', 'rascled']}
KEYS = ['copy_ok', 'copy_stall', 'copy_del', 'copy_phase', 'unknown', 'stepped', 'not_stepped']


def inserts(b, ch='n'):
    out = []
    for i in range(len(b) + 1):
        s = b[:i] + ch + b[i:]
        if s not in out:
            out.append(s)
    return out


def loop_Lk(code):
    li = code.index('l')
    p = li - 1
    while p >= 0 and code[p] != 's':
        p -= 1
    tgt = p + 1 if p >= 0 else 0
    seg = code[tgt:li + 1]
    return sum(1 for ch in seg if ch != 'x'), seg.count('c')


def racld_normal_seeds(lo, hi, step, n):
    out = []
    s = lo
    while len(out) < n and (lo <= s <= hi if step > 0 else hi <= s <= lo):
        d = json.load(open(f'main/mat_d8_mu0p01_s{s:05d}.json'))
        fc = d['final_counts']
        if d['extinct_at'] == -1 and d['ticks_run'] == d['ticks_planned'] and fc and fc[0][0] == 'racld':
            out.append(s)
        s += step
    return out


def pct(xs, q):
    xs = sorted(xs)
    if not xs:
        return float('nan')
    pos = (len(xs) - 1) * q
    lo = math.floor(pos); hi = math.ceil(pos)
    return xs[lo] + (xs[hi] - xs[lo]) * (pos - lo)


def boot(fn, bgs, name):
    rng = random.Random(f'recount24|{name}')
    vals = []
    for _ in range(N_BOOT):
        smp = [bgs[rng.randrange(len(bgs))] for _ in bgs]
        vals.append(fn(smp))
    return pct(vals, 0.025), pct(vals, 0.975)


def main():
    global GRADE
    st = int(sys.argv[sys.argv.index('--set') + 1])
    GRADE = 'meanln' if '--grade-meanln' in sys.argv else 'median'
    print(f'등급 대표값 = {GRADE}')
    problems = []
    # ---- 배경
    if st == 2:
        seeds = racld_normal_seeds(26, 300, 1, 20)
        ibroot, refroot = 'inv15Eib_{}/mu0', 'inv15E_{}/mu0'
    else:
        seeds = json.load(open('_RESULT_pairs20.json'))['seeds']
        mine = sorted(racld_normal_seeds(100, 1, -1, 20))
        if mine != sorted(seeds):
            problems.append(f'세트3 배경: pairs20 seeds {seeds} ≠ 내림차순 재유도 {mine}')
        ibroot, refroot = 'inv24ib_{}/mu0', 'inv20P_{}/mu0'
    print(f'set {st} · 배경 {len(seeds)} {seeds}')
    # ---- N3 실패 쌍
    if st == 2:
        P = json.load(open('_RESULT_pairs15.json'))
        n3 = [tuple(k.split('|')[1:]) for k, v in P.items() if k.startswith('E2|') and v.get('fail')]
        negkeys = {tuple(k.split('|')[1:]) for k in P if k.startswith('E2|')}
    else:
        P = json.load(open('_RESULT_pairs20.json'))
        n3 = [tuple(k.split('|')[1:]) for k, v in P.items() if k.startswith('neg|') and v.get('N3') is False]
        negkeys = {tuple(k.split('|')[1:]) for k in P if k.startswith('neg|')}
        if P['pass'] != BASES[3]:
            problems.append(f"pairs20 pass {P['pass']} ≠ {BASES[3]}")
    print(f'N3 실패 쌍 {len(n3)}: ' + ' · '.join(f'{b}:{a}/{c}' for b, a, c in n3))

    # ---- 원자료
    R, T, A, PH, dead, g0a = {}, {}, {}, {}, {}, [0, 0, 0, 0]  # g0a: 대조 팔, 불일치, 기준없음, 체크섬다름
    codes_of, Lk = {}, {}
    for b in BASES[st]:
        codes = inserts(b)
        codes_of[b] = codes
        for c in codes:
            Lk[c] = loop_Lk(c)
            R[c], T[c], A[c], PH[c], dead[c] = {}, {}, {}, {}, []
        for s in seeds:
            f = os.path.join(ibroot.format(b), f'ib_inv_s{s:05d}_mu0.json')
            d = json.load(open(f))
            chk = [(d['seed'] == s, 'seed'), (d['mode'] == 'inv', 'mode'), (d['mu'] == 0, 'mu'),
                   (d['wt'] == 'racld', 'wt'), (d['ticks'] == 3000, 'ticks'), (d['checksum_match'] is True, 'checksum'),
                   (d['codes'] == codes, 'codes'), (len(d['arms']) == 1 + len(codes), 'arms'),
                   (d['invbud_version'] == '1.0.0', 'invbud_ver'), (d['sim_version'] == '0.3.5', 'sim_ver')]
            for k in ('by_letter', 'density', 'grow_mu', 'ancestor'):
                chk.append((k not in d, f'opt {k}'))
            for ok, nm in chk:
                if not ok:
                    problems.append(f'{f}: {nm}')
            # G0a 기준
            rf = os.path.join(refroot.format(b), f'dms_mat_{b}_s{s:05d}_list.json')
            ref = None
            if os.path.exists(rf):
                ref = json.load(open(rf))
                if ref['grow_checksum'] != d['grow_checksum']:
                    g0a[3] += 1
            refarms = {}
            if ref:
                for a in ref['arms']:
                    refarms[(a['kind'], a['key'])] = a
            for a in d['arms']:
                tot = {k: sum(bn[k][1] for bn in a['bins']) for k in KEYS}
                if a['kind'] == 'ref':
                    key = ('ref', a['key'])
                else:
                    key = ('mut', a['key'])
                # G0a: 끼운 팔 = dms 의 같은 key 팔 · ref 는 key 가 같을 때만
                cand = None
                if ref is None:
                    if a['kind'] != 'ref':
                        g0a[2] += 1
                else:
                    if a['kind'] == 'ref':
                        cand = refarms.get(('ref', a['key']))
                    else:
                        cand = next((v for (kk, kk2), v in refarms.items() if kk != 'ref' and kk2 == a['key']), None)
                        if cand is None:
                            g0a[2] += 1
                    if cand is not None:
                        g0a[0] += 1
                        mine = [list(x[:3]) for x in a['series']]
                        theirs = [list(x[:3]) for x in cand['series']]
                        if mine != theirs or a['stopped'] != cand['stopped']:
                            g0a[1] += 1
                            problems.append(f'G0a 불일치 {b} s{s} {a["key"]}')
                if a['kind'] == 'ref':
                    continue
                c = a['key']
                if tot['copy_del'] != 0:
                    problems.append(f'copy_del≠0 {c} s{s}')
                denom = tot['stepped'] + tot['not_stepped'] + tot['unknown']
                if denom and tot['unknown'] / denom > 0.01:
                    problems.append(f'unknown>1% {c} s{s}')
                ext = a['series'][-1][2] == 0
                if ext != (a['stopped'] != -1):
                    problems.append(f'멸종 판정 두 기준 다름 {c} s{s} m={a["series"][-1][2]} stopped={a["stopped"]}')
                if ext:
                    dead[c].append(s)
                    continue
                if tot['copy_ok'] == 0:
                    problems.append(f'copy_ok 0 {c} s{s}')
                    continue
                att = tot['copy_ok'] + tot['copy_stall'] + tot['copy_del']
                R[c][s] = att / tot['copy_ok']
                T[c][s] = tot['copy_phase'] / tot['copy_ok']
                A[c][s] = att
                PH[c][s] = tot['copy_phase']

    allcodes = [c for b in BASES[st] for c in codes_of[b]]
    pairs = [(b, i, j) for b in BASES[st] for i, j in combinations(codes_of[b], 2) if Lk[i][0] == Lk[j][0]]
    print(f'코드 {len(allcodes)} · 같은-고리 쌍 {len(pairs)} · N3 {len(n3)} (기대 {EXPECT[st]})')
    mypairset = {(b, *sorted((i, j))) for b, i, j in pairs}
    nk = {(b, *sorted((i, j))) for b, i, j in negkeys}
    print(f'같은-고리 쌍 = Δ 파일 음성 쌍 목록: {mypairset == nk} (Δ 파일 {len(nk)})')
    print(f'G0a: 대조 팔 {g0a[0]} · 불일치 {g0a[1]} · 기준 없음 {g0a[2]} · grow_checksum 다른 배경 {g0a[3]}')
    print('멸종 팔: ' + (' · '.join(f'{c} {v}' for c, v in dead.items() if v) or '0'))

    # ---- 판정 가능 조건
    viol = []
    if (len(allcodes), len(pairs), len(n3)) != EXPECT[st]:
        viol.append('개수')
    for c in allcodes:
        if len(dead[c]) > DEAD_MAX:
            viol.append(f'멸종 {c} {len(dead[c])}')
    common = {b: [s for s in seeds if all(s in R[c] for c in codes_of[b])] for b in BASES[st]}
    for b in BASES[st]:
        if len(common[b]) < N_MIN:
            viol.append(f'공통배경 {b} {len(common[b])}')
    def pbg(i, j):
        return [s for s in seeds if s in R[i] and s in R[j]]
    for b, i, j in pairs:
        if len(pbg(i, j)) < N_MIN:
            viol.append(f'쌍배경 {i}/{j}')
    for b, a, c in n3:
        if a not in R or c not in R:
            viol.append(f'N3 코드 없음 {a}/{c}')
        elif len(pbg(a, c)) < N_MIN:
            viol.append(f'N3배경 {a}/{c}')
    print(f'점검 문제 {len(problems)}' + ('' if not problems else ': ' + ' | '.join(problems[:20])))
    print(f'판정 가능 조건 위반 {len(viol)}' + ('' if not viol else ': ' + ' | '.join(viol)))

    # ---- 자기 검산
    sc = {}
    for c in allcodes:
        tr = sum(PH[c].values()) / sum(A[c].values())
        L, k = Lk[c]
        sc[c] = (tr, tr - L / k, L, k)
    lo_c = min(sc, key=lambda c: sc[c][1]); hi_c = max(sc, key=lambda c: sc[c][1])
    for c in sorted(sc, key=lambda c: sc[c][1])[:3]:
        print(f'  하위 {c} T/R−L/k {sc[c][1]:+.5f}')
    self_ok = all(SELF_LO <= v[1] <= SELF_HI for v in sc.values())
    print(f'자기 검산 {"✅" if self_ok else "❌"} {sum(SELF_LO <= v[1] <= SELF_HI for v in sc.values())}/{len(sc)} · '
          f'최소 {sc[lo_c][1]:+.3f}({lo_c}) · 최대 {sc[hi_c][1]:+.3f}({hi_c})')
    k2 = [c for c in allcodes if Lk[c][1] == 2]
    if k2:
        print(f'  k=2 {len(k2)} 코드 T/R−L/k {min(sc[c][1] for c in k2):+.3f} ~ {max(sc[c][1] for c in k2):+.3f}')

    # ---- 함수
    def pair_mean(i, j, bgs):
        return statistics.fmean(math.log(R[i][s] / R[j][s]) for s in bgs)

    def grades(b):
        gs = sorted({Lk[c][0] for c in codes_of[b]})
        return gs, [(gs[t], gs[t + 1]) for t in range(len(gs) - 1) if gs[t + 1] == gs[t] + 1]

    def grade_eff(b, L1, L2, bgs):
        c1 = [c for c in codes_of[b] if Lk[c][0] == L1]
        c2 = [c for c in codes_of[b] if Lk[c][0] == L2]
        if GRADE == 'median':
            return statistics.fmean(math.log(statistics.median(R[c][s] for c in c1)) -
                                    math.log(statistics.median(R[c][s] for c in c2)) for s in bgs)
        return statistics.fmean(statistics.fmean(math.log(R[c][s]) for c in c1) -
                                statistics.fmean(math.log(R[c][s]) for c in c2) for s in bgs)

    # ---- B1
    print('\nB1 이웃 등급 ln R 차(짧은 − 긴) · 하한 > 0')
    b1_all = True
    for b in BASES[st]:
        gs, adj = grades(b)
        bg = common[b]
        for L1, L2 in adj:
            e = grade_eff(b, L1, L2, bg)
            lo, hi = boot(lambda x: grade_eff(b, L1, L2, x), bg, f'B1|{b}|{L1}')
            ok = lo > 0; b1_all &= ok
            print(f'  {b:9s} L{L1}/L{L2} 배경 {len(bg)}: {e:+.3f} [{lo:+.3f}, {hi:+.3f}] {"✅" if ok else "❌"}')
        if not adj:
            print(f'  {b} 이웃 등급 없음 {gs}'); b1_all = False

    # ---- B2 · B3
    vals = []
    for b, i, j in pairs:
        bg = pbg(i, j)
        m = pair_mean(i, j, bg)
        tri, trj = sc[i][0], sc[j][0]
        vals.append((abs(m), m, abs(tri - trj), b, i, j, len(bg)))
    absr = sorted(v[0] for v in vals)
    mx = max(vals, key=lambda v: v[0]); mx3 = max(vals, key=lambda v: v[2])
    b2 = all(v[0] <= TAU for v in vals); b3 = all(v[2] <= TR_PAIR for v in vals)
    lo, hi = boot(lambda x: pair_mean(mx[4], mx[5], x), pbg(mx[4], mx[5]), f'B2max|{mx[4]}|{mx[5]}')
    print(f'\nB2 {"✅" if b2 else "❌"} {sum(v[0] <= TAU for v in vals)}/{len(vals)} · 중앙값 {statistics.median(absr):.4f} · '
          f'90분위 {pct(absr, 0.9):.4f} · 최대 {mx[0]:.4f} ({mx[4]}/{mx[5]} {mx[1]:+.4f} [{lo:+.3f}, {hi:+.3f}])')
    print(f'B3 {"✅" if b3 else "❌"} {sum(v[2] <= TR_PAIR for v in vals)}/{len(vals)} · 최대 {mx3[2]:.4f} ({mx3[4]}/{mx3[5]})')

    # ---- B4
    print('\nB4 Q = max 같은-고리 |ln R 비| ÷ min 이웃 등급 효과 · 상한 < 0.5')
    b4_all = True
    for b in BASES[st]:
        bg = common[b]
        bp = [(i, j) for bb, i, j in pairs if bb == b]
        _, adj = grades(b)
        def Q(x):
            return max(abs(pair_mean(i, j, x)) for i, j in bp) / min(grade_eff(b, L1, L2, x) for L1, L2 in adj)
        q = Q(bg); lo, hi = boot(Q, bg, f'B4|{b}')
        ok = hi < Q_MAX; b4_all &= ok
        print(f'  {b:9s} 배경 {len(bg)} 쌍 {len(bp)}: Q {q:.3f} [{lo:.3f}, {hi:.3f}] {"✅" if ok else "❌"}')

    # ---- B5
    print('\nB5 N3 실패 쌍 ln(T/R)_a − ln(T/R)_c · 구간 ⊂ [−0.04, +0.04]')
    b5_all = True
    def d5(a, c, x):
        return statistics.fmean(math.log(T[a][s] / R[a][s]) - math.log(T[c][s] / R[c][s]) for s in x)
    for b, a, c in n3:
        bg = pbg(a, c)
        m = d5(a, c, bg); lo, hi = boot(lambda x: d5(a, c, x), bg, f'B5|{a}|{c}')
        ok = -D_MAX <= lo and hi <= D_MAX; b5_all &= ok
        lnr = pair_mean(a, c, bg)
        print(f'  {b:9s} {a}/{c} 배경 {len(bg)}: {m:+.4f} [{lo:+.4f}, {hi:+.4f}] {"✅" if ok else "❌"} (ln R 비 {lnr:+.4f})')

    core = b2 and b4_all and b5_all
    print(f'\n종합 세트 {st}: 판정 가능 {"예" if not viol and not problems and self_ok else "아니오"} · '
          f'B2 {"✅" if b2 else "❌"} · B4 {"✅" if b4_all else "❌"} · B5 {"✅" if b5_all else "❌"} → 핵심 {"✅" if core else "❌"} · '
          f'B1 {"✅" if b1_all else "❌"} · B3 {"✅" if b3 else "❌"}')


if __name__ == '__main__':
    main()
