#!/usr/bin/env node
/* 코드 배양접시 — 이긴 코드 해부(가상 DMS) 계측기. 판정은 하지 않는다(dms_analyze.py).
 *
 * 1. 배경: 본실험 접시(조건 · 시드)를 같은 설정으로 50,000틱까지 다시 키운다 → 본실험 기록의 체크섬과 대조.
 * 2. 팔(arm)마다 그 끝 상태를 복제해, 살아 있는 개체의 frac 을 한 코드로 바꿔 끼우고 계통 표지 tag=1 을 붙인 뒤
 *    돌연변이 없이(μ=0) ticks 틱 더 돌려 표지 계통의 개체 수를 every 틱마다 센다.
 *    바꿔 끼울 개체는 배경 시드로만 정한다 → 한 배경 안의 모든 팔이 같은 자리에 끼운다(짝 설계).
 * 3. 팔 종류
 *    ref       WT 를 끼운다 · 뒤 난수 A
 *    neutral   WT 를 끼운다 · 뒤 난수 B1..B10 (가짜 돌연변이 — 잡음 기준)
 *    mut       한 글자 돌연변이(바뀜 · 빠뜨림 · 끼어듦 전부, 같은 결과 코드는 하나로) · 뒤 난수 A
 *    표지 계통이 사라지면(돌연변이가 없어 되살아나지 않는다) 그 팔은 거기서 멈춘다.
 * 4. 계측기 점검(프로세스마다): 작은 접시로 복제 충실도 — 원본과 복제본을 같은 난수로 500틱 돌려 전체 상태 비교.
 *
 * 사용: node dms.js --cond space --seed 3 --wt reasccld [--frac 0.1] [--ticks 3000] [--every 250]
 *                   [--grow 50000] [--arms all|pilot] [--main-dir ../main] [--out dir]
 */
'use strict';
const fs = require('fs');
const path = require('path');
const Sim = require('./sim.js');

const args = {};
for (let i = 2; i < process.argv.length; i += 2) args[process.argv[i].replace(/^--/, '')] = process.argv[i + 1];
const COND = {
  space: { mat: false, energy: false },
  mat: { mat: true, energy: false },
  energy: { mat: false, energy: true }
};
const cond = args.cond;
if (!COND[cond]) { console.error(`알 수 없는 cond: ${cond}`); process.exit(2); }
if (!['0.3.1', '0.3.2', '0.3.3', '0.3.4', '0.3.5'].includes(Sim.VERSION)) { console.error(`sim ${Sim.VERSION} — 0.3.1~0.3.5 가 필요하다`); process.exit(2); }
const seed = parseInt(args.seed, 10);
const WT = args.wt;
const frac = parseFloat(args.frac || '0.1');
const T = parseInt(args.ticks || '3000', 10);
const every = parseInt(args.every || '250', 10);
const grow = parseInt(args.grow || '50000', 10);
const armsMode = args.arms || 'all';
// 해부 2A(09-15 추가): 팔을 돌리는 동안의 복제 오류율. 기본 0 이면 해부 1 과 같은 동작.
const muAssay = parseFloat(args['mu-assay'] || '0');
// 해부 3(09-16 추가): --arms list 면 --codes 로 준 코드만 팔로 돌린다. 기본값(all · pilot)은 그대로.
const codes = (args.codes || '').split(',').map(x => x.trim()).filter(Boolean);
// 해부 6(09-16 추가): --exact 1 이면 기록 줄 끝에 '표지 계통 안에서 끼운 코드 그대로인 개체' · '나머지 중 WT 그대로인 개체' 두 칸을 더한다. 기본값(끔)이면 출력이 전과 바이트 단위로 같다.
const exact = args.exact === '1';
// 해부 7(09-16 추가): --find-first 1 이면 배경은 원래 규칙으로 키우고(체크섬 대조 그대로) **팔을 돌리는 동안만**
// sim v0.3.3 의 재료 먼저 규칙을 쓴다. 기본값(끔)이면 출력이 전과 같다.
const findFirst = args['find-first'] === '1';
if (findFirst && !['0.3.3', '0.3.4', '0.3.5'].includes(Sim.VERSION)) { console.error('--find-first 는 sim 0.3.3~0.3.5 가 필요하다'); process.exit(2); }
// 해부 9(09-17 추가): --remember-die 1 이면 팔을 돌리는 동안만 sim v0.3.4 의 기억하는 주사위 규칙을 쓴다. 기본값(끔)이면 출력이 전과 같다.
const rememberDie = args['remember-die'] === '1';
if (rememberDie && !['0.3.4', '0.3.5'].includes(Sim.VERSION)) { console.error('--remember-die 는 sim 0.3.4/0.3.5 가 필요하다'); process.exit(2); }
// 해부 10(09-18 추가): --cosmic 0 이면 팔 동안만 우주선 돌연변이를 끈다(sim 0.3.5). --ancestor CODE 면 배경을 그 시조로 키운다(본실험 기록이 없어 체크섬 대조는 'no-main').
const cosmicOff = args.cosmic === '0';
if (cosmicOff && Sim.VERSION !== '0.3.5') { console.error('--cosmic 0 은 sim 0.3.5 가 필요하다'); process.exit(2); }
const ancestor = args.ancestor || null;
// 해부 18(09-28 추가): --density D 면 배경을 그 밀도로 키운다(기본 8 이면 출력 전과 같음). 8 이 아니면 main 기록이 없어 체크섬은 'no-main' · 파일 이름 끝 _dD · 출력 density.
const densityV = parseFloat(args.density || '8');
const dtag = densityV === 8 ? '' : '_d' + String(densityV).replace('.', 'p');
// 해부 18(09-28 추가): --grow-mu M 이면 배경을 그 오류율로 키운다(기본 0.01 이면 출력 전과 같음). 0.01 이 아니면 파일 이름 끝 _gmuM · 출력 grow_mu · main 대조 없음.
const growMu = parseFloat(args['grow-mu'] || '0.01');
const gtag = growMu === 0.01 ? '' : '_gmu' + String(growMu).replace('.', 'p');
// 해부 11(09-20 추가): --age0 N --age-var N 이면 **배경과 팔 모두** 그 수명으로 돈다(sim 0.3.5 · 세계 전체 옵션). 기본값(안 줌)이면 출력이 전과 같다.
const age0V = args.age0 ? parseInt(args.age0, 10) : null;
const ageVarV = args['age-var'] ? parseInt(args['age-var'], 10) : null;
if ((age0V !== null || ageVarV !== null) && Sim.VERSION !== '0.3.5') { console.error('--age0 은 sim 0.3.5 가 필요하다'); process.exit(2); }
// 해부 9 B(09-17 추가): --no-mat 1 이면 팔을 돌리는 동안만 재료 규칙을 끈다(w.o.mat = false · w.mat = null → 헛손질 0 · 거름 0).
// 배경은 그대로 재료 세계에서 키운다(체크섬 대조). 재료 보존 점검은 뜻이 없어 참으로 둔다. 기본값(끔)이면 출력이 전과 같다.
const noMat = args['no-mat'] === '1';
// 해부 12(09-20 추가): --kids 1 이면 팔을 돌리는 동안 **표지 계통의 출생을 밖에서 읽어** 부모가 끼운 코드 그대로일 때의 자식 코드 분포를 적는다
// (틱 전에 '자식 칸 → 부모 코드' 를 적어 두고 틱 뒤에 새 id 를 그 칸에서 찾는다 · 읽기만 — 궤적 · 옛 필드 동일).
// 복사 없이 코드가 바뀐 것(우주선)도 따로 센다. --arms nbr --nbr-of CODE 면 CODE 의 한 글자 이웃 전부와 CODE 자신을 끼운다.
const kidsOn = args.kids === '1';
const nbrOf = args['nbr-of'] || null;
if (armsMode === 'nbr' && !nbrOf) { console.error('--arms nbr 에는 --nbr-of CODE 가 필요하다'); process.exit(2); }
const mainDir = args['main-dir'] || path.join(__dirname, 'main');
const outDir = args.out || '.';
const N_NEUTRAL = 10;
const { LET, K, DISH, NB } = Sim;
const enc = s => [...s].map(ch => { const k = LET.indexOf(ch); if (k < 0) throw new Error(`글자 아님: ${ch}`); return k; });
const opts = Object.assign({ seed, density: densityV, mu: growMu, light: 1 }, COND[cond]);
if (ancestor) opts.ancestor = ancestor;
if (age0V !== null) opts.age0 = age0V;
if (ageVarV !== null) opts.ageVar = ageVarV;
const INJ_SEED = (seed * 1000003 + 17) >>> 0;
const POST_A = (seed * 7919 + 1) >>> 0;
const t0 = Date.now();

/* ---------- 돌연변이 목록 ---------- */
function mutantList(wt) {
  const out = new Map();
  const add = (key, label) => {
    if (key === wt || key.length < Sim.MINLEN || key.length > Sim.MAXLEN) return;
    if (!out.has(key)) out.set(key, []);
    out.get(key).push(label);
  };
  for (let p = 0; p < wt.length; p++) for (const ch of LET) if (ch !== wt[p]) add(wt.slice(0, p) + ch + wt.slice(p + 1), `sub:${p}:${ch}`);
  for (let p = 0; p < wt.length; p++) add(wt.slice(0, p) + wt.slice(p + 1), `del:${p}`);
  for (let p = 0; p <= wt.length; p++) for (const ch of LET) add(wt.slice(0, p) + ch + wt.slice(p), `ins:${p}:${ch}`);
  return [...out.entries()];
}

/* ---------- 세계 복제 ---------- */
function cloneWorld(w, rngSeed) {
  if (w.born.length) throw new Error('born 이 비어 있지 않은 시점에 복제하려 했다');
  const orgs = new Map();
  for (const [id, o] of w.orgs) {
    orgs.set(id, Object.assign({}, o, { g: o.g.slice(), child: o.child ? { i: o.child.i, g: o.child.g.slice() } : null }));
  }
  const alive = w.alive.map(o => orgs.get(o.id));
  if (alive.length !== orgs.size || alive.some(o => !o)) throw new Error('alive 와 orgs 가 어긋난다');
  const over = {
    o: Object.assign({}, w.o), rnd: Sim.mulberry32(rngSeed),
    cell: w.cell.slice(), hueAt: w.hueAt.slice(), mat: w.mat ? w.mat.slice() : null,
    orgs, alive, born: [],
    hist: { pop: w.hist.pop.slice(), kinds: w.hist.kinds.slice(), len: w.hist.len.slice() }, summary: null
  };
  for (const [k, v] of Object.entries(w)) {
    if (v !== null && typeof v === 'object' && !(k in over)) throw new Error(`복제에서 빠진 참조 필드: ${k}`);
  }
  return Object.assign({}, w, over);
}

function fullState(w) {
  const orgs = [...w.orgs.values()].sort((a, b) => a.id - b.id)
    .map(o => [o.id, o.i, o.key, o.ip, o.ctx, o.rh, o.face, o.age, o.maxAge, +o.E.toFixed(6), o.bornAt, o.last, o.lastHost, o.tag,
      o.child ? [o.child.i, o.child.g.join('')] : null]);
  return JSON.stringify([w.tick, w.nextId, w.births, w.deaths, w.alive.map(o => o.id), orgs,
    Array.from(w.cell), w.mat ? Array.from(w.mat) : null]);
}

function cloneFidelityTest() {
  const w = Sim.makeWorld(Object.assign({}, opts, { seed: 424242 }));
  for (let t = 0; t < 1500 && w.extinctAt < 0; t++) Sim.tick(w);
  const c = cloneWorld(w, 777);
  w.rnd = Sim.mulberry32(777);
  for (let t = 0; t < 500; t++) { Sim.tick(w); Sim.tick(c); }
  return { pop: w.alive.length, same: fullState(w) === fullState(c), checksum: [Sim.checksum(w), Sim.checksum(c)] };
}

/* ---------- 끼워 넣기 ---------- */
function takeNear(w, i, k, rnd) {
  const m = w.mat;
  if (m[i * K + k] > 0) { m[i * K + k]--; return i; }
  for (let d = 0; d < 4; d++) {
    const t = NB[i * 4 + d];
    if (t >= 0 && m[t * K + k] > 0) { m[t * K + k]--; return t; }
  }
  for (let q = 0; q < 200000; q++) {
    const t = DISH[(rnd() * DISH.length) | 0];
    if (m[t * K + k] > 0) { m[t * K + k]--; return t; }
  }
  return -1;
}

function inject(w, G) {
  const rnd = Sim.mulberry32(INJ_SEED);
  const list = w.alive.slice().sort((a, b) => a.id - b.id);
  const n = list.length, target = Math.round(n * frac);
  for (let q = 0; q < target; q++) {
    const j = q + ((rnd() * (n - q)) | 0);
    const tmp = list[q]; list[q] = list[j]; list[j] = tmp;
  }
  const key = Sim.keyOf(G), hue = Sim.hueOf(key);
  let done = 0, skipped = 0;
  for (let q = 0; q < target; q++) {
    const org = list[q];
    if (org.child) {
      if (w.mat) for (const k of org.child.g) w.mat[org.child.i * K + k]++;
      w.cell[org.child.i] = -1; org.child = null;
    }
    if (w.mat) {
      for (const k of org.g) w.mat[org.i * K + k]++;
      const taken = [];
      let ok = true;
      for (const k of G) { const at = takeNear(w, org.i, k, rnd); if (at < 0) { ok = false; break; } taken.push([at, k]); }
      if (!ok) {
        for (const [at, k] of taken) w.mat[at * K + k]++;
        for (const k of org.g) w.mat[org.i * K + k]--;
        skipped++;
        continue;
      }
    }
    org.g = G.slice(); org.key = key; org.hue = hue; w.hueAt[org.i] = hue;
    org.ip = 0; org.ctx = org.id; org.rh = 0; org.tag = 1;
    done++;
  }
  return { n, target, done, skipped };
}

/* ---------- 팔 하나 ---------- */
function runArm(base, kind, key, labels, postSeed) {
  const ta = Date.now();
  const w = cloneWorld(base, postSeed);
  w.o.mu = muAssay;
  if (findFirst) w.o.findFirst = true;
  if (rememberDie) w.o.rememberDie = true;
  if (noMat) { w.o.mat = false; w.mat = null; }
  if (cosmicOff) w.o.cosmic = false;
  const inj = inject(w, enc(key));
  const consAfterInj = w.mat ? Sim.countMat(w) === w.matTotal0 : true;
  const count = exact
    ? () => {
      let m = 0, mx = 0, rx = 0;
      for (const o of w.alive) {
        if (o.tag === 1) { m++; if (o.key === key) mx++; } else if (o.key === WT) rx++;
      }
      return [w.alive.length, m, mx, rx];
    }
    : () => { let m = 0; for (const o of w.alive) if (o.tag === 1) m++; return [w.alive.length, m]; };
  const series = [[0, ...count()]];
  let stopped = -1;
  const KD = kidsOn ? { b_exact: 0, b_other: 0, unknown: 0, kids: new Map(), cosmic: new Map() } : null;
  const bump = (m, k) => m.set(k, (m.get(k) || 0) + 1);
  for (let t = 1; t <= T; t++) {
    let pre = null, id0 = 0, exactList = null;
    if (kidsOn) {
      pre = new Map(); id0 = w.nextId; exactList = [];
      for (const o of w.alive) if (o.tag === 1) { if (o.child) pre.set(o.child.i, o.key); if (o.key === key) exactList.push(o); }
    }
    Sim.tick(w);
    if (kidsOn) {
      for (const o of w.alive) if (o.id >= id0 && o.tag === 1) {
        const pk = pre.get(o.i);
        if (pk === undefined) KD.unknown++;
        else if (pk === key) { KD.b_exact++; if (o.key !== key) bump(KD.kids, o.key); }
        else KD.b_other++;
      }
      for (const o of exactList) if (!o.dead && o.key !== key) bump(KD.cosmic, o.key);
    }
    if (t % every === 0 || w.extinctAt >= 0) {
      const c = count();
      series.push([t, ...c]);
      if (c[1] === 0 || w.extinctAt >= 0) { stopped = t; break; }
    }
  }
  const consEnd = w.mat ? Sim.countMat(w) === w.matTotal0 : true;
  const ret = { kind, key, labels, post_seed: postSeed, inj, cons_ok: consAfterInj && consEnd, stopped, series, elapsed_ms: Date.now() - ta };
  if (kidsOn) ret.kids = { b_exact: KD.b_exact, b_other: KD.b_other, unknown: KD.unknown, kids: [...KD.kids.entries()].sort((a, b) => b[1] - a[1]), cosmic: [...KD.cosmic.entries()].sort((a, b) => b[1] - a[1]) };
  return ret;
}

/* ---------- 실행 ---------- */
const fidelity = cloneFidelityTest();
const w = Sim.makeWorld(opts);
while (w.tick < grow && w.extinctAt < 0) Sim.tick(w);
const growChecksum = Sim.checksum(w);
let mainChecksum = null, mainTop = null;
const mainPath = path.join(mainDir, `${cond}_d${String(densityV).replace('.', 'p')}_mu0p01_s${String(seed).padStart(5, '0')}.json`);
if (grow === 50000 && growMu === 0.01 && fs.existsSync(mainPath)) {
  const z = JSON.parse(fs.readFileSync(mainPath, 'utf8'));
  mainChecksum = z.checksum; mainTop = z.final_counts.length ? z.final_counts[0][0] : null;
}
const counts = new Map();
for (const o of w.alive) counts.set(o.key, (counts.get(o.key) || 0) + 1);
const top = [...counts.entries()].sort((a, b) => b[1] - a[1] || (a[0] < b[0] ? -1 : 1));
const growMs = Date.now() - t0;

const arms = [];
const muts = mutantList(WT);
if (w.extinctAt < 0) {
  arms.push(runArm(w, 'ref', WT, ['wt'], POST_A));
  for (let b = 1; b <= N_NEUTRAL; b++) arms.push(runArm(w, 'neutral', WT, [`neutral:${b}`], (POST_A + b * 104729) >>> 0));
  const chosen = armsMode === 'pilot'
    ? muts.filter(([, labels]) => labels.includes(`del:${WT.lastIndexOf('d')}`))
    : armsMode === 'list' ? muts.filter(([key]) => codes.includes(key))
    : armsMode === 'nbr' ? [[nbrOf, ['self']], ...mutantList(nbrOf).filter(([key]) => key !== WT)]
    : muts;
  if (armsMode === 'list' && chosen.length !== codes.length) {
    console.error(`--codes 중 ${codes.length - chosen.length}개가 ${WT} 의 한 글자 이웃이 아니다 — 중단`);
    process.exit(2);
  }
  for (const [key, labels] of chosen) arms.push(runArm(w, 'mut', key, labels, POST_A));
}

const out = {
  sim_version: Sim.VERSION, cond, seed, wt: WT, opts, frac, ticks: T, every, grow, arms_mode: armsMode, mu_assay: muAssay, codes,
  inj_seed: INJ_SEED, post_a: POST_A, n_mutants_total: muts.length,
  fidelity, grow_checksum: growChecksum, main_checksum: mainChecksum, checksum_match: mainChecksum === growChecksum,
  main_top: mainTop, grow_top: top.slice(0, 5), grow_pop: w.alive.length, grow_extinct: w.extinctAt,
  grow_ms: growMs, elapsed_ms: Date.now() - t0, arms
};
fs.mkdirSync(outDir, { recursive: true });
if (exact) out.exact_series = ['t', 'pop', 'tag1', 'tag1_exact_injected', 'rest_exact_wt'];
if (findFirst) out.find_first = true;
if (rememberDie) out.remember_die = true;
if (noMat) out.no_mat = true;
if (cosmicOff) out.cosmic_off = true;
if (ancestor) out.ancestor = ancestor;
if (densityV !== 8) out.density = densityV;
if (growMu !== 0.01) out.grow_mu = growMu;
if (kidsOn) out.kids_on = true;
if (nbrOf) out.nbr_of = nbrOf;
if (age0V !== null) out.age0 = age0V;
if (ageVarV !== null) out.age_var = ageVarV;
const tag = `dms_${cond}_${WT}_s${String(seed).padStart(5, '0')}_${armsMode}${nbrOf ? '_' + nbrOf : ''}${dtag}${gtag}`;
fs.writeFileSync(path.join(outDir, `${tag}.json`), JSON.stringify(out));
const ref = arms.find(a => a.kind === 'ref');
console.log([
  `DONE ${tag}`, `v${Sim.VERSION}`, `fidelity ${fidelity.same ? 'OK' : 'FAIL'}`,
  `checksum ${mainChecksum === null ? 'no-main' : (out.checksum_match ? 'MATCH' : 'MISMATCH')}`,
  `top ${top.length ? top[0][0] : '-'}`, `pop ${w.alive.length}`, `arms ${arms.length}`,
  ref ? `ref inj ${ref.inj.done}/${ref.inj.target} skip ${ref.inj.skipped}` : 'no-ref',
  `cons ${arms.every(a => a.cons_ok) ? 'OK' : 'FAIL'}`,
  `grow ${(growMs / 1000).toFixed(0)}s`, `total ${((Date.now() - t0) / 1000).toFixed(0)}s`
].join(' · '));
