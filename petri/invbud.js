#!/usr/bin/env node
/* 코드 배양접시 — 계통별 시간 예산 계측기(침입 팔 · 단일 코드 접시). 판정은 하지 않는다(invbud_analyze.py).
 *
 * 해부 7 ③b: 오류가 없을 때(μ 0) 씨앗이 찬 접시 침입에서 이기는 이유를 **잰다**.
 *   침입 모드: dms.js 와 똑같이 배경을 다시 키워(본실험 체크섬 대조) 살아 있는 개체의 10% 에 코드를 끼우고(표지 1)
 *              3,000틱 돌린다. 팔 = ref(WT 를 끼움) + --codes. 난수 · 끼우기가 dms.js 와 같아 **계통 개체 수 기록이
 *              해부 1 의 같은 팔과 같아야 한다**(회귀 검사).
 *   단일 모드: --mono CODE — 그 코드 하나로 빈 접시(mat 밀도 8)를 키운다. 섞이지 않았을 때의 기준.
 * 칸(bin)마다 **표지별로** 센다: 개체-틱 · 출생(자식의 표지 = 부모의 표지) · 사망 · 복사 성공/헛손질/빠뜨림/할 일 없음 ·
 *   복사 중인 틱(자식이 있고 다 못 베낀 상태) · 자리잡기 성공/막힘/할 일 없음/모름 · 나누기.
 *   글자 하나 얻는 데 걸린 시간 = 복사 중인 틱 / 복사 성공.
 *   사망은 원인별로도 센다: 나이(수명 넘음) · 차례 전(자기 차례가 오기 전에 죽음 = 먹힘) · 그 밖(무작위 · 차례 뒤 먹힘).
 *   칸마다 표지별 나이 합도 적는다(평균 나이 = 나이 합 / 개체-틱).
 * sim.js 는 바꾸지 않는다 — 틱 사이에 읽기만 한다.
 *
 * 사용: node invbud.js --seed 3 [--codes rascld,rsacld,racldx] [--mu 0] [--ticks 3000] [--bin 250] [--main-dir main] [--out dir]
 *       node invbud.js --mono rascld --seed 801 [--ticks 20000] [--bin 500] [--out dir]
 */
'use strict';
const fs = require('fs');
const path = require('path');
const Sim = require('./sim.js');
const INVBUD_VERSION = '1.0.0';
// 해부 14 후속(09-28): --by-letter 1 이면 헛손질을 **요청한 글자 종류별**로 센다(칸마다 stall_by[표지][글자 k]). 읽기만 · 기본값(끔)이면 출력이 전과 같다.

const args = {};
for (let i = 2; i < process.argv.length; i += 2) args[process.argv[i].replace(/^--/, '')] = process.argv[i + 1];
if (!['0.3.2', '0.3.3', '0.3.4', '0.3.5'].includes(Sim.VERSION)) { console.error(`sim ${Sim.VERSION} — 0.3.2~0.3.5 가 필요하다`); process.exit(2); }
const seed = parseInt(args.seed, 10);
const MONO = args.mono || null;
const WT = 'racld';
const codes = (args.codes || 'rascld,rsacld,racldx').split(',').map(x => x.trim()).filter(Boolean);
const mu = parseFloat(args.mu || '0');
const byLetter = args['by-letter'] === '1';
// 해부 18(09-28 추가): --density D (기본 8 이면 출력 전과 같음) · 8 이 아니면 main 대조 없음 · 파일 이름 끝 _dD · 출력 density.
const densityV = parseFloat(args.density || '8');
const dtag = densityV === 8 ? '' : '_d' + String(densityV).replace('.', 'p');
// 해부 18(09-28 추가): --grow-mu M · --ancestor CODE (기본 0.01 · 없음 이면 출력 전과 같음) · 0.01 이 아니면 파일 이름 끝 _gmuM · main 대조 없음.
const growMu = parseFloat(args['grow-mu'] || '0.01');
const ancestorV = args.ancestor || null;
const gtag = growMu === 0.01 ? '' : '_gmu' + String(growMu).replace('.', 'p');
const T = parseInt(args.ticks || (MONO ? '20000' : '3000'), 10);
const BIN = parseInt(args.bin || (MONO ? '500' : '250'), 10);
const every = 250;
const frac = 0.1;
const grow = 50000;
const mainDir = args['main-dir'] || path.join(__dirname, 'main');
const outDir = args.out || '.';
const { LET, K, DISH, NB } = Sim;
const enc = s => [...s].map(ch => { const k = LET.indexOf(ch); if (k < 0) throw new Error(`글자 아님: ${ch}`); return k; });
const opts = { seed, density: densityV, mu: growMu, light: 1, mat: true, energy: false };
if (ancestorV) opts.ancestor = ancestorV;
const INJ_SEED = (seed * 1000003 + 17) >>> 0;
const POST_A = (seed * 7919 + 1) >>> 0;
const t0 = Date.now();
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

/* ---------- 계통별 예산 ---------- */
const KEYS = ['stepped', 'not_stepped', 'unknown', 'births', 'deaths', 'death_age', 'death_prestep', 'death_other', 'age_sum',
  'copy_ok', 'copy_stall', 'copy_del', 'copy_idle', 'copy_phase',
  'alloc_ok', 'alloc_blocked', 'alloc_idle', 'alloc_dead', 'div_with_child', 'div_idle'];
function newBin(t, w) {
  const b = { t, n: [0, 0] };
  for (const o of w.alive) b.n[o.tag === 1 ? 1 : 0]++;
  for (const k of KEYS) b[k] = [0, 0];
  if (byLetter) b.stall_by = [new Array(K).fill(0), new Array(K).fill(0)];
  return b;
}
function run(w, countSeries) {
  const bins = [];
  const count = () => { let m = 0; for (const o of w.alive) if (o.tag === 1) m++; return [w.alive.length, m]; };
  const series = countSeries ? [[0, ...count()]] : null;
  let stopped = -1;
  let age0 = new Int32Array(1 << 15), len0 = new Int32Array(1 << 15), rh0 = new Int32Array(1 << 15), gl0 = new Int32Array(1 << 15);
  let ch0 = new Array(1 << 15);
  for (let t = 1; t <= T; t++) {
    if ((t - 1) % BIN === 0) bins.push(newBin(t - 1, w));
    const B = bins[bins.length - 1];
    const pre = w.alive.slice();
    const n = pre.length;
    if (n > age0.length) {
      let c = age0.length; while (c < n) c *= 2;
      const g = A => { const x = new Int32Array(c); x.set(A); return x; };
      age0 = g(age0); len0 = g(len0); rh0 = g(rh0); gl0 = g(gl0); ch0.length = c;
    }
    for (let q = 0; q < n; q++) {
      const o = pre[q];
      age0[q] = o.age; ch0[q] = o.child; len0[q] = o.child ? o.child.g.length : -1; rh0[q] = o.rh; gl0[q] = o.g.length;
    }
    Sim.tick(w);
    const bt = w.tick - 1;
    for (const o of w.alive) if (o.bornAt === bt) B.births[o.tag === 1 ? 1 : 0]++;
    for (let q = 0; q < n; q++) {
      const o = pre[q], L = o.tag === 1 ? 1 : 0;
      const stepped = o.age === age0[q] + 1;
      B.age_sum[L] += age0[q];
      if (o.dead) {
        B.deaths[L]++;
        if (!stepped) B.death_prestep[L]++;
        else if (o.age > o.maxAge) B.death_age[L]++;
        else B.death_other[L]++;
      }
      if (!stepped) { B.not_stepped[L]++; continue; }
      B.stepped[L]++;
      const c = ch0[q];
      const copying = c && rh0[q] < gl0[q];
      if (copying) B.copy_phase[L]++;
      const host = o.lastHost === o.id ? o : w.orgs.get(o.lastHost);
      if (!host) { B.unknown[L]++; continue; }
      const ins = host.g[o.last];
      if (ins === 4) {
        if (!copying) B.copy_idle[L]++;
        else if (c.g.length > len0[q]) B.copy_ok[L]++;
        else if (o.rh > rh0[q]) B.copy_del[L]++;
        else { B.copy_stall[L]++; if (byLetter) B.stall_by[L][o.g[rh0[q]]]++; }
      } else if (ins === 3) {
        if (c) B.alloc_idle[L]++;
        else if (o.child) B.alloc_ok[L]++;
        else if (o.dead) B.alloc_dead[L]++;
        else B.alloc_blocked[L]++;
      } else if (ins === 6) {
        if (c) B.div_with_child[L]++; else B.div_idle[L]++;
      }
    }
    for (let q = 0; q < n; q++) ch0[q] = null;
    if (countSeries && (t % every === 0 || w.extinctAt >= 0)) {
      const cc = count();
      series.push([t, ...cc]);
      if (cc[1] === 0 || w.extinctAt >= 0) { stopped = t; break; }
    }
    if (!countSeries && w.extinctAt >= 0) { stopped = t; break; }
  }
  return { bins, series, stopped };
}

/* ---------- 실행 ---------- */
fs.mkdirSync(outDir, { recursive: true });
if (MONO) {
  const w = Sim.makeWorld({ seed, density: 8, mu: 0, light: 1, ancestor: MONO, mat: true, energy: false });
  const r = run(w, false);
  const out = { invbud_version: INVBUD_VERSION, sim_version: Sim.VERSION, mode: 'mono', code: MONO, seed, ticks: T, bin: BIN,
    extinct_at: w.extinctAt, pop_end: w.alive.length, checksum: Sim.checksum(w), bins: r.bins, keys: KEYS, elapsed_ms: Date.now() - t0 };
  const tag = `ib_mono_${MONO}_s${String(seed).padStart(5, '0')}`;
  fs.writeFileSync(path.join(outDir, `${tag}.json`), JSON.stringify(out));
  console.log(`DONE ${tag} · v${Sim.VERSION}/ib${INVBUD_VERSION} · pop ${w.alive.length} · ext ${w.extinctAt} · chk ${out.checksum} · ${((Date.now() - t0) / 1000).toFixed(1)}s`);
} else {
  const w = Sim.makeWorld(opts);
  while (w.tick < grow && w.extinctAt < 0) Sim.tick(w);
  const growChecksum = Sim.checksum(w);
  const mainPath = path.join(mainDir, `mat_d${String(densityV).replace('.', 'p')}_mu0p01_s${String(seed).padStart(5, '0')}.json`);
  const main = (growMu === 0.01 && !ancestorV && fs.existsSync(mainPath)) ? JSON.parse(fs.readFileSync(mainPath, 'utf8')) : null;
  const arms = [];
  if (w.extinctAt < 0) {
    for (const [kind, key] of [['ref', WT], ...codes.map(c => ['mut', c])]) {
      const ta = Date.now();
      const a = cloneWorld(w, POST_A);
      a.o.mu = mu;
      const inj = inject(a, enc(key));
      const r = run(a, true);
      arms.push({ kind, key, post_seed: POST_A, inj, stopped: r.stopped, series: r.series, bins: r.bins, elapsed_ms: Date.now() - ta });
    }
  }
  const out = { invbud_version: INVBUD_VERSION, sim_version: Sim.VERSION, mode: 'inv', seed, wt: WT, codes, mu, frac, ticks: T, bin: BIN, by_letter: byLetter || undefined, density: densityV !== 8 ? densityV : undefined, grow_mu: growMu !== 0.01 ? growMu : undefined, ancestor: ancestorV || undefined,
    inj_seed: INJ_SEED, post_a: POST_A, grow_checksum: growChecksum, main_checksum: main ? main.checksum : null,
    checksum_match: !!main && main.checksum === growChecksum, keys: KEYS, arms, elapsed_ms: Date.now() - t0 };
  const tag = `ib_inv_s${String(seed).padStart(5, '0')}_mu${String(mu).replace('.', 'p')}${dtag}${gtag}`;
  fs.writeFileSync(path.join(outDir, `${tag}.json`), JSON.stringify(out));
  console.log(`DONE ${tag} · v${Sim.VERSION}/ib${INVBUD_VERSION} · checksum ${out.checksum_match ? 'MATCH' : 'MISMATCH'} · arms ${arms.length} · ${((Date.now() - t0) / 1000).toFixed(0)}s`);
}
