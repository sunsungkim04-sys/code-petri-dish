#!/usr/bin/env node
/* 코드 배양접시 — 시간 예산 계측기(개체-틱이 어디에 쓰였나). 판정은 하지 않는다(budget_analyze.py).
 *
 * 단일 코드 접시를 μ = 0 으로 돌리며, 틱마다 **각 개체가 실행한 명령과 그 결과**를 밖에서 읽는다.
 * sim.js 는 바꾸지 않는다: 틱 전에 개체마다 (나이 · 자식 객체 · 자식 길이 · 복사 위치 rh)를 적어 두고,
 * 틱 뒤에 org.last(실행한 자리)와 이 기록을 견준다. 읽기만 하므로 궤적 불변(dens 기록과 체크섬 대조).
 *
 *   복사 c : 성공(자식이 한 글자 늘었다) · 헛손질(자식이 있고 다 못 베꼈는데 그대로 — 가까운 칸에 글자가 없었다)
 *            · 빠뜨림(rh 만 늘었다 · μ0 에선 0 이어야 한다) · 할 일 없음(자식 없음 또는 다 베낌)
 *   자리잡기 a : 성공 · 막힘(바라보는 칸이 차 있거나 접시 밖) · 할 일 없음(자식이 이미 있다) · 알 수 없음(그 틱에 죽음)
 *   나누기 d : 자식이 있었다 · 할 일 없음
 *   나머지 명령은 명령 글자별로 센다.
 * v1.1.0: 창을 고정하지 않고 **bin 틱 칸**마다 센다(칸 시작 개체 수를 함께 적는다). 채우는 시기 · 다 찬 시기는
 *         판정 스크립트가 개체 수로 가른다 — 고정 창 [0, 1만)은 밀도마다 채우는 시기가 달라 섞였다(v1.0.0 시험에서 확인).
 *
 * v1.2.0(해부 11): --age0 N --age-var N(sim 0.3.5 · 세계 전체) · 칸마다 죽음 수 · 죽은 나이 합 · 칸 시작의 글자별 자유 재료를 더 적는다.
 *         읽기만 한다 — 옵션을 안 주면 궤적 · 옛 필드가 v1.1.0 과 같다(reg11_check.py).
 *
 * 사용: node budget.js --cond mat --density 8 --code racld --seed 601 [--ticks 20000] [--bin 500] [--out dir]
 */
'use strict';
const fs = require('fs');
const path = require('path');
const Sim = require('./sim.js');
const BUDGET_VERSION = '1.2.0';

const args = {};
for (let i = 2; i < process.argv.length; i += 2) args[process.argv[i].replace(/^--/, '')] = process.argv[i + 1];
if (!['0.3.2', '0.3.3', '0.3.4', '0.3.5'].includes(Sim.VERSION)) { console.error(`sim ${Sim.VERSION} — 0.3.2~0.3.5 가 필요하다`); process.exit(2); }
const COND = { space: { mat: false, energy: false }, mat: { mat: true, energy: false } };
const cond = args.cond;
if (!COND[cond]) { console.error(`알 수 없는 cond: ${cond}`); process.exit(2); }
const code = args.code;
const seed = parseInt(args.seed, 10);
const density = parseFloat(args.density || '8');
const T = parseInt(args.ticks || '20000', 10);
const BIN = parseInt(args.bin || '500', 10);
const outDir = args.out || '.';
const age0V = args.age0 ? parseInt(args.age0, 10) : null;
const ageVarV = args['age-var'] ? parseInt(args['age-var'], 10) : null;
if ((age0V !== null || ageVarV !== null) && Sim.VERSION !== '0.3.5') { console.error('--age0 은 sim 0.3.5 가 필요하다'); process.exit(2); }
for (const ch of code) if (!Sim.LET.includes(ch)) { console.error(`코드에 없는 글자: ${ch}`); process.exit(2); }

const t0 = Date.now();
const wopts = Object.assign({ seed, density, mu: 0, light: 1, ancestor: code }, COND[cond]);
if (age0V !== null) wopts.age0 = age0V;
if (ageVarV !== null) wopts.ageVar = ageVarV;
const w = Sim.makeWorld(wopts);
const KEYS = ['stepped', 'not_stepped', 'unknown',
  'copy_ok', 'copy_stall', 'copy_del', 'copy_idle',
  'alloc_ok', 'alloc_blocked', 'alloc_idle', 'alloc_dead',
  'div_with_child', 'div_idle', 'births'];
const freeNow = () => { if (!w.mat) return null; const f = new Array(Sim.K).fill(0); const m = w.mat; for (let j = 0; j < m.length; j++) f[j % Sim.K] += m[j]; return f; };
const mk = (t, pop) => { const o = { t0: t, pop0: pop, deaths: 0, death_age_sum: 0, free0: freeNow() }; for (const k of KEYS) o[k] = 0; o.by_ins = new Array(Sim.K).fill(0); return o; };
const bins = [];

let n = 0, age0 = new Int32Array(1 << 15), len0 = new Int32Array(1 << 15), rh0 = new Int32Array(1 << 15);
let ch0 = new Array(1 << 15);
function ensure(m) {
  if (m <= age0.length) return;
  let c = age0.length; while (c < m) c *= 2;
  const g = (A) => { const b = new Int32Array(c); b.set(A); return b; };
  age0 = g(age0); len0 = g(len0); rh0 = g(rh0); ch0.length = c;
}
while (w.tick < T && w.extinctAt < 0) {
  if (w.tick % BIN === 0) bins.push(mk(w.tick, w.alive.length));
  const W = bins[bins.length - 1];
  const pre = w.alive.slice();
  n = pre.length; ensure(n);
  for (let q = 0; q < n; q++) {
    const o = pre[q];
    age0[q] = o.age; ch0[q] = o.child; len0[q] = o.child ? o.child.g.length : -1; rh0[q] = o.rh;
  }
  const b0 = w.births;
  Sim.tick(w);
  W.births += w.births - b0;
  for (let q = 0; q < n; q++) if (pre[q].dead) { W.deaths++; W.death_age_sum += pre[q].age; }
  for (let q = 0; q < n; q++) {
    const o = pre[q];
    if (o.age !== age0[q] + 1) { W.not_stepped++; continue; }
    W.stepped++;
    const host = o.lastHost === o.id ? o : w.orgs.get(o.lastHost);
    if (!host) { W.unknown++; continue; }
    const ins = host.g[o.last];
    W.by_ins[ins]++;
    const c = ch0[q];
    if (ins === 4) {
      if (!c || rh0[q] >= o.g.length) W.copy_idle++;
      else if (c.g.length > len0[q]) W.copy_ok++;
      else if (o.rh > rh0[q]) W.copy_del++;
      else W.copy_stall++;
    } else if (ins === 3) {
      if (c) W.alloc_idle++;
      else if (o.child) W.alloc_ok++;
      else if (o.dead) W.alloc_dead++;
      else W.alloc_blocked++;
    } else if (ins === 6) {
      if (c) W.div_with_child++; else W.div_idle++;
    }
  }
  for (let q = 0; q < n; q++) ch0[q] = null;
}
const last = w.alive.length;
const out = {
  budget_version: BUDGET_VERSION, sim_version: Sim.VERSION, opts: wopts, mat_total0: w.matTotal0, cond, density: COND[cond].mat ? density : null, code, seed,
  ticks_planned: T, ticks_run: w.tick, bin: BIN, extinct_at: w.extinctAt, pop_end: last, checksum: Sim.checksum(w),
  letters: Sim.LET, bins, elapsed_ms: Date.now() - t0
};
fs.mkdirSync(outDir, { recursive: true });
const tag = `bud_${cond}${COND[cond].mat ? '_d' + String(density).padStart(2, '0') : ''}_${code}_s${String(seed).padStart(5, '0')}${age0V !== null ? '_a' + age0V : ''}`;
fs.writeFileSync(path.join(outDir, `${tag}.json`), JSON.stringify(out));
const sum = k => bins.reduce((a, b) => a + b[k], 0);
console.log(`DONE ${tag} · v${Sim.VERSION}/b${BUDGET_VERSION} · pop ${last} · 칸 ${bins.length} · 헛손질 ${sum('copy_stall')} · 막힘 ${sum('alloc_blocked')} · 빠뜨림 ${sum('copy_del')} · 모름 ${sum('unknown') + sum('alloc_dead')} · chk ${out.checksum} · ${((Date.now() - t0) / 1000).toFixed(1)}s`);
