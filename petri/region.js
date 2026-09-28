#!/usr/bin/env node
/* 코드 배양접시 — 에너지 접시의 구역별 조성 계측기(해부 4C). 판정은 하지 않는다(region_analyze.py).
 *
 * 본실험 에너지 접시를 같은 설정으로 50,000틱 다시 키운 뒤(본실험 체크섬과 대조) 끝 상태를 구역으로 나눠 센다.
 * 구역은 sim.js 의 빛 얼룩과 같은 자리다 — L1 (0.3N, 0.35N) · L2 (0.7N, 0.65N) · 반지름 0.18N 안이면 그 얼룩, 아니면 어두운 곳.
 * 각 구역: 개체 수 · h 든 비율 · 상위 10 코드 · **무작위 절반 두 벌**(구역 안 변이를 재는 대조).
 * 사용: node region.js --seed 1 [--grow 50000] [--main-dir ../main] [--out dir]
 */
'use strict';
const fs = require('fs');
const path = require('path');
const Sim = require('./sim.js');

const args = {};
for (let i = 2; i < process.argv.length; i += 2) args[process.argv[i].replace(/^--/, '')] = process.argv[i + 1];
const seed = parseInt(args.seed, 10);
const grow = parseInt(args.grow || '50000', 10);
const mainDir = args['main-dir'] || path.join(__dirname, 'main');
const outDir = args.out || '.';
if (!['0.3.2', '0.3.3', '0.3.4', '0.3.5'].includes(Sim.VERSION)) { console.error(`sim ${Sim.VERSION} — 0.3.2~0.3.5 가 필요하다`); process.exit(2); }
const N = Sim.N;
const L1 = [0.3 * N, 0.35 * N], L2 = [0.7 * N, 0.65 * N], RAD = 0.18 * N;
const t0 = Date.now();
const w = Sim.makeWorld({ seed, density: 8, mu: 0.01, light: 1, mat: false, energy: true });
while (w.tick < grow && w.extinctAt < 0) Sim.tick(w);
const chk = Sim.checksum(w);
let mainChk = null;
const mainPath = path.join(mainDir, `energy_d8_mu0p01_s${String(seed).padStart(5, '0')}.json`);
if (fs.existsSync(mainPath)) mainChk = JSON.parse(fs.readFileSync(mainPath, 'utf8')).checksum;

const rnd = Sim.mulberry32((seed * 31337 + 7) >>> 0);
const regions = { L1: [], L2: [], dark: [] };
for (const o of w.alive) {
  const x = o.i % N, y = (o.i / N) | 0;
  const d1 = Math.hypot(x - L1[0], y - L1[1]), d2 = Math.hypot(x - L2[0], y - L2[1]);
  (d1 <= RAD ? regions.L1 : d2 <= RAD ? regions.L2 : regions.dark).push(o.key);
}
function comp(keys) {
  const c = new Map();
  for (const k of keys) c.set(k, (c.get(k) || 0) + 1);
  return [...c.entries()].sort((a, b) => b[1] - a[1] || (a[0] < b[0] ? -1 : 1));
}
const out = { sim_version: Sim.VERSION, seed, grow, checksum: chk, main_checksum: mainChk, checksum_match: mainChk === chk,
  pop: w.alive.length, extinct_at: w.extinctAt, regions: {}, elapsed_ms: 0 };
for (const [name, keys] of Object.entries(regions)) {
  const halves = [[], []];
  for (const k of keys) halves[rnd() < 0.5 ? 0 : 1].push(k);
  out.regions[name] = {
    n: keys.length,
    h_share: keys.length ? keys.filter(k => k.includes('h')).length / keys.length : 0,
    kinds: new Set(keys).size,
    top: comp(keys).slice(0, 10),
    half_a: comp(halves[0]).slice(0, 40),
    half_b: comp(halves[1]).slice(0, 40),
    full: comp(keys).slice(0, 200)
  };
}
out.elapsed_ms = Date.now() - t0;
fs.mkdirSync(outDir, { recursive: true });
fs.writeFileSync(path.join(outDir, `region_energy_s${String(seed).padStart(5, '0')}.json`), JSON.stringify(out));
console.log(`DONE region_energy_s${String(seed).padStart(5, '0')} · checksum ${mainChk === null ? 'no-main' : (out.checksum_match ? 'MATCH' : 'MISMATCH')} · pop ${w.alive.length} · ` +
  Object.entries(out.regions).map(([k, v]) => `${k} ${v.n}(h ${(v.h_share * 100).toFixed(0)}%)`).join(' · ') + ` · ${((Date.now() - t0) / 1000).toFixed(0)}s`);
