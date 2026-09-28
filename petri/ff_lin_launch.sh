#!/usr/bin/env bash
# 해부 7 ② 공급 · ③a 교체 경로 — 계보 추적을 두 규칙으로. 사전등록 해부7-사전등록-2026-09-16.md §2 · §3
#   규칙 {원래(v0.3.2 와 같음) · 재료 먼저} × 코드 {racld(지켜볼 코드 rascld) · rascld(racld) · rsacld(racld)}
#   × μ 0.1% · 0.5% · 1% × 시드 701~730 · mat 밀도 8 · 20,000틱 → 2 × 3 × 3 × 30 = 540 접시
# 마른 실행: PETRI_SEEDS="691 692" PETRI_ROOT=smoke7lin ./ff_lin_launch.sh
set -u
export NVM_DIR="$HOME/.nvm"; . "$NVM_DIR/nvm.sh"
cd "$HOME/petri" || exit 1
grep -qE "VERSION = '0\.3\.[345]'" sim.js || { echo "sim.js 가 v0.3.3~0.3.5 가 아니다 — 발사 중단"; exit 1; }
grep -qE "LINEAGE_VERSION = '1\.[23]\.0'" lineage.js || { echo "lineage.js 가 v1.2.0/1.3.0 이 아니다 — 발사 중단"; exit 1; }
SEEDS=${PETRI_SEEDS:-$(seq 701 730)}
ROOT=${PETRI_ROOT:-lin7}
L=${ROOT}.log; : > "$L"
J=${ROOT}_jobs.txt; : > "$J"
for ff in 0 1; do
  sub=$([ "$ff" = 1 ] && echo ff || echo orig)
  for mu in 0.001 0.005 0.01; do
    for s in $SEEDS; do
      echo "racld rascld $mu $s $ff $ROOT/$sub" >> "$J"
      echo "rascld racld $mu $s $ff $ROOT/$sub" >> "$J"
      echo "rsacld racld $mu $s $ff $ROOT/$sub" >> "$J"
    done
  done
done
echo "작업 $(wc -l < "$J") · $(date '+%F %T')" >> "$L"
t0=$(date +%s)
xargs -P 60 -n 6 sh -c 'node lineage.js --code "$0" --watch "$1" --mu "$2" --seed "$3" --find-first "$4" --out "$5"' < "$J" >> "$L" 2>&1
echo "===== 끝 $(( $(date +%s) - t0 ))초 =====" >> "$L"
