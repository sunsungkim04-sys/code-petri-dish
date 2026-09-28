#!/usr/bin/env bash
# 해부 8 C — 부풀림을 주사위 몫 R 과 거름 몫 F 로. 사전등록 해부8-사전등록-2026-09-17.md §4
#   규칙 {원래 · 재료 먼저} × 코드 {racld · rascld · rsacld} × 밀도 {4 · 8 · 32} × μ {0.1% · 1%} × 시드 901~920 · 20,000틱
#   → 2 × 3 × 3 × 2 × 20 = 720 접시 (lineage.js v1.3.0 --spec 1). 판정은 spec_analyze.py.
# 마른 실행: PETRI_SEEDS="9301 9302" PETRI_ROOT=smoke8c ./spec8_launch.sh
set -u
export NVM_DIR="$HOME/.nvm"; . "$NVM_DIR/nvm.sh"
cd "$HOME/petri" || exit 1
grep -qE "VERSION = '0\.3\.[345]'" sim.js || { echo "sim.js 가 v0.3.3~0.3.5 가 아니다 — 발사 중단"; exit 1; }
grep -qE "LINEAGE_VERSION = '1\.[345]\.0'" lineage.js || { echo "lineage.js 가 v1.3.0~1.5.0 이 아니다 — 발사 중단"; exit 1; }
SEEDS=${PETRI_SEEDS:-$(seq 901 920)}
ROOT=${PETRI_ROOT:-spec8}
L=${ROOT}.log; : > "$L"
J=${ROOT}_jobs.txt; : > "$J"
for ff in 0 1; do
  sub=$([ "$ff" = 1 ] && echo ff || echo orig)
  for d in 4 8 32; do
    for mu in 0.001 0.01; do
      for s in $SEEDS; do
        for code in racld rascld rsacld; do
          echo "$code $mu $s $ff $d $ROOT/$sub/d$d" >> "$J"
        done
      done
    done
  done
done
[ "$(awk 'NF!=6' "$J" | wc -l)" = 0 ] || { echo "작업 줄의 칸 수가 6 이 아니다 — 발사 중단"; exit 1; }
echo "작업 $(wc -l < "$J") · $(date '+%F %T')" >> "$L"
t0=$(date +%s)
xargs -P 40 -n 6 sh -c 'node lineage.js --code "$0" --mu "$1" --seed "$2" --find-first "$3" --density "$4" --spec 1 --out "$5"' < "$J" >> "$L" 2>&1
echo "===== 끝 $(( $(date +%s) - t0 ))초 =====" >> "$L"
