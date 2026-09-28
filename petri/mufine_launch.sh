#!/usr/bin/env bash
# 해부 5B — 촘촘한 μ 사다리. 사전등록 해부5-사전등록-2026-09-16.md §3
#   mat 밀도 8 · racld 대 rascld · μ 0·0.001·0.002·0.003·0.004·0.005·0.0075·0.01 · 시드 501~560 · 20,000틱
#   → 2 코드 × 8 μ × 60 시드 = 960 접시. 판정은 mufine_analyze.py.
# 마른 실행: PETRI_SEEDS="501 502" PETRI_ROOT=smoke5mu ./mufine_launch.sh  (기본값은 사전등록 값 그대로)
set -u
export NVM_DIR="$HOME/.nvm"; . "$NVM_DIR/nvm.sh"
cd "$HOME/petri" || exit 1
grep -qE "VERSION = '0\.3\.[2345]'" sim.js || { echo "sim.js 가 v0.3.2~0.3.5 가 아니다 — 발사 중단"; exit 1; }
grep -q "MONO_VERSION = '1.1.0'" mono.js || { echo "mono.js 가 v1.1.0 이 아니다 — 발사 중단"; exit 1; }
SEEDS=${PETRI_SEEDS:-$(seq 501 560)}
ROOT=${PETRI_ROOT:-mufine}
L=${ROOT}.log; : > "$L"
J=${ROOT}_jobs.txt; : > "$J"
for code in racld rascld; do
  for mu in 0 0.001 0.002 0.003 0.004 0.005 0.0075 0.01; do
    for s in $SEEDS; do echo "$code $mu $s $ROOT" >> "$J"; done
  done
done
echo "작업 $(wc -l < "$J") · $(date '+%F %T')" >> "$L"
t0=$(date +%s)
xargs -P 60 -n 4 sh -c 'node mono.js --cond mat --density 8 --code "$0" --mu "$1" --seed "$2" --ticks 20000 --out "$3"' < "$J" >> "$L" 2>&1
echo "===== 끝 $(( $(date +%s) - t0 ))초 =====" >> "$L"
