#!/usr/bin/env bash
# 해부 6A — 계보 추적. 사전등록 해부6-사전등록-2026-09-16.md §2
#   d08: racld · rascld · rsacld  /  d32: racld · rascld   ×  μ 0.1% · 0.5% · 1%  ×  시드 601~630  ×  20,000틱
#   → (3 + 2) × 3 × 30 = 450 접시. 판정은 lineage_analyze.py.
# 마른 실행: PETRI_SEEDS="591 592" PETRI_ROOT=smoke6lin ./lin_launch.sh
set -u
export NVM_DIR="$HOME/.nvm"; . "$NVM_DIR/nvm.sh"
cd "$HOME/petri" || exit 1
grep -qE "VERSION = '0\.3\.[2345]'" sim.js || { echo "sim.js 가 v0.3.2~0.3.5 가 아니다 — 발사 중단"; exit 1; }
grep -qE "LINEAGE_VERSION = '1\.[123]\.0'" lineage.js || { echo "lineage.js 가 v1.1.0~1.3.0 이 아니다 — 발사 중단"; exit 1; }
SEEDS=${PETRI_SEEDS:-$(seq 601 630)}
ROOT=${PETRI_ROOT:-lin6}
L=${ROOT}.log; : > "$L"
J=${ROOT}_jobs.txt; : > "$J"
for mu in 0.001 0.005 0.01; do
  for s in $SEEDS; do
    for code in racld rascld rsacld; do echo "$code 8 $mu $s $ROOT/d08" >> "$J"; done
    for code in racld rascld; do echo "$code 32 $mu $s $ROOT/d32" >> "$J"; done
  done
done
echo "작업 $(wc -l < "$J") · $(date '+%F %T')" >> "$L"
t0=$(date +%s)
xargs -P 60 -n 5 sh -c 'node lineage.js --code "$0" --density "$1" --mu "$2" --seed "$3" --out "$4"' < "$J" >> "$L" 2>&1
echo "===== 끝 $(( $(date +%s) - t0 ))초 =====" >> "$L"
