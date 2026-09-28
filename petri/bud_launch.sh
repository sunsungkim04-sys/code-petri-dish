#!/usr/bin/env bash
# 해부 6B — 시간 예산. 사전등록 해부6-사전등록-2026-09-16.md §3
#   칸 mat 밀도 4·8·14·32 + space  ×  racld · rascld · rsacld  ×  μ 0  ×  시드 601~630  ×  20,000틱
#   → 5 × 3 × 30 = 450 접시. 판정은 budget_analyze.py.
#   점검: 첫 시드의 밀도 8 racld 를 mono.js 로도 돌려 체크섬 대조(budget.js 가 궤적을 안 바꾸나)
# 마른 실행: PETRI_SEEDS="591 592" PETRI_ROOT=smoke6bud ./bud_launch.sh
set -u
export NVM_DIR="$HOME/.nvm"; . "$NVM_DIR/nvm.sh"
cd "$HOME/petri" || exit 1
grep -qE "VERSION = '0\.3\.[2345]'" sim.js || { echo "sim.js 가 v0.3.2~0.3.5 가 아니다 — 발사 중단"; exit 1; }
grep -q "BUDGET_VERSION = '1.1.0'" budget.js || { echo "budget.js 가 v1.1.0 이 아니다 — 발사 중단"; exit 1; }
SEEDS=${PETRI_SEEDS:-$(seq 601 630)}
ROOT=${PETRI_ROOT:-bud6}
L=${ROOT}.log; : > "$L"
J=${ROOT}_jobs.txt; : > "$J"
for s in $SEEDS; do
  for code in racld rascld rsacld; do
    for d in 4 8 14 32; do echo "mat $d $code $s $ROOT" >> "$J"; done
    echo "space 8 $code $s $ROOT" >> "$J"
  done
done
echo "작업 $(wc -l < "$J") · $(date '+%F %T')" >> "$L"
t0=$(date +%s)
xargs -P 60 -n 5 sh -c 'node budget.js --cond "$0" --density "$1" --code "$2" --seed "$3" --out "$4"' < "$J" >> "$L" 2>&1
S1=$(echo $SEEDS | awk '{print $1}')
node mono.js --cond mat --density 8 --code racld --seed "$S1" --mu 0 --ticks 20000 --out "${ROOT}_mono" >> "$L" 2>&1
echo "===== 끝 $(( $(date +%s) - t0 ))초 =====" >> "$L"
