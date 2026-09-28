#!/usr/bin/env bash
# 해부 5A — 재료 밀도 사다리. 사전등록 해부5-사전등록-2026-09-16.md §2
#   칸 = mat 밀도 4·8·14·20·32 + space(재료 규칙 없음)
#   코드 = acld(4) racld(5) rascld(6) rsacld(6) rascled(7)  ·  μ 0  ·  시드 501~530  ·  20,000틱
#   → 6 칸 × 5 코드 × 30 시드 = 900 접시. 판정은 dens_analyze.py.
# 마른 실행: PETRI_SEEDS="501 502" PETRI_ROOT=smoke5 ./dens_launch.sh  (기본값은 사전등록 값 그대로)
set -u
export NVM_DIR="$HOME/.nvm"; . "$NVM_DIR/nvm.sh"
cd "$HOME/petri" || exit 1
grep -qE "VERSION = '0\.3\.[2345]'" sim.js || { echo "sim.js 가 v0.3.2~0.3.5 가 아니다 — 발사 중단"; exit 1; }
grep -q "MONO_VERSION = '1.1.0'" mono.js || { echo "mono.js 가 v1.1.0 이 아니다 — 발사 중단"; exit 1; }
SEEDS=${PETRI_SEEDS:-$(seq 501 530)}
ROOT=${PETRI_ROOT:-dens}
L=${ROOT}.log; : > "$L"
J=${ROOT}_jobs.txt; : > "$J"
for code in acld racld rascld rsacld rascled; do
  for d in 4 8 14 20 32; do
    for s in $SEEDS; do echo "mat $d $code $s $ROOT/mat_d$(printf %02d "$d")" >> "$J"; done
  done
  for s in $SEEDS; do echo "space 8 $code $s $ROOT/space" >> "$J"; done
done
echo "작업 $(wc -l < "$J") · $(date '+%F %T')" >> "$L"
t0=$(date +%s)
xargs -P 60 -n 5 sh -c 'node mono.js --cond "$0" --density "$1" --code "$2" --seed "$3" --mu 0 --ticks 20000 --out "$4"' < "$J" >> "$L" 2>&1
# 점검 ⑤ — 같은 시드 두 번 (첫 시드로)
S1=$(echo $SEEDS | awk '{print $1}')
for code in racld rascld; do
  node mono.js --cond mat --density 8 --code "$code" --seed "$S1" --mu 0 --ticks 20000 --out "${ROOT}_twice" >> "$L" 2>&1
done
echo "===== 끝 $(( $(date +%s) - t0 ))초 =====" >> "$L"
