#!/usr/bin/env bash
# 파일럿 — 계측기 · 시간척도 · 재료 밀도 확인용. 되감기 질문의 판정에는 쓰지 않는다.
# 본실험 시드(1~100)와 겹치지 않게 9001~9020 을 쓴다.
# 조합: space·energy(밀도 무관, 4 로 표기) + mat·both × 밀도 4·6·8·10·14  →  시드당 12개 × 20 = 240
set -u
export NVM_DIR="$HOME/.nvm"; . "$NVM_DIR/nvm.sh"
cd "$HOME/petri" || exit 1
J=pilot_jobs.txt; : > "$J"
for s in $(seq 9001 9020); do
  for c in space energy; do echo "$c 4 $s" >> "$J"; done
  for d in 4 6 8 10 14; do for c in mat both; do echo "$c $d $s" >> "$J"; done; done
done
L=pilot.log; : > "$L"
echo "작업 $(wc -l < "$J")" >> "$L"
t0=$(date +%s)
xargs -P 60 -n 3 sh -c 'node replay.js --cond "$0" --density "$1" --seed "$2" --ticks 20000 --out pilot' < "$J" >> "$L" 2>&1
echo "===== 끝 $(( $(date +%s) - t0 ))초 =====" >> "$L"
