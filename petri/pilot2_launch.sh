#!/usr/bin/env bash
# 파일럿 2 — sim v0.3.0(무늬 글자 x 추가). 판정에 쓰지 않는다.
#   ① 결정론 재확인(space 시드 9001 3,000틱 두 번)
#   ② 재료 밀도 8 생존 재확인(글자 종류가 10 → 11 로 늘어 글자당 재료가 줄었다) · 밀도 6 · 10 은 감도 확인
#   ③ 50,000틱 소요 시간 · 무늬 x 가 중립 대조로 쓸 만한지(pilot2_check.py)
# 본실험 시드(1~100)와 겹치지 않게 9001~9020 을 쓴다. 조합: space · energy · mat 6/8/10  →  시드당 5개 × 20 = 100
set -u
export NVM_DIR="$HOME/.nvm"; . "$NVM_DIR/nvm.sh"
cd "$HOME/petri" || exit 1
L=pilot2.log; : > "$L"
rm -rf smoke2 pilot2
node replay.js --cond space --density 8 --seed 9001 --ticks 3000 --out smoke2/a >> "$L" 2>&1
node replay.js --cond space --density 8 --seed 9001 --ticks 3000 --out smoke2/b >> "$L" 2>&1
J=pilot2_jobs.txt; : > "$J"
for s in $(seq 9001 9020); do
  for c in space energy; do echo "$c 8 $s" >> "$J"; done
  for d in 6 8 10; do echo "mat $d $s" >> "$J"; done
done
echo "작업 $(wc -l < "$J")" >> "$L"
t0=$(date +%s)
xargs -P 60 -n 3 sh -c 'node replay.js --cond "$0" --density "$1" --seed "$2" --ticks 50000 --out pilot2' < "$J" >> "$L" 2>&1
echo "===== 끝 $(( $(date +%s) - t0 ))초 =====" >> "$L"
