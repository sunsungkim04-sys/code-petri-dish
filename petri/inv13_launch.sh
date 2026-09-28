#!/usr/bin/env bash
# 해부 13 — 둘째 코드 가족에서 짝지은 고리 사다리. 사전등록 해부13-사전등록-2026-09-27.md §3
#   밑 코드 {acld · rascled} × μ {0 · 0.003} × mat 배경 3~22 → 80 작업
#   계측기는 해부 12 것 그대로(dms.js --arms nbr --nbr-of CODE) — 새 패치 없음
#   판정: pairs13_analyze.py
# 마른 실행: PETRI_NBG=1 PETRI_BG_FROM=23 PETRI_MUS="0" PETRI_ROOT=pilot13 ./inv13_launch.sh
set -u
export NVM_DIR="$HOME/.nvm"; . "$NVM_DIR/nvm.sh"
cd "${PETRI_HOME:-$HOME/petri}" || exit 1
grep -q "VERSION = '0.3.5'" sim.js || { echo "sim.js 가 v0.3.5 가 아니다 — 발사 중단"; exit 1; }
grep -q "nbrOf" dms.js || { echo "dms.js 에 --nbr-of 가 없다 — 발사 중단"; exit 1; }
NBG=${PETRI_NBG:-20}
MUS=${PETRI_MUS:-"0 0.003"}
NBASES=${PETRI_NBASES:-"acld rascled"}
ROOT=${PETRI_ROOT:-inv13}
L=${ROOT}.log; : > "$L"
J=${ROOT}_jobs.txt; : > "$J"
BGS=$(NBG=$NBG python3 - <<'PYEOF'
import json, os
nbg = int(os.environ["NBG"]); first = int(os.environ.get("PETRI_BG_FROM", "1")); picked = []
for s in range(first, 101):
    z = json.load(open(f"main/mat_d8_mu0p01_s{s:05d}.json"))
    if z["extinct_at"] < 0 and z["final_counts"] and z["final_counts"][0][0] == "racld":
        picked.append(s)
    if len(picked) == nbg:
        break
print(" ".join(map(str, picked)))
PYEOF
)
for b in $NBASES; do for mu in $MUS; do tag=mu$(echo "$mu" | sed 's/\./p/'); for s in $BGS; do
  echo "$mu $s $b ${ROOT}nbr_$b/$tag" >> "$J"
done; done; done
[ "$(awk 'NF!=4' "$J" | wc -l)" = 0 ] || { echo "작업 줄의 칸 수가 4 가 아니다 — 발사 중단"; exit 1; }
echo "배경 $BGS" >> "$L"
echo "작업 $(wc -l < "$J") · 밑 코드 $NBASES · μ $MUS · $(date '+%F %T')" >> "$L"
t0=$(date +%s)
xargs -P 60 -L 1 sh -c 'node dms.js --cond mat --wt racld --seed "$1" --mu-assay "$0" --exact 1 --arms nbr --nbr-of "$2" --out "$3"' < "$J" >> "$L" 2>&1
echo "===== 끝 $(( $(date +%s) - t0 ))초 =====" >> "$L"
