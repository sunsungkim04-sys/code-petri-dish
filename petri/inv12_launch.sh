#!/usr/bin/env bash
# 해부 12 — 해석 모형의 입력을 잰다. 사전등록 해부12-사전등록-2026-09-20.md §2~§3
#   P  이웃 DMS:   dms --arms nbr --nbr-of CODE · μ 0 · CODE ∈ {racld · rascld · rsacld · racldx} × mat 배경 3~22 → 80 (팔 110~140 · 작업당 30~40분)
#   Q  자식 세기:  dms --arms list --codes rascld,rsacld,racldx --kids 1 · μ {0 · 0.1 · 0.3 · 0.5 · 1%} × 같은 배경 → 100
#                  (팔 기록은 해부 6C inv6 과 같아야 한다 — 판정기가 대조한다)
#   판정: model12_analyze.py
# 마른 실행: PETRI_NBG=1 PETRI_BG_FROM=23 PETRI_MUS="0 0.003" PETRI_ROOT=pilot12 ./inv12_launch.sh
set -u
export NVM_DIR="$HOME/.nvm"; . "$NVM_DIR/nvm.sh"
cd "${PETRI_HOME:-$HOME/petri}" || exit 1
grep -q "VERSION = '0.3.5'" sim.js || { echo "sim.js 가 v0.3.5 가 아니다 — 발사 중단"; exit 1; }
grep -q "kidsOn" dms.js || { echo "dms.js 에 --kids 가 없다 — 발사 중단"; exit 1; }
NBG=${PETRI_NBG:-20}
MUS=${PETRI_MUS:-"0 0.001 0.003 0.005 0.01"}
NCODES=${PETRI_NCODES:-"racld rascld rsacld racldx"}
PARTS=${PETRI_PARTS:-"P Q"}
ROOT=${PETRI_ROOT:-inv12}
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
for p in $PARTS; do case $p in
  P) for code in $NCODES; do for s in $BGS; do echo "P 0 $s $code ${ROOT}nbr_$code/mu0" >> "$J"; done; done;;
  Q) for mu in $MUS; do tag=mu$(echo "$mu" | sed 's/\./p/'); for s in $BGS; do echo "Q $mu $s - ${ROOT}kids/$tag" >> "$J"; done; done;;
esac; done
[ "$(awk 'NF!=5' "$J" | wc -l)" = 0 ] || { echo "작업 줄의 칸 수가 5 가 아니다 — 발사 중단"; exit 1; }
echo "배경 $BGS" >> "$L"
echo "작업 $(wc -l < "$J") (P $(grep -c '^P ' "$J") · Q $(grep -c '^Q ' "$J")) · $(date '+%F %T')" >> "$L"
t0=$(date +%s)
xargs -P 60 -L 1 sh -c '
  if [ "$0" = P ]; then F="--arms nbr --nbr-of $3"; else F="--arms list --codes rascld,rsacld,racldx --kids 1"; fi
  node dms.js --cond mat --wt racld --seed "$2" --mu-assay "$1" --exact 1 $F --out "$4"' < "$J" >> "$L" 2>&1
echo "===== 끝 $(( $(date +%s) - t0 ))초 =====" >> "$L"
