#!/usr/bin/env bash
# 해부 6C — 촘촘한 침입 사다리. 사전등록 해부6-사전등록-2026-09-16.md §4
#   배경: 해부 1 · 2A 와 같은 mat 20 접시(본실험 시드 1~100 중 끝 우세 racld 인 앞 20개)
#   팔 = ref + neutral 10 + rascld · rsacld · racldx  ·  μ 0·0.1·0.2·0.3·0.4·0.5·0.75·1%  ·  --exact 1
#   → 8 × 20 = 160 작업. μ 0 · 1% 는 해부 1 · 2A 기록과 궤적이 같아야 한다(회귀 검사). 판정은 inv_analyze.py.
# 마른 실행: PETRI_NBG=1 PETRI_MUS="0 0.01" PETRI_ROOT=smoke6inv ./inv_launch.sh
set -u
export NVM_DIR="$HOME/.nvm"; . "$NVM_DIR/nvm.sh"
cd "$HOME/petri" || exit 1
grep -qE "VERSION = '0\.3\.[2345]'" sim.js || { echo "sim.js 가 v0.3.2~0.3.5 가 아니다 — 발사 중단"; exit 1; }
grep -q "const exact = args.exact === '1';" dms.js || { echo "dms.js 에 --exact 가 없다 — 발사 중단"; exit 1; }
NBG=${PETRI_NBG:-20}
MUS=${PETRI_MUS:-0 0.001 0.002 0.003 0.004 0.005 0.0075 0.01}
ROOT=${PETRI_ROOT:-inv6}
L=${ROOT}.log; : > "$L"
J=${ROOT}_jobs.txt; : > "$J"
BGS=$(NBG=$NBG python3 - <<'PYEOF'
import json, os
nbg = int(os.environ["NBG"]); picked = []
for s in range(1, 101):
    z = json.load(open(f"main/mat_d8_mu0p01_s{s:05d}.json"))
    if z["extinct_at"] < 0 and z["final_counts"] and z["final_counts"][0][0] == "racld":
        picked.append(s)
    if len(picked) == nbg:
        break
print(" ".join(map(str, picked)))
PYEOF
)
for mu in $MUS; do
  tag=mu$(echo "$mu" | sed 's/^0$/0/; s/\./p/')
  for s in $BGS; do echo "$mu $s $ROOT/$tag" >> "$J"; done
done
echo "배경 $BGS" >> "$L"
echo "작업 $(wc -l < "$J") · $(date '+%F %T')" >> "$L"
t0=$(date +%s)
xargs -P 60 -n 3 sh -c 'node dms.js --cond mat --wt racld --seed "$1" --arms list --codes rascld,rsacld,racldx --mu-assay "$0" --exact 1 --out "$2"' < "$J" >> "$L" 2>&1
echo "===== 끝 $(( $(date +%s) - t0 ))초 =====" >> "$L"
