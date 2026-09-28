#!/usr/bin/env bash
# 해부 7 ② 인과 — 재료 먼저 규칙으로 침입 팔을 돌린다. 사전등록 해부7-사전등록-2026-09-16.md §2
#   배경: 해부 6C 와 같은 mat 20 접시(배경은 원래 규칙으로 키움 · 본실험 체크섬 대조)
#   팔 = ref + neutral 10 + rascld · rsacld · racldx · --exact 1 --find-first 1 · μ 0 · 0.1 · 0.3 · 0.5 · 1%
#   → 5 × 20 = 100 작업. 원래 규칙 쪽은 해부 6C(inv6/)의 같은 배경 · 같은 μ 를 짝으로 쓴다.
# 마른 실행: PETRI_NBG=1 PETRI_MUS="0 0.01" PETRI_ROOT=smoke7inv ./ff_inv_launch.sh
set -u
export NVM_DIR="$HOME/.nvm"; . "$NVM_DIR/nvm.sh"
cd "$HOME/petri" || exit 1
grep -qE "VERSION = '0\.3\.[345]'" sim.js || { echo "sim.js 가 v0.3.3~0.3.5 가 아니다 — 발사 중단"; exit 1; }
grep -q "const findFirst = args\['find-first'\] === '1';" dms.js || { echo "dms.js 에 --find-first 가 없다 — 발사 중단"; exit 1; }
NBG=${PETRI_NBG:-20}
MUS=${PETRI_MUS:-0 0.001 0.003 0.005 0.01}
ROOT=${PETRI_ROOT:-inv7ff}
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
xargs -P 60 -n 3 sh -c 'node dms.js --cond mat --wt racld --seed "$1" --arms list --codes rascld,rsacld,racldx --mu-assay "$0" --exact 1 --find-first 1 --out "$2"' < "$J" >> "$L" 2>&1
echo "===== 끝 $(( $(date +%s) - t0 ))초 =====" >> "$L"
