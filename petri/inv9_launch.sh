#!/usr/bin/env bash
# 해부 9 — 기억하는 주사위 세계(A) · 재료 끔 세계(B)의 침입 μ 사다리. 사전등록 해부9-사전등록-2026-09-17.md §2 · §3
#   해부 6C · 7② 와 같은 mat 배경 20(본실험 시드 중 끝 우세 racld 인 첫 20 = 3~22) × μ {0 · .001 · .003 · .005 · .01}
#   × 팔 ref + neutral 10 + rascld · rsacld · racldx(--arms list · --exact 1) × 세계 {mem: --remember-die 1 · nm: --no-mat 1} → 200 작업
#   판정은 inv9_analyze.py (원래 = inv6/ · 재료 먼저 = inv7ff/ 기록과 짝).
# 마른 실행: PETRI_NBG=1 PETRI_MUS="0 0.01" PETRI_ROOT=smoke9 ./inv9_launch.sh
set -u
export NVM_DIR="$HOME/.nvm"; . "$NVM_DIR/nvm.sh"
cd "$HOME/petri" || exit 1
grep -qE "VERSION = '0\.3\.[45]'" sim.js || { echo "sim.js 가 v0.3.4/0.3.5 가 아니다 — 발사 중단"; exit 1; }
grep -q "remember-die" dms.js || { echo "dms.js 에 --remember-die 가 없다 — 발사 중단"; exit 1; }
grep -q "no-mat" dms.js || { echo "dms.js 에 --no-mat 가 없다 — 발사 중단"; exit 1; }
NBG=${PETRI_NBG:-20}
MUS=${PETRI_MUS:-"0 0.001 0.003 0.005 0.01"}
ROOT=${PETRI_ROOT:-inv9}
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
for world in mem nm; do
  for mu in $MUS; do
    tag=mu$(echo "$mu" | sed 's/\./p/')
    for s in $BGS; do echo "$mu $s $world ${ROOT}${world}/$tag" >> "$J"; done
  done
done
[ "$(awk 'NF!=4' "$J" | wc -l)" = 0 ] || { echo "작업 줄의 칸 수가 4 가 아니다 — 발사 중단"; exit 1; }
echo "배경 $BGS" >> "$L"
echo "작업 $(wc -l < "$J") · $(date '+%F %T')" >> "$L"
t0=$(date +%s)
xargs -P 60 -n 4 sh -c 'if [ "$2" = mem ]; then F="--remember-die 1"; else F="--no-mat 1"; fi; node dms.js --cond mat --wt racld --seed "$1" --arms list --codes rascld,rsacld,racldx --mu-assay "$0" --exact 1 $F --out "$3"' < "$J" >> "$L" 2>&1
echo "===== 끝 $(( $(date +%s) - t0 ))초 =====" >> "$L"
