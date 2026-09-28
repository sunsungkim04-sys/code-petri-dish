#!/usr/bin/env bash
# 해부 7 ③b — μ 0 침입에서 씨앗이 이기는 이유. 사전등록 해부7-사전등록-2026-09-16.md §3
#   침입: 해부 1 과 같은 mat 20 배경 × 팔 ref · rascld · rsacld · racldx × μ 0 × 3,000틱 (계통 기록 = 해부 1 과 같아야 함)
#   단일: racld · rascld · rsacld × 시드 801~820 × mat 밀도 8 × 20,000틱
#   → 20 + 60 = 80 작업. 판정은 invbud_analyze.py.
# 마른 실행: PETRI_NBG=1 PETRI_SEEDS="801" PETRI_ROOT=smoke7ib ./ib_launch.sh
set -u
export NVM_DIR="$HOME/.nvm"; . "$NVM_DIR/nvm.sh"
cd "$HOME/petri" || exit 1
grep -qE "VERSION = '0\.3\.[2345]'" sim.js || { echo "sim.js 가 v0.3.2~0.3.5 가 아니다 — 발사 중단"; exit 1; }
grep -q "INVBUD_VERSION = '1.0.0'" invbud.js || { echo "invbud.js 가 v1.0.0 이 아니다 — 발사 중단"; exit 1; }
NBG=${PETRI_NBG:-20}
SEEDS=${PETRI_SEEDS:-$(seq 801 820)}
ROOT=${PETRI_ROOT:-ib7}
L=${ROOT}.log; : > "$L"
J1=${ROOT}_inv_jobs.txt; : > "$J1"
J2=${ROOT}_mono_jobs.txt; : > "$J2"
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
for s in $BGS; do echo "$s $ROOT/inv" >> "$J1"; done
for code in racld rascld rsacld; do for s in $SEEDS; do echo "$code $s $ROOT/mono" >> "$J2"; done; done
echo "배경 $BGS" >> "$L"
echo "작업 $(( $(wc -l < "$J1") + $(wc -l < "$J2") )) · $(date '+%F %T')" >> "$L"
t0=$(date +%s)
cat "$J1" | xargs -P 20 -n 2 sh -c 'node invbud.js --seed "$0" --mu 0 --out "$1"' >> "$L" 2>&1 &
cat "$J2" | xargs -P 40 -n 3 sh -c 'node invbud.js --mono "$0" --seed "$1" --out "$2"' >> "$L" 2>&1 &
wait
echo "===== 끝 $(( $(date +%s) - t0 ))초 =====" >> "$L"
