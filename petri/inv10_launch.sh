#!/usr/bin/env bash
# 해부 10 — 규칙 의존 사다리(우주선 끔 · 확산 · 수명) + 둘째 시조. 사전등록 해부10-사전등록-2026-09-18.md §2~§5
#   C  우주선 끔 침입:  mat 배경 3~22 × μ 5 × {c0: 원래 규칙 + --cosmic 0 · ffc0: --find-first 1 --cosmic 0} → 200 (dms.js · 짝 = inv6 · inv7ff)
#   G  둘째 시조 침입:  --ancestor racld 로 키운 배경(시드 1101~1120) × μ 5 · 원래 규칙 → 100 (본실험 기록 없음 → 체크섬 no-main)
#   D  확산 사다리:     lineage --spec · 밀도 8 · μ 1% · diffuse {0.05 · 0.15 · 0.5} × 코드 3 × 시드 1201~1220 → 180
#   E  수명 사다리:     lineage --spec · 밀도 8 · μ {0.1% · 1%} · age {300 · 600 · 1200}(age0 = ageVar) × 코드 3 × 시드 1301~1320 → 360
#   판정: inv10_analyze.py(C · G) · lin10_analyze.py(D · E)
# 마른 실행: PETRI_NBG=1 PETRI_MUS="0 0.01" PETRI_SEEDS="1201" PETRI_ROOT=smoke10 ./inv10_launch.sh
set -u
export NVM_DIR="$HOME/.nvm"; . "$NVM_DIR/nvm.sh"
cd "$HOME/petri" || exit 1
grep -q "VERSION = '0.3.5'" sim.js || { echo "sim.js 가 v0.3.5 가 아니다 — 발사 중단"; exit 1; }
grep -q "LINEAGE_VERSION = '1.5.0'" lineage.js || { echo "lineage.js 가 v1.5.0 이 아니다 — 발사 중단"; exit 1; }
grep -q "cosmicOff" dms.js || { echo "dms.js 에 --cosmic 이 없다 — 발사 중단"; exit 1; }
NBG=${PETRI_NBG:-20}
MUS=${PETRI_MUS:-"0 0.001 0.003 0.005 0.01"}
SEEDS=${PETRI_SEEDS:-$(seq 1201 1220)}
ASEEDS=${PETRI_ASEEDS:-$(seq 1301 1320)}
GSEEDS=${PETRI_GSEEDS:-$(seq 1101 1120)}
ROOT=${PETRI_ROOT:-inv10}
L=${ROOT}.log; : > "$L"
J1=${ROOT}_dms_jobs.txt; : > "$J1"
J2=${ROOT}_lin_jobs.txt; : > "$J2"
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
  tag=mu$(echo "$mu" | sed 's/\./p/')
  for s in $BGS; do
    echo "$mu $s c0 ${ROOT}c0/$tag" >> "$J1"
    echo "$mu $s ffc0 ${ROOT}ffc0/$tag" >> "$J1"
  done
  n=0; for s in $GSEEDS; do n=$((n+1)); [ $n -le $NBG ] && echo "$mu $s anc ${ROOT}anc/$tag" >> "$J1"; done
done
for s in $SEEDS; do for code in racld rascld rsacld; do for df in 0.05 0.15 0.5; do
  echo "$code 0.01 $s df $df ${ROOT}df" >> "$J2"
done; done; done
for s in $ASEEDS; do for code in racld rascld rsacld; do for mu in 0.001 0.01; do for age in 300 600 1200; do
  echo "$code $mu $s age $age ${ROOT}age" >> "$J2"
done; done; done; done
[ "$(awk 'NF!=4' "$J1" | wc -l)" = 0 ] || { echo "dms 작업 줄의 칸 수가 4 가 아니다 — 발사 중단"; exit 1; }
[ "$(awk 'NF!=6' "$J2" | wc -l)" = 0 ] || { echo "lineage 작업 줄의 칸 수가 6 이 아니다 — 발사 중단"; exit 1; }
echo "배경 $BGS" >> "$L"
echo "작업 dms $(wc -l < "$J1") · lineage $(wc -l < "$J2") · $(date '+%F %T')" >> "$L"
t0=$(date +%s)
xargs -P 60 -n 4 sh -c 'case "$2" in c0) F="--cosmic 0";; ffc0) F="--find-first 1 --cosmic 0";; anc) F="--ancestor racld";; esac; node dms.js --cond mat --wt racld --seed "$1" --arms list --codes rascld,rsacld,racldx --mu-assay "$0" --exact 1 $F --out "$3"' < "$J1" >> "$L" 2>&1
xargs -P 60 -n 6 sh -c 'if [ "$3" = df ]; then F="--diffuse $4"; else F="--age0 $4 --age-var $4"; fi; node lineage.js --code "$0" --mu "$1" --seed "$2" --find-first 0 --spec 1 $F --out "$5"' < "$J2" >> "$L" 2>&1
echo "===== 끝 $(( $(date +%s) - t0 ))초 =====" >> "$L"
