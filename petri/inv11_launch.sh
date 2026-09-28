#!/usr/bin/env bash
# 해부 11 — 수명 사다리의 시간 예산 · 밀도 개입 · 침입 Δ 의 수명 의존 · 고리 길이 사다리. 사전등록 해부11-사전등록-2026-09-20.md §2~§5
#   A  시간 예산:       budget.js v1.2.0 · μ 0 · 수명 {300 · 600 · 1200}(age0 = ageVar) × 밀도 {8 · DHI} × 코드 3 × 시드 1401~1420 → 360
#   A′ 부풀림 직접:     lineage --spec · μ 0.1% · 수명 {300 · 1200} × 밀도 {8 · DHI} × 코드 3 × 시드 1421~1440 → 240 (밀도별 폴더 — 파일 꼬리에 밀도가 없다)
#   B  침입 Δ 수명:     dms --ancestor racld --age0 A --age-var A · 수명 {600 · 1200} × μ 5 × 배경 시드 1101~1120 → 200
#                       배경 · 팔 길이를 수명 배수 f = A/300 로 늘린다(grow 5만·f · ticks 3000·f · every 250·f). 수명 300 짝 = 해부 10 G(inv10anc)
#   C  고리 길이 사다리: dms · mat 배경 3~22 · μ 0 · 코드 rascld(고리 2) · rsacld(3) · racldx(4) · racldn(4 + 한 틱) · nracld(5) · racldr(4 + 묶인 글자) → 20
#   A · A′ 은 접시 길이도 수명 배수 f 로 늘린다(파일럿: 수명 1200 은 2만 틱에 정상 상태에 못 닿는다) — ticks 2만·f · 칸 500·f · 창 시작 5000·f
#   판정: an11_analyze.py
# 마른 실행: PETRI_NBG=1 PETRI_MUS="0 0.01" PETRI_SEEDS="1491" PETRI_LSEEDS="1491" PETRI_GSEEDS="1191" PETRI_ROOT=smoke11 ./inv11_launch.sh
set -u
export NVM_DIR="$HOME/.nvm"; . "$NVM_DIR/nvm.sh"
cd "${PETRI_HOME:-$HOME/petri}" || exit 1
grep -q "VERSION = '0.3.5'" sim.js || { echo "sim.js 가 v0.3.5 가 아니다 — 발사 중단"; exit 1; }
grep -q "LINEAGE_VERSION = '1.5.0'" lineage.js || { echo "lineage.js 가 v1.5.0 이 아니다 — 발사 중단"; exit 1; }
grep -q "BUDGET_VERSION = '1.2.0'" budget.js || { echo "budget.js 가 v1.2.0 이 아니다 — 발사 중단"; exit 1; }
grep -q "age0V" dms.js || { echo "dms.js 에 --age0 이 없다 — 발사 중단"; exit 1; }
NBG=${PETRI_NBG:-20}
MUS=${PETRI_MUS:-"0 0.001 0.003 0.005 0.01"}
SEEDS=${PETRI_SEEDS:-$(seq 1401 1420)}
LSEEDS=${PETRI_LSEEDS:-$(seq 1421 1440)}
GSEEDS=${PETRI_GSEEDS:-$(seq 1101 1120)}
AGES=${PETRI_AGES:-"300 600 1200"}
BAGES=${PETRI_BAGES:-"600 1200"}
DENS=${PETRI_DENS:-"8 14"}
PARTS=${PETRI_PARTS:-"A L B C"}
ROOT=${PETRI_ROOT:-inv11}
case "$DENS" in *UNSET*) echo "밀도 개입 값(DHI)이 정해지지 않았다 — 발사 중단"; exit 1;; esac
L=${ROOT}.log; : > "$L"
J1=${ROOT}_dms_jobs.txt; : > "$J1"
J2=${ROOT}_lin_jobs.txt; : > "$J2"
J3=${ROOT}_bud_jobs.txt; : > "$J3"
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
  A) for s in $SEEDS; do for code in racld rascld rsacld; do for d in $DENS; do for age in $AGES; do
       echo "$code $s $d $age ${ROOT}bud" >> "$J3"
     done; done; done; done;;
  L) for s in $LSEEDS; do for code in racld rascld rsacld; do for d in $DENS; do for age in 300 1200; do
       echo "$code 0.001 $s $d $age ${ROOT}lin_d$d" >> "$J2"
     done; done; done; done;;
  B) for mu in $MUS; do tag=mu$(echo "$mu" | sed 's/\./p/'); for age in $BAGES; do
       n=0; for s in $GSEEDS; do n=$((n+1)); [ $n -le $NBG ] && echo "$mu $s age$age ${ROOT}age$age/$tag" >> "$J1"; done
     done; done;;
  C) for s in $BGS; do echo "0 $s loop ${ROOT}loop/mu0" >> "$J1"; done;;
esac; done
[ "$(awk 'NF!=4' "$J1" | wc -l)" = 0 ] || { echo "dms 작업 줄의 칸 수가 4 가 아니다 — 발사 중단"; exit 1; }
[ "$(awk 'NF!=6' "$J2" | wc -l)" = 0 ] || { echo "lineage 작업 줄의 칸 수가 6 이 아니다 — 발사 중단"; exit 1; }
[ "$(awk 'NF!=5' "$J3" | wc -l)" = 0 ] || { echo "budget 작업 줄의 칸 수가 5 가 아니다 — 발사 중단"; exit 1; }
echo "배경 $BGS" >> "$L"
echo "작업 dms $(wc -l < "$J1") · lineage $(wc -l < "$J2") · budget $(wc -l < "$J3") · $(date '+%F %T')" >> "$L"
t0=$(date +%s)
# 긴 작업(B · 수명 1200 은 작업당 20분 남짓)을 먼저 채운다 — 세 묶음을 한 줄에 세워 60 코어를 나눈다
{
  { grep ' age1200 ' "$J1"; grep ' age600 ' "$J1"; grep -v ' age1200 \| age600 ' "$J1"; } | sed 's/^/D /'
  sed 's/^/L /' "$J2"
  sed 's/^/U /' "$J3"
} > "${ROOT}_all_jobs.txt"
xargs -P 60 -L 1 sh -c '
  kind="$0"
  if [ "$kind" = D ]; then
    case "$3" in
      loop) F="--codes rascld,rsacld,racldx,racldn,nracld,racldr";;
      age*) A=${3#age}; f=$((A/300)); F="--codes rascld,rsacld,racldx --ancestor racld --age0 $A --age-var $A --grow $((50000*f)) --ticks $((3000*f)) --every $((250*f))";;
    esac
    node dms.js --cond mat --wt racld --seed "$2" --arms list --mu-assay "$1" --exact 1 $F --out "$4"
  elif [ "$kind" = L ]; then
    f=$(($5/300)); node lineage.js --code "$1" --mu "$2" --seed "$3" --find-first 0 --spec 1 --density "$4" --age0 "$5" --age-var "$5" --ticks $((20000*f)) --win0 $((5000*f)) --out "$6"
  else
    f=$(($4/300)); node budget.js --cond mat --code "$1" --seed "$2" --density "$3" --age0 "$4" --age-var "$4" --ticks $((20000*f)) --bin $((500*f)) --out "$5"
  fi' < "${ROOT}_all_jobs.txt" >> "$L" 2>&1
echo "===== 끝 $(( $(date +%s) - t0 ))초 =====" >> "$L"
