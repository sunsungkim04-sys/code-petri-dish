#!/bin/sh
# mini27_launch.sh — 해부 27 발사기 (로컬 맥 또는 lab101 · python3 + numpy 만 필요)
#   판정:   PETRI_ROOT=mini27_out ./mini27_launch.sh                         (시드 27001~27024 · ρ 8 · 32)
#   파일럿: PETRI_ROOT=pilot27 PETRI_SEEDS="27901" ./mini27_launch.sh
#   재현:   PETRI_ROOT=mini27_repro PETRI_SEEDS="27001" ./mini27_launch.sh   (판정 뒤 · 바이트 대조)
#   G2a:    PETRI_ROOT=ident27 ./mini27_launch.sh                            (시드 26001 · glob 만 · 원 사다리 → mini26_out 과 바이트 대조)
# 분리 실행(로컬): PETRI_PY=/usr/bin/python3 PETRI_ROOT=mini27_out PETRI_PAR=8 nohup caffeinate -i ./mini27_launch.sh > /dev/null 2>&1 & echo $! > mini27_out.pid
# v1.1(검토 S3): 스레드 1 고정(lab101 규칙 · 로컬도 같게) · ρ 32 먼저(긴 작업 먼저 → 꼬리 줄임) · 같은 ROOT 다시 돌리면 옛 로그를 .prev-시각 으로 보존
# 한 작업 = (ρ, 시드): 기하 셋(glob · loc · place)의 성장 6,000틱 + 팔(glob 180 · loc 180 · place 54) × 3,000틱.
set -eu
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
PY=${PETRI_PY:-python3}
ROOT=${PETRI_ROOT:?PETRI_ROOT 필요}
RHOS=${PETRI_RHOS:-"32 8"}
PAR=${PETRI_PAR:-8}
EXTRA=""
case "$ROOT" in
  ident*) SEEDS=${PETRI_SEEDS:-26001}; EXTRA="--geoms glob --mus 0,0.001,0.002,0.003,0.004,0.005,0.0075,0.01" ;;
  *)      SEEDS=${PETRI_SEEDS:-$(seq 27001 27024)} ;;
esac
if [ "$PAR" -gt 24 ]; then echo "병렬 $PAR > 24 — 거부"; exit 1; fi
for s in $SEEDS; do
  case "$ROOT" in
    pilot*) [ "$s" -ge 27901 ] && [ "$s" -le 27910 ] || { echo "파일럿 시드는 27901~27910 만: $s"; exit 1; } ;;
    ident*) [ "$s" = 26001 ] || { echo "G2a 시드는 26001 만: $s"; exit 1; } ;;
    *)      [ "$s" -ge 27001 ] && [ "$s" -le 27030 ] || { echo "판정 시드는 27001~27030 만: $s"; exit 1; } ;;
  esac
done
mkdir -p "$ROOT"
L="$ROOT.log"
J="$ROOT.jobs"
[ -f "$L" ] && mv "$L" "$L.prev-$(date '+%H%M%S')"
: > "$J"
for r in $RHOS; do for s in $SEEDS; do echo "$r $s" >> "$J"; done; done
echo "===== 시작 $(date '+%F %T') · 작업 $(wc -l < "$J") · 병렬 $PAR · $($PY -c 'import numpy,sys;print("python",sys.version.split()[0],"numpy",numpy.__version__)')" > "$L"
shasum -a 256 mini27.py 2>/dev/null >> "$L" || sha256sum mini27.py >> "$L"
xargs -P "$PAR" -n 2 sh -c "$PY mini27.py job --rho \"\$0\" --seed \"\$1\" --out \"$ROOT\" $EXTRA || echo \"ERROR 작업 \$0 \$1\"" < "$J" >> "$L" 2>&1
echo "===== 끝 $(date '+%F %T') · DONE $(grep -c 'qualified=' "$L") · ERROR $(grep -c '^ERROR' "$L")" >> "$L"
