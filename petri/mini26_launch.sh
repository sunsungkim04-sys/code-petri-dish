#!/bin/sh
# mini26_launch.sh — 해부 26 발사기 (로컬 맥 또는 lab101 · python3 + numpy 만 필요 · node 불필요)
#   판정:   PETRI_ROOT=mini26_out ./mini26_launch.sh
#   파일럿: PETRI_ROOT=pilot26 PETRI_SEEDS="26901" ./mini26_launch.sh
#   재현:   PETRI_ROOT=mini26_repro PETRI_SEEDS="26001" ./mini26_launch.sh   (판정 뒤 · 두 ρ 의 같은 시드를 다시 돌려 바이트 대조)
# 한 작업 = (ρ, 시드): 성장 6,000틱 + 규칙 2 × μ 8 × 팔 6(끼운 1 + 기준 5) × 3,000틱.
set -eu
cd "$(dirname "$0")"
PY=${PETRI_PY:-python3}
ROOT=${PETRI_ROOT:?PETRI_ROOT 필요}
SEEDS=${PETRI_SEEDS:-$(seq 26001 26024)}
RHOS=${PETRI_RHOS:-"8 32"}
PAR=${PETRI_PAR:-8}
if [ "$PAR" -gt 24 ]; then echo "병렬 $PAR > 24 — 거부"; exit 1; fi
for s in $SEEDS; do
  case "$ROOT" in
    pilot*) [ "$s" -ge 26901 ] && [ "$s" -le 26910 ] || { echo "파일럿 시드는 26901~26910 만: $s"; exit 1; } ;;
    *)      [ "$s" -ge 26001 ] && [ "$s" -le 26030 ] || { echo "판정 시드는 26001~26030 만: $s"; exit 1; } ;;
  esac
done
mkdir -p "$ROOT"
L="$ROOT.log"
J="$ROOT.jobs"
: > "$J"
for r in $RHOS; do for s in $SEEDS; do echo "$r $s" >> "$J"; done; done
echo "===== 시작 $(date '+%F %T') · 작업 $(wc -l < "$J") · 병렬 $PAR · $($PY -c 'import numpy,sys;print("python",sys.version.split()[0],"numpy",numpy.__version__)')" > "$L"
shasum -a 256 mini26.py 2>/dev/null >> "$L" || sha256sum mini26.py >> "$L"
xargs -P "$PAR" -n 2 sh -c "$PY mini26.py job --rho \"\$0\" --seed \"\$1\" --out \"$ROOT\" || echo \"ERROR 작업 \$0 \$1\"" < "$J" >> "$L" 2>&1
echo "===== 끝 $(date '+%F %T') · DONE $(grep -c 'qualified=' "$L") · ERROR $(grep -c '^ERROR' "$L")" >> "$L"
