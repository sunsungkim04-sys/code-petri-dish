#!/bin/sh
# mini31_launch.sh — 해부 31 발사기 (사전등록 `해부31-사전등록-2026-10-06.md` §3)
#   파일럿: PETRI_ROOT=pilot31 PETRI_SEEDS="31901 31902" ./mini31_launch.sh
#   판정:   PETRI_ROOT=mini31_out ./mini31_launch.sh                        (시드 31001~31024)
#   재현:   PETRI_ROOT=mini31_repro PETRI_SEEDS="31001" ./mini31_launch.sh  (G2 · 바이트 대조)
#   G4:     PETRI_ROOT=ident31 ./mini31_launch.sh                           (mini26 을 같은 시드·같은 사다리로)
# 분리 실행: PETRI_ROOT=mini31_out PETRI_PAR=8 nohup ./mini31_launch.sh > /dev/null 2>&1 & echo $! > mini31_out.pid
#
# 🔑 동결된 `mini27_launch.sh` 를 고치지 않고 새로 썼다. 그쪽은 시드(27001~27030)와
#    사다리(ident 만 --mus 를 준다)를 코드에 박고 거부한다. 해부 31 은 시드 31001~ 이고
#    사다리가 0…4% 라 그 검사에 걸린다. 구조·스레드 고정·로그 형식은 그대로 따랐다.
#
# 해부 27 과 다른 점 셋
#   ① μ 사다리 = 원 사다리 × 4.0 → 0 … 4%  (사전등록 §3 · 예측 2.14% 를 밟지 않는다)
#   ② 기하가 밀도마다 다르다 — ρ8 은 loc+glob, ρ32 는 loc 만. place 는 안 돌린다
#   ③ G4 는 mini26.py 를 같은 시드·같은 사다리로 돌려 glob 과 바이트 대조한다
#
# 한 작업 = (ρ, 시드, 기하). 작업 줄이 3 칸이므로 xargs 도 -n 3 이다.
#   🚩 09-16 에 7 칸을 -n 6 으로 읽어 결과 11 개가 엉뚱한 곳에 떨어졌다. 칸 수와 -n 을 꼭 맞출 것.
set -eu
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
PY=${PETRI_PY:-python3}
ROOT=${PETRI_ROOT:?PETRI_ROOT 필요}
PAR=${PETRI_PAR:-8}

MUS31="0,0.004,0.008,0.012,0.016,0.020,0.030,0.040"

case "$ROOT" in
  ident*) SEEDS=${PETRI_SEEDS:-31001}; SCRIPT=mini26.py ;;
  pilot*) SEEDS=${PETRI_SEEDS:?파일럿은 시드를 명시할 것}; SCRIPT=mini27.py ;;
  *)      SEEDS=${PETRI_SEEDS:-$(seq 31001 31024)}; SCRIPT=mini27.py ;;
esac

for s in $SEEDS; do
  case "$ROOT" in
    pilot*) [ "$s" -ge 31901 ] && [ "$s" -le 31910 ] || { echo "파일럿 시드는 31901~31910 만: $s"; exit 1; } ;;
    ident*) [ "$s" = 31001 ] || { echo "G4 시드는 31001 만: $s"; exit 1; } ;;
    *)      [ "$s" -ge 31001 ] && [ "$s" -le 31030 ] || { echo "판정 시드는 31001~31030 만: $s"; exit 1; } ;;
  esac
done

if [ "$PAR" -gt 24 ]; then echo "병렬 $PAR > 24 — 거부"; exit 1; fi
[ -f "$SCRIPT" ] || { echo "$SCRIPT 가 없다 — 올릴 것"; exit 1; }

# 🔑 판 고정 — 해부 27 판정과 같은 python·numpy 여야 한다.
#    10-06 에 기본 python3(3.8.10 · numpy 1.22.0)으로 파일럿이 돌았다. 바이트 동일성 가드(G2 · G4)는
#    numpy 판을 건너 성립하지 않으므로, 틀린 판으로 돌면 가드가 조용히 거짓이 된다.
#    lab101 에서 맞는 것: PETRI_PY=$HOME/miniforge3/bin/python3
NEED_NP=${PETRI_NUMPY:-2.4.3}
GOT=$($PY -c 'import numpy;print(numpy.__version__)' 2>/dev/null || echo none)
if [ "$GOT" != "$NEED_NP" ]; then
  echo "numpy $GOT ≠ $NEED_NP ($PY) — 거부."
  echo "  해부 27 판정은 python 3.13.12 · numpy 2.4.3 이었다. 바이트 동일성 가드가 깨진다."
  echo "  lab101: PETRI_PY=\$HOME/miniforge3/bin/python3 ./mini31_launch.sh"
  echo "  (일부러 다른 판으로 돌리려면 PETRI_NUMPY=$GOT 를 명시할 것 — 그 사실을 기록에 남긴다)"
  exit 1
fi

mkdir -p "$ROOT"
L="$ROOT.log"
J="$ROOT.jobs"
[ -f "$L" ] && mv "$L" "$L.prev-$(date '+%H%M%S')"
: > "$J"
# ρ32 를 먼저(긴 작업 먼저 → 꼬리 줄임). 기하는 밀도마다 다르다.
for s in $SEEDS; do
  case "$ROOT" in
    ident*) echo "8 $s -"          >> "$J" ;;   # mini26 은 전역 풀뿐 — 기하 인자가 없다
    *)      echo "32 $s loc"       >> "$J"
            echo "8 $s loc,glob"   >> "$J" ;;
  esac
done

echo "===== 시작 $(date '+%F %T') · 작업 $(wc -l < "$J") · 병렬 $PAR · $($PY -c 'import numpy,sys;print("python",sys.version.split()[0],"numpy",numpy.__version__)')" > "$L"
echo "사다리 $MUS31 · 스크립트 $SCRIPT" >> "$L"
{ shasum -a 256 "$SCRIPT" 2>/dev/null || sha256sum "$SCRIPT"; } >> "$L"

if [ "$SCRIPT" = mini26.py ]; then
  # G4 — mini26 은 --geoms 를 받지 않는다(전역 풀만)
  xargs -P "$PAR" -n 3 sh -c \
    "$PY mini26.py job --rho \"\$0\" --seed \"\$1\" --mus $MUS31 --out \"$ROOT\" || echo \"ERROR 작업 \$0 \$1\"" \
    < "$J" >> "$L" 2>&1
else
  xargs -P "$PAR" -n 3 sh -c \
    "$PY mini27.py job --rho \"\$0\" --seed \"\$1\" --geoms \"\$2\" --mus $MUS31 --out \"$ROOT\" || echo \"ERROR 작업 \$0 \$1 \$2\"" \
    < "$J" >> "$L" 2>&1
fi

echo "===== 끝 $(date '+%F %T') · DONE $(grep -c 'qualified=' "$L") · ERROR $(grep -c '^ERROR' "$L")" >> "$L"
