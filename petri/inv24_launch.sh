#!/usr/bin/env bash
# 해부 24 — 직접 계수(invbud R · T)를 셋째 배경 세트(해부 20 짝 배경 100 … 76)에서. 사전등록 해부24-사전등록-2026-10-01.md §2 · §8
#   계측기 = invbud.js v1.0.0 그대로(해부 14 IB 줄과 같은 명령 · 새 패치 없음):
#     node invbud.js --seed s --mu 0 --codes <밑 코드 B 의 n 삽입 전부> --out <ROOT>ib_<B>/mu0
#   둘째 세트(배경 26~50 의 20 · inv15Eib_*)는 이미 있다 → 새로 돌리지 않는다(사전등록된 재분석 · pairs24_analyze.py --set 2).
#   부분(PETRI_PARTS):
#     IB   본 작업. 밑 코드 = _RESULT_screen20.json 의 pass(= 해부 20 짝 단계 밑 코드 8) · 배경 = _RESULT_pairs20.json 의 seeds(20) → 160 작업
#     D    (파일럿 전용) 같은 배경 · 밑 코드의 dms.js 목록 모드 μ 0 — 해부 20 짝 단계와 같은 명령. G0a(팔 기록 대조) 경로 연습용.
#          본 세트에는 해부 20 의 inv20P_* 가 이미 있으므로 D 를 돌리지 않는다.
#     REG  (파일럿 전용) 회귀: 지금 invbud.js 로 기존 출력 두 개를 다시 만든다 — 배경 3 · racld(inv14ib_racld) · 배경 26 · rascled(inv15Eib_rascled)
#          → pairs24_analyze.py --reg 가 팔 전부(bins · series · stopped · inj)를 바이트 대조(G0r). 판정 통계량은 계산하지 않는다.
# 본 발사(검토 · 동결 뒤에만):  PETRI_ROOT=inv24 PETRI_PAR=12 setsid nohup ./inv24_launch.sh > inv24_launcher.out 2>&1 < /dev/null &
# 마른 실행(작업 목록만 · 아무것도 안 돌림):  PETRI_DRYRUN=1 ./inv24_launch.sh
# 파일럿(겹치지 않는 배경 53 · 54):  PETRI_ROOT=pilot24 PETRI_SEEDS="53 54" PETRI_BASES="racld rascled reasccld" PETRI_PARTS="IB D REG" PETRI_PAR=12 ./inv24_launch.sh
# 재개(검토 개정 M3 · 사전등록 §9b — rc≠0 · 출력 누락 · 무인 재부팅 뒤에만):  PETRI_ROOT=inv24 PETRI_RESUME=1 PETRI_PAR=12 setsid nohup ./inv24_launch.sh >> inv24_launcher.out 2>&1 < /dev/null &
#   → 출력 파일이 없거나 깨진(JSON 아님 · 배경/팔 수 다름) 작업만 같은 명령 · 같은 폴더로 다시 돌린다(결정론적). 목록은 inv24_resume_<n>.txt 로 남고 inv24.log 에 이어 적는다.
#   재부팅 뒤라면 사람이 새 uptime -s 를 PETRI_BOOT_NOW="…" 로 준다(재개 모드에서만 · 실제 uptime -s 와 같아야 함 · 로그에 남김).
# 검토 개정(10-01 · 독립 검토 · 결과 전): M1 비-inv24 루트는 PETRI_SEEDS 필수 · 'pilot' 으로 시작 · M3 재개 모드 · S1 다른 세션 node 수 + PAR ≤ 36 · S2 본 루트는 동결 · 원자료 manifest 검사
set -u
export NVM_DIR="$HOME/.nvm"; . "$NVM_DIR/nvm.sh"
cd "${PETRI_HOME:-$HOME/petri}" || exit 1
BOOT_EXPECT="2026-09-26 23:34:08"
RESUME=${PETRI_RESUME:-0}
if [ "$RESUME" = 1 ] && [ -n "${PETRI_BOOT_NOW:-}" ]; then
  [ "$PETRI_BOOT_NOW" = "$(uptime -s)" ] || { echo "PETRI_BOOT_NOW($PETRI_BOOT_NOW) ≠ 지금 uptime -s $(uptime -s) — 재개 중단"; exit 1; }
  echo "재개 · 재부팅 뒤 — 사람이 준 새 uptime -s $PETRI_BOOT_NOW (원래 $BOOT_EXPECT)"; BOOT_EXPECT=$PETRI_BOOT_NOW
fi
[ "$(uptime -s)" = "$BOOT_EXPECT" ] || { echo "uptime -s $(uptime -s) ≠ $BOOT_EXPECT — 무인 재부팅 의심 · 발사 중단(사람에게 보고)"; exit 1; }
grep -q "VERSION = '0.3.5'" sim.js || { echo "sim.js 가 v0.3.5 가 아니다 — 발사 중단"; exit 1; }
grep -q "INVBUD_VERSION = '1.0.0'" invbud.js || { echo "invbud.js 가 v1.0.0 이 아니다 — 발사 중단"; exit 1; }
IB_SHA=28bfc9c11f27dcdd3e911fa34a3c8d066e2118bc93d08b4a732e6cea783e9246   # 10-01 서버 = 볼트 · 해부 14 후속(--by-letter) · 18(--density 등) 뒤 · 기본값 출력은 G0r 로 확인
[ "$(sha256sum invbud.js | cut -d' ' -f1)" = "$IB_SHA" ] || { echo "invbud.js sha256 이 사전등록 값($IB_SHA)과 다르다 — 발사 중단"; exit 1; }
node --version | grep -q '^v20' || { echo "node 가 v20 이 아니다($(node --version)) — nvm 확인 · 발사 중단"; exit 1; }
PAR=${PETRI_PAR:-12}
[ "$PAR" -ge 1 ] && [ "$PAR" -le 12 ] || { echo "병렬 $PAR — 해부 24 코어 예산은 1~12 · 발사 중단"; exit 1; }
ROOT=${PETRI_ROOT:-inv24}
PARTS=${PETRI_PARTS:-IB}
DRY=${PETRI_DRYRUN:-0}
# M1 — 본 루트가 아니면 배경을 반드시 따로 준다(빈 PETRI_SEEDS 로 본 배경 76~100 을 다른 루트에 미리 만드는 길을 닫음) · 루트 이름은 pilot 으로 시작
if [ "$ROOT" != inv24 ]; then
  [ -n "${PETRI_SEEDS:-}" ] || { echo "루트 $ROOT(≠ inv24) 에는 PETRI_SEEDS 가 필수 — 발사 중단"; exit 1; }
  case "$ROOT" in pilot*) ;; *) echo "루트 $ROOT — 본 루트 inv24 가 아니면 'pilot' 으로 시작해야 한다 — 발사 중단"; exit 1;; esac
fi
# S2 — 본 루트는 동결 파일과 원자료 manifest 가 맞아야 한다
if [ "$ROOT" = inv24 ] && [ "$DRY" != 1 ]; then
  sha256sum -c --quiet dms24_frozen.sha256 || { echo "dms24_frozen.sha256 검사 실패 — 발사 중단"; exit 1; }
  sha256sum -c --quiet dms24_rawdata.sha256 || { echo "dms24_rawdata.sha256(원자료 manifest) 검사 실패 — 발사 중단"; exit 1; }
fi
# S1 — 다른 세션의 node invbud/dms/sim 수 + 이 발사 병렬 ≤ 36(해부 24 ≤ 12 · 해부 25 ≤ 24)
OTHER=$(pgrep -c -f '(^|/)node (invbud|dms|sim)\.js' || true); OTHER=${OTHER:-0}
[ $(( OTHER + PAR )) -le 36 ] || { echo "다른 세션 node $OTHER + 병렬 $PAR > 36 — 코어 예산 초과 · 발사 중단(나중에 다시)"; exit 1; }

# 본 세트 배경 · 밑 코드: 해부 20 결과 파일에서 읽고, 100 에서 내려가며 다시 세어 같은지 확인한다
MAIN_SEEDS=$(python3 - <<'PYEOF'
import json
z = json.load(open("_RESULT_pairs20.json"))
picked = []
for s in range(100, 52, -1):
    m = json.load(open("main/mat_d8_mu0p01_s%05d.json" % s))
    if m["extinct_at"] < 0 and m["final_counts"] and m["final_counts"][0][0] == "racld":
        picked.append(s)
assert sorted(z["seeds"]) == sorted(picked[:20]) and len(z["seeds"]) == 20 and not z["pilot"], "pairs20 배경 불일치"
print(" ".join(map(str, sorted(z["seeds"]))))
PYEOF
) || { echo "_RESULT_pairs20.json 의 배경이 '100 에서 내려가며 20' 과 다르다 — 발사 중단"; exit 1; }
MAIN_BASES=$(python3 -c 'import json; s = json.load(open("_RESULT_screen20.json")); p = json.load(open("_RESULT_pairs20.json")); assert s["complete"] and s["pass"] == p["pass"]; print(" ".join(s["pass"]))') \
  || { echo "_RESULT_screen20.json pass ≠ _RESULT_pairs20.json pass — 발사 중단"; exit 1; }
# 겹치면 안 되는 배경: 첫 세트(3~22) · 둘째 세트(26 부터 racld 정상 20) · 셋째 세트(위)
USED=$(python3 - <<'PYEOF'
import json
def pick(first, n):
    out = []
    for s in range(first, 101):
        m = json.load(open("main/mat_d8_mu0p01_s%05d.json" % s))
        if m["extinct_at"] < 0 and m["final_counts"] and m["final_counts"][0][0] == "racld":
            out.append(s)
        if len(out) == n: break
    return out
print(" ".join(map(str, pick(1, 20) + pick(26, 20))))
PYEOF
)

if [ -n "${PETRI_SEEDS:-}" ]; then
  [ "$ROOT" != inv24 ] || { echo "PETRI_SEEDS 를 주면 PETRI_ROOT 는 inv24 가 아니어야 한다(파일럿 루트) — 발사 중단"; exit 1; }
  SEEDS=$PETRI_SEEDS
  for s in $SEEDS; do for u in $MAIN_SEEDS $USED; do
    [ "$s" != "$u" ] || { echo "파일럿 배경 $s 가 본 세트 배경과 겹친다 — 발사 중단"; exit 1; }
  done; done
else
  SEEDS=$MAIN_SEEDS
  [ "$(echo $SEEDS | wc -w)" = 20 ] || { echo "배경 수 ≠ 20 — 발사 중단"; exit 1; }
fi
BASES=${PETRI_BASES:-$MAIN_BASES}
if [ "$ROOT" = inv24 ]; then
  [ "$BASES" = "$MAIN_BASES" ] || { echo "본 루트 inv24 에서는 밑 코드를 바꿀 수 없다 — 발사 중단"; exit 1; }
  [ "$PARTS" = IB ] || { echo "본 루트 inv24 는 IB 만 — 발사 중단"; exit 1; }
fi

codes_of() {  # $1 = 밑 코드 → n 삽입 코드 목록(쉼표) — pairs13/14/20 판정기의 loop_ticks · inserts 와 같은 함수
  B="$1" python3 - <<'PYEOF'
import os
def loop_ticks(code):
    if "l" not in code or "c" not in code: return None
    li = code.index("l"); si = code.rfind("s", 0, li); start = si + 1 if si >= 0 else 0
    if start > li: return None
    seg = code[start:li + 1]
    return sum(1 for ch in seg if ch != "x") if "c" in seg else None
def inserts(base, ch):
    out = {}
    for p in range(len(base) + 1):
        k = base[:p] + ch + base[p:]; L = loop_ticks(k)
        if L is not None: out.setdefault(k, L)
    return out
print(",".join(inserts(os.environ["B"], "n")))
PYEOF
}
CODEMAP=$(for b in $BASES; do echo "$b:$(codes_of "$b")"; done)
L=${ROOT}.log
J=${ROOT}_jobs.txt
if [ "$RESUME" = 1 ]; then
  [ "$PARTS" = IB ] && [ -e "$L" ] || { echo "재개: 부분 IB · 기존 $L 이 있어야 한다(파일럿 루트에서도 같은 경로 — 재개 시험용) — 중단"; exit 1; }
elif [ "$DRY" != 1 ]; then
  for p in $PARTS; do case $p in
    IB) for b in $BASES; do [ ! -e "${ROOT}ib_$b" ] || { echo "${ROOT}ib_$b 가 이미 있다 — 덮어쓰지 않는다 · 발사 중단"; exit 1; }; done;;
    D)  for b in $BASES; do [ ! -e "${ROOT}D_$b" ] || { echo "${ROOT}D_$b 가 이미 있다 — 발사 중단"; exit 1; }; done;;
    REG) [ ! -e "${ROOT}reg" ] || { echo "${ROOT}reg 가 이미 있다 — 발사 중단"; exit 1; };;
    *) echo "모르는 부분 $p — 발사 중단"; exit 1;;
  esac; done
  [ "$ROOT" != inv24 ] || [ ! -e "$L" ] || { echo "$L 이 이미 있다 — 본 발사는 한 번만 · 발사 중단"; exit 1; }
fi
: > "$J"
for b in $BASES; do
  codes=$(echo "$CODEMAP" | awk -F: -v b="$b" '$1==b{print $2}')
  [ -n "$codes" ] || { echo "밑 코드 $b 의 코드 목록이 비었다 — 발사 중단"; exit 1; }
  for p in $PARTS; do case $p in
    IB) for s in $SEEDS; do echo "IB $s $b $codes ${ROOT}ib_$b/mu0" >> "$J"; done;;
    D)  for s in $SEEDS; do echo "D $s $b $codes ${ROOT}D_$b/mu0" >> "$J"; done;;
  esac; done
done
case " $PARTS " in *" REG "*)
  echo "IB 3 racld $(codes_of racld) ${ROOT}reg/s3_racld" >> "$J"
  echo "IB 26 rascled $(codes_of rascled) ${ROOT}reg/s26_rascled" >> "$J";;
esac
[ "$(awk 'NF!=5' "$J" | wc -l)" = 0 ] || { echo "작업 줄의 칸 수가 5 가 아니다 — 발사 중단"; exit 1; }
if [ "$DRY" = 1 ]; then
  echo "마른 실행 — 루트 $ROOT · 부분 $PARTS · 배경 $SEEDS · 밑 코드 $BASES · 병렬 $PAR"
  echo "작업 $(wc -l < "$J") (IB $(grep -c '^IB ' "$J") · D $(grep -c '^D ' "$J"))"
  cut -c1-160 "$J"
  exit 0
fi
if [ "$RESUME" = 1 ]; then   # M3 — 출력이 없거나 깨진 작업만 남긴다
  n=1; while [ -e "${ROOT}_resume_$n.txt" ]; do n=$((n + 1)); done; RJ=${ROOT}_resume_$n.txt
  python3 - "$J" > "$RJ" <<'PYEOF' || { echo "재개 목록 만들기 실패 — 중단"; exit 1; }
import json, os, sys
for line in open(sys.argv[1]):
    kind, s, b, codes, out = line.split()
    f = os.path.join(out, "ib_inv_s%05d_mu0.json" % int(s))
    try:
        z = json.load(open(f)); ok = z["seed"] == int(s) and len(z["arms"]) == 1 + len(codes.split(","))
    except Exception:
        ok = False
    if not ok: sys.stdout.write(line)
PYEOF
  k=$(wc -l < "$RJ")
  echo "===== 재실행 $n · 작업 $k · 병렬 $PAR · 다른 세션 node $OTHER · $(date '+%F %T') · uptime -s $(uptime -s) · 목록 $RJ =====" >> "$L"
  [ "$k" -gt 0 ] || { echo "재개: 다시 돌릴 작업 없음" >> "$L"; echo "재개: 다시 돌릴 작업 없음"; exit 0; }
  J=$RJ
else
: > "$L"
echo "해부 24 · 루트 $ROOT · 부분 $PARTS · 배경 $SEEDS · 밑 코드 $BASES" >> "$L"
echo "코드 목록: $CODEMAP" | tr '\n' ' ' >> "$L"; echo >> "$L"
echo "작업 $(wc -l < "$J") · 병렬 $PAR · 다른 세션 node $OTHER · $(date '+%F %T') · uptime -s $(uptime -s) · node $(node --version) · invbud $(sha256sum invbud.js | cut -c1-12) · sim $(sha256sum sim.js | cut -c1-12) · dms $(sha256sum dms.js | cut -c1-12)" >> "$L"
fi
t0=$(date +%s)
xargs -P "$PAR" -L 1 sh -c '
  if [ "$0" = IB ]; then node invbud.js --seed "$1" --mu 0 --codes "$3" --out "$4";
  else node dms.js --cond mat --wt "$2" --seed "$1" --mu-assay 0 --exact 1 --arms list --codes "$3" --out "$4"; fi
  echo "rc $? $0 $1 $2"' < "$J" >> "$L" 2>&1
echo "===== 끝 $(( $(date +%s) - t0 ))초 · $(date '+%F %T') · uptime -s $(uptime -s) =====" >> "$L"
