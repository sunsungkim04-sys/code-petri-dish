#!/usr/bin/env bash
# 해부 25 — 둘째 밀도에서 짝 · 교차 사다리. 사전등록 해부25-사전등록-2026-10-01.md §1 · §2
#   배경 성장(grow25_launch.sh)이 끝난 뒤에 돈다. 배경 = <GROW> 의 밀도 D 기록 중 '멸종 없음 · 끝 우세 racld' 인 시드를
#   시드 순으로 앞 20 개(밀도 8 의 자격 규칙과 같다 — dms_launch.sh · inv15_launch.sh · inv21_launch.sh).
#   계측기 · 명령은 기존 것 그대로, --density D --main-dir <GROW> 만 붙는다(dms.js 가 같은 배경을 다시 키워 체크섬 대조):
#     P   짝     node dms.js --cond mat --wt B --seed s --mu-assay 0 --exact 1 --arms list --codes <B 의 n 삽입 전부>      (해부 15 E · 20 P)
#               B = racld rascld rsacld racldx acld rascled rcald crald reasccld                                → 9 × 20 = 180
#     IB  계수   node invbud.js --seed s --mu 0 --codes <racld 의 n 삽입 전부>  (해부 15 Eib — 조작 확인 서술)              → 20
#     L   사다리 node dms.js --cond mat --wt racld --seed s --arms list --codes rascld,rsacld,racldx --mu-assay μ --exact 1
#               μ 0 · 0.1 · 0.2 · 0.3 · 0.4 · 0.5 · 0.75 · 1%   (해부 6C · 21 과 같은 μ 점)                         → 8 × 20 = 160
#   판정: pairs25_analyze.py
# 본 발사:  PETRI_ROOT=inv25 setsid nohup ./inv25_launch.sh > inv25.out 2>&1 < /dev/null &
# 파일럿:   PETRI_DENS=<D> PETRI_GROW=pilot25g PETRI_SEEDS="a b" PETRI_ROOT=pilot25 ./inv25_launch.sh
set -u
export NVM_DIR="$HOME/.nvm"; . "$NVM_DIR/nvm.sh"
cd "${PETRI_HOME:-$HOME/petri}" || exit 1
node --version | grep -q '^v20\.' || { echo "node 가 v20 이 아니다($(node --version)) — nvm 확인 · 발사 중단"; exit 1; }
# 계측기는 새 패치 없이 그대로(사전등록 §1) — 해시로 묶는다(dms20 · dms21 동결과 같은 sim.js · dms.js · 해부 15 의 invbud.js)
sha256sum -c --quiet - <<'EOF' || { echo "계측기 해시가 동결과 다르다 — 발사 중단"; exit 1; }
b459bb4a9933a38bbc3e06e3681aeb977dfa0a7cf680dd4d375d744415cf06b2  sim.js
9c28666a30d4ca265d4f6120331166f9a3c96b5cedafedec5cc002711c3199fd  dms.js
28bfc9c11f27dcdd3e911fa34a3c8d066e2118bc93d08b4a732e6cea783e9246  invbud.js
EOF
grep -q "VERSION = '0.3.5'" sim.js || { echo "sim.js 가 v0.3.5 가 아니다 — 발사 중단"; exit 1; }
grep -q "INVBUD_VERSION = '1.0.0'" invbud.js || { echo "invbud.js 가 v1.0.0 이 아니다 — 발사 중단"; exit 1; }
UP=$(uptime -s)
[ "$UP" = "${PETRI_UPTIME:-2026-09-26 23:34:08}" ] || { echo "uptime -s $UP ≠ 2026-09-26 23:34:08 — 무인 재부팅 가능성 · 발사 중단(보고)"; exit 1; }

# 🔒 사전등록 §2 의 밀도 규칙으로 파일럿 뒤 정한 밀도(Δ 를 보기 전 · 자격 비율만 보고). 본 발사는 이 값만 받는다.
DENS_FIXED=16
ROOT=${PETRI_ROOT:-inv25}
GROW=${PETRI_GROW:-main25}
DENS=${PETRI_DENS:-$DENS_FIXED}
PAR=${PETRI_PAR:-24}
STAGES=${PETRI_STAGES:-"P IB L"}
[ "$PAR" -le 24 ] || { echo "병렬 $PAR > 코어 예산 24 — 발사 중단"; exit 1; }
case $DENS in 12|16) ;; *) echo "밀도 '$DENS' — 12 또는 16 이어야 한다(DENS_FIXED 가 아직 안 정해졌나?) · 발사 중단"; exit 1;; esac
if [ "$ROOT" = inv25 ]; then   # 본 발사: 손잡이를 닫는다
  [ "$DENS" = "$DENS_FIXED" ] || { echo "본 발사 밀도 $DENS ≠ 고정 $DENS_FIXED — 발사 중단"; exit 1; }
  [ "$GROW" = main25 ] && [ -z "${PETRI_SEEDS:-}" ] && [ "$STAGES" = "P IB L" ] || { echo "본 발사에서 GROW · SEEDS · STAGES 를 바꿀 수 없다 — 발사 중단"; exit 1; }
fi
MUS="0 0.001 0.002 0.003 0.004 0.005 0.0075 0.01"
PBASES="racld rascld rsacld racldx acld rascled rcald crald reasccld"

# 배경 목록 — 자격 규칙(밀도 8 과 같다) · 성장 기록이 완결이어야 한다
SEEDS=$(GROW="$GROW" DENS="$DENS" PSEEDS="${PETRI_SEEDS:-}" python3 - <<'PYEOF'
import json, os, sys
g, d = os.environ["GROW"], os.environ["DENS"]
cands = [int(x) for x in os.environ["PSEEDS"].split()] if os.environ["PSEEDS"] else list(range(2501, 2551))
picked = []
for s in cands:
    f = os.path.join(g, "mat_d%s_mu0p01_s%05d.json" % (d, s))
    if not os.path.exists(f):
        print("성장 기록 없음 %s" % f, file=sys.stderr); sys.exit(1)
    z = json.load(open(f))
    # 검토 개정 R9: 판정기 G0(expected_seeds)와 같은 검사 — 후보 전부의 재료 보존 위반이 0 이어야 한다(65분 계산 뒤가 아니라 지금 멈춘다)
    if z["conservation_violations"]:
        print("성장 재료 보존 위반 %s (%s) — 판정기 G0 가 멈출 자료" % (f, z["conservation_violations"]), file=sys.stderr); sys.exit(1)
    if z["ticks_planned"] != 50000 or abs(z["opts"]["mu"] - 0.01) > 1e-12 or float(z["opts"]["density"]) != float(d) or z["cond"] != "mat":
        print("성장 설정이 다르다 %s" % f, file=sys.stderr); sys.exit(1)
    if z["extinct_at"] < 0 and z["final_counts"] and z["final_counts"][0][0] == "racld":
        picked.append(s)
need = len(cands) if os.environ["PSEEDS"] else 20
if len(picked) < need:
    print("자격 배경 %d < %d — 판정 불가(사전등록 §2) · 발사 중단" % (len(picked), need), file=sys.stderr); sys.exit(1)
print(" ".join(map(str, picked[:need])))
PYEOF
) || { echo "배경 목록을 못 정했다 — 발사 중단"; exit 1; }

CODEMAP=$(BASES="$PBASES" python3 - <<'PYEOF'
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
for b in os.environ["BASES"].split():
    print(b + ":" + ",".join(inserts(b, "n")))
PYEOF
)
mt() { echo "mu$(echo "$1" | sed 's/\./p/')"; }
for st in $STAGES; do
  case $st in
    P)  for b in $PBASES; do
          codes=$(echo "$CODEMAP" | awk -F: -v b="$b" '$1==b{print $2}')
          [ -n "$codes" ] || { echo "밑 코드 $b 의 코드 목록이 비었다 — 발사 중단"; exit 1; }
          [ -e "${ROOT}P_$b" ] && { echo "${ROOT}P_$b 가 이미 있다(섞지 않는다) — 발사 중단"; exit 1; }
        done;;
    IB) [ -e "${ROOT}ib_racld" ] && { echo "${ROOT}ib_racld 가 이미 있다 — 발사 중단"; exit 1; };;
    L)  [ -e "${ROOT}L" ] && { echo "${ROOT}L 이 이미 있다 — 발사 중단"; exit 1; };;
    *)  echo "단계 $st — P · IB · L 만"; exit 1;;
  esac
done
# 코어 예산(검토 개정 R7): 다른 node(해부 24 등) + 이 발사의 병렬 ≤ 36. 넘으면 1분마다 다시 보고, PETRI_BUDGET_WAIT 초(기본 3600) 뒤에도 넘으면 멈춘다.
#   성장(grow25)이 끝난 뒤 이어서 돌므로 그쪽 node 는 이미 없다. 판정과 무관한 손잡이(기다리는 시간)만 있다 — 예산 36 은 고정.
BW=${PETRI_BUDGET_WAIT:-3600}; waited=0
while :; do
  OTHER=$(pgrep -c -x node || true); OTHER=${OTHER:-0}
  [ $(( OTHER + PAR )) -le 36 ] && break
  [ "$waited" -ge "$BW" ] && { echo "다른 node $OTHER + 병렬 $PAR > 36 — ${waited}초 기다려도 넘는다 · 발사 중단"; exit 1; }
  echo "다른 node $OTHER + 병렬 $PAR > 36 — 60초 뒤 다시 본다($(date '+%F %T'))"; sleep 60; waited=$(( waited + 60 ))
done
L=${ROOT}.log; : > "$L"
J=${ROOT}_jobs.txt; : > "$J"
for st in $STAGES; do
  case $st in
    P)  for b in $PBASES; do codes=$(echo "$CODEMAP" | awk -F: -v b="$b" '$1==b{print $2}')
          for s in $SEEDS; do echo "D 0 $s $b $codes ${ROOT}P_$b/mu0" >> "$J"; done; done;;
    IB) codes=$(echo "$CODEMAP" | awk -F: '$1=="racld"{print $2}')
        for s in $SEEDS; do echo "IB 0 $s racld $codes ${ROOT}ib_racld/mu0" >> "$J"; done;;
    L)  for mu in $MUS; do for s in $SEEDS; do echo "D $mu $s racld rascld,rsacld,racldx ${ROOT}L/$(mt "$mu")" >> "$J"; done; done;;
  esac
done
[ "$(awk 'NF!=6' "$J" | wc -l)" = 0 ] || { echo "작업 줄의 칸 수가 6 이 아니다 — 발사 중단"; exit 1; }
echo "해부 25 · 단계 $STAGES · 밀도 $DENS · 성장 $GROW · 배경 $SEEDS" >> "$L"
echo "코드 목록: $CODEMAP" | tr '\n' ' ' >> "$L"; echo >> "$L"
echo "작업 $(wc -l < "$J") (P $(grep -c '^D 0 .* '"${ROOT}"'P_' "$J") · IB $(grep -c '^IB ' "$J") · L $(grep -c " ${ROOT}L/" "$J")) · 병렬 $PAR · 다른 node $OTHER(합 $(( OTHER + PAR )) ≤ 36 · 기다림 ${waited}초) · $(date '+%F %T') · uptime -s $UP · node $(node --version)" >> "$L"
t0=$(date +%s)
DD="$DENS" GG="$GROW" xargs -P "$PAR" -L 1 sh -c '
  if [ "$0" = IB ]; then node invbud.js --density "$DD" --main-dir "$GG" --seed "$2" --mu "$1" --codes "$4" --out "$5" || echo "ERROR 작업 $*";
  else node dms.js --cond mat --density "$DD" --main-dir "$GG" --wt "$3" --seed "$2" --mu-assay "$1" --exact 1 --arms list --codes "$4" --out "$5" || echo "ERROR 작업 $*"; fi' < "$J" >> "$L" 2>&1
echo "===== 끝 $(( $(date +%s) - t0 ))초 · $(date '+%F %T') · uptime -s $(uptime -s) =====" >> "$L"
