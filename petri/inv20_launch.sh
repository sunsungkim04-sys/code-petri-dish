#!/usr/bin/env bash
# 해부 20 — 가족 밖 밑 코드의 짝지은 고리 효과. 사전등록 해부20-사전등록-2026-10-01.md §2 · §3
#   두 단계, 순차. 둘 다 같은 계측기(해부 15 E 와 같은 명령):
#     node dms.js --cond mat --wt B --seed s --mu-assay 0 --exact 1 --arms list --codes <B 의 n 삽입 전부>
#   S (선별)  밑 코드 = 후보(_RESULT_cand20.json 의 picked) + 기존 여섯(규칙의 양성 · 음성 대조) · 선별 배경 10 → 9 × 10 = 90 작업
#             판정: screen20_analyze.py → _RESULT_screen20.json (통과 밑 코드 목록 · Δ 는 읽지 않는다)
#   P (짝)    밑 코드 = _RESULT_screen20.json 의 pass(규칙이 고른 것만) · 짝 배경 20 → 통과 수 × 20 작업
#             판정: pairs20_analyze.py
#   배경 = main/mat_d8_mu0p01 중 '멸종 안 함 · 끝 우세 racld' 인 시드를 **100 에서 내려가며** 센다(53 부터는 기존 실험이 안 썼다):
#             짝 배경 = 처음 20 (100 … 76) · 선별 배경 = 다음 10 (74 … 64). 53~63 은 남겨 둔다.
# 실행:     PETRI_STAGE=S PETRI_ROOT=inv20 setsid nohup ./inv20_launch.sh > inv20S.out 2>&1 &
#           (선별 판정 뒤) PETRI_STAGE=P PETRI_ROOT=inv20 setsid nohup ./inv20_launch.sh > inv20P.out 2>&1 &
# 마른 실행: PETRI_STAGE=S PETRI_SEEDS="24 25" PETRI_ROOT=pilot20 ./inv20_launch.sh
#           PETRI_STAGE=P PETRI_SEEDS="24 25" PETRI_ROOT=pilot20 PETRI_SCREEN=_RESULT_pilot20_screen.json ./inv20_launch.sh
set -u
export NVM_DIR="$HOME/.nvm"; . "$NVM_DIR/nvm.sh"
cd "${PETRI_HOME:-$HOME/petri}" || exit 1
grep -q "VERSION = '0.3.5'" sim.js || { echo "sim.js 가 v0.3.5 가 아니다 — 발사 중단"; exit 1; }
node --version | grep -q '^v20' || { echo "node 가 v20 이 아니다($(node --version)) — nvm 확인 · 발사 중단"; exit 1; }
STAGE=${PETRI_STAGE:?PETRI_STAGE=S 또는 P}
ROOT=${PETRI_ROOT:-inv20}
PAR=${PETRI_PAR:-36}
[ "$PAR" -le 36 ] || { echo "병렬 $PAR > 코어 예산 36 — 발사 중단"; exit 1; }
OLD6="racld rascld rsacld racldx acld rascled"
L=${ROOT}${STAGE}.log; : > "$L"
J=${ROOT}${STAGE}_jobs.txt; : > "$J"

pick_desc() {  # $1 = 건너뛸 개수 · $2 = 개수 — 100 에서 내려가며 racld 정상 배경
  SKIP=$1 NBG=$2 python3 - <<'PYEOF'
import json, os
skip, nbg = int(os.environ["SKIP"]), int(os.environ["NBG"]); picked = []
for s in range(100, 52, -1):
    z = json.load(open(f"main/mat_d8_mu0p01_s{s:05d}.json"))
    if z["extinct_at"] < 0 and z["final_counts"] and z["final_counts"][0][0] == "racld":
        picked.append(s)
print(" ".join(map(str, picked[skip:skip + nbg])))
PYEOF
}

case $STAGE in
  S) BASES=${PETRI_BASES:-"$(python3 -c 'import json; print(" ".join(json.load(open("_RESULT_cand20.json"))["picked"]))') $OLD6"}
     SEEDS=${PETRI_SEEDS:-$(pick_desc 20 10)}
     NEED=10;;
  P) SCREEN=${PETRI_SCREEN:-_RESULT_screen20.json}
     [ -f "$SCREEN" ] || { echo "$SCREEN 없음 — 선별 판정(screen20_analyze.py)을 먼저 · 발사 중단"; exit 1; }
     # 검토 개정 R1(10-01): 미완결 선별과 '통과 없음' 을 가른다(전에는 assert 실패도 '통과 없음 · exit 0' 으로 보였다)
     SCREEN="$SCREEN" python3 -c 'import json, os, sys; z = json.load(open(os.environ["SCREEN"])); sys.exit(0 if z.get("complete") else 1)' \
       || { echo "$SCREEN 이 완결되지 않았다(complete ≠ true · 점검 실패) — 발사 중단"; exit 1; }
     # 검토 개정 R2(10-01): 통과한 후보가 없으면 멈춘다 — 대조만으로 돌릴지는 사람이 정한다(사전등록 §5 '선별 전멸')
     NCP=$(SCREEN="$SCREEN" python3 -c 'import json, os; print(len(json.load(open(os.environ["SCREEN"]))["candidates_pass"]))')
     [ "$NCP" -ge 1 ] || [ "${PETRI_ALLOW_CTRL_ONLY:-0}" = 1 ] || { echo "통과한 후보 0 — 선별 전멸(§5) · 대조만 돌리려면 사람이 PETRI_ALLOW_CTRL_ONLY=1 · 발사 중단"; exit 1; }
     BASES=${PETRI_BASES:-$(SCREEN="$SCREEN" python3 -c 'import json, os; z = json.load(open(os.environ["SCREEN"])); print(" ".join(z["pass"]))')}
     [ -n "$BASES" ] || { echo "선별 통과 밑 코드가 없다 — 짝 단계 없음 · 발사 중단"; exit 1; }
     SEEDS=${PETRI_SEEDS:-$(pick_desc 0 20)}
     NEED=20;;
  *) echo "PETRI_STAGE 는 S 또는 P"; exit 1;;
esac
[ -n "${PETRI_SEEDS:-}" ] || [ "$(echo $SEEDS | wc -w)" = "$NEED" ] || { echo "배경 수 $(echo $SEEDS | wc -w) ≠ $NEED — 발사 중단"; exit 1; }

CODEMAP=$(BASES="$BASES" python3 - <<'PYEOF'
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
for b in $BASES; do
  codes=$(echo "$CODEMAP" | awk -F: -v b="$b" '$1==b{print $2}')
  [ -n "$codes" ] || { echo "밑 코드 $b 의 코드 목록이 비었다 — 발사 중단"; exit 1; }
  for s in $SEEDS; do echo "$s $b $codes ${ROOT}${STAGE}_$b/mu0" >> "$J"; done
done
[ "$(awk 'NF!=4' "$J" | wc -l)" = 0 ] || { echo "작업 줄의 칸 수가 4 가 아니다 — 발사 중단"; exit 1; }
echo "단계 $STAGE · 배경 $SEEDS · 밑 코드 $BASES" >> "$L"
echo "코드 목록: $CODEMAP" | tr '\n' ' ' >> "$L"; echo >> "$L"
echo "작업 $(wc -l < "$J") · 병렬 $PAR · $(date '+%F %T') · uptime -s $(uptime -s) · node $(node --version)" >> "$L"
t0=$(date +%s)
xargs -P "$PAR" -L 1 sh -c 'node dms.js --cond mat --wt "$1" --seed "$0" --mu-assay 0 --exact 1 --arms list --codes "$2" --out "$3"' < "$J" >> "$L" 2>&1
echo "===== 끝 $(( $(date +%s) - t0 ))초 · $(date '+%F %T') =====" >> "$L"
