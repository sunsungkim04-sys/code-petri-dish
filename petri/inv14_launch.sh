#!/usr/bin/env bash
# 해부 14 — 짝 코드 전수의 글자당 주사위 횟수 직접 계수. 사전등록 해부14-사전등록-2026-09-27.md §1~§2
#   IB  invbud.js --codes <밑 코드 B 의 n 삽입 전부> --mu 0                         × 밑 코드 6 × mat 배경 3~22 → 120
#   Q   dms.js --wt B --arms list --codes <같은 목록> --kids 1 --exact 1 --mu-assay 0.003  × 같은 → 120
#       (--arms list 는 코드가 --wt 의 한 글자 이웃이어야 한다 → --wt 를 밑 코드로. ref 팔은 밑 코드를 끼우지만 U 계산엔 안 쓴다 · 배경 성장은 --wt 와 무관)
#   계측기는 해부 7(invbud.js v1.0.0) · 해부 12(dms.js --kids) 것 그대로 — 새 패치 없음
#   판정: pairs14_analyze.py
# 마른 실행(파일럿): PETRI_NBG=2 PETRI_BG_FROM=23 PETRI_NBASES="racld rascled" PETRI_ROOT=pilot14 ./inv14_launch.sh
set -u
export NVM_DIR="$HOME/.nvm"; . "$NVM_DIR/nvm.sh"
cd "${PETRI_HOME:-$HOME/petri}" || exit 1
grep -q "VERSION = '0.3.5'" sim.js || { echo "sim.js 가 v0.3.5 가 아니다 — 발사 중단"; exit 1; }
grep -q "INVBUD_VERSION = '1.0.0'" invbud.js || { echo "invbud.js 가 v1.0.0 이 아니다 — 발사 중단"; exit 1; }
grep -q "kidsOn" dms.js || { echo "dms.js 에 --kids 가 없다 — 발사 중단"; exit 1; }
NBG=${PETRI_NBG:-20}
NBASES=${PETRI_NBASES:-"racld rascld rsacld racldx acld rascled"}
PARTS=${PETRI_PARTS:-"IB Q"}
MUQ=${PETRI_MUQ:-0.003}
ROOT=${PETRI_ROOT:-inv14}
L=${ROOT}.log; : > "$L"
J=${ROOT}_jobs.txt; : > "$J"
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
# 밑 코드마다 n 삽입 목록 — pairs13_analyze.py 의 loop_ticks · inserts 와 같은 함수(고리 없는 삽입은 뺀다)
CODEMAP=$(NBASES="$NBASES" python3 - <<'PYEOF'
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
for b in os.environ["NBASES"].split():
    print(b + ":" + ",".join(inserts(b, "n")))
PYEOF
)
for b in $NBASES; do
  codes=$(echo "$CODEMAP" | awk -F: -v b="$b" '$1==b{print $2}')
  [ -n "$codes" ] || { echo "밑 코드 $b 의 코드 목록이 비었다 — 발사 중단"; exit 1; }
  for p in $PARTS; do case $p in
    IB) for s in $BGS; do echo "IB 0 $s $b $codes ${ROOT}ib_$b/mu0" >> "$J"; done;;
    Q)  tag=mu$(echo "$MUQ" | sed 's/\./p/'); for s in $BGS; do echo "Q $MUQ $s $b $codes ${ROOT}kids_$b/$tag" >> "$J"; done;;
  esac; done
done
[ "$(awk 'NF!=6' "$J" | wc -l)" = 0 ] || { echo "작업 줄의 칸 수가 6 이 아니다 — 발사 중단"; exit 1; }
echo "배경 $BGS" >> "$L"
echo "코드 목록: $CODEMAP" | tr '\n' ' ' >> "$L"; echo >> "$L"
echo "작업 $(wc -l < "$J") (IB $(grep -c '^IB ' "$J") · Q $(grep -c '^Q ' "$J")) · μ_Q $MUQ · $(date '+%F %T') · uptime -s $(uptime -s)" >> "$L"
t0=$(date +%s)
xargs -P 60 -L 1 sh -c '
  if [ "$0" = IB ]; then node invbud.js --seed "$2" --mu "$1" --codes "$4" --out "$5";
  else node dms.js --cond mat --wt "$3" --seed "$2" --mu-assay "$1" --exact 1 --arms list --codes "$4" --kids 1 --out "$5"; fi' < "$J" >> "$L" 2>&1
echo "===== 끝 $(( $(date +%s) - t0 ))초 · $(date '+%F %T') =====" >> "$L"
