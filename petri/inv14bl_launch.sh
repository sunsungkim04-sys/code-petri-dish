#!/usr/bin/env bash
# 해부 14 후속(③ · 09-28) — 헛손질을 요청 글자 종류별로 센다. 해부 14 IB 와 같은 40 코드 · 배경 3~22 · μ 0 · invbud.js --by-letter 1 → 120
#   회귀: 팔 기록 [t, n, m] 이 inv14ib_*/mu0 의 같은 팔과 같아야 한다(플래그는 읽기만) · 칸마다 stall_by 합 = copy_stall
# 마른 실행: PETRI_NBG=1 PETRI_BG_FROM=23 PETRI_NBASES="rascled" PETRI_ROOT=pilot14bl2 ./inv14bl_launch.sh
set -u
export NVM_DIR="$HOME/.nvm"; . "$NVM_DIR/nvm.sh"
cd "${PETRI_HOME:-$HOME/petri}" || exit 1
grep -q "byLetter" invbud.js || { echo "invbud.js 에 --by-letter 가 없다 — 발사 중단"; exit 1; }
NBG=${PETRI_NBG:-20}
NBASES=${PETRI_NBASES:-"racld rascld rsacld racldx acld rascled"}
ROOT=${PETRI_ROOT:-inv14bl}
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
  for s in $BGS; do echo "$s $codes ${ROOT}_$b/mu0" >> "$J"; done
done
echo "배경 $BGS · 작업 $(wc -l < "$J") · $(date '+%F %T') · uptime -s $(uptime -s)" >> "$L"
t0=$(date +%s)
xargs -P 60 -n 3 sh -c 'node invbud.js --seed "$0" --mu 0 --codes "$1" --by-letter 1 --out "$2"' < "$J" >> "$L" 2>&1
echo "===== 끝 $(( $(date +%s) - t0 ))초 · $(date '+%F %T') =====" >> "$L"
