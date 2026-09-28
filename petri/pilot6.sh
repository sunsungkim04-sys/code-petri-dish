#!/usr/bin/env bash
# 해부 6 설계용 탐색 — 판정 아님. 시드 591~593 (판정 시드 601~ 과 겹치지 않는다)
set -u
export NVM_DIR="$HOME/.nvm"; . "$NVM_DIR/nvm.sh"
cd "$HOME/petri" || exit 1
J=pilot6_jobs.txt; : > "$J"
for code in racld rascld; do
  for mu in 0.001 0.005 0.01; do for s in 591 592 593; do echo "node lineage.js --code $code --seed $s --mu $mu --out pilot6/lin" >> "$J"; done; done
  for d in 4 8 14 32; do for s in 591 592 593; do echo "node budget.js --cond mat --density $d --code $code --seed $s --out pilot6/bud" >> "$J"; done; done
  for s in 591 592 593; do echo "node budget.js --cond space --code $code --seed $s --out pilot6/bud" >> "$J"; done
done
for code in rsacld; do for d in 8; do for s in 591 592 593; do echo "node budget.js --cond mat --density $d --code $code --seed $s --out pilot6/bud" >> "$J"; done; done
  for s in 591 592 593; do echo "node budget.js --cond space --code $code --seed $s --out pilot6/bud" >> "$J"; done; done
tr '\n' '\0' < "$J" | xargs -0 -P 60 -I{} sh -c '{}' > pilot6.log 2>&1
echo "===== 끝 =====" >> pilot6.log
