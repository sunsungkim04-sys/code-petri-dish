#!/usr/bin/env bash
# 해부 5 본 발사 — A(밀도 사다리) 먼저, 끝나면 B(촘촘한 μ). 코어를 나눠 쓰지 않으려고 차례로 돈다.
cd "$HOME/petri" || exit 1
./dens_launch.sh
./mufine_launch.sh
echo "===== 해부5 전부 끝 $(date '+%F %T') =====" >> dens.log
