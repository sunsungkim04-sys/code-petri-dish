#!/usr/bin/env bash
# 해부 8 본 발사 — B(나이별 인구학) · A(재료 먼저 되감기) · C(베끼기 결과) 차례로
cd "$HOME/petri" || exit 1
./lh8_launch.sh
./rw8_launch.sh
./spec8_launch.sh
echo "===== 해부8 끝 $(date "+%F %T") =====" > run8.done
