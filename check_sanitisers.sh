#!/usr/bin/env bash
# AddressSanitizer + UndefinedBehaviorSanitizer build; runs the RK4 self-test and one trajectory (Section 7.6).
set -e
cd "$(dirname "$0")"
g++ -std=c++20 -O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer -pthread code_annotated.cpp -o /tmp/engine_san
tmp=$(mktemp -d)
/tmp/engine_san --selftest input.txt
/tmp/engine_san input.txt "$tmp/traj.txt" "$tmp/err.txt" > /dev/null
echo "sanitisers: PASS (no ASan/UBSan reports on self-test and a trajectory run)"
