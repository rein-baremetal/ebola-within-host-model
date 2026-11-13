#!/usr/bin/env bash
# Screening output must be byte-identical for any thread count (Table 8, Table 10, Section 7.6).
# Builds the engine, runs the reference screen on 1 and 4 threads, compares SHA-256 of the sample CSV.
set -e
cd "$(dirname "$0")"
g++ -std=c++20 -O2 -pthread code_annotated.cpp -o /tmp/engine_det
tmp=$(mktemp -d)
for n in 1 4; do
  sed "s|^threads *=.*|threads = $n|; s|^prefix *=.*|prefix = $tmp/t$n|" screen.txt > "$tmp/screen_$n.txt"
  /tmp/engine_det --screen "$tmp/screen_$n.txt" input.txt "$tmp/err_$n.txt" > /dev/null
done
a=$(sha256sum "$tmp/t1_samples.csv" | cut -d' ' -f1); b=$(sha256sum "$tmp/t4_samples.csv" | cut -d' ' -f1)
echo "1 thread : $a"; echo "4 threads: $b"
[ "$a" = "$b" ] && echo "determinism: PASS (identical SHA-256)" || { echo "determinism: FAIL"; exit 1; }
