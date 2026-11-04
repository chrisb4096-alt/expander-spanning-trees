#!/usr/bin/env bash
# Run every check in this folder and save its output under out/. Stops at the first failure.
# Needs g++ and uv. Single-threaded BLAS keeps the floating-point output reproducible.
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p out
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
g++ -O2 -o enum_trees enum_trees.cpp
uv run --no-project -q constants.py            | tee out/constants.txt
uv run --no-project -q check_forest.py         | tee out/check_forest.txt
uv run --no-project -q check_tensorization.py  | tee out/check_tensorization.txt
: > out/check_chain.txt
for g in prism_K4 prism_C6 prism_K33 C10_square grid_4x4; do
  uv run --no-project -q check_chain.py --graph "$g" --cuts 3 --search 2 --seed 1 | tee -a out/check_chain.txt
done
uv run --no-project -q independent_check.py   | tee out/independent_check.txt
echo "ALL CHECKS PASSED"
