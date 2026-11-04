#!/usr/bin/env bash
# Run every check in this folder and save its output under out/. Stops at the first failure.
# Needs g++ and uv. Single-threaded BLAS keeps the floating-point output reproducible.
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p out
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
g++ -O2 -o enum_trees enum_trees.cpp
uv run --no-project -q constants.py            | tee out/constants.txt
uv run --no-project -q check_forest.py         | tee out/check_forest.txt
uv run --no-project -q check_tensorization.py  | tee out/check_tensorization.txt
: > out/check_chain.txt
for g in prism_K4 prism_C6 prism_K33 C10_square grid_4x4; do
  uv run --no-project -q check_chain.py --graph "$g" --cuts 3 --search 2 --seed 1 | tee -a out/check_chain.txt
done
uv run --no-project -q independent_check.py   | tee out/independent_check.txt
# vertex-touch argument (2-connected expanders)
uv run --no-project -q check_blocks.py          | tee out/check_blocks.txt
: > out/check_vertex_exact.txt
for g in K33 cube petersen heawood mobius_kantor pappus rr3_14 rr3_16 rr4_12; do
  uv run --no-project -q check_vertex_exact.py --graph "$g" --cuts 2 --search 1 --steps 100 --seed 1 | tee -a out/check_vertex_exact.txt
done
: > out/check_vertex_medium.txt
for g in tutte_coxeter foster tutte12 rr3_200 rr3_600; do
  uv run --no-project -q check_vertex_medium.py --graph "$g" --sep 3 --samples 20 --seed 2 | tee -a out/check_vertex_medium.txt
done
uv run --no-project -q check_deletion.py        | tee out/check_deletion.txt
uv run --no-project -q constants_touch.py       | tee out/constants_touch.txt
# complete graph and matroids
uv run --no-project -q check_complete.py        | tee out/check_complete.txt
uv run --no-project -q check_edge_pairs.py      | tee out/check_edge_pairs.txt
uv run --no-project -q check_matroid.py         | tee out/check_matroid.txt
# circuits with subtraction
uv run --no-project -q check_signed.py          | tee out/check_signed.txt
uv run --no-project -q check_signed_restricted.py | tee out/check_signed_restricted.txt
uv run --no-project -q check_signed_approx.py   | tee out/check_signed_approx.txt
uv run --no-project -q check_signed_packing.py  | tee out/check_signed_packing.txt
echo "ALL CHECKS PASSED"
