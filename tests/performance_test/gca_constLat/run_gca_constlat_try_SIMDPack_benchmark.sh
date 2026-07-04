#!/usr/bin/env bash
set -euo pipefail


#salloc -N 1 -q interactive -t 04:00:00 -C cpu -A m2637 --ntasks=1 --cpus-per-task=16
REPO_ROOT="/global/homes/h/hyvchen/AccuSphGeom"
SCRIPT_DIR="${REPO_ROOT}/tests/performance_test/gca_constLat"
BUILD_DIR="${REPO_ROOT}/build_benchmark"
TARGET="benchmark_gca_constlat_try_SIMDPack"

# Conservative larger try-API benchmark size.
DATA_SIZE="${DATA_SIZE:-100000000}"
NUM_TESTS="${NUM_TESTS:-50}"
NUM_REPEATS="${NUM_REPEATS:-7}"

# OpenMP thread count for the benchmark.
# This is intentionally called THREADS, not CORES.
THREADS="${THREADS:-16}"

OUTPUT_DIR="${SCRIPT_DIR}/output"
CSV_PATH="${OUTPUT_DIR}/gca_constlat_SIMDPack_timing_summary.csv"
REPEATS_CSV_PATH="${OUTPUT_DIR}/gca_constlat_SIMDPack_timing_repeats.csv"
PLOT_PATH="${OUTPUT_DIR}/gca_constlat_SIMDPack_timing_try.png"

cd "${REPO_ROOT}"

echo "==> Cleaning build directory"
rm -rf "${BUILD_DIR}"

echo "==> Configuring CMake Release build"
cmake -S . -B "${BUILD_DIR}" \
  -DCMAKE_BUILD_TYPE=Release \
  -DACCUSPHGEOM_BUILD_TESTS=ON \
  -DACCUSPHGEOM_ENABLE_EIGEN_ADAPTERS=ON \
  -DACCUSPHGEOM_PREDICATES_USE_FLOAT=OFF

echo "==> Building ${TARGET}"
cmake --build "${BUILD_DIR}" --target "${TARGET}" -j 8

echo "==> Running benchmark"
echo "    DATA_SIZE=${DATA_SIZE}"
echo "    NUM_TESTS=${NUM_TESTS}"
echo "    NUM_REPEATS=${NUM_REPEATS}"
echo "    THREADS=${THREADS}"

export OMP_NUM_THREADS="${THREADS}"
export OMP_PROC_BIND=close
export OMP_PLACES=cores

srun -n 1 -c "${THREADS}" --cpu-bind=cores \
  "${BUILD_DIR}/${TARGET}" "${DATA_SIZE}" "${NUM_TESTS}" "${NUM_REPEATS}"

echo "==> Plotting benchmark results"
module load python
python3 "${SCRIPT_DIR}/plot_gca_constlat_try_SIMDPack.py" \
  "GCA Constant-Latitude Try-API Performance" \
  --csv "${CSV_PATH}" \
  --output "${PLOT_PATH}"

echo "==> Done"
echo "    Summary CSV: ${CSV_PATH}"
echo "    Repeats CSV: ${REPEATS_CSV_PATH}"
echo "    Plot: ${PLOT_PATH}"
