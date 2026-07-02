#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="/global/homes/h/hyvchen/AccuSphGeom"
BUILD_DIR="${REPO_ROOT}/build_benchmark"
TARGET="${BUILD_DIR}/benchmark_gca_constlat_SIMDPack"
OUTPUT_DIR="${REPO_ROOT}/tests/performance_test/output"
SATURATION_DIR="${OUTPUT_DIR}/saturation_runs"
SUMMARY_CSV="${OUTPUT_DIR}/gca_constlat_SIMDPack_timing_summary.csv"
REPEATS_CSV="${OUTPUT_DIR}/gca_constlat_SIMDPack_timing_repeats.csv"
COMBINED_CSV="${OUTPUT_DIR}/gca_constlat_saturation_sweep.csv"
INSPECT_SCRIPT="${REPO_ROOT}/tests/performance_test/inspect_perlmutter_compute_unit.sh"

CORES="${CORES:-16}"
NUM_REPEATS="${NUM_REPEATS:-7}"
SRUN_JOBID="${SRUN_JOBID:-}"

DATA_SIZES=(
  100000
  200000
  500000
  1000000
  2000000
  5000000
  10000000
  20000000
  50000000
)

num_tests_for_size() {
  local data_size="$1"
  if (( data_size <= 500000 )); then
    echo 100
  elif (( data_size <= 5000000 )); then
    echo 50
  else
    echo 20
  fi
}

count_cpus_in_list() {
  local cpu_list="$1"
  local total=0
  local chunk start end

  IFS=',' read -r -a parts <<< "${cpu_list}"
  for chunk in "${parts[@]}"; do
    if [[ -z "${chunk}" ]]; then
      continue
    fi
    if [[ "${chunk}" == *-* ]]; then
      start="${chunk%-*}"
      end="${chunk#*-}"
      total=$((total + end - start + 1))
    else
      total=$((total + 1))
    fi
  done
  echo "${total}"
}

append_summary_to_combined() {
  local data_size="$1"
  local summary_path="$2"
  local combined_path="$3"

  awk -F',' -v OFS=',' -v data_size="${data_size}" '
    NR == 1 {
      for (i = 1; i <= NF; ++i) {
        gsub(/^[[:space:]]+|[[:space:]]+$/, "", $i)
        if ($i == "method") method_col = i
        else if (($i == "threadsNum") || ($i == "threads_num")) threads_col = i
        else if (($i == "vec_width") || ($i == "vecWidth")) width_col = i
        else if ($i == "min_time") min_col = i
        else if ($i == "median_time") median_col = i
        else if ($i == "mean_time") mean_col = i
      }

      if (!(method_col && threads_col && width_col && min_col && median_col && mean_col)) {
        print "Unexpected summary CSV header in " FILENAME > "/dev/stderr"
        exit 2
      }
      next
    }

    NF > 0 {
      print data_size, $method_col, $threads_col, $width_col, $min_col, $median_col, $mean_col
    }
  ' "${summary_path}" >> "${combined_path}"
}

mkdir -p "${OUTPUT_DIR}" "${SATURATION_DIR}"

cpus_allowed_list="$(grep Cpus_allowed_list /proc/self/status | awk '{print $2}' || true)"
visible_cpus=""
if [[ -n "${cpus_allowed_list}" ]]; then
  visible_cpus="$(count_cpus_in_list "${cpus_allowed_list}")"
fi

srun_cores="${CORES}"
if [[ -n "${visible_cpus}" ]] && (( srun_cores > visible_cpus )); then
  srun_cores="${visible_cpus}"
fi

declare -a srun_prefix
if [[ -n "${SRUN_JOBID}" ]]; then
  srun_prefix=(srun --jobid="${SRUN_JOBID}")
else
  srun_prefix=(srun)
fi

cd "${REPO_ROOT}"

echo "==> Cleaning build directory"
rm -rf "${BUILD_DIR}"

echo "==> Configuring CMake Release build"
cmake -S . -B "${BUILD_DIR}" \
  -DCMAKE_BUILD_TYPE=Release \
  -DACCUSPHGEOM_BUILD_TESTS=ON \
  -DACCUSPHGEOM_PREDICATES_USE_FLOAT=OFF

echo "==> Building benchmark_gca_constlat_SIMDPack"
cmake --build "${BUILD_DIR}" --target benchmark_gca_constlat_SIMDPack -j 8

if [[ "${SKIP_INSPECT:-0}" != "1" ]]; then
  echo "==> Inspecting compute unit"
  "${INSPECT_SCRIPT}"
fi

export OMP_PROC_BIND=close
export OMP_PLACES=cores

echo "data_size,method,threadsNum,vec_width,min_time,median_time,mean_time" > "${COMBINED_CSV}"

for data_size in "${DATA_SIZES[@]}"; do
  if (( data_size % 8 != 0 )); then
    echo "data_size ${data_size} is not divisible by 8" >&2
    exit 1
  fi

  num_tests="$(num_tests_for_size "${data_size}")"

  echo "==> Running saturation point"
  echo "    DATA_SIZE=${data_size}"
  echo "    NUM_TESTS=${num_tests}"
  echo "    NUM_REPEATS=${NUM_REPEATS}"
  echo "    CORES=${CORES}"
  echo "    SRUN_CORES=${srun_cores}"
  if [[ -n "${SRUN_JOBID}" ]]; then
    echo "    SRUN_JOBID=${SRUN_JOBID}"
  fi

  "${srun_prefix[@]}" -n 1 -c "${srun_cores}" --cpu-bind=cores \
    "${TARGET}" "${data_size}" "${num_tests}" "${NUM_REPEATS}"

  cp "${SUMMARY_CSV}" "${SATURATION_DIR}/gca_constlat_summary_data_${data_size}.csv"
  cp "${REPEATS_CSV}" "${SATURATION_DIR}/gca_constlat_repeats_data_${data_size}.csv"

  append_summary_to_combined "${data_size}" "${SUMMARY_CSV}" "${COMBINED_CSV}"
done

echo
echo "Combined saturation CSV: ${COMBINED_CSV}"
echo "Per-data-size summaries: ${SATURATION_DIR}"
