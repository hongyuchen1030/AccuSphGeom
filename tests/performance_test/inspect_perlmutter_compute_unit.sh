#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="/global/homes/h/hyvchen/AccuSphGeom"
OUTPUT_DIR="${REPO_ROOT}/tests/performance_test/output"
OUTPUT_PATH="${OUTPUT_DIR}/perlmutter_compute_unit_info.txt"

mkdir -p "${OUTPUT_DIR}"

cpus_allowed_list="$(grep Cpus_allowed_list /proc/self/status | awk '{print $2}' || true)"

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

visible_logical_cpus=""
if [[ -n "${cpus_allowed_list}" ]]; then
  visible_logical_cpus="$(count_cpus_in_list "${cpus_allowed_list}")"
fi

unit_core=1
unit_small="unresolved"
if [[ -n "${visible_logical_cpus}" ]]; then
  if (( visible_logical_cpus >= 4 )); then
    unit_small=4
  elif (( visible_logical_cpus >= 2 )); then
    unit_small=2
  elif (( visible_logical_cpus >= 1 )); then
    unit_small=1
  fi
fi

unit_requested="${CORES:-16}"
if [[ -n "${visible_logical_cpus}" ]] && (( unit_requested > visible_logical_cpus )); then
  unit_requested="${visible_logical_cpus}"
fi

{
  echo "Perlmutter compute unit inspection"
  echo "================================="
  echo

  echo "[basic]"
  echo "hostname: $(hostname)"
  echo "date: $(date -Is)"
  echo "SLURM_JOB_ID: ${SLURM_JOB_ID:-}"
  echo "SLURM_JOB_NUM_NODES: ${SLURM_JOB_NUM_NODES:-}"
  echo "SLURM_CPUS_ON_NODE: ${SLURM_CPUS_ON_NODE:-}"
  echo "SLURM_CPUS_PER_TASK: ${SLURM_CPUS_PER_TASK:-}"
  echo "SLURM_JOB_CPUS_PER_NODE: ${SLURM_JOB_CPUS_PER_NODE:-}"
  echo "SLURM_CPU_BIND: ${SLURM_CPU_BIND:-}"
  echo "OMP_NUM_THREADS: ${OMP_NUM_THREADS:-}"
  echo "OMP_PLACES: ${OMP_PLACES:-}"
  echo "OMP_PROC_BIND: ${OMP_PROC_BIND:-}"
  echo

  echo "[derived summary]"
  echo "visible_logical_cpus: ${visible_logical_cpus:-unknown}"
  echo "unit_core: ${unit_core}"
  echo "unit_small: ${unit_small}"
  echo "unit_requested: ${unit_requested}"
  echo

  echo "[nproc]"
  nproc
  echo

  echo "[proc status cpuset]"
  grep Cpus_allowed_list /proc/self/status || true
  echo

  echo "[compiler version]"
  if command -v c++ >/dev/null 2>&1; then
    c++ --version || true
  elif command -v g++ >/dev/null 2>&1; then
    g++ --version || true
  else
    echo "No c++ or g++ found"
  fi
  echo

  echo "[cmake version]"
  cmake --version || true
  echo

  echo "[lscpu]"
  lscpu || true
  echo

  echo "[lscpu -e]"
  lscpu -e || true
  echo

  echo "[numactl --hardware]"
  if command -v numactl >/dev/null 2>&1; then
    numactl --hardware || true
  else
    echo "numactl not available"
  fi
  echo

  echo "[hwloc-ls]"
  if command -v hwloc-ls >/dev/null 2>&1; then
    hwloc-ls || true
  else
    echo "hwloc-ls not available"
  fi
  echo

  echo "[lstopo]"
  if command -v lstopo >/dev/null 2>&1; then
    lstopo --of console || true
  else
    echo "lstopo not available"
  fi
  echo
} | tee "${OUTPUT_PATH}"

echo
echo "Wrote compute-unit inspection to ${OUTPUT_PATH}"
