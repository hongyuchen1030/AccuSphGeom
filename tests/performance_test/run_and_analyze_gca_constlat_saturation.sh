#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="/global/homes/h/hyvchen/AccuSphGeom"
INSPECT_SCRIPT="${REPO_ROOT}/tests/performance_test/inspect_perlmutter_compute_unit.sh"
SWEEP_SCRIPT="${REPO_ROOT}/tests/performance_test/run_gca_constlat_saturation_sweep.sh"
ANALYZE_SCRIPT="${REPO_ROOT}/tests/performance_test/analyze_gca_constlat_saturation.py"

REPORT_PATH="${REPO_ROOT}/tests/performance_test/output/gca_constlat_saturation_report.txt"
SWEEP_PLOT="${REPO_ROOT}/tests/performance_test/output/gca_constlat_saturation_sweep.png"
RATIO_PLOT="${REPO_ROOT}/tests/performance_test/output/gca_constlat_saturation_best_ratio.png"

"${INSPECT_SCRIPT}"
SKIP_INSPECT=1 "${SWEEP_SCRIPT}"
python3 "${ANALYZE_SCRIPT}"

echo
echo "Saturation report:"
echo "  ${REPORT_PATH}"
echo
echo "Plots:"
echo "  ${SWEEP_PLOT}"
echo "  ${RATIO_PLOT}"
