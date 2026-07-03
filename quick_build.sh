#!/usr/bin/env bash
set -euo pipefail

# From project root
ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BUILD_DIR="${ROOT_DIR}/build"
BUILD_TESTS="OFF"
ENABLE_EIGEN_ADAPTERS="OFF"

usage() {
  cat <<'EOF'
Usage: ./quick_build.sh [options]

Options:
  --tests                  Configure with -DACCUSPHGEOM_BUILD_TESTS=ON
  --no-tests               Configure with -DACCUSPHGEOM_BUILD_TESTS=OFF
  --eigen-adapters         Configure with -DACCUSPHGEOM_ENABLE_EIGEN_ADAPTERS=ON
  --no-eigen-adapters      Configure with -DACCUSPHGEOM_ENABLE_EIGEN_ADAPTERS=OFF
  -h, --help               Show this help message
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --tests)
      BUILD_TESTS="ON"
      ;;
    --no-tests)
      BUILD_TESTS="OFF"
      ;;
    --eigen-adapters)
      ENABLE_EIGEN_ADAPTERS="ON"
      ;;
    --no-eigen-adapters)
      ENABLE_EIGEN_ADAPTERS="OFF"
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown option: $1" >&2
      usage >&2
      exit 1
      ;;
  esac
  shift
done

mkdir -p "${BUILD_DIR}"
cd "${BUILD_DIR}"

# Prefer Ninja if available; otherwise fall back to Makefiles.
GENERATOR="Unix Makefiles"
if command -v ninja >/dev/null 2>&1; then
  GENERATOR="Ninja"
fi

cmake -G "${GENERATOR}" \
  -DCMAKE_BUILD_TYPE=Release \
  -DACCUSPHGEOM_BUILD_TESTS="${BUILD_TESTS}" \
  -DACCUSPHGEOM_ENABLE_EIGEN_ADAPTERS="${ENABLE_EIGEN_ADAPTERS}" \
  -DACCUSPHGEOM_PREDICATES_USE_FLOAT=OFF \
  "${ROOT_DIR}"

cmake --build . -j
