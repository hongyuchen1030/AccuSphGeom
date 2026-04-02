#!/usr/bin/env bash
set -euo pipefail

# From project root
ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BUILD_DIR="${ROOT_DIR}/build"

mkdir -p "${BUILD_DIR}"
cd "${BUILD_DIR}"

# Prefer Ninja if available; otherwise fall back to Makefiles.
GENERATOR="Unix Makefiles"
if command -v ninja >/dev/null 2>&1; then
  GENERATOR="Ninja"
fi

cmake -G "${GENERATOR}" \
  -DCMAKE_BUILD_TYPE=Release \
  "${ROOT_DIR}"

cmake --build . -j
