#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
EXADIS_ROOT="${EXADIS_ROOT:-$ROOT_DIR/core/exadis}"
BUILD_DIR="${BUILD_DIR:-$EXADIS_ROOT/build-audit}"
PYTHON_BIN="${PYTHON_BIN:-python3}"
BUILD_JOBS="${BUILD_JOBS:-4}"
EVENT_AUDIT="${EVENT_AUDIT:-ON}"

if [[ "$EVENT_AUDIT" != "ON" && "$EVENT_AUDIT" != "OFF" ]]; then
  echo "EVENT_AUDIT must be ON or OFF" >&2
  exit 2
fi

git -C "$EXADIS_ROOT" submodule update --init kokkos python/pybind11

cmake_args=(
  -S "$EXADIS_ROOT"
  -B "$BUILD_DIR"
  -DEXADIS_ENABLE_EVENT_AUDIT="$EVENT_AUDIT"
  -DEXADIS_PYTHON_BINDING=ON
  -DEXADIS_BUILD_EXAMPLES=OFF
  -DEXADIS_FFT=ON
  -DKokkos_ENABLE_OPENMP=ON
  -DKokkos_ENABLE_SERIAL=ON
  -DPYTHON_EXECUTABLE="$PYTHON_BIN"
  -DPYEXADIS_OUTPUT_DIR="$BUILD_DIR/python"
)

if [[ -n "${CMAKE_CXX_COMPILER:-}" ]]; then
  cmake_args+=("-DCMAKE_CXX_COMPILER=$CMAKE_CXX_COMPILER")
fi
if [[ -n "${FFTW_INC_DIR:-}" ]]; then
  cmake_args+=("-DFFTW_INC_DIR=$FFTW_INC_DIR")
fi
if [[ -n "${FFTW_LIB_DIR:-}" ]]; then
  cmake_args+=("-DFFTW_LIB_DIR=$FFTW_LIB_DIR")
fi

cmake "${cmake_args[@]}"
cmake --build "$BUILD_DIR" --target pyexadis -j"$BUILD_JOBS"

echo "built pyexadis with EXADIS_ENABLE_EVENT_AUDIT=$EVENT_AUDIT in $BUILD_DIR/python"
