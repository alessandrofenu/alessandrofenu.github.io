#!/bin/bash
# Builds ../ho.js + ../ho.wasm from homology_optimized.cpp (needs emsdk, see README.md).
set -e
cd "$(dirname "$0")"
GMP=${GMP:-$HOME/emsdk-build/gmp-wasm}
em++ -O3 -std=c++17 -fwasm-exceptions ho_wasm.cpp -Icompat -I"$GMP/include" "$GMP/lib/libgmpxx.a" "$GMP/lib/libgmp.a" \
    -sMODULARIZE=1 -sEXPORT_NAME=createHO -sENVIRONMENT=worker,node \
    -sALLOW_MEMORY_GROWTH=1 -sINITIAL_MEMORY=64MB -sMAXIMUM_MEMORY=4GB -sSTACK_SIZE=16MB \
    -sEXPORTED_FUNCTIONS=_ho_compute -sEXPORTED_RUNTIME_METHODS=ccall,UTF8ToString \
    -o ../ho.js
