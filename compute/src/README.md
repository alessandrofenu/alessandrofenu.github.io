# Compute! — sources

Everything served by `compute.html` is built from this folder.

| file | role |
|---|---|
| `homology_optimized.cpp` | the engine for the torus and punctured torus (copy of `idrissi-rocailucio/recontre/homology_optimized.cpp`) |
| `ho_wasm.cpp` | WebAssembly entry point `ho_compute(k, pu)`: same pipeline and same JSON as `./homology_optimized k [pu] --json`, one thread |
| `build_wasm.sh` | builds `../ho.js` + `../ho.wasm` |
| `compat/bits/stdc++.h` | stand-in for GCC's `<bits/stdc++.h>` in the wasm build |
| `build_data.py` | builds `../data.json` (the precomputed tables) |
| `../py/` | `colcx.py` and its dependencies, run by Pyodide for the plane and annulus |
| `../worker.js` | the background thread that runs either engine |

## Rebuilding the WebAssembly engine

1. Emscripten: `git clone https://github.com/emscripten-core/emsdk ~/emsdk && cd ~/emsdk && ./emsdk install latest && ./emsdk activate latest && source ./emsdk_env.sh`
2. GMP for wasm (once), from the gmp-6.3.0 sources:
   `emconfigure ./configure --host=none --disable-assembly --enable-cxx --disable-shared --prefix=$HOME/emsdk-build/gmp-wasm CC_FOR_BUILD=gcc CFLAGS=-O3 CXXFLAGS="-O3 -fwasm-exceptions" && emmake make -j8 && emmake make install`
   (`CC_FOR_BUILD=gcc` is needed: GMP generates its tables with build-time programs, and if emcc compiles them the build fails in `mpn/comb_tables.c`.)
3. `./build_wasm.sh` (links with `em++`; `compat/bits/stdc++.h` replaces the GCC header that Emscripten's libc++ lacks, so `homology_optimized.cpp` stays an unmodified copy)

The wasm engine was checked against the native binary: identical output for k <= 22, torus and punctured torus.

## Rebuilding the tables

`python3 build_data.py plane.jsonl annulus.jsonl ../data.json`, where each line of the `.jsonl` files is
`{"k": k, "bm": {n: [rank, {p: [valuations]}]}}` as printed by a loop over `colcx.bm_homology(k, surface, primes <= k)`.
The torus tables are read from `uconf-torus/tables/homology-torus-k35.json`.

## Estimates shown on the page

`SPACES[...].bench` in `compute.html` holds single-core native measurements
(`HO_THREADS=1 ./homology_optimized --benchmark 10 33 --ref 0 [pu]`, and CPython timings of `colcx.py`);
`OVERHEAD` holds the measured browser/native factors.
