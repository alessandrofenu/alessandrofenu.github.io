// Background thread for compute.html: runs one computation and posts the Borel-Moore homology back.
//   {engine: 'ho', k, pu}       -> homology_optimized.cpp compiled to WebAssembly (ho.js / ho.wasm)
//   {engine: 'py', k, surface}  -> colcx.py (plane, annulus) under Pyodide
const PYODIDE = 'https://cdn.jsdelivr.net/pyodide/v0.29.5/full/';
const PY_FILES = ['colcx.py', 'napolitano.py', 'exact.py', 'exact_snf.py', 'reduction.py'];

function fail(msg) { postMessage({ type: 'error', message: String(msg) }); }

async function runHO(k, pu) {
    postMessage({ type: 'status', message: 'Loading the WebAssembly engine…' });
    importScripts('ho.js');
    const M = await createHO({ onAbort: what => fail('The engine stopped: ' + what + ' (most likely the browser refused more memory).') });
    postMessage({ type: 'status', message: 'Computing…', started: true });
    const out = JSON.parse(M.ccall('ho_compute', 'string', ['number', 'number'], [k, pu ? 1 : 0]));
    if (out.error) return fail(out.error);
    postMessage({ type: 'done', format: 'bm-invariant', k, bm: out.bm });
}

async function runPy(k, surface) {
    postMessage({ type: 'status', message: 'Downloading Python (Pyodide, about 15 MB, only the first time)…' });
    importScripts(PYODIDE + 'pyodide.js');
    const py = await loadPyodide({ indexURL: PYODIDE });
    await py.loadPackage('numpy');
    for (const f of PY_FILES) {
        const src = await (await fetch('py/' + f)).text();
        py.FS.writeFile(f, src);
    }
    postMessage({ type: 'status', message: 'Computing…', started: true });
    const res = py.runPython(`
import sys, json
sys.path.insert(0, '.')
from colcx import bm_homology
k = ${k}
primes = [p for p in range(2, max(k, 2) + 1) if all(p % q for q in range(2, p))]
h = bm_homology(k, '${surface}', primes)
json.dumps({str(n): [b, {str(p): v for p, v in t.items()}] for n, (b, t) in h.items() if b or t})
`);
    postMessage({ type: 'done', format: 'bm-primary', k, bm: JSON.parse(res) });
}

onmessage = async e => {
    const { engine, k, pu, surface } = e.data;
    try {
        if (engine === 'ho') await runHO(k, pu);
        else await runPy(k, surface);
    } catch (err) {
        fail(err && err.message ? err.message : err);
    }
};
