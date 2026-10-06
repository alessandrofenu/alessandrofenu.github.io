// WebAssembly entry point for homology_optimized.cpp, used by compute.html (single-threaded).
// ho_compute(k, pu) runs the default pipeline of the command-line tool and returns the same
// JSON as `./homology_optimized k [pu] --json`, or {"error": "..."}.
#define main ho_cli_main
#include "homology_optimized.cpp"
#undef main
#include <emscripten/emscripten.h>

static string OUT;

extern "C" EMSCRIPTEN_KEEPALIVE const char* ho_compute(int k, int pu) {
    NTHREADS = 1;
    binom(66, 0);
    try {
        RES.build(k); RES.verify_identities(k, true);
        RESP.build(k); RESP.verify_identities(k, true);
        HomologyResult H = homology_padic(k, pu, false);
        ostringstream o;
        o << "{\"k\": " << k << ", \"punctured\": " << (pu ? "true" : "false") << ", \"bm\": {";
        set<int> ns; for (auto& [n, f] : H.free) ns.insert(n); for (auto& [n, t] : H.tors) ns.insert(n);
        bool first = true;
        for (int n : ns) {
            o << (first ? "" : ", ") << "\"" << n << "\": [" << (H.free.count(n) ? H.free.at(n) : 0) << ", ["; first = false;
            if (H.tors.count(n)) { bool f2 = true; for (auto& x : H.tors.at(n)) { o << (f2 ? "" : ", ") << "\"" << x.get_str() << "\""; f2 = false; } }
            o << "]]";
        }
        o << "}}";
        OUT = o.str();
    } catch (exception& e) {
        OUT = string("{\"error\": \"") + e.what() + "\"}";
    }
    return OUT.c_str();
}
