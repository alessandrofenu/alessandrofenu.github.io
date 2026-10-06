// =============================================================================================
// homology_optimized.cpp
//
// Exact integral (Borel-Moore and ordinary) homology of the unordered configuration spaces
//     B_k(Sigma_1)            (closed torus)          and
//     B_k(Sigma_1 \ {pt})     (punctured torus)
// computed from Napolitano's Borel-Moore cell complex, as in homology.py / reduction.py /
// napolitano.py, but through a mathematically equivalent complex that is exponentially smaller.
//
// Build:
//     g++ -O2 -std=c++17 -pthread -o homology_optimized homology_optimized.cpp -lgmpxx -lgmp
//
// Usage:
//     ./homology_optimized K [pu] [--bm] [--json]  exact homology for weight K (closed; 'pu' = punctured)
//                                                  (--snf-cert: use the Section 5 pipeline instead of Section 6b)
//     ./homology_optimized K1 K2 [pu]           every weight in K1..K2
//     ./homology_optimized --verify [KMAX] [--ref KREF] [pu|both]   (defaults: KMAX 26, KREF 13)
//     ./homology_optimized --benchmark K1 K2 [--ref KREF] [pu]
//     ./homology_optimized --reference K [pu]   the original algorithm on the full cell complex
//
// ---------------------------------------------------------------------------------------------
// THE MATHEMATICS (full proofs in the accompanying report, REPORT.md)
//
// (0) The complex of napolitano.py.  A cell (e,v;(m_1,v_1),...,(m_d,v_d)) is read as the
//     Hochschild chain  a(e,v) (x) [a(m_1,v_1)|...|a(m_d,v_d)]  of the graded algebra
//         A = Z{a(m,v) : m >= 0, v in {0,1}},   a(m,v) a(m',v') = P(m,m') a(m+m', v+v')
//     (zero if v+v' = 2), |a(m,v)| = m, weight m+v, with the derivation
//         D a(m,0) = -(1+(-1)^m) a(m-1,1),   D a(m,1) = 0.
//     As a graded algebra A = Gamma[x_2] (x) Lambda[x_1] (x) Z[y]/y^2 with a(2i+e,v) = g_i x_1^e y^v.
//     The five face types of napolitano.boundary are EXACTLY (d_int - d_ext) on
//     A (x)_{A^e} B(A,A,A), B the normalised two-sided bar resolution (a cell-by-cell identity,
//     checked in --verify); this is isomorphic to the standard Hochschild complex by
//     (-1)^{#letters}.  The punctured complex is (A/yA) (x)_{A^e} B(A,A,A).
//
// (1) Any semifree resolution R -> A of A as a dg A^e-module computes the same homology:
//     H(M (x)_{A^e} B) = H(M (x)_{A^e} R) for M = A, A/yA.  We use
//         R_0 = K_Gamma (x) K_x (x) K_y,
//     K_Gamma = two-sided bar resolution of Gamma[x_2]   (generators: words E in g_1,g_2,...),
//     K_x     = Koszul resolution of Lambda[x_1]         (generators u^[a]),
//     K_y     = 2-periodic resolution of Z[y]/y^2       (generators e_r),
//     a resolution of A_0 = (A, D=0) by Kunneth, and perturb its differential by D_R so that
//     R_D = (R_0, d_0 + D_R) is a semifree resolution of (A, D).  D_R is constructed by the
//     homological perturbation recursion and then VERIFIED on every generator:
//         (d_0 + D_R)^2 = 0,   eps o (d_0 + D_R) = D o eps,   D_R raises the pin number.
//     These identities plus the (finite, in each weight) pin filtration make R_D a
//     resolution; so M (x)_{A^e} R_D is chain homotopy equivalent to Napolitano's complex.
//     Size in weight 26: 139 144 cells (closed) instead of 3 389 154 437 772.
//
// (2) Every face of Napolitano's complex preserves Q = sum floor(m/2) + (number of pins):
//     the complex splits as a direct sum over Q (checked in --verify), and so does R_D.
//
// (2b) In weights <= k, after inverting the primes <= k/2, Gamma[x_2] may be replaced by the polynomial
//     algebra Z[x_2]:  phi: A' = Z[x_2](x)Lambda[x_1](x)Z[y]/y^2 -> A,  x_2^i |-> i! gamma_i,  is a map of
//     dg-algebras and an isomorphism there (checked on the whole basis in --verify).  Z[x_2] is polynomial, so
//     its Koszul resolution has ONE generator and the whole Hochschild complex of A' has O(k^2) cells (2 452 at
//     k = 35).  That tiny complex gives the exact rank over Q of every Q-block and degree, and the exact torsion
//     at every prime > k/2.  The ranks then determine the ranks of the big complex by the recursion
//     r_{n+1} = c_n - b_n - r_n, and each prime p <= k/2 is handled on the big complex by chain-level
//     elimination over Z/p^N (Bockstein stages), where nothing grows.  This is the default pipeline; the one of
//     (3) below is kept as an independent cross-check (--snf-cert).
//
// (3) Each Q-block is reduced by exact Gaussian elimination on unit (+-1) incidences
//     (algebraic Morse reduction -- the same operation as reduction.morse_reduce), then
//     every boundary matrix gets an exact Smith normal form:
//       * fraction-free elimination over Z[1/N_S] with S-unit pivots gives the exact rank
//         over Q and proves that no prime outside S divides an elementary divisor;
//       * for each p in S, elimination over Z/p^N with minimal-valuation pivots gives the
//         exact p-adic valuations (exact because the rank over Q is known).
//     All arithmetic is overflow-checked: int64, then int128, then GMP integers, chosen per block / per
//     matrix; the p-adic stage of (2b) needs only int64/int128.  Independent problems run in parallel
//     (threads; HO_THREADS=n to override).  Tables for every k <= 35 are in tables/ (see REPORT.md).
// =============================================================================================
#include <bits/stdc++.h>
#include <gmpxx.h>
#include <sys/resource.h>
#include <sys/wait.h>
#include <unistd.h>
using namespace std;
typedef long long ll;
typedef unsigned long long ull;
typedef __int128 i128;

// =============================================================================================
// Section 1.  Exact integer arithmetic: checked int64, GMP fallback
// =============================================================================================
static int MIN_TIER = 0;                               // testing only: 0 = int64, 1 = int128, 2 = GMP from the start
static bool FORCE_PLOCAL_MPZ = false;                  // testing only: use the GMP p-local routine (low precision first)
struct OverflowError : std::exception { const char* what() const noexcept override { return "integer overflow in an int64/int128 fast path"; } };
static inline ll zadd(ll a, ll b) { ll r; if (__builtin_add_overflow(a, b, &r)) throw OverflowError(); return r; }
static inline ll zsub(ll a, ll b) { ll r; if (__builtin_sub_overflow(a, b, &r)) throw OverflowError(); return r; }
static inline ll zmul(ll a, ll b) { ll r; if (__builtin_mul_overflow(a, b, &r)) throw OverflowError(); return r; }
static inline ll zneg(ll a) { if (a == LLONG_MIN) throw OverflowError(); return -a; }
static inline ll zabs(ll a) { return a < 0 ? zneg(a) : a; }
static inline ll zdivexact(ll a, ll b) { return a / b; }
static inline bool zunit(ll a) { return a == 1 || a == -1; }
static inline ll zgcd(ll a, ll b) { a = zabs(a); b = zabs(b); while (b) { ll t = a % b; a = b; b = t; } return a; }
static inline int zbits(ll a) { ull x = a < 0 ? (ull)(-(a + 1)) + 1 : (ull)a; return x ? 64 - __builtin_clzll(x) : 0; }
static inline bool zdivisible(ll a, ll p) { return a % p == 0; }
static inline ll zdivp(ll a, ll p) { return a / p; }
static inline mpz_class to_mpz(ll a) { mpz_class r; mpz_set_si(r.get_mpz_t(), a); return r; }

static const i128 I128_MAX = (i128)(((unsigned __int128)1 << 127) - 1), I128_MIN = -I128_MAX - 1;
static inline i128 zadd(i128 a, i128 b) { i128 r; if (__builtin_add_overflow(a, b, &r)) throw OverflowError(); return r; }
static inline i128 zsub(i128 a, i128 b) { i128 r; if (__builtin_sub_overflow(a, b, &r)) throw OverflowError(); return r; }
static inline i128 zmul(i128 a, i128 b) { i128 r; if (__builtin_mul_overflow(a, b, &r)) throw OverflowError(); return r; }
static inline i128 zneg(i128 a) { if (a == I128_MIN) throw OverflowError(); return -a; }
static inline i128 zabs(i128 a) { return a < 0 ? zneg(a) : a; }
static inline i128 zdivexact(i128 a, i128 b) { return a / b; }
static inline bool zunit(i128 a) { return a == 1 || a == -1; }
static inline i128 zgcd(i128 a, i128 b) {
    a = zabs(a); b = zabs(b);
    while (b) { if ((a >> 63) == 0 && (b >> 63) == 0) return zgcd((ll)a, (ll)b); i128 t = a % b; a = b; b = t; }
    return a;
}
static inline int zbits(i128 a) { unsigned __int128 x = a < 0 ? (unsigned __int128)(-(a + 1)) + 1 : (unsigned __int128)a; int b = 0; while (x) { x >>= 1; b++; } return b; }
static inline bool zdivisible(i128 a, ll p) { return a % p == 0; }
static inline i128 zdivp(i128 a, ll p) { return a / p; }
static inline mpz_class to_mpz(i128 a) {
    bool neg = a < 0; unsigned __int128 u = neg ? (unsigned __int128)(-(a + 1)) + 1 : (unsigned __int128)a;
    mpz_class hi, lo; mpz_set_ui(hi.get_mpz_t(), (unsigned long)(u >> 64)); mpz_set_ui(lo.get_mpz_t(), (unsigned long)(u & ~0ull));
    mpz_class r = (hi << 64) + lo; return neg ? mpz_class(-r) : r;
}
static inline ll mod_ll(ll x, ll m) { ll r = x % m; return r < 0 ? r + m : r; }
static inline ll mod_ll(i128 x, ll m) { i128 r = x % m; return (ll)(r < 0 ? r + m : r); }
static inline ll mod_ll(const mpz_class& x, ll m) { return (ll)mpz_fdiv_ui(x.get_mpz_t(), (unsigned long)m); }
static inline bool fits_ll(const mpz_class& x) { return mpz_fits_slong_p(x.get_mpz_t()); }

static inline mpz_class zadd(const mpz_class& a, const mpz_class& b) { return a + b; }
static inline mpz_class zsub(const mpz_class& a, const mpz_class& b) { return a - b; }
static inline mpz_class zmul(const mpz_class& a, const mpz_class& b) { return a * b; }
static inline mpz_class zneg(const mpz_class& a) { return -a; }
static inline mpz_class zabs(const mpz_class& a) { return abs(a); }
static inline mpz_class zdivexact(const mpz_class& a, const mpz_class& b) { mpz_class r; mpz_divexact(r.get_mpz_t(), a.get_mpz_t(), b.get_mpz_t()); return r; }
static inline bool zunit(const mpz_class& a) { return a == 1 || a == -1; }
static inline mpz_class zgcd(const mpz_class& a, const mpz_class& b) { mpz_class r; mpz_gcd(r.get_mpz_t(), a.get_mpz_t(), b.get_mpz_t()); return r; }
static inline int zbits(const mpz_class& a) { return a == 0 ? 0 : (int)mpz_sizeinbase(a.get_mpz_t(), 2); }
static inline bool zdivisible(const mpz_class& a, ll p) { return mpz_divisible_ui_p(a.get_mpz_t(), (unsigned long)p) != 0; }
static inline mpz_class zdivp(const mpz_class& a, ll p) { mpz_class r; mpz_divexact_ui(r.get_mpz_t(), a.get_mpz_t(), (unsigned long)p); return r; }
static inline mpz_class to_mpz(const mpz_class& a) { return a; }
template <class Z> static inline Z zfrom(ll a);
template <> inline ll zfrom<ll>(ll a) { return a; }
template <> inline mpz_class zfrom<mpz_class>(ll a) { return to_mpz(a); }
template <> inline i128 zfrom<i128>(ll a) { return (i128)a; }

static ll binom(int n, int r) {
    if (r < 0 || r > n) return 0;
    static vector<vector<ll>> T;
    if ((int)T.size() <= n) {
        int old = T.size(); T.resize(n + 1);
        for (int i = old; i <= n; i++) {
            T[i].assign(i + 1, 1);
            for (int j = 1; j < i; j++) T[i][j] = zadd(T[i - 1][j - 1], T[i - 1][j]);
        }
    }
    return T[n][r];
}
// the parity binomial P(a,b) of [IR26b, 6.1] = the Gaussian binomial [a+b choose a] at q = -1
static inline ll Ppar(int a, int b) { if ((a & 1) && (b & 1)) return 0; return binom((a + b) / 2, a / 2); }

static bool is_prime_ll(ull n) {           // deterministic Miller-Rabin for 64-bit integers
    if (n < 2) return false;
    for (ull p : {2ull, 3ull, 5ull, 7ull, 11ull, 13ull, 17ull, 19ull, 23ull, 29ull, 31ull, 37ull})
        if (n % p == 0) return n == p;
    ull d = n - 1; int s = 0; while (!(d & 1)) { d >>= 1; s++; }
    auto mulmod = [&](ull a, ull b) { return (ull)((unsigned __int128)a * b % n); };
    auto powmod = [&](ull a, ull e) { ull r = 1; while (e) { if (e & 1) r = mulmod(r, a); a = mulmod(a, a); e >>= 1; } return r; };
    for (ull a : {2ull, 3ull, 5ull, 7ull, 11ull, 13ull, 17ull, 19ull, 23ull, 29ull, 31ull, 37ull}) {
        ull x = powmod(a, d);
        if (x == 1 || x == n - 1) continue;
        bool comp = true;
        for (int i = 1; i < s; i++) { x = mulmod(x, x); if (x == n - 1) { comp = false; break; } }
        if (comp) return false;
    }
    return true;
}
// exact prime factorisation (distinct primes) of a non-zero integer; throws if it cannot be done
static vector<ll> prime_factors(mpz_class x) {
    x = abs(x);
    vector<ll> out;
    for (ll p = 2; p <= 1000000 && x > 1; p += (p == 2 ? 1 : 2)) {
        if (mpz_divisible_ui_p(x.get_mpz_t(), p)) {
            out.push_back(p);
            while (mpz_divisible_ui_p(x.get_mpz_t(), p)) mpz_divexact_ui(x.get_mpz_t(), x.get_mpz_t(), p);
        }
        if (mpz_class((long)p) * (long)p > x) break;
    }
    if (x > 1) {
        if (!mpz_fits_ulong_p(x.get_mpz_t()) || !is_prime_ll(mpz_get_ui(x.get_mpz_t())))
            throw runtime_error("prime_factors: cannot certify the factorisation of " + x.get_str());
        out.push_back((ll)mpz_get_ui(x.get_mpz_t()));
    }
    return out;
}

static double now_sec() { return chrono::duration<double>(chrono::steady_clock::now().time_since_epoch()).count(); }
static long peak_rss_kb() { struct rusage r; getrusage(RUSAGE_SELF, &r); return r.ru_maxrss; }

// =============================================================================================
// Section 2.  The algebra A = Gamma[x_2] (x) Lambda[x_1] (x) Z[y]/y^2 and A^e = A (x) A
// =============================================================================================
// basis element a(m,v) <-> id = 2m+v.  degree m, weight m+v, pins v, Q = floor(m/2)+v.
static inline int A_id(int m, int v) { return 2 * m + v; }
static inline int A_m(int id) { return id >> 1; }
static inline int A_v(int id) { return id & 1; }
struct Coef { ll c; int id; };
// poly = false: A  = Gamma[x_2] (x) Lambda[x_1] (x) Z[y]/y^2, a(2i+e,v) = gamma_i x_1^e y^v   (Napolitano's algebra)
// poly = true : A' = Z[x_2]    (x) Lambda[x_1] (x) Z[y]/y^2, a(2i+e,v) = x_2^i  x_1^e y^v
//               with D(x_2^i) = -2i x_2^(i-1) x_1 y.  x_2^i -> i! gamma_i is a dg-algebra map A' -> A,
//               an isomorphism in weights <= k once the primes <= k/2 are inverted (Section 6b).
static inline Coef mulA(int x, int y, bool poly = false) {       // a(m,v) a(m',v') = P(m,m') a(m+m',v+v')
    if (A_v(x) + A_v(y) > 1) return {0, -1};
    ll c = poly ? (((A_m(x) & 1) && (A_m(y) & 1)) ? 0 : 1) : Ppar(A_m(x), A_m(y));
    if (!c) return {0, -1};
    return {c, A_id(A_m(x) + A_m(y), A_v(x) + A_v(y))};
}
static inline Coef DA(int x, bool poly = false) {                // D a(m,0) = -(1+(-1)^m) a(m-1,1)   (A')  -m a(m-1,1)
    int m = A_m(x);
    if (A_v(x) || (m & 1) || m == 0) return {0, -1};
    return {poly ? -(ll)m : -2, A_id(m - 1, 1)};
}
// A^e element a (x) b; product (a(x)b)(c(x)d) = (-1)^{|b||c|} ac (x) bd;  D^e = D(x)1 + 1(x)D (Koszul)
struct Lam { int a, b; };
static inline int degE(Lam l) { return A_m(l.a) + A_m(l.b); }
static inline int pinsE(Lam l) { return A_v(l.a) + A_v(l.b); }

// =============================================================================================
// Section 3.  The small resolution R_D
// =============================================================================================
// A generator is (E, a, r): E = (i_1..i_s), i_j >= 1, a word in the divided powers g_i = x_2^[i]
// (bar resolution of Gamma), a = exponent of u (Koszul resolution of Lambda[x_1]), r = index of
// e_r (periodic resolution of Z[y]/y^2).  degree sum(2i_j+1)+2a+r, weight 2 sum i_j + a + r,
// pins r, Q = sum i_j + r.
struct Gen { vector<int> E; int a, r; int deg, wt, Q, sumE; };
struct Term { int a, b, g; ll c; };                               // c * (a(x)b) . generator g
typedef vector<Term> Elem;

static inline ull tkey(const Term& t) { return ((ull)(unsigned)t.g << 20) | ((ull)t.a << 10) | (ull)t.b; }
static void normalize(Elem& e) {
    sort(e.begin(), e.end(), [](const Term& x, const Term& y) { return tkey(x) < tkey(y); });
    size_t w = 0;
    for (size_t i = 0; i < e.size();) {
        size_t j = i; ll s = 0;
        while (j < e.size() && tkey(e[j]) == tkey(e[i])) { s = zadd(s, e[j].c); j++; }
        if (s) { e[w] = e[i]; e[w].c = s; w++; }
        i = j;
    }
    e.resize(w);
}

struct Resolution {
    bool poly = false;                           // false: R_D for A (bar resolution of Gamma); true: tiny R'_D for A' (Koszul in x_2)
    int K = -1;                                  // generators of weight <= K are built
    vector<Gen> gens;
    unordered_map<ull, int> gid;
    vector<Elem> d0, DR;                         // d_0(g) and D_R(g) as elements of R
    vector<vector<int>> byWeight;                // generator ids by weight
    int unitGen = -1;
    double t_build = 0;

    static ull gkey(const vector<int>& E, int a, int r) {
        ull mask = 0; int j = 0;
        for (size_t t = 0; t < E.size(); t++) { j += E[t]; if (t + 1 < E.size()) mask |= 1ull << (j - 1); }
        return ((ull)j << 44) | (mask << 12) | ((ull)a << 6) | (ull)r;
    }
    int find(const vector<int>& E, int a, int r) const {
        auto it = gid.find(gkey(E, a, r));
        if (it == gid.end()) throw runtime_error("internal: generator not found");
        return it->second;
    }
    int add_gen(const vector<int>& E, int a, int r) {
        Gen g; g.E = E; g.a = a; g.r = r; g.sumE = 0; g.deg = 2 * a + r;
        for (int i : E) { g.sumE += i; g.deg += 2 * i + 1; }
        g.wt = 2 * g.sumE + a + r; g.Q = g.sumE + r;
        int id = gens.size(); gens.push_back(g); gid[gkey(E, a, r)] = id;
        return id;
    }

    // d_0 on a generator: bar differential of Gamma, Koszul differential, periodic differential
    Elem d0_gen(int gi) const {
        const Gen& g = gens[gi];
        Elem out; int s = g.E.size();
        const int ONE = A_id(0, 0);
        if (s) {
            vector<int> tail(g.E.begin() + 1, g.E.end());
            out.push_back({A_id(2 * g.E[0], 0), ONE, find(tail, g.a, g.r), 1});
            for (int j = 1; j < s; j++) {
                vector<int> F(g.E.begin(), g.E.begin() + j - 1);
                F.push_back(g.E[j - 1] + g.E[j]);
                F.insert(F.end(), g.E.begin() + j + 1, g.E.end());
                ll c = binom(g.E[j - 1] + g.E[j], g.E[j - 1]);
                out.push_back({ONE, ONE, find(F, g.a, g.r), (j & 1) ? -c : c});
            }
            vector<int> head(g.E.begin(), g.E.end() - 1);
            out.push_back({ONE, A_id(2 * g.E[s - 1], 0), find(head, g.a, g.r), (s & 1) ? -1 : 1});
        }
        if (g.a) {
            int h = find(g.E, g.a - 1, g.r);
            out.push_back({A_id(1, 0), ONE, h, 1});
            out.push_back({ONE, A_id(1, 0), h, -1});
        }
        if (g.r) {
            int h = find(g.E, g.a, g.r - 1);
            ll sg = (s & 1) ? -1 : 1;
            out.push_back({A_id(0, 1), ONE, h, sg});
            out.push_back({ONE, A_id(0, 1), h, (g.r & 1) ? -sg : sg});
        }
        normalize(out);
        return out;
    }
    // lam * elem, with the sign (-1)^{|lam|} of the dg-module rule d(lam g) = D^e(lam) g + (-1)^{|lam|} lam dg
    void lmul_into(Lam l, const Elem& el, ll c, Elem& out) const {
        for (const Term& t : el) {
            ll s = ((A_m(l.b) * A_m(t.a)) & 1) ? -1 : 1;
            Coef c1 = mulA(l.a, t.a, poly); if (!c1.c) continue;
            Coef c2 = mulA(l.b, t.b, poly); if (!c2.c) continue;
            out.push_back({c1.id, c2.id, t.g, zmul(zmul(zmul(s, c), zmul(c1.c, c2.c)), t.c)});
        }
    }
    Elem apply_d0(const Elem& el) const {
        Elem out;
        for (const Term& t : el) lmul_into({t.a, t.b}, d0[t.g], (degE({t.a, t.b}) & 1) ? -t.c : t.c, out);
        normalize(out); return out;
    }
    Elem apply_DR(const Elem& el) const {
        Elem out;
        for (const Term& t : el) {
            Coef da = DA(t.a, poly);
            if (da.c) out.push_back({da.id, t.b, t.g, zmul(t.c, da.c)});
            Coef db = DA(t.b, poly);
            if (db.c) out.push_back({t.a, db.id, t.g, zmul(t.c, (A_m(t.a) & 1) ? -db.c : db.c)});
            lmul_into({t.a, t.b}, DR[t.g], (degE({t.a, t.b}) & 1) ? -t.c : t.c, out);
        }
        normalize(out); return out;
    }
    Elem apply_dT(const Elem& el) const {
        Elem a = apply_d0(el), b = apply_DR(el);
        a.insert(a.end(), b.begin(), b.end()); normalize(a); return a;
    }
    Elem eps_map(const Elem& el) const {                          // augmentation R -> A, as {A-id: coef}
        Elem out;
        for (const Term& t : el) if (t.g == unitGen) {
            Coef c = mulA(t.a, t.b, poly);
            if (c.c) out.push_back({c.id, 0, 0, zmul(c.c, t.c)});
        }
        normalize(out); return out;
    }
    // Z-linear contraction h_0 of R_0 = K_Gamma (x) K_x (x) K_y (the standard tensor-trick contraction
    // built from h_Gamma(a'[E]a'') = [a'|E]a'', the Koszul contraction in z = x'-x'', and the
    // periodic contraction in the basis {1, y'', t_{r+1} = y'+(-1)^{r+1}y'', y'y''}).  It is only used
    // to FIND D_R; correctness rests on the identities verified afterwards, not on h_0.
    void h0_term(const Term& t, Elem& out) const {
        int m = A_m(t.a), v = A_v(t.a), m2 = A_m(t.b), v2 = A_v(t.b);
        const Gen& g = gens[t.g];
        int i = m / 2, ep = m % 2, i2 = m2 / 2, ep2 = m2 % 2;
        ll sg1 = ((ep + ep2) & 1) ? -1 : 1;
        if (!poly && i >= 1) {
            vector<int> E2; E2.push_back(i); E2.insert(E2.end(), g.E.begin(), g.E.end());
            out.push_back({A_id(ep, v), A_id(2 * i2 + ep2, v2), find(E2, g.a, g.r), zmul(sg1, t.c)});
        }
        if (poly && i >= 1 && g.E.empty()) {      // Koszul: h(x'^i x''^i2) = sum_{j<i} x'^j x''^(i-1-j+i2) v,  v = [x_2]
            int h = find({1}, g.a, g.r);
            for (int j = 0; j < i; j++) out.push_back({A_id(2 * j + ep, v), A_id(2 * (i - 1 - j + i2) + ep2, v2), h, zmul(sg1, t.c)});
        }
        if (!g.E.empty()) return;
        ll cg = poly ? 1 : binom(i + i2, i);
        if (ep == 1) {
            int h = find({}, g.a + 1, g.r);
            if (ep2 == 0) out.push_back({A_id(0, v), A_id(2 * (i + i2), v2), h, zmul(cg, t.c)});
            else out.push_back({A_id(0, v), A_id(2 * (i + i2) + 1, v2), h, zmul(cg, t.c)});
        }
        if (g.a) return;
        if (ep && ep2) return;
        int xe = ep + ep2;
        ll sg = (xe & 1) ? -1 : 1;
        if (v == 1) {
            int h = find({}, 0, g.r + 1);
            out.push_back({A_id(0, 0), A_id(2 * (i + i2) + xe, v2), h, zmul(zmul(sg, cg), t.c)});
        }
    }
    Elem apply_h0(const Elem& el) const {
        Elem out; for (const Term& t : el) h0_term(t, out); normalize(out); return out;
    }

    void gens_of_weight(int w) {
        // compositions E of j, then a + r = w - 2j   (poly: E = () or (1) only -- the Koszul generator v = [x_2])
        for (int j = 0; 2 * j <= w && (!poly || j <= 1); j++) {
            vector<vector<int>> comps;
            if (j == 0) comps.push_back({});
            else for (ull mask = 0; mask < (1ull << (j - 1)); mask++) {
                vector<int> E; int cur = 1;
                for (int t = 0; t < j - 1; t++) { if (mask >> t & 1) { E.push_back(cur); cur = 1; } else cur++; }
                E.push_back(cur); comps.push_back(E);
            }
            for (auto& E : comps) for (int a = 0; a <= w - 2 * j; a++) {
                int id = add_gen(E, a, w - 2 * j - a);
                byWeight[w].push_back(id);
            }
        }
    }

    void build(int Kmax) {
        if (Kmax <= K) return;
        double t0 = now_sec();
        byWeight.resize(Kmax + 1);
        for (int w = K + 1; w <= Kmax; w++) gens_of_weight(w);
        if (unitGen < 0) unitGen = find({}, 0, 0);
        d0.resize(gens.size()); DR.resize(gens.size());
        for (int w = K + 1; w <= Kmax; w++) {
            vector<int> ids = byWeight[w];
            stable_sort(ids.begin(), ids.end(), [&](int x, int y) { return gens[x].deg < gens[y].deg; });
            for (int gi : ids) d0[gi] = d0_gen(gi);
            for (int gi : ids) {
                if (gi == unitGen) { DR[gi].clear(); continue; }
                // D_R(g) = h_inf(Y), Y = -D_R(d_0 g), h_inf = sum_n h_0 (-D_R h_0)^n  (perturbation lemma)
                Elem Y = apply_DR(d0[gi]);
                for (auto& t : Y) t.c = zneg(t.c);
                Elem X, W = Y;
                int guard = 0;
                while (!W.empty()) {
                    Elem Xn = apply_h0(W);
                    X.insert(X.end(), Xn.begin(), Xn.end());
                    W = apply_DR(Xn);
                    for (auto& t : W) t.c = zneg(t.c);
                    if (++guard > 4 * Kmax + 10) throw runtime_error("perturbation series did not terminate");
                }
                normalize(X);
                DR[gi] = X;
            }
        }
        K = Kmax;
        t_build += now_sec() - t0;
    }

    // The identities that make R_D a semifree resolution of (A,D) (see header, (1)).
    // Returns the number of generators checked; throws on the first failure.
    long verify_identities(int Kmax, bool quiet = false) const {
        long n = 0;
        for (int w = 0; w <= Kmax; w++) for (int gi : byWeight[w]) {
            const Gen& g = gens[gi];
            Elem dT = d0[gi]; dT.insert(dT.end(), DR[gi].begin(), DR[gi].end()); normalize(dT);
            if (!apply_dT(dT).empty()) throw runtime_error("R_D: (d_0 + D_R)^2 != 0 on a generator");
            if (!apply_d0(d0[gi]).empty()) throw runtime_error("R_0: d_0^2 != 0 on a generator");
            if (gi != unitGen && !eps_map(dT).empty()) throw runtime_error("R_D: eps o d_T != D o eps");
            for (const Term& t : dT) {
                const Gen& h = gens[t.g];
                Lam l{t.a, t.b};
                if (degE(l) + h.deg != g.deg - 1) throw runtime_error("R_D: d_T not of degree -1");
                if (A_m(t.a) + A_v(t.a) + A_m(t.b) + A_v(t.b) + h.wt != g.wt) throw runtime_error("R_D: weight");
                int q = A_m(t.a) / 2 + A_v(t.a) + A_m(t.b) / 2 + A_v(t.b) + h.Q;
                if (q != g.Q) throw runtime_error("R_D: Q-grading violated");
                if (h.wt >= g.wt && !(h.wt == g.wt && h.deg < g.deg)) throw runtime_error("R_D: not semifree-ordered");
            }
            for (const Term& t : DR[gi])
                if (pinsE({t.a, t.b}) + gens[t.g].r <= g.r) throw runtime_error("R_D: D_R does not raise pins");
            n++;
        }
        if (!quiet) fprintf(stderr, "    [resolution] identities verified on %ld generators (weight <= %d)\n", n, Kmax);
        return n;
    }
};

// =============================================================================================
// Section 4.  Based chain complexes, exact unit (Morse) reduction, exact Smith normal form
// =============================================================================================
struct BlockComplex {                            // a based complex over Z (one Q-block)
    vector<int> deg;                             // degree of each cell
    vector<vector<pair<int, ll>>> bnd;           // boundary, sorted by face index
};

// ---- exact Gaussian elimination on +-1 incidences (chain-level; a chain homotopy equivalence) ----
template <class Z>
struct ChainReducer {
    int n;
    vector<int> deg;
    vector<vector<pair<int, Z>>> d;
    vector<vector<int>> up;                      // lazy co-boundary lists (validated on use)
    vector<int> upcnt;                           // exact co-boundary sizes
    vector<char> alive;
    map<int, long> pivots;                       // number of cancelled pairs (t,s) with deg t = n
    int maxbits = 0;

    explicit ChainReducer(const BlockComplex& B) {
        n = B.deg.size(); deg = B.deg; d.resize(n); up.resize(n); upcnt.assign(n, 0); alive.assign(n, 1);
        for (int c = 0; c < n; c++) {
            d[c].reserve(B.bnd[c].size());
            for (auto& fc : B.bnd[c]) { d[c].push_back({fc.first, zfrom<Z>(fc.second)}); up[fc.first].push_back(c); upcnt[fc.first]++; }
        }
    }
    const Z* coef(int t, int s) const {
        auto& v = d[t];
        auto it = lower_bound(v.begin(), v.end(), s, [](const pair<int, Z>& a, int b) { return a.first < b; });
        return (it != v.end() && it->first == s) ? &it->second : nullptr;
    }
    typedef tuple<ll, int, int> HE;
    priority_queue<HE, vector<HE>, greater<HE>> heap;
    void push(int t) {
        if (!alive[t]) return;
        ll best = LLONG_MAX; int bs = -1;
        ll lt = (ll)d[t].size() - 1;
        for (auto& fc : d[t]) if (zunit(fc.second)) {
            ll cost = lt * (ll)(upcnt[fc.first] - 1);
            if (cost < best) { best = cost; bs = fc.first; }
        }
        if (bs >= 0) heap.push({best, t, bs});
    }
    void run() {
        for (int t = 0; t < n; t++) push(t);
        vector<int> mark(n, 0); int stamp = 0;
        while (!heap.empty()) {
            auto [k0, t, s] = heap.top(); heap.pop();
            if (!alive[t]) continue;
            const Z* cp = (alive[s] ? coef(t, s) : nullptr);
            if (!cp || !zunit(*cp)) { push(t); continue; }
            ll cost = ((ll)d[t].size() - 1) * (ll)(upcnt[s] - 1);
            if (cost > k0) { push(t); continue; }
            Z eps = *cp;
            vector<pair<int, Z>> dt = d[t];
            // d(tp) <- d(tp) - (<d tp, s> / eps) d(t)  for every other tp with s in d(tp)
            vector<int> touched; stamp++;
            for (int tp : up[s]) {
                if (tp == t || !alive[tp] || mark[tp] == stamp) continue;
                mark[tp] = stamp;
                const Z* c = coef(tp, s);
                if (!c) continue;
                touched.push_back(tp);
                Z lam = zmul(*c, eps);                 // eps = +-1, so 1/eps = eps
                vector<pair<int, Z>> nw; nw.reserve(d[tp].size() + dt.size());
                auto& old = d[tp];
                size_t i = 0, j = 0;
                while (i < old.size() || j < dt.size()) {
                    if (j == dt.size() || (i < old.size() && old[i].first < dt[j].first)) { nw.push_back(old[i]); i++; }
                    else if (i == old.size() || dt[j].first < old[i].first) {
                        Z v = zneg(zmul(lam, dt[j].second));
                        nw.push_back({dt[j].first, v}); up[dt[j].first].push_back(tp); upcnt[dt[j].first]++;
                        maxbits = max(maxbits, zbits(v)); j++;
                    } else {
                        Z v = zsub(old[i].second, zmul(lam, dt[j].second));
                        if (v != 0) { nw.push_back({old[i].first, v}); maxbits = max(maxbits, zbits(v)); }
                        else upcnt[old[i].first]--;
                        i++; j++;
                    }
                }
                d[tp].swap(nw);
            }
            // drop t (and its row in d_{n+1}), and s (now only in d(t))
            for (int y : up[t]) if (alive[y]) {
                auto& v = d[y];
                auto it = lower_bound(v.begin(), v.end(), t, [](const pair<int, Z>& a, int b) { return a.first < b; });
                if (it != v.end() && it->first == t) v.erase(it);
            }
            for (auto& fc : d[t]) upcnt[fc.first]--;
            d[t].clear(); d[t].shrink_to_fit(); alive[t] = 0;
            for (int y : up[s]) if (alive[y] && coef(y, s)) throw runtime_error("internal: unit reduction left s in a boundary");
            for (auto& fc : d[s]) upcnt[fc.first]--;
            d[s].clear(); d[s].shrink_to_fit(); alive[s] = 0;
            up[t].clear(); up[t].shrink_to_fit(); up[s].clear(); up[s].shrink_to_fit();
            pivots[deg[t]]++;
            for (int tp : touched) push(tp);
        }
    }
};

// ---- exact Smith normal form of a sparse integer matrix (given by rows) ----
struct SNFResult { int rank = 0; vector<ll> S; map<ll, vector<int>> val; int maxbits = 0; };
// The elimination of Section 5 starts from this set of primes (they are treated as units).  For the
// tiny complex of Section 6b it is set to every prime <= k/2, i.e. the elimination runs over Lambda_k,
// and SNF_ONLY_NEW then asks only for the p-parts at the primes outside it.
static vector<ll> SNF_S0 = {2, 3, 5, 7};
static bool SNF_ONLY_NEW = false;

// (a) rank over Q and a prime set S with coker free over Z[1/prod S]; every step is invertible over Z[1/prod S]
template <class Z>
static int rank_and_primes(vector<vector<pair<int, Z>>> rows, int ncols, vector<ll>& S, int& maxbits) {
    auto spart_div = [&](Z x) { for (ll p : S) while (x != 0 && zdivisible(x, p)) x = zdivp(x, p); return x; };
    auto is_sunit = [&](const Z& x) { Z r = spart_div(x); return zunit(r); };
    vector<vector<int>> colrows(ncols); vector<int> colcnt(ncols, 0);
    vector<int> act;
    for (int i = 0; i < (int)rows.size(); i++) if (!rows[i].empty()) {
        act.push_back(i);
        for (auto& e : rows[i]) { colrows[e.first].push_back(i); colcnt[e.first]++; }
    }
    vector<char> live(rows.size(), 0); for (int i : act) live[i] = 1;
    int rank = 0;
    vector<int> mark(rows.size(), 0); int stamp = 0;
    while (true) {
        // prune empty rows
        { size_t w = 0; for (int i : act) if (live[i] && !rows[i].empty()) act[w++] = i; else live[i] = 0; act.resize(w); }
        if (act.empty()) break;
        // pivot: an S-unit entry; prefer +-1, then Markowitz, then small |x|
        // pivot: an S-unit entry; +-1 first, then least Markowitz fill-in, then small |x|
        int pi = -1, pj = -1; tuple<int, ll, int> best{2, LLONG_MAX, INT_MAX};
        for (int i : act) {
            ll lr = (ll)rows[i].size() - 1;
            for (auto& e : rows[i]) {
                tuple<int, ll, int> cost{zunit(e.second) ? 0 : 1, lr * (ll)(colcnt[e.first] - 1), zbits(e.second)};
                if (cost < best && is_sunit(e.second)) { best = cost; pi = i; pj = e.first; }
            }
            if (get<0>(best) == 0 && get<1>(best) == 0) break;
        }
        if (pi < 0) {                                         // no S-unit entry: enlarge S
            Z mn = 0; bool first = true;
            for (int i : act) for (auto& e : rows[i]) { Z a = zabs(e.second); if (first || a < mn) { mn = a; first = false; } }
            for (ll p : prime_factors(to_mpz(mn))) if (find(S.begin(), S.end(), p) == S.end()) S.push_back(p);
            sort(S.begin(), S.end());
            continue;
        }
        vector<pair<int, Z>> prow; prow.swap(rows[pi]); live[pi] = 0;
        for (auto& e : prow) colcnt[e.first]--;
        Z c = prow[lower_bound(prow.begin(), prow.end(), pj, [](const pair<int, Z>& a, int b) { return a.first < b; }) - prow.begin()].second;
        stamp++;
        for (int i : colrows[pj]) {
            if (!live[i] || mark[i] == stamp) continue;
            mark[i] = stamp;
            auto& r = rows[i];
            auto it = lower_bound(r.begin(), r.end(), pj, [](const pair<int, Z>& a, int b) { return a.first < b; });
            if (it == r.end() || it->first != pj) continue;
            Z b = it->second;
            Z g = zgcd(c, b); Z cm = zdivexact(c, g), bm = zdivexact(b, g);   // cm is an S-unit
            vector<pair<int, Z>> nw; nw.reserve(r.size() + prow.size());
            size_t x = 0, y = 0;
            while (x < r.size() || y < prow.size()) {
                if (y == prow.size() || (x < r.size() && r[x].first < prow[y].first)) { nw.push_back({r[x].first, zmul(cm, r[x].second)}); x++; }
                else if (x == r.size() || prow[y].first < r[x].first) { nw.push_back({prow[y].first, zneg(zmul(bm, prow[y].second))}); y++; }
                else { Z v = zsub(zmul(cm, r[x].second), zmul(bm, prow[y].second)); if (v != 0) nw.push_back({r[x].first, v}); x++; y++; }
            }
            // remove the S-part of the content (invertible over Z[1/prod S])
            if (!nw.empty()) {
                Z gg = 0; for (auto& e : nw) { gg = zgcd(gg, e.second); if (zunit(gg)) break; }
                if (!zunit(gg)) {
                    Z rest = spart_div(gg); Z sp = zdivexact(gg, rest);
                    if (!zunit(sp)) for (auto& e : nw) e.second = zdivexact(e.second, sp);
                }
            }
            // update column bookkeeping
            for (auto& e : r) colcnt[e.first]--;
            for (auto& e : nw) { colcnt[e.first]++; colrows[e.first].push_back(i); maxbits = max(maxbits, zbits(e.second)); }
            r.swap(nw);
        }
        rank++;
    }
    return rank;
}

// (b) p-adic valuations of the non-zero elementary divisors, by elimination over Z/p^N.
//     Exact when the rank over Q is known: the SNF of (A mod p^N) is (SNF A) mod p^N.
template <class Z>
static bool plocal_ll(const vector<vector<pair<int, Z>>>& rows_in, int ncols, ll p, int rank, int N, vector<int>& vals) {
    // entries reduced into [0, p^N), p^N < 2^62; products in 128 bits
    ll pN = 1; for (int i = 0; i < N; i++) pN *= p;
    vector<vector<pair<int, ll>>> rows(rows_in.size());
    for (size_t i = 0; i < rows_in.size(); i++)
        for (auto& e : rows_in[i]) { ll x = mod_ll(e.second, pN); if (x) rows[i].push_back({e.first, x}); }
    vector<vector<int>> colrows(ncols); vector<int> colcnt(ncols, 0);
    vector<char> live(rows.size(), 0); vector<int> act;
    for (int i = 0; i < (int)rows.size(); i++) if (!rows[i].empty()) { live[i] = 1; act.push_back(i); for (auto& e : rows[i]) { colrows[e.first].push_back(i); colcnt[e.first]++; } }
    auto val = [&](ll x) { int v = 0; while (x % p == 0) { x /= p; v++; } return v; };
    vector<int> mark(rows.size(), 0); int stamp = 0;
    vals.clear();
    while ((int)vals.size() < rank) {
        { size_t w = 0; for (int i : act) if (live[i] && !rows[i].empty()) act[w++] = i; else live[i] = 0; act.resize(w); }
        int pi = -1, pj = -1; pair<int, ll> best{INT_MAX, LLONG_MAX};
        for (int i : act) {
            ll lr = (ll)rows[i].size() - 1;
            for (auto& e : rows[i]) {
                pair<int, ll> cost{val(e.second), lr * (ll)(colcnt[e.first] - 1)};
                if (cost < best) { best = cost; pi = i; pj = e.first; }
            }
            if (best.first == 0 && best.second == 0) break;
        }
        if (pi < 0) return false;                               // active block = 0 mod p^N
        int v = best.first;
        vector<pair<int, ll>> prow; prow.swap(rows[pi]); live[pi] = 0;
        for (auto& e : prow) colcnt[e.first]--;
        ll pv = 1; for (int i = 0; i < v; i++) pv *= p;
        ll u = 0; for (auto& e : prow) if (e.first == pj) u = e.second / pv;
        mpz_class U = to_mpz(u), M = to_mpz(pN), Ui; mpz_invert(Ui.get_mpz_t(), U.get_mpz_t(), M.get_mpz_t());
        ll uinv = Ui.get_si();
        stamp++;
        for (int i : colrows[pj]) {
            if (!live[i] || mark[i] == stamp) continue;
            mark[i] = stamp;
            auto& r = rows[i];
            auto it = lower_bound(r.begin(), r.end(), pj, [](const pair<int, ll>& a, int b) { return a.first < b; });
            if (it == r.end() || it->first != pj) continue;
            ll lam = (ll)((i128)(it->second / pv) * uinv % pN);
            vector<pair<int, ll>> nw; nw.reserve(r.size() + prow.size());
            size_t x = 0, y = 0;
            while (x < r.size() || y < prow.size()) {
                if (y == prow.size() || (x < r.size() && r[x].first < prow[y].first)) { nw.push_back(r[x]); x++; }
                else {
                    ll base = (x < r.size() && r[x].first == prow[y].first) ? r[x].second : 0;
                    ll w = (ll)(((i128)base - (i128)lam * prow[y].second) % pN); if (w < 0) w += pN;
                    if (w) nw.push_back({prow[y].first, w});
                    if (x < r.size() && r[x].first == prow[y].first) x++;
                    y++;
                }
            }
            for (auto& e : r) colcnt[e.first]--;
            for (auto& e : nw) { colcnt[e.first]++; colrows[e.first].push_back(i); }
            r.swap(nw);
        }
        vals.push_back(v);
    }
    sort(vals.begin(), vals.end());
    return true;
}
static bool plocal_mpz(const vector<vector<pair<int, mpz_class>>>& rows_in, int ncols, ll p, int rank, int N, vector<int>& vals) {
    mpz_class pN; mpz_ui_pow_ui(pN.get_mpz_t(), p, N);
    vector<map<int, mpz_class>> rows(rows_in.size());
    for (size_t i = 0; i < rows_in.size(); i++) for (auto& e : rows_in[i]) { mpz_class x; mpz_fdiv_r(x.get_mpz_t(), e.second.get_mpz_t(), pN.get_mpz_t()); if (x != 0) rows[i][e.first] = x; }
    vector<char> live(rows.size(), 1);
    auto val = [&](mpz_class x) { int v = 0; while (mpz_divisible_ui_p(x.get_mpz_t(), p)) { mpz_divexact_ui(x.get_mpz_t(), x.get_mpz_t(), p); v++; } return v; };
    vals.clear();
    while ((int)vals.size() < rank) {
        int pi = -1, pj = -1, bv = INT_MAX;
        for (size_t i = 0; i < rows.size(); i++) if (live[i]) for (auto& e : rows[i]) { int v = val(e.second); if (v < bv) { bv = v; pi = i; pj = e.first; } }
        if (pi < 0) return false;
        auto prow = rows[pi]; live[pi] = 0; rows[pi].clear();
        mpz_class pv; mpz_ui_pow_ui(pv.get_mpz_t(), p, bv);
        mpz_class u = prow[pj] / pv, ui; mpz_invert(ui.get_mpz_t(), u.get_mpz_t(), pN.get_mpz_t());
        for (size_t i = 0; i < rows.size(); i++) if (live[i] && rows[i].count(pj)) {
            mpz_class lam = (rows[i][pj] / pv) * ui % pN;
            for (auto& e : prow) {
                mpz_class w = rows[i][e.first] - lam * e.second; mpz_fdiv_r(w.get_mpz_t(), w.get_mpz_t(), pN.get_mpz_t());
                if (w != 0) rows[i][e.first] = w; else rows[i].erase(e.first);
            }
        }
        vals.push_back(bv);
    }
    sort(vals.begin(), vals.end());
    return true;
}

template <class Z>
static SNFResult snf_exact(const vector<vector<pair<int, Z>>>& rows, int ncols) {
    SNFResult R;
    R.S = SNF_S0;
    vector<ll> S0 = SNF_S0;
    R.rank = rank_and_primes<Z>(rows, ncols, R.S, R.maxbits);
    for (ll p : R.S) {
        if (SNF_ONLY_NEW && find(S0.begin(), S0.end(), p) != S0.end()) { R.val[p] = vector<int>(R.rank, 0); continue; }
        vector<int> v;
        int N = 0; { ll x = 1; while (x <= (1ll << 62) / p) { x *= p; N++; } }            // p^N < 2^62
        if (FORCE_PLOCAL_MPZ || !plocal_ll<Z>(rows, ncols, p, R.rank, N, v)) {
            vector<vector<pair<int, mpz_class>>> rm(rows.size());
            for (size_t i = 0; i < rows.size(); i++) for (auto& e : rows[i]) rm[i].push_back({e.first, to_mpz(e.second)});
            int NN = FORCE_PLOCAL_MPZ ? 4 : 2 * N;
            while (!plocal_mpz(rm, ncols, p, R.rank, NN, v)) NN *= 2;
        }
        R.val[p] = v;
    }
    return R;
}

// the elementary divisors from the p-local valuations: d_t = prod_p p^{val_p[t]}
static vector<mpz_class> divisors_from(const SNFResult& R) {
    vector<mpz_class> out(R.rank, 1);
    for (auto& [p, v] : R.val) for (int t = 0; t < R.rank; t++) { mpz_class q; mpz_ui_pow_ui(q.get_mpz_t(), p, v[t]); out[t] *= q; }
    return out;
}

// =============================================================================================
// Section 5.  Homology of a based complex split into blocks
// =============================================================================================
struct DegData { long cells = 0, cells_red = 0, rank = 0, pivots = 0; vector<mpz_class> divs; };
struct BlockResult { map<int, DegData> deg; int maxbits_morse = 0, maxbits_snf = 0; vector<ll> S; bool used_gmp = false; };

static double T_MORSE = 0, T_SNF = 0;                // wall-clock time of the two phases
static int NTHREADS = max(1, (int)thread::hardware_concurrency());

template <class F> static void parallel_for(int n, F f) {
    int nt = min(NTHREADS, n);
    if (nt <= 1) { for (int i = 0; i < n; i++) f(i); return; }
    atomic<int> next(0); exception_ptr err = nullptr; mutex mu; vector<thread> th;
    for (int t = 0; t < nt; t++) th.emplace_back([&] {
        for (;;) { int i = next++; if (i >= n) break;
            try { f(i); } catch (...) { lock_guard<mutex> g(mu); if (!err) err = current_exception(); } } });
    for (auto& x : th) x.join();
    if (err) rethrow_exception(err);
}

// exact chain-level unit reduction in the cheapest arithmetic that does not overflow
struct Reduced { vector<int> alive; vector<vector<pair<int, ll>>> bl; vector<vector<pair<int, mpz_class>>> bz; bool big = false; map<int, long> piv; int maxbits = 0; int tier = 0; };
template <class Z> static void chain_reduce_into(const BlockComplex& B, Reduced& R) {
    ChainReducer<Z> CR(B); CR.run();
    R.maxbits = CR.maxbits; R.piv = CR.pivots; R.alive.clear();
    for (int c = 0; c < CR.n; c++) if (CR.alive[c]) R.alive.push_back(c);
    if constexpr (is_same<Z, ll>::value) { R.big = false; R.bl.assign(CR.n, {}); for (int c : R.alive) R.bl[c] = CR.d[c]; }
    else { R.big = true; R.bz.assign(CR.n, {}); for (int c : R.alive) for (auto& e : CR.d[c]) R.bz[c].push_back({e.first, to_mpz(e.second)}); }
}
static Reduced chain_reduce(const BlockComplex& B) {
    Reduced R;
    if (MIN_TIER <= 0) try { chain_reduce_into<ll>(B, R); R.tier = 0; return R; } catch (OverflowError&) {}
    if (MIN_TIER <= 1) try { chain_reduce_into<i128>(B, R); R.tier = 1; return R; } catch (OverflowError&) {}
    chain_reduce_into<mpz_class>(B, R); R.tier = 2; return R;
}
static bool mpz_to_i128(const mpz_class& x, i128& out) {
    if (mpz_sizeinbase(x.get_mpz_t(), 2) > 125) return false;
    mpz_class a = abs(x), hi = a >> 64, lo = a - (hi << 64);
    unsigned __int128 u = ((unsigned __int128)mpz_get_ui(hi.get_mpz_t()) << 64) | (unsigned __int128)mpz_get_ui(lo.get_mpz_t());
    out = x < 0 ? -(i128)u : (i128)u; return true;
}
// one boundary matrix: exact SNF, int64 -> int128 -> GMP
struct SnfTask { int blk, n, ncols; bool big; vector<vector<pair<int, ll>>> rl; vector<vector<pair<int, mpz_class>>> rz; SNFResult R; int tier = 0; };
static void run_snf_task(SnfTask& T) {
    if (!T.big && MIN_TIER <= 0) { try { T.R = snf_exact<ll>(T.rl, T.ncols); T.tier = 0; return; } catch (OverflowError&) {} }
    if (MIN_TIER <= 1) try {
        vector<vector<pair<int, i128>>> r(T.big ? T.rz.size() : T.rl.size());
        bool ok = true;
        for (size_t i = 0; i < r.size() && ok; i++) {
            if (!T.big) for (auto& e : T.rl[i]) r[i].push_back({e.first, (i128)e.second});
            else for (auto& e : T.rz[i]) { i128 v; if (!mpz_to_i128(e.second, v)) { ok = false; break; } r[i].push_back({e.first, v}); }
        }
        if (ok) { T.R = snf_exact<i128>(r, T.ncols); T.tier = 1; return; }
    } catch (OverflowError&) {}
    vector<vector<pair<int, mpz_class>>> r(T.big ? T.rz.size() : T.rl.size());
    for (size_t i = 0; i < r.size(); i++) {
        if (!T.big) for (auto& e : T.rl[i]) r[i].push_back({e.first, to_mpz(e.second)});
        else r[i] = T.rz[i];
    }
    T.R = snf_exact<mpz_class>(r, T.ncols); T.tier = 2;
}

// reduce every block, then the exact SNF of every boundary matrix of every block, in parallel
static vector<BlockResult> process_blocks(const vector<const BlockComplex*>& BS) {
    int nb = BS.size();
    vector<BlockResult> res(nb);
    vector<Reduced> red(nb);
    double t0 = now_sec();
    parallel_for(nb, [&](int b) { red[b] = chain_reduce(*BS[b]); });
    T_MORSE += now_sec() - t0;
    double t1 = now_sec();
    vector<SnfTask> tasks;
    vector<map<int, vector<int>>> bydeg(nb);
    for (int b = 0; b < nb; b++) {
        const BlockComplex& B = *BS[b];
        BlockResult& R = res[b];
        for (int c = 0; c < (int)B.deg.size(); c++) R.deg[B.deg[c]].cells++;
        for (auto& [n, p] : red[b].piv) R.deg[n].pivots = p;
        R.maxbits_morse = red[b].maxbits; R.used_gmp = (red[b].tier == 2);
        for (int c : red[b].alive) bydeg[b][B.deg[c]].push_back(c);
        vector<int> loc(B.deg.size(), -1);
        for (auto& [n, cs] : bydeg[b]) for (size_t i = 0; i < cs.size(); i++) loc[cs[i]] = i;
        for (auto& [n, cs] : bydeg[b]) {
            R.deg[n].cells_red = cs.size();
            SnfTask T; T.blk = b; T.n = n; T.ncols = cs.size(); T.big = red[b].big;
            int nrows = bydeg[b].count(n - 1) ? bydeg[b][n - 1].size() : 0;
            if (!T.big) { T.rl.resize(nrows); for (size_t j = 0; j < cs.size(); j++) for (auto& fe : red[b].bl[cs[j]]) T.rl[loc[fe.first]].push_back({(int)j, fe.second}); }
            else { T.rz.resize(nrows); for (size_t j = 0; j < cs.size(); j++) for (auto& fe : red[b].bz[cs[j]]) T.rz[loc[fe.first]].push_back({(int)j, fe.second}); }
            tasks.push_back(std::move(T));
        }
    }
    red.clear();
    // largest matrices first
    vector<int> order(tasks.size()); iota(order.begin(), order.end(), 0);
    auto sz = [&](const SnfTask& T) { size_t z = 0; for (auto& r : T.rl) z += r.size(); for (auto& r : T.rz) z += r.size(); return z; };
    vector<size_t> tsz(tasks.size()); for (size_t i = 0; i < tasks.size(); i++) tsz[i] = sz(tasks[i]);
    sort(order.begin(), order.end(), [&](int a, int b) { return tsz[a] > tsz[b]; });
    parallel_for(order.size(), [&](int i) { run_snf_task(tasks[order[i]]); });
    vector<set<ll>> S(nb);
    for (auto& T : tasks) {
        BlockResult& R = res[T.blk];
        R.maxbits_snf = max(R.maxbits_snf, T.R.maxbits);
        R.used_gmp |= (T.tier == 2);
        for (ll p : T.R.S) S[T.blk].insert(p);
        R.deg[T.n].rank = T.R.rank;
        R.deg[T.n].divs = divisors_from(T.R);
    }
    for (int b = 0; b < nb; b++) res[b].S.assign(S[b].begin(), S[b].end());
    T_SNF += now_sec() - t1;
    return res;
}
static BlockResult reduce_and_snf(const BlockComplex& B) { return process_blocks({&B})[0]; }

struct HomologyResult {                              // Borel-Moore homology, by degree
    map<int, long> free;
    map<int, vector<mpz_class>> tors;                // invariant factors > 1 (normalised)
    map<int, map<int, pair<long, vector<mpz_class>>>> byQ;   // Q -> n -> (free, divisors>1)
    long cells = 0, cells_red = 0;
    int maxbits_morse = 0, maxbits_snf = 0;
    set<ll> S; bool used_gmp = false;
    double t_assemble = 0, t_morse_snf = 0;
    map<int, long> ncells_by_deg;
    vector<pair<int, BlockResult>> blocks;           // per-Q block data (for the verification)
};

static vector<mpz_class> normalize_invariants(vector<mpz_class> xs) {   // SNF of a diagonal matrix
    vector<mpz_class> v; for (auto& x : xs) if (x > 1) v.push_back(x);
    // collect prime powers, then rebuild the divisibility chain
    map<ll, vector<int>> pv;
    for (auto& x : v) {
        mpz_class y = x;
        for (ll p : prime_factors(y)) { int e = 0; while (mpz_divisible_ui_p(y.get_mpz_t(), p)) { mpz_divexact_ui(y.get_mpz_t(), y.get_mpz_t(), p); e++; } pv[p].push_back(e); }
    }
    size_t L = 0; for (auto& [p, es] : pv) L = max(L, es.size());
    vector<mpz_class> out(L, 1);
    for (auto& [p, es] : pv) {
        sort(es.begin(), es.end());
        for (size_t t = 0; t < es.size(); t++) { mpz_class q; mpz_ui_pow_ui(q.get_mpz_t(), p, es[t]); out[L - es.size() + t] *= q; }
    }
    return out;
}

static void accumulate_block(HomologyResult& H, int Q, const BlockResult& br) {
    for (auto& [n, dd] : br.deg) {
        H.cells += dd.cells; H.cells_red += dd.cells_red; H.ncells_by_deg[n] += dd.cells;
        long r_n = dd.rank, r_n1 = br.deg.count(n + 1) ? br.deg.at(n + 1).rank : 0;
        long fr = dd.cells_red - r_n - r_n1;
        H.free[n] += fr;
        vector<mpz_class> t;
        if (br.deg.count(n + 1)) for (auto& x : br.deg.at(n + 1).divs) if (x > 1) t.push_back(x);
        auto& slot = H.byQ[Q][n]; slot.first = fr; slot.second = t;
        auto& T = H.tors[n]; T.insert(T.end(), t.begin(), t.end());
    }
    H.maxbits_morse = max(H.maxbits_morse, br.maxbits_morse);
    H.maxbits_snf = max(H.maxbits_snf, br.maxbits_snf);
    for (ll p : br.S) H.S.insert(p);
    H.used_gmp |= br.used_gmp;
    H.blocks.push_back({Q, br});
}
static void finish(HomologyResult& H) {
    for (auto& [n, t] : H.tors) t = normalize_invariants(t);
    for (auto it = H.free.begin(); it != H.free.end();) {
        if (it->second == 0 && (!H.tors.count(it->first) || H.tors[it->first].empty())) { H.tors.erase(it->first); it = H.free.erase(it); }
        else ++it;
    }
    for (auto it = H.tors.begin(); it != H.tors.end();) { if (it->second.empty()) it = H.tors.erase(it); else ++it; }
}

// =============================================================================================
// Section 6.  The new pipeline: M (x)_{A^e} R_D in weight k, split by Q
// =============================================================================================
static Resolution RES;

// cells: (A-basis element m, generator g) with wt(m)+wt(g) = k; punctured: m unpinned
// Q(cell) = floor(m/2) + pins(m) + Q(g); every face keeps Q, so the complex splits over Q.
static void cells_by_Q(const Resolution& R, int k, bool pu, map<int, vector<pair<int, int>>>& out) {
    for (int w0 = 0; w0 <= k; w0++) {
        vector<int> ms{A_id(w0, 0)};
        if (w0 >= 1 && !pu) ms.push_back(A_id(w0 - 1, 1));
        for (int aid : ms) for (int g : R.byWeight[k - w0])
            out[A_m(aid) / 2 + A_v(aid) + R.gens[g].Q].push_back({aid, g});
    }
}
// the block of one Q, as a based complex over Z
static void assemble_block(const Resolution& R, bool pu, int Q, const vector<pair<int, int>>& cells, BlockComplex& B) {
    unordered_map<ull, int> index;
    index.reserve(cells.size() * 2);
    for (size_t i = 0; i < cells.size(); i++) index[((ull)(unsigned)cells[i].second << 12) | (ull)cells[i].first] = i;
    auto idx = [&](int aid, int g) {
        auto it = index.find(((ull)(unsigned)g << 12) | (ull)aid);
        if (it == index.end()) throw runtime_error("internal: face outside its Q-block");
        return it->second;
    };
    B.deg.resize(cells.size()); B.bnd.assign(cells.size(), {});
    for (size_t i = 0; i < cells.size(); i++) B.deg[i] = A_m(cells[i].first) + R.gens[cells[i].second].deg;
    vector<pair<int, ll>> out;
    for (size_t i = 0; i < cells.size(); i++) {
        int aid = cells[i].first, g = cells[i].second;
        out.clear();
        Coef dm = DA(aid, R.poly);
        if (dm.c && !(pu && A_v(dm.id))) out.push_back({idx(dm.id, g), dm.c});
        ll sg = (A_m(aid) & 1) ? -1 : 1;
        for (const Elem* src : {&R.d0[g], &R.DR[g]})
            for (const Term& t : *src) {
                Coef c1 = mulA(aid, t.a, R.poly); if (!c1.c) continue;
                Coef c2 = mulA(c1.id, t.b, R.poly); if (!c2.c) continue;
                if (pu && A_v(c2.id)) continue;
                if (A_m(c2.id) / 2 + A_v(c2.id) + R.gens[t.g].Q != Q) throw runtime_error("internal: Q-grading violated");
                out.push_back({idx(c2.id, t.g), zmul(zmul(sg, t.c), zmul(c1.c, c2.c))});
            }
        sort(out.begin(), out.end());
        auto& b = B.bnd[i];
        for (auto& e : out) {
            if (!b.empty() && b.back().first == e.first) b.back().second = zadd(b.back().second, e.second);
            else b.push_back(e);
        }
        b.erase(remove_if(b.begin(), b.end(), [](const pair<int, ll>& e) { return e.second == 0; }), b.end());
    }
}

// the S-unit / p-local pipeline of Section 5 (exact, self-certifying; the reference method)
static HomologyResult homology_small(int k, bool pu, bool keep_blocks = false, map<int, BlockComplex>* blocks_out = nullptr,
                                     Resolution* Rp = nullptr) {
    Resolution& R = Rp ? *Rp : RES;
    R.build(k);
    HomologyResult H;
    double t0 = now_sec();
    map<int, vector<pair<int, int>>> byQ;
    cells_by_Q(R, k, pu, byQ);
    map<int, BlockComplex> blocks;
    for (auto& [Q, cl] : byQ) assemble_block(R, pu, Q, cl, blocks[Q]);
    H.t_assemble = now_sec() - t0;
    double t1 = now_sec();
    vector<const BlockComplex*> BS; vector<int> Qs;
    for (auto& [Q, B] : blocks) { BS.push_back(&B); Qs.push_back(Q); }
    vector<BlockResult> brs = process_blocks(BS);
    for (size_t b = 0; b < BS.size(); b++) accumulate_block(H, Qs[b], brs[b]);
    H.t_morse_snf = now_sec() - t1;
    finish(H);
    if (blocks_out) *blocks_out = std::move(blocks);
    return H;
}

// =============================================================================================
// Section 6b.  The tiny complex over Z[1/(primes <= k/2)], and the p-adic pipeline it drives
// =============================================================================================
// A' -> A, x_2^i |-> i! gamma_i, is a map of dg-algebras; in weights <= k it becomes an isomorphism
// after inverting every prime <= k/2 (then i! is a unit for every i <= k/2).  Hochschild homology
// commutes with that flat base change, so with Lambda_k = Z[1/p : p <= k/2]:
//     H_*(C(A))_k  (x) Lambda_k  =  H_*(C(A'))_k (x) Lambda_k.
// C(A') is computed from the Koszul resolution of Z[x_2] and has O(k^2) cells, so it gives
//   * the exact rational Betti numbers of every Q-block and degree, and
//   * the exact torsion at every prime > k/2,
// with no elimination over Z on the big complex at all.  The primes p <= k/2 are then handled on
// the big complex by chain-level elimination over Z/p^N (bounded arithmetic, no coefficient growth),
// which is exact once the rank over Q is known.
static Resolution RESP = [] { Resolution r; r.poly = true; return r; }();   // the tiny resolution R'_D for A'
static vector<ll> primes_upto(ll n) { vector<ll> v; for (ll q = 2; q <= n; q++) { bool ok = true; for (ll d = 2; d * d <= q; d++) if (q % d == 0) { ok = false; break; } if (ok) v.push_back(q); } return v; }

// ---- chain-level reduction over Z/p^N: cancel unit incidences, then divide by p and repeat ----
// Stage j cancels the pairs whose incidence has p-adic valuation exactly j, i.e. the elementary
// divisors p^j of the boundary matrices.  Everything stays in [0, p^N).
struct ModPReducer {
    ll p, pN; int N;
    int n; vector<int> deg;
    vector<vector<pair<int, ll>>> d;
    vector<vector<int>> up; vector<int> upcnt; vector<char> alive;
    map<pair<int, int>, long> counts;                      // (degree of tau, stage) -> cancellations
    ModPReducer(const vector<int>& deg_, const vector<vector<pair<int, ll>>>& bnd, ll p_, int N_) : p(p_), N(N_) {
        pN = 1; for (int i = 0; i < N; i++) pN *= p;
        n = deg_.size(); deg = deg_; d.resize(n); up.resize(n); upcnt.assign(n, 0); alive.assign(n, 1);
        for (int c = 0; c < n; c++) for (auto& e : bnd[c]) {
            ll x = e.second % pN; if (x < 0) x += pN;
            if (!x) continue;
            d[c].push_back({e.first, x}); up[e.first].push_back(c); upcnt[e.first]++;
        }
    }
    static ll inv_mod(ll a, ll m) {                        // a invertible mod m = p^N
        i128 t = 0, nt = 1, r = m, nr = a % m;
        while (nr) { i128 q = r / nr; i128 x = t - q * nt; t = nt; nt = x; x = r - q * nr; r = nr; nr = x; }
        if (r > 1) throw runtime_error("internal: non-invertible pivot mod p^N");
        if (t < 0) t += m;
        return (ll)t;
    }
    const ll* coef(int t, int s) const {
        auto& v = d[t];
        auto it = lower_bound(v.begin(), v.end(), s, [](const pair<int, ll>& a, int b) { return a.first < b; });
        return (it != v.end() && it->first == s) ? &it->second : nullptr;
    }
    typedef tuple<ll, int, int> HE;
    priority_queue<HE, vector<HE>, greater<HE>> heap;
    void push(int t) {
        if (!alive[t]) return;
        ll best = LLONG_MAX; int bs = -1, lt = (int)d[t].size() - 1;
        for (auto& fc : d[t]) if (fc.second % p) {
            ll cost = (ll)lt * (ll)(upcnt[fc.first] - 1);
            if (cost < best) { best = cost; bs = fc.first; }
        }
        if (bs >= 0) heap.push({best, t, bs});
    }
    // one stage: cancel every unit incidence.  Returns the number of cancellations.
    long stage(int j) {
        long done = 0;
        heap = decltype(heap)();
        for (int t = 0; t < n; t++) push(t);
        vector<int> mark(n, 0); int stamp = 0;
        while (!heap.empty()) {
            auto [k0, t, s] = heap.top(); heap.pop();
            if (!alive[t]) continue;
            const ll* cp = (alive[s] ? coef(t, s) : nullptr);
            if (!cp || *cp % p == 0) { push(t); continue; }
            ll cost = ((ll)d[t].size() - 1) * (ll)(upcnt[s] - 1);
            if (cost > k0) { push(t); continue; }
            ll epsinv = inv_mod(*cp, pN);
            vector<pair<int, ll>> dt = d[t];
            vector<int> touched; stamp++;
            for (int tp : up[s]) {
                if (tp == t || !alive[tp] || mark[tp] == stamp) continue;
                mark[tp] = stamp;
                const ll* c = coef(tp, s);
                if (!c) continue;
                touched.push_back(tp);
                ll lam = (ll)((i128)(*c) * epsinv % pN);
                auto& old = d[tp];
                vector<pair<int, ll>> nw; nw.reserve(old.size() + dt.size());
                size_t i = 0, jj = 0;
                while (i < old.size() || jj < dt.size()) {
                    if (jj == dt.size() || (i < old.size() && old[i].first < dt[jj].first)) { nw.push_back(old[i]); i++; }
                    else {
                        ll base = (i < old.size() && old[i].first == dt[jj].first) ? old[i].second : 0;
                        ll v = (ll)(((i128)base - (i128)lam * dt[jj].second) % pN); if (v < 0) v += pN;
                        bool had = (i < old.size() && old[i].first == dt[jj].first);
                        if (v) { nw.push_back({dt[jj].first, v}); if (!had) { up[dt[jj].first].push_back(tp); upcnt[dt[jj].first]++; } }
                        else if (had) upcnt[old[i].first]--;
                        if (had) i++;
                        jj++;
                    }
                }
                old.swap(nw);
            }
            for (int y : up[t]) if (alive[y]) {
                auto& v = d[y];
                auto it = lower_bound(v.begin(), v.end(), t, [](const pair<int, ll>& a, int b) { return a.first < b; });
                if (it != v.end() && it->first == t) v.erase(it);
            }
            for (auto& fc : d[t]) upcnt[fc.first]--;
            d[t].clear(); d[t].shrink_to_fit(); alive[t] = 0;
            for (auto& fc : d[s]) upcnt[fc.first]--;
            d[s].clear(); d[s].shrink_to_fit(); alive[s] = 0;
            up[t].clear(); up[t].shrink_to_fit(); up[s].clear(); up[s].shrink_to_fit();
            counts[{deg[t], j}]++; done++;
            for (int tp : touched) push(tp);
        }
        return done;
    }
    // run stages until every boundary matrix has reached its known rank over Q
    void run(const map<int, long>& rank) {
        for (int j = 0; j < N; j++) {
            stage(j);
            bool complete = true;
            for (auto& [nn, r] : rank) {
                long got = 0;
                for (auto& [key, c] : counts) if (key.first == nn) got += c;
                if (got != r) { complete = false; break; }
            }
            if (complete) {
                for (int c = 0; c < n; c++) if (alive[c] && !d[c].empty())
                    throw runtime_error("internal: p-adic reduction: residual differential after full rank");
                return;
            }
            if (j + 1 >= N) throw runtime_error("p-adic reduction: p^N precision exhausted");
            // every remaining entry must now be divisible by p; divide and continue one valuation up
            pN /= p;
            for (int c = 0; c < n; c++) if (alive[c]) for (auto& e : d[c]) {
                if (e.second % p) throw runtime_error("internal: p-adic reduction: unit left after a stage");
                e.second /= p;
                if (e.second >= pN) e.second %= pN;
            }
            for (int c = 0; c < n; c++) if (alive[c]) {
                auto& v = d[c];
                for (auto& e : v) if (!e.second) upcnt[e.first]--;
                v.erase(remove_if(v.begin(), v.end(), [](const pair<int, ll>& e) { return e.second == 0; }), v.end());
            }
        }
        throw runtime_error("p-adic reduction: p^N precision exhausted");
    }
};

// ---- the p-adic pipeline: exact, and with bounded arithmetic on the big complex ----
static long D2_CHECKED = 0;
static HomologyResult homology_padic(int k, bool pu, bool verbose = false, map<int, BlockComplex>* blocks_out = nullptr, bool check_d2 = false) {
    if (!RESP.poly) throw runtime_error("internal: RESP must be the polynomial-model resolution");
    RES.build(k); RESP.build(k);
    HomologyResult H;
    double t0 = now_sec();
    // (1) the tiny complex: rational Betti numbers per (Q, degree) and the torsion at primes > k/2
    vector<ll> P = primes_upto(k / 2);                       // the primes that A' -> A inverts
    vector<ll> saveS = SNF_S0; bool saveO = SNF_ONLY_NEW;
    SNF_S0 = P; SNF_ONLY_NEW = true;                         // the tiny complex is computed over Lambda_k
    if (P.empty()) SNF_ONLY_NEW = false;                     // k <= 3: Lambda_k = Z, the tiny complex gives everything
    HomologyResult T = homology_small(k, pu, false, nullptr, &RESP);
    SNF_S0 = saveS; SNF_ONLY_NEW = saveO;
    double t_tiny = now_sec() - t0; T_SNF += t_tiny;
    // (2) the big complex, one Q-block at a time
    double t1 = now_sec();
    map<int, vector<pair<int, int>>> byQ;
    cells_by_Q(RES, k, pu, byQ);
    vector<int> order;
    for (auto& [Q, cl] : byQ) order.push_back(Q);
    sort(order.begin(), order.end(), [&](int a, int b) { return byQ[a].size() > byQ[b].size(); });
    double t_asm = 0, t_morse = 0, t_padic = 0;
    for (int Q : order) {
        double ta = now_sec();
        BlockComplex B;
        assemble_block(RES, pu, Q, byQ[Q], B);
        vector<pair<int, int>>().swap(byQ[Q]);
        if (blocks_out) (*blocks_out)[Q] = B;
        if (check_d2) {                                       // d^2 = 0 on this block, before any reduction
            parallel_for(B.deg.size(), [&](int i) {
                map<int, ll> acc;
                for (auto& [f, c] : B.bnd[i]) for (auto& [g, c2] : B.bnd[f]) acc[g] = zadd(acc[g], zmul(c, c2));
                for (auto& [g, x] : acc) if (x) throw runtime_error("d^2 != 0 on the small complex");
            });
            D2_CHECKED += B.deg.size();
        }
        t_asm += now_sec() - ta;
        BlockResult R;
        for (int c = 0; c < (int)B.deg.size(); c++) R.deg[B.deg[c]].cells++;
        double tb = now_sec();
        Reduced red = chain_reduce(B);                        // exact unit-pivot reduction over Z
        vector<vector<pair<int, ll>>>().swap(B.bnd);
        R.maxbits_morse = red.maxbits; R.used_gmp = (red.tier == 2);
        for (auto& [n, c] : red.piv) R.deg[n].pivots = c;
        t_morse += now_sec() - tb; T_MORSE += now_sec() - tb;
        // renumber the survivors, degree by degree
        map<int, vector<int>> bydeg;
        for (int c : red.alive) bydeg[B.deg[c]].push_back(c);
        vector<int> loc(B.deg.size(), -1), rdeg;
        vector<vector<pair<int, ll>>> rb;
        for (auto& [n, cs] : bydeg) for (int c : cs) { loc[c] = rdeg.size(); rdeg.push_back(n); rb.push_back({}); }
        for (int c : red.alive) {
            auto& dst = rb[loc[c]];
            if (red.big) { for (auto& e : red.bz[c]) { if (!fits_ll(e.second)) throw runtime_error("p-adic pipeline: coefficient too large for the mod p^N path"); dst.push_back({loc[e.first], e.second.get_si()}); } }
            else for (auto& e : red.bl[c]) dst.push_back({loc[e.first], e.second});
            sort(dst.begin(), dst.end());
        }
        red.bl.clear(); red.bl.shrink_to_fit(); red.bz.clear(); red.bz.shrink_to_fit();
        // (3) ranks over Q from the tiny complex:  c_n - r_n - r_{n+1} = b_n
        map<int, long> cnt, rank;
        for (int nn : rdeg) cnt[nn]++;
        for (auto& [n, c] : cnt) R.deg[n].cells_red = c;
        int nmin = cnt.empty() ? 0 : cnt.begin()->first, nmax = cnt.empty() ? -1 : cnt.rbegin()->first;
        long prev = 0;
        for (int nn = nmin; nn <= nmax + 1; nn++) {
            rank[nn] = prev;
            long b = 0;
            if (T.byQ.count(Q) && T.byQ.at(Q).count(nn)) b = T.byQ.at(Q).at(nn).first;
            long c = cnt.count(nn) ? cnt[nn] : 0;
            prev = c - b - prev;                              // r_{n+1}
            if (prev < 0) throw runtime_error("rank recursion went negative: the tiny complex and the big one disagree");
        }
        if (prev != 0) throw runtime_error("rank recursion did not close at the top degree");
        for (auto& [n, r] : rank) if (cnt.count(n)) R.deg[n].rank = r;
        // (4) p-adic chain reduction for every prime <= k/2, in parallel
        vector<map<pair<int, int>, long>> cts(P.size());
        double tc = now_sec();
        parallel_for(P.size(), [&](int i) {
            ll q = P[i];
            int N = 0; { ll x = 1; while (x <= (1ll << 62) / q) { x *= q; N++; } }
            map<int, long> tgt;
            for (auto& [n, r] : rank) if (r) tgt[n] = r;
            ModPReducer M(rdeg, rb, q, N);
            M.run(tgt);
            cts[i] = M.counts;
        });
        t_padic += now_sec() - tc; T_SNF += now_sec() - tc;
        for (size_t i = 0; i < P.size(); i++) for (auto& [key, c] : cts[i]) if (key.second > 0) {
            mpz_class pj; mpz_ui_pow_ui(pj.get_mpz_t(), P[i], key.second);
            for (long t = 0; t < c; t++) R.deg[key.first].divs.push_back(pj);
        }
        // (5) the torsion at primes > k/2, from the tiny complex (H_n there feeds d_{n+1} here)
        if (T.byQ.count(Q)) for (auto& [n, fr_t] : T.byQ.at(Q)) for (auto& x : fr_t.second) {
            mpz_class y = x;
            for (ll q : prime_factors(y)) {
                if (find(P.begin(), P.end(), q) != P.end()) continue;
                int e = 0; mpz_class z = x;
                while (mpz_divisible_ui_p(z.get_mpz_t(), q)) { mpz_divexact_ui(z.get_mpz_t(), z.get_mpz_t(), q); e++; }
                mpz_class pe; mpz_ui_pow_ui(pe.get_mpz_t(), q, e);
                R.deg[n + 1].divs.push_back(pe);
            }
        }
        for (ll q : P) R.S.push_back(q);
        accumulate_block(H, Q, R);
        // the free ranks must reproduce the tiny complex exactly
        for (auto& [n, dd] : R.deg) {
            long b = (T.byQ.count(Q) && T.byQ.at(Q).count(n)) ? T.byQ.at(Q).at(n).first : 0;
            long r1 = R.deg.count(n + 1) ? R.deg[n + 1].rank : 0;
            if (dd.cells_red - dd.rank - r1 != b) throw runtime_error("free rank does not match the tiny complex");
        }
    }
    H.t_assemble = t_asm + t_tiny;
    H.t_morse_snf = now_sec() - t1;
    finish(H);
    if (verbose) {
        long mx = 0, mt = 0; for (auto& e : RESP.DR) { mt = max(mt, (long)e.size()); for (auto& t : e) mx = max(mx, (long)llabs(t.c)); }
        fprintf(stderr, "    [tiny] %ld cells, %ld gens, D_R max terms %ld max|c| %ld, morse %d bits, snf %d bits\n",
                T.cells, (long)RESP.gens.size(), mt, mx, T.maxbits_morse, T.maxbits_snf);
    }
    if (verbose) fprintf(stderr, "    [k=%d%s] tiny %.2fs, assemble %.2fs, unit reduction %.2fs, p-adic %.2fs (primes <= %d)\n",
                         k, pu ? " pu" : "", t_tiny, t_asm, t_morse, t_padic, k / 2);
    return H;
}

// =============================================================================================
// Section 7.  The reference: Napolitano's full complex (napolitano.py) and the original algorithm
// =============================================================================================
// A cell (e,v,F), F = ((m_1,v_1),...).  The fibre word F of weight N >= 1 is encoded as a base-3
// string of length N: each letter of weight w and pin v becomes the symbol v (start) followed by
// w-1 symbols 2 (continue); the first symbol is 0 or 1, so there are 2*3^(N-1) words.
struct NapCell { int e, v; vector<pair<int, int>> F; };
struct NapIndex {
    int k; bool pu;
    vector<ll> off;                               // offset of (v,e): index 2*e+v
    vector<ll> pow3;
    ll total = 0;
    NapIndex(int k_, bool pu_) : k(k_), pu(pu_) {
        pow3.assign(k + 2, 1); for (int i = 1; i <= k + 1; i++) pow3[i] = pow3[i - 1] * 3;
        off.assign(2 * k + 4, -1);
        for (int v = 0; v <= (pu ? 0 : 1); v++) for (int e = 0; e + v <= k; e++) {
            int N = k - e - v; off[2 * e + v] = total; total += (N == 0) ? 1 : 2 * pow3[N - 1];
        }
    }
    ll encode(const NapCell& c) const {
        int N = k - c.e - c.v;
        ll x = 0;
        for (auto& [m, vv] : c.F) { int w = m + vv; x = x * 3 + vv; for (int t = 1; t < w; t++) x = x * 3 + 2; }
        (void)N;
        return off[2 * c.e + c.v] + x;
    }
    NapCell decode(ll idx) const {
        NapCell c;
        // find (e,v)
        int best = -1;
        for (int t = 0; t < (int)off.size(); t++) if (off[t] >= 0 && off[t] <= idx && (best < 0 || off[t] > off[best])) best = t;
        c.e = best / 2; c.v = best % 2;
        ll x = idx - off[best];
        int N = k - c.e - c.v;
        vector<int> dig(N);
        for (int t = N - 1; t >= 0; t--) { dig[t] = x % 3; x /= 3; }
        for (int t = 0; t < N; t++) {
            if (dig[t] < 2) c.F.push_back({1 - dig[t], dig[t]});     // new letter: weight 1 -> (1,0) or (0,1)
            else c.F.back().first++;                                   // continue: one more off-pole point
        }
        return c;
    }
    static int degree(const NapCell& c) { int n = c.e; for (auto& f : c.F) n += f.first + 1; return n; }
    static int Qof(const NapCell& c) { int q = c.e / 2 + c.v; for (auto& f : c.F) q += f.first / 2 + f.second; return q; }
};

// napolitano.boundary, verbatim: the five face types with their signed incidence numbers
static vector<pair<NapCell, ll>> nap_boundary(const NapCell& c) {
    int e = c.e, v = c.v, d = c.F.size();
    const auto& F = c.F;
    vector<pair<NapCell, ll>> out;
    auto sgn = [](ll x) { return (x & 1) ? -1 : 1; };
    // (1) two adjacent fibres merge, unless both are pinned
    int pref = 0;                                          // sum of m over F[:i]
    for (int i = 1; i < d; i++) {
        pref += F[i - 1].first;
        auto [mi, vi] = F[i - 1]; auto [mj, vj] = F[i];
        if (vi == 1 && vj == 1) continue;
        ll co = sgn(e + i - 1 + pref) * Ppar(mi, mj);
        if (!co) continue;
        NapCell f{e, v, {}};
        f.F.assign(F.begin(), F.begin() + i - 1); f.F.push_back({mi + mj, max(vi, vj)}); f.F.insert(f.F.end(), F.begin() + i + 1, F.end());
        out.push_back({f, co});
    }
    // (2) an off-pole point of a fibre reaches that fibre's north pole
    pref = 0;
    for (int i = 1; i <= d; i++) {
        auto [mi, vi] = F[i - 1];
        if (vi == 0 && mi >= 1) {
            ll co = sgn(e + i - 1 + pref) * (1 + sgn(mi));
            if (co) { NapCell f = c; f.F[i - 1] = {mi - 1, 1}; out.push_back({f, co}); }
        }
        pref += mi;
    }
    // (3) a base-circle point reaches the base north pole
    if (v == 0 && e >= 1) { ll co = -(1 + sgn(e)); if (co) { NapCell f = c; f.e = e - 1; f.v = 1; out.push_back({f, co}); } }
    // (4) the leftmost fibre slides onto the base circle
    if (d >= 1) {
        auto [m1, v1] = F[0];
        if (!(v == 1 && v1 == 1)) {
            ll co = -sgn(e) * Ppar(e, m1);
            if (co) { NapCell f{e + m1, max(v, v1), vector<pair<int, int>>(F.begin() + 1, F.end())}; out.push_back({f, co}); }
        }
    }
    // (5) the rightmost fibre slides onto the base circle
    if (d >= 1) {
        auto [md, vd] = F[d - 1];
        if (!(v == 1 && vd == 1)) {
            int S = d - 1; for (int i = 0; i < d - 1; i++) S += F[i].first;
            ll co = sgn(e + (ll)S * (1 + md)) * Ppar(e, md);
            if (co) { NapCell f{e + md, max(v, vd), vector<pair<int, int>>(F.begin(), F.end() - 1)}; out.push_back({f, co}); }
        }
    }
    return out;
}

// the reference homology: build each Q-block of the FULL complex, check d^2 = 0 and the Q-grading on
// every cell, reduce by unit pivots (as reduction.morse_reduce), exact SNF.
struct RefStats { long cells = 0; long d2_checked = 0; double t = 0; long rss_kb = 0; };
static HomologyResult homology_reference(int k, bool pu, RefStats* st = nullptr) {
    double t0 = now_sec();
    NapIndex I(k, pu);
    ll T = I.total;
    vector<uint8_t> Qc(T);
    vector<int> deg(T);
    for (ll x = 0; x < T; x++) { NapCell c = I.decode(x); Qc[x] = NapIndex::Qof(c); deg[x] = NapIndex::degree(c); }
    int Qmax = 0; for (ll x = 0; x < T; x++) Qmax = max(Qmax, (int)Qc[x]);
    HomologyResult H;
    vector<int> loc(T, -1);
    long d2 = 0;
    for (int Q = 0; Q <= Qmax; Q++) {
        vector<ll> ids; for (ll x = 0; x < T; x++) if (Qc[x] == Q) ids.push_back(x);
        if (ids.empty()) continue;
        for (size_t i = 0; i < ids.size(); i++) loc[ids[i]] = i;
        BlockComplex B; B.deg.resize(ids.size()); B.bnd.resize(ids.size());
        for (size_t i = 0; i < ids.size(); i++) {
            NapCell c = I.decode(ids[i]); B.deg[i] = deg[ids[i]];
            map<int, ll> acc;
            for (auto& [f, co] : nap_boundary(c)) {
                if (pu && f.v == 1) continue;                    // the punctured complex: quotient by v = 1
                ll fi = I.encode(f);
                if (Qc[fi] != Q) throw runtime_error("reference: Q-grading violated by a face");
                if (deg[fi] != B.deg[i] - 1) throw runtime_error("reference: face of wrong degree");
                acc[loc[fi]] = zadd(acc[loc[fi]], co);
            }
            for (auto& [j, co] : acc) if (co) B.bnd[i].push_back({j, co});
        }
        // d^2 = 0 on every cell of the block (napolitano.check_d2)
        for (size_t i = 0; i < ids.size(); i++) {
            map<int, ll> acc;
            for (auto& [f, co] : B.bnd[i]) for (auto& [g, co2] : B.bnd[f]) acc[g] = zadd(acc[g], zmul(co, co2));
            for (auto& [g, x] : acc) if (x) throw runtime_error("reference: d^2 != 0 on the full complex");
            d2++;
        }
        accumulate_block(H, Q, reduce_and_snf(B));
        for (ll x : ids) loc[x] = -1;
    }
    finish(H);
    if (st) { st->cells = T; st->d2_checked = d2; st->t = now_sec() - t0; st->rss_kb = peak_rss_kb(); }
    return H;
}

// =============================================================================================
// Section 8.  Output
// =============================================================================================
static string divs_str(const vector<mpz_class>& t) {
    map<string, int> cnt; vector<string> order;
    for (auto& x : t) { string s = x.get_str(); if (!cnt.count(s)) order.push_back(s); cnt[s]++; }
    string out;
    for (auto& s : order) out += " + Z/" + s + (cnt[s] > 1 ? "^" + to_string(cnt[s]) : "");
    return out;
}
// ordinary homology via [IR26b, Lemma 6.1]: rank H_i = rank H^BM_{2k-i}, Tors H_i = Tors H^BM_{2k-i-1}
struct Ordinary { map<int, long> free; map<int, map<ll, int>> pcount; map<int, vector<mpz_class>> tors; };
static Ordinary to_ordinary(int k, const HomologyResult& H) {
    Ordinary O;
    for (auto& [n, f] : H.free) if (f) O.free[2 * k - n] = f;
    for (auto& [n, t] : H.tors) {
        int i = 2 * k - n - 1;
        O.tors[i] = t;
        for (auto& x : t) { mpz_class y = x; for (ll p : prime_factors(y)) O.pcount[i][p]++; }
    }
    return O;
}
static void print_result(int k, bool pu, const HomologyResult& H, bool bm) {
    const char* S = pu ? "Sigma_1^o" : "Sigma_1";
    if (bm) {
        printf("H^BM_n(B_%d(%s);Z):\n", k, S);
        set<int> ns; for (auto& [n, f] : H.free) ns.insert(n); for (auto& [n, t] : H.tors) ns.insert(n);
        for (int n : ns) printf("  H^BM_%-2d = Z^%ld%s\n", n, H.free.count(n) ? H.free.at(n) : 0, H.tors.count(n) ? divs_str(H.tors.at(n)).c_str() : "");
    }
    Ordinary O = to_ordinary(k, H);
    printf("H_i(B_%d(%s);Z), exact:\n", k, S);
    set<int> is; for (auto& [i, f] : O.free) is.insert(i); for (auto& [i, t] : O.tors) is.insert(i);
    for (int i : is) {
        string s = "Z^" + to_string(O.free.count(i) ? O.free.at(i) : 0);
        if (O.pcount.count(i)) for (auto& [p, c] : O.pcount.at(i)) s += " + (Z/" + to_string(p) + ")" + (c > 1 ? "^" + to_string(c) : "");
        printf("  H_%-2d = %-60s", i, s.c_str());
        if (O.tors.count(i)) printf("  [inv. factors:%s]", divs_str(O.tors.at(i)).c_str());
        printf("\n");
    }
}

// =============================================================================================
// Section 9.  Independent checks used by --verify
// =============================================================================================
static int N_PASS = 0, N_FAIL = 0;
static vector<string> FAILED;
static void report(const string& name, bool ok, const string& detail = "") {
    (ok ? N_PASS : N_FAIL)++;
    if (!ok) FAILED.push_back(name);
    printf("  [%s] %s%s%s\n", ok ? "PASS" : "FAIL", name.c_str(), detail.empty() ? "" : "  ", detail.c_str());
    fflush(stdout);
}

// ---- rank over F_q of a sparse integer matrix: plain column-by-column Gaussian elimination mod q.
//      Deliberately independent of the unit reduction and of the S-unit / p-local code. ----
static long rank_mod_q(const vector<vector<pair<int, ll>>>& cols, ull q) {
    auto mulm = [&](ull a, ull b) { return (ull)((unsigned __int128)a * b % q); };
    auto powm = [&](ull a, ull e) { ull r = 1; while (e) { if (e & 1) r = mulm(r, a); a = mulm(a, a); e >>= 1; } return r; };
    // pivot rows: pivrow[r] = normalised sparse row with leading entry 1 at column r (rows = cells of degree n-1)
    unordered_map<int, vector<pair<int, ull>>> piv;          // leading row index -> reduced vector (sorted)
    long rank = 0;
    vector<pair<int, ull>> v, w;
    for (auto& c : cols) {
        v.clear();
        for (auto& e : c) { ll x = e.second % (ll)q; if (x < 0) x += q; if (x) v.push_back({e.first, (ull)x}); }
        sort(v.begin(), v.end());
        // reduce v against the pivot vectors, smallest leading index first
        while (!v.empty()) {
            auto it = piv.find(v[0].first);
            if (it == piv.end()) break;
            const auto& pv = it->second; ull lam = v[0].second;
            w.clear(); size_t i = 0, j = 0;
            while (i < v.size() || j < pv.size()) {
                if (j == pv.size() || (i < v.size() && v[i].first < pv[j].first)) { w.push_back(v[i]); i++; }
                else if (i == v.size() || pv[j].first < v[i].first) { w.push_back({pv[j].first, (q - mulm(lam, pv[j].second)) % q}); j++; }
                else { ull x = (v[i].second + q - mulm(lam, pv[j].second)) % q; if (x) w.push_back({v[i].first, x}); i++; j++; }
            }
            v.swap(w);
        }
        if (v.empty()) continue;
        ull inv = powm(v[0].second, q - 2);
        for (auto& e : v) e.second = mulm(e.second, inv);
        piv[v[0].first] = v; rank++;
    }
    return rank;
}

// ---- cone of the augmentation eps: R_w -> A_w, for R_0 (D = 0) or R_D; acyclic <=> quasi-iso ----
static bool cone_acyclic(int w, bool perturbed, long* size = nullptr, Resolution* Rp = nullptr) {
    Resolution& RR = Rp ? *Rp : RES;
    RR.build(w);
    vector<tuple<int, int, int>> rb;                      // (a, b, g) basis of R in weight w
    for (int wl = 0; wl <= w; wl++) for (int w1 = 0; w1 <= wl; w1++) {
        vector<int> as{A_id(w1, 0)}; if (w1) as.push_back(A_id(w1 - 1, 1));
        vector<int> bs{A_id(wl - w1, 0)}; if (wl - w1) bs.push_back(A_id(wl - w1 - 1, 1));
        for (int a : as) for (int b : bs) for (int g : RR.byWeight[w - wl]) rb.push_back({a, b, g});
    }
    vector<int> ab{A_id(w, 0)}; if (w) ab.push_back(A_id(w - 1, 1));
    map<tuple<int, int, int>, int> idx; for (size_t i = 0; i < rb.size(); i++) idx[rb[i]] = i;
    int nR = rb.size();
    BlockComplex B; B.deg.resize(nR + ab.size()); B.bnd.resize(nR + ab.size());
    for (int i = 0; i < nR; i++) { auto [a, b, g] = rb[i]; B.deg[i] = degE({a, b}) + RR.gens[g].deg + 1; }
    for (size_t t = 0; t < ab.size(); t++) B.deg[nR + t] = A_m(ab[t]);
    for (int i = 0; i < nR; i++) {
        auto [a, b, g] = rb[i];
        Elem el{{a, b, g, 1}};
        Elem dd = perturbed ? RR.apply_dT(el) : RR.apply_d0(el);
        map<int, ll> acc;
        for (auto& t : dd) acc[idx.at({t.a, t.b, t.g})] -= t.c;
        for (auto& t : RR.eps_map(el)) { int j = nR + (t.a == ab[0] ? 0 : 1); acc[j] += t.c; }
        for (auto& [j, c] : acc) if (c) B.bnd[i].push_back({j, c});
    }
    if (perturbed) for (size_t t = 0; t < ab.size(); t++) { Coef c = DA(ab[t], RR.poly); if (c.c) B.bnd[nR + t].push_back({nR + (c.id == ab[0] ? 0 : 1), c.c}); }
    for (size_t i = 0; i < B.deg.size(); i++) {             // d^2 = 0 on the cone
        map<int, ll> acc; for (auto& [f, c] : B.bnd[i]) for (auto& [g, c2] : B.bnd[f]) acc[g] += c * c2;
        for (auto& [g, x] : acc) if (x) return false;
    }
    if (size) *size = B.deg.size();
    BlockResult br = reduce_and_snf(B);
    for (auto& [n, dd] : br.deg) {
        long r1 = br.deg.count(n + 1) ? br.deg[n + 1].rank : 0;
        if (dd.cells_red - dd.rank - r1 != 0) return false;
        for (auto& x : dd.divs) if (x != 1) return false;
    }
    return true;
}

// ---- the standard two-sided-bar Hochschild differential on A (x)_{A^e} B(A,A,A), cell basis ----
static map<ll, ll> std_hochschild(const NapCell& c, int sign_ext, const NapIndex& I) {
    int e = c.e, v = c.v, n = c.F.size();
    int m = A_id(e, v);
    map<ll, ll> out;
    auto put = [&](int aid, const vector<pair<int, int>>& F, ll co) { if (!co) return; NapCell f{A_m(aid), A_v(aid), F}; out[I.encode(f)] += co; };
    auto lid = [&](pair<int, int> f) { return A_id(f.first, f.second); };
    auto lpair = [&](int id) { return make_pair(A_m(id), A_v(id)); };
    Coef dm = DA(m); if (dm.c) put(dm.id, c.F, dm.c);
    vector<int> eps(n + 1); eps[0] = e; for (int i = 1; i <= n; i++) eps[i] = eps[i - 1] + c.F[i - 1].first + 1;
    ll sg = (e & 1) ? -1 : 1;
    for (int i = 1; i <= n; i++) {                            // d(s a) = -s(Da)
        Coef da = DA(lid(c.F[i - 1])); if (!da.c) continue;
        auto F = c.F; F[i - 1] = lpair(da.id);
        put(m, F, ((eps[i - 1] & 1) ? 1 : -1) * da.c);
    }
    if (n) {
        Coef c1 = mulA(m, lid(c.F[0]));
        if (c1.c) put(c1.id, vector<pair<int, int>>(c.F.begin() + 1, c.F.end()), sign_ext * sg * c1.c);
        for (int i = 1; i < n; i++) {
            Coef cm = mulA(lid(c.F[i - 1]), lid(c.F[i])); if (!cm.c) continue;
            vector<pair<int, int>> F(c.F.begin(), c.F.begin() + i - 1); F.push_back(lpair(cm.id)); F.insert(F.end(), c.F.begin() + i + 1, c.F.end());
            put(m, F, sign_ext * sg * (((eps[i] - e) & 1) ? -1 : 1) * cm.c);
        }
        int w = eps[n - 1] - e; int an = lid(c.F[n - 1]);
        Coef cr = mulA(m, an);
        if (cr.c) put(cr.id, vector<pair<int, int>>(c.F.begin(), c.F.end() - 1), sign_ext * sg * ((w & 1) ? 1 : -1) * (((ll)A_m(an) * w & 1) ? -1 : 1) * cr.c);
    }
    for (auto it = out.begin(); it != out.end();) if (!it->second) it = out.erase(it); else ++it;
    return out;
}

// ---- phi: A' -> A,  x_2^i x_1^e y^v |-> i! gamma_i x_1^e y^v, is a map of dg-algebras, and an
//      isomorphism once the primes <= W/2 are inverted (i! is then a unit).  Checked on the basis. ----
static bool phi_is_dga_map(int W, string& why) {
    auto fac = [](int i) { mpz_class f = 1; for (int t = 2; t <= i; t++) f *= t; return f; };
    auto phi = [&](int id) { return fac(A_m(id) / 2); };                    // the scalar of phi on a(m,v)
    for (int m = 0; m <= W; m++) for (int v = 0; v <= 1; v++) {
        if (m + v > W) continue;
        int x = A_id(m, v);
        Coef dp = DA(x, true), dg = DA(x, false);                           // phi(D' a) = D(phi a)
        if ((dp.c == 0) != (dg.c == 0)) { why = "D support"; return false; }
        if (dp.c) {
            if (dp.id != dg.id) { why = "D target"; return false; }
            if (to_mpz(dp.c) * phi(dp.id) != phi(x) * to_mpz(dg.c)) { why = "D"; return false; }
        }
        for (int m2 = 0; m2 <= W - m - v; m2++) for (int v2 = 0; v2 <= 1; v2++) {
            if (m + v + m2 + v2 > W) continue;
            int y = A_id(m2, v2);
            Coef cp = mulA(x, y, true), cg = mulA(x, y, false);             // phi(ab) = phi(a) phi(b)
            if ((cp.c == 0) != (cg.c == 0)) { why = "product support"; return false; }
            if (!cp.c) continue;
            if (phi(x) * phi(y) * to_mpz(cg.c) != phi(cp.id) * to_mpz(cp.c)) { why = "product"; return false; }
        }
    }
    return true;
}

// ---- Knudsen's rational Chevalley-Eilenberg model (port of ce_model.py) ----
struct CESurface { vector<int> th; map<pair<int, int>, int> cup; };
static CESurface CE_CLOSED() { CESurface S; S.th = {0, 1, 1, 2}; for (int i = 0; i < 4; i++) { S.cup[{0, i}] = i; S.cup[{i, 0}] = i; } S.cup[{1, 2}] = 3; S.cup[{2, 1}] = 3; return S; }
static CESurface CE_PUNCT() { CESurface S; S.th = {1, 1, 2}; S.cup[{0, 1}] = 2; S.cup[{1, 0}] = 2; return S; }
static map<int, map<int, long>> ce_betti(const CESurface& S, int K) {
    struct G { int deg, wt; bool w; int th; };
    vector<G> gens; int T = S.th.size();
    for (int i = 0; i < T; i++) gens.push_back({2 - S.th[i], 1, false, i});
    for (int i = 0; i < T; i++) gens.push_back({3 - S.th[i], 2, true, i});
    auto widx = [&](int th) { return T + th; };
    map<int, map<int, long>> R; R[0][0] = 1;
    for (int k = 1; k <= K; k++) {
        vector<vector<int>> mons; vector<int> cur;
        function<void(int, int)> rec = [&](int j, int w) {
            if (w == 0) { mons.push_back(cur); return; }
            if (j == (int)gens.size()) return;
            rec(j + 1, w);
            int top = (gens[j].deg & 1) ? 1 : w / gens[j].wt;
            for (int r = 1; r <= top && r * gens[j].wt <= w; r++) { for (int t = 0; t < r; t++) cur.push_back(j); rec(j + 1, w - r * gens[j].wt); for (int t = 0; t < r; t++) cur.pop_back(); }
        };
        rec(0, k);
        auto dg = [&](const vector<int>& m) { int d = 0; for (int j : m) d += gens[j].deg; return d; };
        map<int, vector<vector<int>>> byd; for (auto& m : mons) { auto s = m; sort(s.begin(), s.end()); byd[dg(s)].push_back(s); }
        map<int, long> rank;
        for (auto& [d, src] : byd) {
            map<vector<int>, int> tidx; if (byd.count(d - 1)) for (size_t i = 0; i < byd[d - 1].size(); i++) tidx[byd[d - 1][i]] = i;
            vector<vector<pair<int, ll>>> rows(tidx.size());
            for (size_t j = 0; j < src.size(); j++) {
                const auto& lst = src[j];
                map<vector<int>, ll> out;
                for (size_t a = 0; a < lst.size(); a++) for (size_t b = a + 1; b < lst.size(); b++) {
                    int ia = lst[a], ib = lst[b];
                    if (gens[ia].w || gens[ib].w) continue;
                    auto it = S.cup.find({gens[ia].th, gens[ib].th}); if (it == S.cup.end()) continue;
                    ll sgn = 1;
                    for (size_t t = 0; t < a; t++) if ((gens[lst[t]].deg * gens[ia].deg) & 1) sgn = -sgn;
                    for (size_t t = 0; t < b; t++) if (t != a && ((gens[lst[t]].deg * gens[ib].deg) & 1)) sgn = -sgn;
                    vector<int> arr{widx(it->second)};
                    for (size_t t = 0; t < lst.size(); t++) if (t != a && t != b) arr.push_back(lst[t]);
                    for (size_t r = 0; r < arr.size(); r++) for (size_t q = 0; q + 1 < arr.size(); q++)
                        if (arr[q] > arr[q + 1]) { if ((gens[arr[q]].deg & 1) && (gens[arr[q + 1]].deg & 1)) sgn = -sgn; swap(arr[q], arr[q + 1]); }
                    bool rep = false; for (size_t t = 0; t + 1 < arr.size(); t++) if (arr[t] == arr[t + 1] && (gens[arr[t]].deg & 1)) rep = true;
                    if (rep) continue;
                    out[arr] += sgn;
                }
                for (auto& [m, c] : out) if (c) rows[tidx.at(m)].push_back({(int)j, c});
            }
            for (auto& r : rows) sort(r.begin(), r.end());
            vector<ll> Sp{2}; int mb = 0;
            rank[d] = rows.empty() ? 0 : rank_and_primes<ll>(rows, src.size(), Sp, mb);
        }
        for (auto& [d, v] : byd) { long b = v.size() - rank[d] - (rank.count(d + 1) ? rank[d + 1] : 0); if (b) R[k][d] = b; }
    }
    return R;
}

// ---- Bianchi-Stavrou's closed formula for H_*(B_k(Sigma_{1,1});F_p) [BS24, Cor. C.6] (port of bianchi_stavrou.py) ----
static map<int, map<int, long>> bs_mod_p_betti(int p, int K, const map<int, map<int, long>>& rat) {
    vector<tuple<int, int, bool>> gens;
    for (ll q = p; 2 * q <= K; q *= p) gens.push_back({2 * q, 2 * q - 1, true});          // alpha_i, exterior
    for (ll q = p; 2 * q <= K; q *= p) gens.push_back({2 * q, 2 * q - 2, false});         // beta_{i-1}
    map<pair<int, int>, long> A; A[{0, 0}] = 1;
    for (auto [w, dg, ext] : gens) {
        map<pair<int, int>, long> B;
        for (auto& [wd, c] : A) for (int r = 0; wd.first + r * w <= K; r++) { B[{wd.first + r * w, wd.second + r * dg}] += c; if (ext && r >= 1) break; }
        A = B;
    }
    map<int, map<int, long>> out;
    for (int k = 0; k <= K; k++) for (auto& [wd, c] : A) if (wd.first <= k && rat.count(k - wd.first))
        for (auto& [i, b] : rat.at(k - wd.first)) out[k][i + wd.second] += c * b;
    return out;
}

// ---- number of cells of Napolitano's complex by degree (no enumeration), for the Euler characteristic ----
static map<int, ll> nap_cell_counts(int k, bool pu) {
    vector<map<int, ll>> f(k + 1); f[0][0] = 1;          // fibre words: weight w, degree n
    for (int w = 1; w <= k; w++) for (int wl = 1; wl <= w; wl++) for (auto& [n, c] : f[w - wl]) { f[w][n + wl + 1] += c; f[w][n + wl] += c; }
    map<int, ll> out;
    for (int v = 0; v <= (pu ? 0 : 1); v++) for (int e = 0; e + v <= k; e++) for (auto& [n, c] : f[k - e - v]) out[n + e] += c;
    return out;
}

// ---- published / known values (external data, used only for comparison) ----
// [IR26b]/write-up Table for k = 14 (closed torus), Borel-Moore degrees: n -> (free, {divisor: count})
static const vector<tuple<int, long, vector<pair<int, int>>>> PUB_K14_BM = {
    {13, 6, {{2, 8}}}, {14, 19, {{2, 26}, {6, 2}, {30, 1}}}, {15, 25, {{2, 48}, {6, 6}, {30, 3}, {210, 1}}},
    {16, 23, {{2, 62}, {6, 8}, {30, 5}}}, {17, 21, {{2, 65}, {6, 9}, {30, 3}}}, {18, 19, {{2, 58}, {6, 9}, {30, 2}}},
    {19, 17, {{2, 47}, {6, 7}, {30, 1}}}, {20, 15, {{2, 36}, {6, 5}}}, {21, 13, {{2, 25}, {6, 3}}}, {22, 11, {{2, 16}, {6, 2}}},
    {23, 9, {{2, 10}, {6, 1}}}, {24, 7, {{2, 6}}}, {25, 5, {{2, 3}}}, {26, 3, {{2, 1}}}, {27, 2, {}}, {28, 1, {}}};
// write-up, weight 15 in full (closed torus), ordinary degrees: i -> (free, {p: count})
static const vector<tuple<int, long, vector<pair<int, int>>>> PUB_K15_ORD = {
    {0, 1, {}}, {1, 2, {{2, 1}}}, {2, 3, {{2, 3}}}, {3, 5, {{2, 6}}}, {4, 7, {{2, 11}}}, {5, 9, {{2, 18}, {3, 2}}},
    {6, 11, {{2, 28}, {3, 3}}}, {7, 13, {{2, 41}, {3, 5}}}, {8, 15, {{2, 57}, {3, 7}, {5, 1}}}, {9, 17, {{2, 74}, {3, 10}, {5, 2}}},
    {10, 19, {{2, 88}, {3, 14}, {5, 3}}}, {11, 21, {{2, 95}, {3, 17}, {5, 5}}}, {12, 23, {{2, 89}, {3, 17}, {5, 7}, {7, 1}}},
    {13, 25, {{2, 66}, {3, 13}, {5, 7}, {7, 2}}}, {14, 27, {{2, 32}, {3, 5}, {5, 3}, {7, 1}}}, {15, 22, {{2, 7}}}, {16, 8, {}}};

static bool same_homology(const HomologyResult& a, const HomologyResult& b) {
    if (a.free != b.free) return false;
    if (a.tors.size() != b.tors.size()) return false;
    for (auto& [n, t] : a.tors) { if (!b.tors.count(n) || b.tors.at(n) != t) return false; }
    return true;
}
static bool same_byQ(const HomologyResult& a, const HomologyResult& b, string& why) {
    set<pair<int, int>> keys;
    for (auto& [Q, m] : a.byQ) for (auto& [n, x] : m) if (x.first || !x.second.empty()) keys.insert({Q, n});
    for (auto& [Q, m] : b.byQ) for (auto& [n, x] : m) if (x.first || !x.second.empty()) keys.insert({Q, n});
    for (auto [Q, n] : keys) {
        pair<long, vector<mpz_class>> xa{0, {}}, xb{0, {}};
        if (a.byQ.count(Q) && a.byQ.at(Q).count(n)) xa = a.byQ.at(Q).at(n);
        if (b.byQ.count(Q) && b.byQ.at(Q).count(n)) xb = b.byQ.at(Q).at(n);
        auto na = normalize_invariants(xa.second), nb = normalize_invariants(xb.second);
        if (xa.first != xb.first || na != nb) { why = "Q=" + to_string(Q) + " n=" + to_string(n); return false; }
    }
    return true;
}

// =============================================================================================
// Section 10.  Drivers
// =============================================================================================
static void run_verify(int KMAX, int KREF, bool do_cl, bool do_pu, int KCERT) {
    printf("=== VERIFICATION (weights <= %d; full-complex reference <= %d; both pipelines compared <= %d) ===\n", KMAX, KREF, KCERT);
    fflush(stdout);
    double t0 = now_sec();
    // V1 -- the identities that make R_D and the tiny R'_D semifree resolutions
    printf("V1  the resolutions R_D (for A) and R'_D (for A', Section 6b)\n");
    RES.build(KMAX); RESP.build(KMAX);
    bool ok = true; long ng = 0; string why;
    try { ng = RES.verify_identities(KMAX, true); } catch (exception& e) { ok = false; why = e.what(); }
    report("R_D: (d_0+D_R)^2 = 0, d_0^2 = 0, eps d = D eps, D_R raises pins, Q/degree/weight homogeneous, semifree order",
           ok, ok ? to_string(ng) + " generators, weight <= " + to_string(KMAX) : why);
    ok = true; long ngp = 0;
    try { ngp = RESP.verify_identities(KMAX, true); } catch (exception& e) { ok = false; why = e.what(); }
    report("R'_D (Koszul in x_2): the same identities", ok, ok ? to_string(ngp) + " generators" : why);
    { string w2; bool okp = phi_is_dga_map(KMAX, w2);
      report("phi: A' -> A, x_2^i |-> i! gamma_i, is a map of dg-algebras on the whole basis (weight <= " + to_string(KMAX) + ")", okp, okp ? "" : w2); }
    // V2 -- R_0 -> A_0 and R_D -> A are quasi-isomorphisms (cone acyclic), small weights
    printf("V2  R_0 and R_D are resolutions (cone of the augmentation is acyclic)\n");
    int W2 = min(KMAX, 10);
    { bool a = true, b = true, a2 = true, b2 = true; long sz = 0;
      for (int w = 0; w <= W2; w++) {
          a = a && cone_acyclic(w, false); b = b && cone_acyclic(w, true, &sz);
          a2 = a2 && cone_acyclic(w, false, nullptr, &RESP); b2 = b2 && cone_acyclic(w, true, nullptr, &RESP);
      }
      report("Cone(R_0 -> A_0) acyclic, every weight <= " + to_string(W2), a);
      report("Cone(R_D -> (A,D)) acyclic, every weight <= " + to_string(W2), b, "(" + to_string(sz) + " basis elements at the top weight)");
      report("Cone(R'_0 -> A'_0) and Cone(R'_D -> (A',D)) acyclic, every weight <= " + to_string(W2), a2 && b2); }
    // V3 -- napolitano.boundary IS the Hochschild differential
    printf("V3  Napolitano's complex is the Hochschild complex of (A,D)\n");
    int W3 = min(KMAX, 10);
    { long tot = 0, b1 = 0, b2 = 0, bq = 0;
      for (int k = 1; k <= W3; k++) { NapIndex I(k, false);
        for (ll x = 0; x < I.total; x++) {
            NapCell c = I.decode(x); if (I.encode(c) != x) b1++;
            map<ll, ll> nb; for (auto& [f, co] : nap_boundary(c)) { nb[I.encode(f)] += co; if (NapIndex::Qof(f) != NapIndex::Qof(c)) bq++; }
            for (auto it = nb.begin(); it != nb.end();) if (!it->second) it = nb.erase(it); else ++it;
            if (std_hochschild(c, -1, I) != nb) b1++;
            auto tw = std_hochschild(c, +1, I); map<ll, ll> t2;
            for (auto& [f, co] : tw) { NapCell fc = I.decode(f); t2[f] = ((c.F.size() + fc.F.size()) & 1) ? -co : co; }
            if (t2 != nb) b2++;
            tot++; } }
      report("boundary(cell) == (d_int - d_ext)(cell) on A (x)_{A^e} B(A,A,A), every cell, k <= " + to_string(W3), b1 == 0, to_string(tot) + " cells");
      report("boundary == phi (d_int + d_ext) phi^-1, phi = (-1)^{#letters} (standard Hochschild complex)", b2 == 0);
      report("every face preserves Q = sum floor(m/2) + #pins, k <= " + to_string(W3), bq == 0); }
    // V4 -- comparison with the reference implementation (full complex, original algorithm)
    printf("V4  new pipeline vs reference implementation of the original algorithm (full Napolitano complex)\n");
    printf("     %3s %4s %12s %9s | %8s %9s | %s\n", "k", "surf", "full cells", "ref time", "small", "new time", "result");
    for (int k = 1; k <= KREF; k++) for (int pu = 0; pu <= 1; pu++) {
        if ((pu && !do_pu) || (!pu && !do_cl)) continue;
        RefStats st; HomologyResult R;
        bool ok = true; string why;
        try { R = homology_reference(k, pu, &st); } catch (exception& e) { ok = false; why = e.what(); }
        double ta = now_sec(); HomologyResult N = homology_small(k, pu); double tn = now_sec() - ta;
        bool same = ok && same_homology(R, N);
        string wq; bool sq = ok && same_byQ(R, N, wq);
        printf("     %3d %4s %12lld %8.2fs | %8ld %8.3fs | %s\n", k, pu ? "pu" : "cl", (long long)st.cells, st.t, N.cells, tn,
               same && sq ? "identical (every Q-block, every degree)" : ("DIFFERENT " + wq + " " + why).c_str());
        report("k=" + to_string(k) + (pu ? " punctured" : " closed") + ": d^2=0 on all " + to_string(st.d2_checked) + " cells; exact H^BM_* and every (Q,n) summand identical", same && sq);
    }
    // V5..V8 -- every weight up to KMAX
    printf("V5-V8  checks on every weight k <= %d (both surfaces)\n", KMAX);
    auto ceC = ce_betti(CE_CLOSED(), KMAX), ceP = ce_betti(CE_PUNCT(), KMAX);
    map<int, map<int, map<int, long>>> bs; for (int p : {3, 5, 7}) bs[p] = bs_mod_p_betti(p, KMAX, ceP);
    const vector<ull> Qp = {2, 3, 5, 7, 11, 13, 2305843009213693951ull};
    map<int, HomologyResult> closedRes;
    for (int k = 1; k <= KMAX; k++) for (int pu = 0; pu <= 1; pu++) {
        if ((pu && !do_pu) || (!pu && !do_cl)) continue;
        map<int, BlockComplex> blocks;
        double ta = now_sec();
        bool cert = (k <= KCERT);
        HomologyResult H, H2;
        bool methods_agree = true;
        if (cert) {
            H = homology_small(k, pu, true, &blocks);                  // Section 5 pipeline (self-certifying)
            H2 = homology_padic(k, pu);                                // Section 6b pipeline (tiny complex + p-adic)
            methods_agree = same_homology(H, H2);
            string w3; methods_agree = methods_agree && same_byQ(H, H2, w3);
        } else {
            D2_CHECKED = 0;
            H = homology_padic(k, pu, false, nullptr, true);           // Section 6b only (the certified one is too slow)
        }
        string tag = "k=" + to_string(k) + (pu ? " pu" : " cl");
        // d^2 = 0 on the small complex
        bool d2 = true;
        if (!cert) d2 = (D2_CHECKED == H.cells);                       // checked inside the pipeline, block by block
        for (auto& [Q, B] : blocks) for (size_t i = 0; i < B.deg.size(); i++) {
            map<int, ll> acc; for (auto& [f, c] : B.bnd[i]) for (auto& [g, c2] : B.bnd[f]) acc[g] = zadd(acc[g], zmul(c, c2));
            for (auto& [g, x] : acc) if (x) d2 = false; }
        // Euler characteristics: original complex (counted), small complex, topology
        map<int, ll> cc = nap_cell_counts(k, pu); ll chi_o = 0, chi_s = 0, chi_h = 0;
        for (auto& [n, c] : cc) chi_o += (n & 1) ? -c : c;
        for (auto& [n, c] : H.ncells_by_deg) chi_s += (n & 1) ? -c : c;
        for (auto& [n, f] : H.free) chi_h += (n & 1) ? -f : f;
        ll chi_top = pu ? ((k & 1) ? -1 : 1) : 0;                 // chi(B_k(S)) = binom(chi(S), k)
        // SNF divisibility chains and the rank census over F_q on the UNREDUCED boundary matrices:
        // rank_{F_q}(d_n) must equal (#unit pivots at n) + rank_Q(reduced d_n) - #{elementary divisors divisible by q}
        bool chain = true; atomic<bool> census(true); string cwhy; mutex mu;   // census: only where blocks are kept
        struct CT { int Q, n; const DegData* dd; vector<vector<pair<int, ll>>> cols; };
        vector<CT> cts;
        for (auto& [Q, br] : H.blocks) {
            if (!cert) break;                                          // blocks are not kept beyond KCERT
            const BlockComplex& B = blocks.at(Q);
            map<int, vector<int>> bydeg; for (size_t c = 0; c < B.deg.size(); c++) bydeg[B.deg[c]].push_back(c);
            vector<int> loc(B.deg.size()); for (auto& [n, cs] : bydeg) for (size_t i = 0; i < cs.size(); i++) loc[cs[i]] = i;
            for (auto& [n, dd] : br.deg) {
                for (size_t t = 0; t + 1 < dd.divs.size(); t++) if (dd.divs[t + 1] % dd.divs[t] != 0) chain = false;
                CT ct{Q, n, &dd, {}};
                for (int c : bydeg[n]) { ct.cols.push_back({}); for (auto& [f, x] : B.bnd[c]) ct.cols.back().push_back({loc[f], x}); }
                cts.push_back(std::move(ct));
            }
        }
        if (!cert) cts.clear();
        parallel_for(cts.size() * Qp.size(), [&](int t) {
            const CT& ct = cts[t / Qp.size()]; ull q = Qp[t % Qp.size()];
            long rq = rank_mod_q(ct.cols, q), nd = 0;
            for (auto& x : ct.dd->divs) { mpz_class r; mpz_mod(r.get_mpz_t(), x.get_mpz_t(), to_mpz((ll)q).get_mpz_t()); if (r == 0) nd++; }
            if (rq != ct.dd->pivots + ct.dd->rank - nd) { census = false; lock_guard<mutex> g(mu); cwhy = "Q=" + to_string(ct.Q) + " n=" + to_string(ct.n) + " q=" + to_string(q); }
        });
        // rational Betti numbers vs Knudsen; mod-p Betti numbers vs Bianchi-Stavrou (punctured)
        Ordinary O = to_ordinary(k, H);
        map<int, long> b; for (auto& [i, f] : O.free) if (f) b[i] = f;
        map<int, long> ce; for (auto& [i, f] : (pu ? ceP : ceC)[k]) if (f) ce[i] = f;
        bool knud = (b == ce);
        bool bsok = true;
        if (pu) for (int p : {3, 5, 7}) {
            map<int, long> got;
            set<int> is; for (auto& [i, f] : O.free) is.insert(i); for (auto& [i, t] : O.pcount) { is.insert(i); is.insert(i + 1); }
            for (int i : is) { long v = (O.free.count(i) ? O.free[i] : 0);
                if (O.pcount.count(i) && O.pcount[i].count(p)) v += O.pcount[i][p];
                if (O.pcount.count(i - 1) && O.pcount[i - 1].count(p)) v += O.pcount[i - 1][p];
                if (v) got[i] = v; }
            map<int, long> want; for (auto& [i, v] : bs[p][k]) if (v) want[i] = v;
            if (got != want) bsok = false;
        }
        char det[256]; snprintf(det, sizeof det, "cells %ld -> %ld, S=%s, %.2fs", H.cells, H.cells_red,
            [&] { string s; for (ll p : H.S) s += (s.empty() ? "" : ",") + to_string(p); return s; }().c_str(), now_sec() - ta);
        bool all = d2 && chi_o == chi_s && chi_s == chi_h && chi_h == chi_top && chain && census && knud && bsok && methods_agree;
        report(tag + ": d^2=0" + (d2 ? "" : "(NO)") + ", chi " + to_string(chi_o) + "=" + to_string(chi_s) + "=" + to_string(chi_h) + "=binom(chi,k)" +
               (chi_h == chi_top ? "" : "(NO)") + ", SNF chains" + (chain ? "" : "(NO)") + (cert ? string(", F_q-rank census") + (census ? "" : "(NO:" + cwhy + ")") : string("")) +
               ", Knudsen" + (knud ? "" : "(NO)") + (pu ? string(", Bianchi-Stavrou p=3,5,7") + (bsok ? "" : "(NO)") : string("")) +
               (cert ? string(", both pipelines agree") + (methods_agree ? "" : "(NO)") : string(", Section 6b only")), all, det);
        if (!pu) closedRes[k] = H;
    }
    printf("V9  known values\n");
    if (do_cl && closedRes.count(6)) {
        auto& H = closedRes[6]; bool ok = H.free.count(7) && H.free[7] == 9 && H.tors.count(7) && H.tors[7] == vector<mpz_class>(8, 2);
        report("Napolitano [Nap03, Table 2]: H^BM_7(B_6(Sigma_1);Z) = Z^9 + (Z/2)^8", ok);
    }
    if (do_cl) {
        bool h1 = true, thA = true;
        for (auto& [k, H] : closedRes) {
            Ordinary O = to_ordinary(k, H);
            if (k >= 2) { if (!(O.free[1] == 2 && O.tors[1] == vector<mpz_class>{2})) h1 = false; }
            for (int p : {3, 5, 7, 11, 13}) if (k <= 2 * p - 1) for (auto& [i, pc] : O.pcount) if (pc.count(p)) thA = false;
        }
        report("H_1(B_k(Sigma_1);Z) = Z^2 + Z/2 for 2 <= k <= " + to_string(KMAX), h1);
        report("[IR26b, Thm A]/Chen-Zhang: no p-torsion when k <= 2p-1, p = 3,5,7,11,13, k <= " + to_string(KMAX), thA);
        if (closedRes.count(13)) {
            Ordinary O = to_ordinary(13, closedRes[13]);
            vector<mpz_class> want(44, 2); for (int i = 0; i < 7; i++) want.push_back(6); want.push_back(30);
            report("README (exact_snf.py 13): H_8(B_13) = Z^15 + Z/2^44 + Z/6^7 + Z/30", O.free[8] == 15 && O.tors[8] == want);
        }
        if (closedRes.count(14)) {
            auto& H = closedRes[14]; bool ok = true;
            for (auto& [n, f, dv] : PUB_K14_BM) {
                vector<mpz_class> want; for (auto [d, c] : dv) for (int i = 0; i < c; i++) want.push_back(d);
                if ((H.free.count(n) ? H.free[n] : 0) != f || (H.tors.count(n) ? H.tors[n] : vector<mpz_class>{}) != want) ok = false;
            }
            long tot = 0; for (auto& [n, f] : H.free) tot += f; (void)tot;
            report("published table, k = 14 (computed on a larger machine with exact_snf.py): all of H^BM_*(B_14)", ok);
        }
        if (closedRes.count(15)) {
            Ordinary O = to_ordinary(15, closedRes[15]); bool ok = true;
            for (auto& [i, f, pc] : PUB_K15_ORD) {
                if ((O.free.count(i) ? O.free[i] : 0) != f) ok = false;
                map<ll, int> want; for (auto [p, c] : pc) want[p] = c;
                map<ll, int> got; if (O.pcount.count(i)) for (auto& [p, c] : O.pcount[i]) got[p] = c;
                if (got != want) ok = false;
            }
            report("published table, k = 15 (computed on a larger machine with exact_snf.py): all of H_*(B_15)", ok);
        }
    }
    printf("\n%d passed, %d failed%s   (%.1f s)\n", N_PASS, N_FAIL, N_FAIL ? "" : "", now_sec() - t0);
    for (auto& f : FAILED) printf("  FAILED: %s\n", f.c_str());
}

static void run_benchmark(int K1, int K2, int KREF, bool pu) {
    printf("=== BENCHMARK (%s torus) ===\n", pu ? "punctured" : "closed");
    printf("%3s | %16s %9s %9s | %9s %9s %9s %9s %9s | %8s | %s\n", "k", "Napolitano cells", "small", "reduced",
           "t_resol", "t_assem", "t_morse", "t_snf", "total", "peakRSS", "H_* summary");
    for (int k = K1; k <= K2; k++) {
        for (int ref = 0; ref <= 1; ref++) {
            if (ref && k > KREF) continue;
            int fd[2]; if (pipe(fd)) throw runtime_error("pipe");
            pid_t pid = fork();
            if (pid == 0) {
                close(fd[0]);
                double t0 = now_sec(); T_MORSE = T_SNF = 0;
                HomologyResult H; RefStats st;
                double tr = 0;
                if (ref) H = homology_reference(k, pu, &st);
                else { RES.build(k); RESP.build(k); tr = RES.t_build + RESP.t_build; H = homology_padic(k, pu); }
                double tot = now_sec() - t0;
                long fr = 0, tc = 0; set<string> dv;
                for (auto& [n, f] : H.free) fr += f;
                for (auto& [n, t] : H.tors) { tc += t.size(); for (auto& x : t) dv.insert(x.get_str()); }
                string dvs; for (auto& s : dv) dvs += (dvs.empty() ? "" : ",") + s;
                ll full = 0; for (auto& [n, c] : nap_cell_counts(k, pu)) full += c;
                char line[512];
                snprintf(line, sizeof line, "%3d | %16lld %9s %9ld | %9.3f %9.3f %9.3f %9.3f %9.3f | %6.0fMB | %s rank %ld, %ld cyclic torsion summands {%s}\n",
                         k, (long long)full, ref ? "(ref)" : to_string(H.cells).c_str(), H.cells_red, tr, ref ? st.t - T_MORSE - T_SNF : H.t_assemble,
                         T_MORSE, T_SNF, tot, peak_rss_kb() / 1024.0, ref ? "REF" : "NEW", fr, tc, dvs.c_str());
                if (write(fd[1], line, strlen(line)) < 0) _exit(1);
                _exit(0);
            }
            close(fd[1]);
            string out; char buf[1024]; ssize_t r; while ((r = read(fd[0], buf, sizeof buf)) > 0) out.append(buf, r);
            close(fd[0]);
            int status; struct rusage ru; wait4(pid, &status, 0, &ru);
            if (out.empty()) printf("%3d | (child failed, status %d)\n", k, status); else printf("%s", out.c_str());
            fflush(stdout);
        }
    }
}

int main(int argc, char** argv) {
    vector<string> a(argv + 1, argv + argc);
    binom(66, 0);                                            // fill the table before any threads start (C(66,33) < 2^63)
    if (getenv("HO_THREADS")) NTHREADS = max(1, atoi(getenv("HO_THREADS")));
    if (getenv("HO_MIN_TIER")) MIN_TIER = atoi(getenv("HO_MIN_TIER"));
    if (getenv("HO_PLOCAL_MPZ")) FORCE_PLOCAL_MPZ = true;
    auto has = [&](const string& s) { return find(a.begin(), a.end(), s) != a.end(); };
    auto optint = [&](const string& s, int def) { for (size_t i = 0; i + 1 < a.size(); i++) if (a[i] == s) return atoi(a[i + 1].c_str()); return def; };
    bool pu = has("pu") || has("--punctured");
    try {
        if (has("--verify")) {
            int KMAX = 26; for (size_t i = 0; i < a.size(); i++) if (a[i] == "--verify" && i + 1 < a.size() && isdigit(a[i + 1][0])) KMAX = atoi(a[i + 1].c_str());
            int KREF = optint("--ref", 13);
            int KCERT = optint("--cert", 26);
            bool both = has("both") || !pu;
            run_verify(KMAX, KREF, both || !pu, both || pu, KCERT);
            return N_FAIL ? 1 : 0;
        }
        if (has("--benchmark")) {
            int K1 = 1, K2 = 20; for (size_t i = 0; i + 2 < a.size() + 1; i++) if (a[i] == "--benchmark") { if (i + 1 < a.size()) K1 = atoi(a[i + 1].c_str()); if (i + 2 < a.size()) K2 = atoi(a[i + 2].c_str()); }
            run_benchmark(K1, K2, optint("--ref", 12), pu);
            return 0;
        }
        if (has("--reference")) {
            int k = optint("--reference", 8);
            RefStats st; HomologyResult H = homology_reference(k, pu, &st);
            print_result(k, pu, H, true);
            printf("[reference] %lld cells, d^2 = 0 checked on every cell, %ld cells after unit reduction, %.2fs, peak RSS %ld MB\n",
                   (long long)st.cells, H.cells_red, st.t, peak_rss_kb() / 1024);
            return 0;
        }
        vector<int> ks; for (auto& s : a) if (!s.empty() && isdigit(s[0])) ks.push_back(atoi(s.c_str()));
        if (ks.empty()) { fprintf(stderr, "usage: %s K [K2] [pu] [--bm] [--json] [--snf-cert] | --verify [KMAX] [--ref KREF] | --benchmark K1 K2 [--ref KREF] | --reference K\n", argv[0]); return 2; }
        int K1 = ks[0], K2 = ks.size() > 1 ? ks[1] : ks[0];
        bool old_method = has("--snf-cert");                       // Section 5 pipeline (self-certifying prime set)
        for (int k = K1; k <= K2; k++) {
            double t0 = now_sec();
            RES.build(k); RES.verify_identities(k, true);          // the proof obligations of Section 3, every run
            if (!old_method) { RESP.build(k); RESP.verify_identities(k, true); }
            HomologyResult H = old_method ? homology_small(k, pu) : homology_padic(k, pu, has("-v"));
            if (has("--json")) {                                   // machine-readable: {n: [free, [invariant factors > 1]]}
                printf("{\"k\": %d, \"punctured\": %s, \"bm\": {", k, pu ? "true" : "false");
                set<int> ns; for (auto& [n, f] : H.free) ns.insert(n); for (auto& [n, t] : H.tors) ns.insert(n);
                bool first = true;
                for (int n : ns) {
                    printf("%s\"%d\": [%ld, [", first ? "" : ", ", n, H.free.count(n) ? H.free.at(n) : 0); first = false;
                    if (H.tors.count(n)) { bool f2 = true; for (auto& x : H.tors.at(n)) { printf("%s%s", f2 ? "" : ", ", x.get_str().c_str()); f2 = false; } }
                    printf("]]");
                }
                printf("}}\n"); fflush(stdout);
                continue;
            }
            print_result(k, pu, H, has("--bm"));
            string S; for (ll p : H.S) S += (S.empty() ? "" : ",") + to_string(p);
            printf("[k=%d] Napolitano complex: %lld cells; equivalent small complex: %ld cells in %zu Q-blocks -> %ld after unit reduction; "
                   "torsion primes certified complete: {%s}; max coefficient %d bits%s; %.2fs, peak RSS %ld MB\n\n",
                   k, [&] { ll t = 0; for (auto& [n, c] : nap_cell_counts(k, pu)) t += c; return (long long)t; }(), H.cells, H.byQ.size(), H.cells_red,
                   S.c_str(), max(H.maxbits_morse, H.maxbits_snf), H.used_gmp ? " (GMP used)" : "", now_sec() - t0, peak_rss_kb() / 1024);
            fflush(stdout);
        }
    } catch (exception& e) {
        fprintf(stderr, "ERROR: %s\n", e.what());
        return 3;
    }
    return 0;
}
