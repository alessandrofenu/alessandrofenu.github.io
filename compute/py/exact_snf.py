r"""
exact_snf.py -- exact integral Smith normal form, one prime at a time.

`reduction.elementary_divisors` is a textbook Smith normal form over Z.  It is correct, but
its intermediate entries are uncontrolled: the gcd-refinement loop ping-pongs between rows
and columns and the entries can blow up far beyond the (tiny) final divisors.  On the
reduced complexes here that makes it a lottery -- `morse_reduce` is non-deterministic, and
whether weight 12 finishes in six seconds or runs for an hour depends on which reduced
complex came out.  That is no basis for a published table.

Smith normal form localises: if d_1 | ... | d_r are the elementary divisors of A, then

        d_i  =  prod_p  p^{e^p_i},        e^p_1 <= ... <= e^p_r

where the e^p_i are the valuations of the elementary divisors of A over the local ring Z_p.
So the integral answer can be assembled from one computation per prime, and each of those is
well behaved: over Z_p one pivots on an entry of *minimal p-valuation*, every multiplier
b/a then lies in Z_p, and the elimination finishes in exactly r steps with no gcd
refinement at all.  Carried out modulo p^N the entries never exceed p^N, so nothing grows.

`snf_valuations` does that, and `integral_bm` / `integral_ordinary` assemble the integral
homology from it.  `torsion_primes_excluded` is the companion check that no prime outside
the set considered can contribute: a prime q divides some elementary divisor of A if and
only if rank_{F_q} A < rank_Q A, and both of those are computed exactly.

Validated against `reduction.elementary_divisors` on every boundary matrix of both
complexes for k <= 11 (test V14 of `verify.py`).
"""
from napolitano import complex_of
from reduction import morse_reduce, matrices, rank_mod
from exact import rank_Q

__all__ = ["snf_valuations", "integral_bm", "integral_ordinary",
           "torsion_primes_excluded", "PrecisionError"]


class PrecisionError(Exception):
    """p^N was not enough precision; retry with a larger N"""


def _primes_upto(n):
    """the primes <= n, by sieve (keeps this module numpy-only, like the rest)"""
    sieve = bytearray([1]) * (n + 1)
    sieve[0:2] = b"\x00\x00"
    for i in range(2, int(n ** 0.5) + 1):
        if sieve[i]:
            sieve[i * i::i] = bytearray(len(sieve[i * i::i]))
    return [i for i in range(n + 1) if sieve[i]]


def _val(x, p, cap):
    """p-adic valuation of x, capped at `cap` (x = 0 reads as the cap)"""
    if x == 0:
        return cap
    v = 0
    while x % p == 0:
        x //= p
        v += 1
        if v >= cap:
            return cap
    return v


def snf_valuations(rows, p, rank, N=200):
    r"""The p-adic valuations e^p_1 <= ... <= e^p_rank of the elementary divisors of an
    integer matrix, computed exactly by Smith normal form over Z_p modulo p^N.

    `rank` is the rank over Q, computed independently and exactly; it is how many pivots
    must be found.  Pivoting on an entry of minimal p-valuation v makes every multiplier
    b/a lie in Z_p (every active entry has valuation >= v), so each step is an exact
    unimodular change of basis over Z_p.  Each step costs v digits of precision, which is
    tracked in `prec`; running out raises PrecisionError rather than returning a guess.
    """
    pN = p ** N
    M = [[x % pN for x in row] for row in rows]
    m = len(M)
    n = len(M[0]) if m else 0
    arows, acols = list(range(m)), list(range(n))
    vals, prec = [], N
    while len(vals) < rank:
        if prec <= 0:
            raise PrecisionError(f"p={p}: out of precision after {vals}")
        best = None
        for i in arows:                         # find an entry of minimal p-valuation
            Mi = M[i]
            for j in acols:
                if Mi[j]:
                    v = _val(Mi[j], p, prec)
                    if v < prec and (best is None or v < best[0]):
                        best = (v, i, j)
                        if v == 0:
                            break
            if best is not None and best[0] == 0:
                break
        if best is None:                        # active block vanishes mod p^prec
            raise PrecisionError(f"p={p}: active block = 0 mod p^{prec} with "
                                 f"{rank - len(vals)} pivots still to find")
        v, pi, pj = best
        pv = p ** v
        uinv = pow(M[pi][pj] // pv, -1, p ** (prec - v))
        Mp = M[pi]
        for i in arows:                         # clear the pivot column
            if i == pi or M[i][pj] == 0:
                continue
            lam = (M[i][pj] // pv) * uinv % pN
            Mi = M[i]
            for j in acols:
                Mi[j] = (Mi[j] - lam * Mp[j]) % pN
        for j in acols:                         # clear the pivot row
            if j == pj or Mp[j] == 0:
                continue
            lam = (Mp[j] // pv) * uinv % pN
            for i in arows:
                M[i][j] = (M[i][j] - lam * M[i][pj]) % pN
        vals.append(v)
        arows.remove(pi)
        acols.remove(pj)
        prec -= v
    return sorted(vals)


def _divisors(A, rank, primes):
    """the elementary divisors of A, assembled from its p-local valuations"""
    if rank == 0:
        return []
    rows = A.tolist()
    val = {p: snf_valuations(rows, p, rank) for p in primes}
    out = []
    for i in range(rank):
        d = 1
        for p in primes:
            d *= p ** val[p][i]
        out.append(d)
    return out


def integral_bm(k, primes=(2, 3, 5, 7), punctured=False, verbose=False):
    r"""Exact integral Borel-Moore homology: {n: (free rank, [elementary divisors > 1])}.

    Ranks are exact over Q (Bareiss), torsion exact p-locally for each p in `primes`.  Use
    `torsion_primes_excluded` to certify that no prime outside `primes` contributes.
    """
    C, bnd = complex_of(k, punctured)
    rc, rd = morse_reduce(C, bnd, verbose=verbose)
    M = matrices(rc, rd)
    rk = {n: (rank_Q(M[n].tolist()) if M[n].size else 0) for n in M}
    out = {}
    for n in sorted(rc):
        A = M.get(n + 1)
        ed = _divisors(A, rk.get(n + 1, 0), primes) if A is not None and A.size else []
        out[n] = (len(rc[n]) - rk.get(n, 0) - rk.get(n + 1, 0), [d for d in ed if d > 1])
    return out


def integral_ordinary(k, primes=(2, 3, 5, 7), punctured=False, verbose=False):
    r"""Exact integral H_i(B_k(S);Z) as {i: (free rank, [elementary divisors > 1])}.

    By [IR26b, Lemma 6.1], rank H_i = rank H^BM_{2k-i} and Tors H_i = Tors H^BM_{2k-i-1}.
    """
    bm = integral_bm(k, primes, punctured, verbose)
    free = {2 * k - n: fr for n, (fr, _) in bm.items()}
    tors = {2 * k - n - 1: t for n, (_, t) in bm.items() if t}
    return {i: (free.get(i, 0), tors.get(i, [])) for i in sorted(set(free) | set(tors))}


def torsion_primes_excluded(k, primes, qmax=100, punctured=False):
    r"""Check directly that no prime q <= qmax outside `primes` divides any elementary
    divisor of any boundary matrix: q does so iff rank_{F_q} < rank_Q.  Returns the list of
    offending (degree, q), empty if `primes` is complete up to qmax.
    """
    C, bnd = complex_of(k, punctured)
    rc, rd = morse_reduce(C, bnd)
    M = matrices(rc, rd)
    bad = []
    for n in sorted(M):
        A = M[n]
        if not A.size:
            continue
        r = rank_Q(A.tolist())
        for q in _primes_upto(qmax):
            if q in primes:
                continue
            if rank_mod(A, q) < r:
                bad.append((n, q))
    return bad


if __name__ == "__main__":
    import sys
    k = int(sys.argv[1])
    pu = "pu" in sys.argv
    S = "Sigma_1^o" if pu else "Sigma_1"
    primes = (2, 3, 5, 7)
    H = integral_ordinary(k, primes, pu, verbose=True)
    print(f"H_i(B_{k}({S});Z), exact:")
    for i, (fr, t) in sorted(H.items()):
        if not fr and not t:
            continue
        s = f"  H_{i:2d} = " + (f"Z^{fr}" if fr else "0")
        if t:
            from collections import Counter
            s += " + " + " + ".join(f"Z/{d}" + (f"^{c}" if c > 1 else "")
                                    for d, c in sorted(Counter(t).items()))
        print(s)
    bad = torsion_primes_excluded(k, primes, 100, pu)
    print(f"  primes <= 100 outside {primes} contributing torsion: {bad if bad else 'none'}")
