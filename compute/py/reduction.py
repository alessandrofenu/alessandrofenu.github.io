r"""
reduction.py -- exact linear algebra over Z and over F_p for based chain complexes.

  * `morse_reduce`  : algebraic Morse reduction over Z.  Repeatedly cancels a pair
        (tau, sigma) with <d tau, sigma> = +-1, updating the boundaries of the other
        cells.  This is Gaussian elimination on the based complex, hence an exact chain
        homotopy equivalence: it changes neither the homology nor its torsion.
  * `elementary_divisors` : Smith normal form of an integer matrix (divisors only).
  * `rank_mod`      : rank of an integer matrix over F_p (dense, numpy).  Taking p a
        large prime is used as a proxy for the rank over Q.
"""
import heapq
from collections import defaultdict
import numpy as np

__all__ = ["morse_reduce", "elementary_divisors", "rank_mod", "matrices"]


def morse_reduce(cells_by_deg, bnd, verbose=False):
    r"""Algebraic Morse (Gaussian) reduction of a based chain complex over Z.

    Repeatedly picks a pair (t, s) with <d t, s> = eps = +-1 and cancels it.  Writing
    C_n = <t> + C_n' and C_{n-1} = <s> + C_{n-1}' and d_n = [[eps, beta],[gamma, delta]],
    the Gaussian-elimination lemma for complexes says that the complex is chain homotopy
    equivalent to the one with

        C'_n = C_n',   C'_{n-1} = C_{n-1}',
        d'_n     = delta - gamma eps^{-1} beta,
        d'_{n+1} = (drop the <t> component of d_{n+1}),
        d'_{n-1} = (restrict d_{n-1} to C_{n-1}').

    So homology and torsion are unchanged.  Pivots are chosen greedily by the Markowitz
    fill-in estimate (|d t| - 1)(|up s| - 1) using a lazy heap.

    Returns (cells_by_degree, boundary_dict) of the reduced complex.
    """
    deg, d, up = {}, {}, defaultdict(dict)
    for n, cs in cells_by_deg.items():
        for c in cs:
            deg[c] = n
    # `bnd` rebuilds each face as a fresh tuple, so storing those directly would keep one
    # duplicate cell object per non-zero entry.  `canon` maps a face to the single canonical
    # cell object, which the dicts then share.  Keys compare by equality, so this changes
    # nothing but the memory: ~2.4 kB/cell -> ~1.1 kB/cell, which is what puts k = 13 inside
    # physical RAM.
    canon = {c: c for c in deg}
    for n, cs in cells_by_deg.items():
        for c in cs:
            dd = {canon[f]: co for f, co in bnd(c).items() if co and f in canon}
            d[c] = dd
            for f, co in dd.items():
                up[f][c] = co
    del canon
    n0 = len(d)

    def key(t, s):
        return (len(d[t]) - 1) * (len(up[s]) - 1)

    def push(t):
        dd = d.get(t)
        if not dd:
            return
        for s, co in dd.items():
            if abs(co) == 1:
                heapq.heappush(heap, (key(t, s), id(t), t, s))
                return

    heap = []
    for t in list(d):
        push(t)
    while heap:
        k0, _, t, s = heapq.heappop(heap)
        dd = d.get(t)
        if dd is None or s not in dd or abs(dd[s]) != 1:
            if dd is not None:
                push(t)
            continue
        if key(t, s) > k0:                      # stale priority; re-insert with the true one
            heapq.heappush(heap, (key(t, s), id(t), t, s))
            continue
        eps = dd[s]
        dtau = dict(dd)
        touched = [tp for tp in up[s] if tp != t]
        for tp in touched:                      # d(tp) <- d(tp) - (<d tp, s>/eps) d(t)
            lam = up[s][tp] // eps
            dtp = d[tp]
            for f, cf in dtau.items():
                nv = dtp.get(f, 0) - lam * cf
                if nv == 0:
                    if f in dtp:
                        del dtp[f]
                        del up[f][tp]
                else:
                    dtp[f] = nv
                    up[f][tp] = nv
        for y in list(up[t]):                   # drop the <t> component of d_{n+1}
            del d[y][t]
        del up[t]
        for f in list(d[t]):                    # remove t as a source
            del up[f][t]
        del d[t]
        for y in list(up[s]):                   # remove s as a target (now unused)
            del d[y][s]
        del up[s]
        if s in d:                              # remove s as a source
            for f in list(d[s]):
                del up[f][s]
            del d[s]
        for tp in touched:
            if tp in d:
                push(tp)
    out = defaultdict(list)
    for c in d:
        out[deg[c]].append(c)
    if verbose:
        print(f"    morse reduction: {n0} -> {len(d)} cells")
    return dict(out), d


def elementary_divisors(Min):
    """the non-zero elementary divisors of an integer matrix (list of lists)"""
    M = [row[:] for row in Min]
    rows = len(M)
    cols = len(M[0]) if rows else 0
    res, r, c = [], 0, 0
    while r < rows and c < cols:
        piv, best = None, None
        for i in range(r, rows):
            for j in range(c, cols):
                if M[i][j] != 0 and (best is None or abs(M[i][j]) < best):
                    best, piv = abs(M[i][j]), (i, j)
        if piv is None:
            break
        i, j = piv
        M[r], M[i] = M[i], M[r]
        for row in M:
            row[c], row[j] = row[j], row[c]
        while True:
            changed = False
            for i in range(r + 1, rows):
                if M[i][c]:
                    q = M[i][c] // M[r][c]
                    for j in range(c, cols):
                        M[i][j] -= q * M[r][j]
                    if M[i][c]:
                        M[r], M[i] = M[i], M[r]
                        changed = True
            for j in range(c + 1, cols):
                if M[r][j]:
                    q = M[r][j] // M[r][c]
                    for i in range(r, rows):
                        M[i][j] -= q * M[i][c]
                    if M[r][j]:
                        for i in range(r, rows):
                            M[i][c], M[i][j] = M[i][j], M[i][c]
                        changed = True
            if not changed:
                break
        a, ok = M[r][c], True
        for i in range(r + 1, rows):
            for j in range(c + 1, cols):
                if M[i][j] % a:
                    for jj in range(c, cols):
                        M[r][jj] += M[i][jj]
                    ok = False
                    break
            if not ok:
                break
        if not ok:
            continue
        res.append(abs(a))
        r += 1
        c += 1
    return res


def rank_mod(A, p):
    """rank over F_p of an integer numpy matrix"""
    if A.size == 0:
        return 0
    B = (A % p).astype(np.int64)
    r, rows, cols = 0, *B.shape
    for c in range(cols):
        nz = np.nonzero(B[r:, c])[0]
        if nz.size == 0:
            continue
        i = r + nz[0]
        B[[r, i]] = B[[i, r]]
        B[r] = (B[r] * pow(int(B[r, c]), p - 2, p)) % p
        col = B[r + 1:, c].copy()
        nzr = np.nonzero(col)[0]
        if nzr.size:
            B[r + 1:][nzr] = (B[r + 1:][nzr] - np.outer(col[nzr], B[r])) % p
        r += 1
        if r == rows:
            break
    return r


def matrices(cells_by_deg, d):
    """the boundary matrices of a based complex, as {degree: numpy array}"""
    idx = {n: {c: i for i, c in enumerate(cs)} for n, cs in cells_by_deg.items()}
    M = {}
    for n in sorted(cells_by_deg):
        src, tgt = cells_by_deg.get(n, []), cells_by_deg.get(n - 1, [])
        A = np.zeros((len(tgt), len(src)), dtype=np.int64)
        ti = idx.get(n - 1, {})
        for j, c in enumerate(src):
            for f, co in d[c].items():
                if f in ti:
                    A[ti[f], j] += co
        M[n] = A
    return M
