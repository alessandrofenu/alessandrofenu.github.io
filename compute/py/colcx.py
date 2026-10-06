# Column (Fox-Neuwirth / Napolitano) complexes for B_k(H x V), H,V in {R, S^1}.
#   V = S^1 : column algebra A   = Lambda[x] Gamma[z] Z[y]/y^2, D g_n(z) = -2xy g_{n-1}(z)
#   V = R   : column algebra A_0 = A/y (D = 0)
#   H = S^1 : normalized Hochschild complex; H = R : reduced bar complex.
# Cells a(e,v)[a(m1,v1)|...]; a(m,v) = g_{m//2}(z) x^{m%2} y^v, degree m, weight m+v.
# Boundary = Napolitano convention d_int - d_ext (checked in ../../simpletorsion-checks/hoch.py).
import sys
sys.path.insert(0, '/home/fenu/idrissi-rocailucio/uconf-torus')
from napolitano import P
from exact import rank_Q
from exact_snf import snf_valuations
import numpy as np

def mulA(a, b, ymax=1):
    (m, v), (m2, v2) = a, b
    if v + v2 > ymax: return None
    c = P(m, m2)
    if c == 0: return None
    return c, (m + m2, v + v2)

def DA(a):
    m, v = a
    if v == 0 and m >= 1 and m % 2 == 0:
        return [(-2, (m - 1, 1))]
    return []

def hoch_d(a0, letters):
    out = {}
    def add(k, c):
        if c: out[k] = out.get(k, 0) + c
    e = a0[0]; n = len(letters)
    eps = [e]
    for a in letters: eps.append(eps[-1] + a[0] + 1)
    di, de = {}, {}
    for c, b in DA(a0): di[(b, tuple(letters))] = di.get((b, tuple(letters)), 0) + c
    for i, ai in enumerate(letters, start=1):
        for c, b in DA(ai):
            L = list(letters); L[i-1] = b
            k = (a0, tuple(L)); di[k] = di.get(k, 0) - (-1) ** eps[i-1] * c
    s0 = (-1) ** e
    if n >= 1:
        r = mulA(a0, letters[0])
        if r: k = (r[1], tuple(letters[1:])); de[k] = de.get(k, 0) + s0 * r[0]
        for i in range(1, n):
            r = mulA(letters[i-1], letters[i])
            if r:
                L = list(letters[:i-1]) + [r[1]] + list(letters[i+1:])
                k = (a0, tuple(L)); de[k] = de.get(k, 0) + s0 * (-1) ** (eps[i] - e) * r[0]
        r = mulA(a0, letters[-1])
        if r:
            sg = -(-1) ** ((eps[n-1] - e) * (letters[-1][0] + 1))
            k = (r[1], tuple(letters[:-1])); de[k] = de.get(k, 0) + s0 * sg * r[0]
    for k in set(di) | set(de):
        add(k, di.get(k, 0) - de.get(k, 0))
    return out

def letters_of_weight(w, vs):
    """all tuples of letters (m,v) != (0,0), v in vs, total weight w"""
    if w == 0:
        yield (); return
    for v in vs:
        for m in range(0, w - v + 1):
            if (m, v) == (0, 0): continue
            wt = m + v
            if wt > w: continue
            for rest in letters_of_weight(w - wt, vs):
                yield ((m, v),) + rest

def cells(k, surface):
    vs = (0, 1) if surface in ('torus', 'punct') else (0,)
    out = []
    if surface == 'plane':
        for L in letters_of_weight(k, (0,)):
            out.append(((0, 0), L))
        return out
    v0s = (0, 1) if surface == 'torus' else (0,)
    for v0 in v0s:
        for e in range(0, k - v0 + 1):
            for L in letters_of_weight(k - e - v0, vs):
                out.append(((e, v0), L))
    return out

def degree(c):
    a0, L = c
    return a0[0] + sum(m + 1 for m, v in L)

def boundary(c, surface):
    a0, L = c
    d = hoch_d(a0, list(L))
    res = {}
    for (b0, M), coef in d.items():
        if surface == 'plane' and b0 != (0, 0): continue
        if surface in ('plane', 'annulus') and (b0[1] or any(v for m, v in M)): continue
        if surface == 'punct' and b0[1]: continue
        res[(b0, M)] = res.get((b0, M), 0) + coef
    return {k: v for k, v in res.items() if v}

def bm_homology(k, surface, primes=(2, 3, 5, 7)):
    """{n: (rank, {p: [valuations>0]})} Borel-Moore homology of B_k(surface)"""
    cs = cells(k, surface)
    by = {}
    for c in cs: by.setdefault(degree(c), []).append(c)
    idx = {n: {c: i for i, c in enumerate(L)} for n, L in by.items()}
    mats, ranks = {}, {}
    for n, L in by.items():
        if n - 1 not in by: ranks[n] = 0; continue
        M = np.zeros((len(by[n - 1]), len(L)), dtype=object)
        for j, c in enumerate(L):
            for t, coef in boundary(c, surface).items():
                M[idx[n - 1][t], j] += coef
        mats[n] = M
        ranks[n] = rank_Q(M) if M.size else 0
    out = {}
    for n in sorted(by):
        r_out = ranks.get(n, 0); r_in = ranks.get(n + 1, 0)
        b = len(by[n]) - r_out - r_in
        tors = {}
        if n + 1 in mats and r_in:
            rows = mats[n + 1].tolist()
            for p in primes:
                vals = [v for v in snf_valuations(rows, p, r_in) if v > 0]
                if vals: tors[p] = vals
        out[n] = (b, tors)
    return out

if __name__ == '__main__':
    surf = sys.argv[1]; K = int(sys.argv[2])
    for k in range(1, K + 1):
        h = bm_homology(k, surf)
        print(k, {n: v for n, v in h.items() if v[0] or v[1]})
