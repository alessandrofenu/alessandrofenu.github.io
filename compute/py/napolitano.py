r"""
napolitano.py -- Napolitano's Borel-Moore cellular chain complex for the unordered
configuration spaces of the closed torus, in the normalisation of

    N. Idrissi and V. Roca i Lucio, "Odd-primary torsion in the homology of unordered
    configurations on the torus" (2026), Sections 1.4 and 6.1,

which in turn follows F. Napolitano, "On the cohomology of configuration spaces on
surfaces", J. LMS 68 (2003), Sections 2.3-2.4.

A cell is indexed by
        (e, v; (m_1,v_1), ..., (m_d,v_d)),      v, v_i in {0,1},  m_i >= 0,  m_i+v_i >= 1,
where, writing Sigma_1 = S^1 u (S^1 x R):
    e   = number of points on the base circle away from its north pole,
    v   = 1 iff the base north pole is occupied,
    m_i = number of off-pole points in the i-th fibre of the cylinder,
    v_i = 1 iff the i-th fibre contains its north pole ("pinned").
Weight and Borel-Moore degree are
        k = e + v + sum_i (m_i + v_i),          n = e + sum_i (m_i + 1) = e + d + sum_i m_i.

The five types of codimension-one faces and their signed incidence numbers are those
displayed in [IR26b, Section 6.1].  We do not rederive them; `check_d2` verifies d^2 = 0,
and `verify.py` checks the resulting homology against Napolitano's published value.

The punctured complex C^pu(k) is the quotient by the subcomplex of cells with v = 1;
it computes the Borel-Moore homology of B_k(Sigma_1 \ {*}) ~ B_k(Sigma_1^o).
r"""
from functools import lru_cache
from math import comb

__all__ = ["P", "cells", "boundary", "complex_of", "check_d2"]


def P(a, b):
    """the parity-binomial incidence number P(a,b) of [IR26b, Section 6.1]"""
    if a % 2 == 1 and b % 2 == 1:
        return 0
    return comb((a + b) // 2, a // 2)


@lru_cache(maxsize=None)
def _fibre_tuples(total):
    """all ordered tuples ((m_1,v_1),...) with m_i+v_i >= 1 and sum (m_i+v_i) = total"""
    if total == 0:
        return ((),)
    res = []
    for m in range(total + 1):
        for v in (0, 1):
            if not (1 <= m + v <= total):
                continue
            for tail in _fibre_tuples(total - m - v):
                res.append(((m, v),) + tail)
    return tuple(res)


def cells(k):
    """all cells of weight k, as {Borel-Moore degree: [cells]}"""
    out = {}
    for v in (0, 1):
        for e in range(k - v + 1):
            for fibres in _fibre_tuples(k - v - e):
                n = e + len(fibres) + sum(m for m, _ in fibres)
                out.setdefault(n, []).append((e, v, fibres))
    return out


def boundary(cell):
    """the Borel-Moore boundary of a cell of the closed complex, as {face: coefficient}"""
    e, v, F = cell
    d = len(F)
    out = {}

    def add(c, coef):
        if coef:
            out[c] = out.get(c, 0) + coef

    # (1) two adjacent fibres merge, unless both are pinned
    for i in range(1, d):
        (mi, vi), (mj, vj) = F[i - 1], F[i]
        if vi == 1 and vj == 1:
            continue
        coef = (-1) ** (e + i - 1 + sum(m for m, _ in F[:i])) * P(mi, mj)
        add((e, v, F[:i - 1] + ((mi + mj, max(vi, vj)),) + F[i + 1:]), coef)

    # (2) an off-pole point of a fibre reaches that fibre's north pole
    for i in range(1, d + 1):
        mi, vi = F[i - 1]
        if vi == 0 and mi >= 1:
            coef = (-1) ** (e + i - 1 + sum(m for m, _ in F[:i - 1])) * (1 + (-1) ** mi)
            add((e, v, F[:i - 1] + ((mi - 1, 1),) + F[i:]), coef)

    # (3) a base-circle point reaches the base north pole
    if v == 0 and e >= 1:
        add((e - 1, 1, F), -(1 + (-1) ** e))

    # (4) the leftmost fibre slides onto the base circle
    if d >= 1:
        m1, v1 = F[0]
        if not (v == 1 and v1 == 1):
            add((e + m1, max(v, v1), F[1:]), -(-1) ** e * P(e, m1))

    # (5) the rightmost fibre slides onto the base circle
    if d >= 1:
        md, vd = F[-1]
        if not (v == 1 and vd == 1):
            S = (d - 1) + sum(m for m, _ in F[:-1])
            add((e + md, max(v, vd), F[:-1]), (-1) ** (e + S * (1 + md)) * P(e, md))

    return {c: co for c, co in out.items() if co}


def complex_of(k, punctured=False):
    """returns (cells_by_degree, boundary_function) for C^cl(k) or C^pu(k)"""
    C = cells(k)
    if not punctured:
        return C, boundary
    C = {n: [c for c in C[n] if c[1] == 0] for n in C}
    C = {n: v for n, v in C.items() if v}
    keep = {c for v in C.values() for c in v}
    return C, (lambda c: {f: co for f, co in boundary(c).items() if f in keep})


def check_d2(k, punctured=False):
    """returns None if d^2 = 0 on every cell of weight k, else a witness"""
    C, bnd = complex_of(k, punctured)
    for n in sorted(C):
        for c in C[n]:
            acc = {}
            for f, co in bnd(c).items():
                for g, co2 in bnd(f).items():
                    acc[g] = acc.get(g, 0) + co * co2
            for g, val in acc.items():
                if val != 0:
                    return (c, g, val)
    return None
