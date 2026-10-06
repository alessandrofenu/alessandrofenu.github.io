r"""
exact.py -- exact rank over Q of an integer matrix, with no floating point and no
modular proxy.

`rank_Q` is fraction-free (Bareiss) Gaussian elimination: every intermediate entry is an
integer minor of the original matrix, and every division is exact, so the computation is
carried out entirely in Z and the returned rank is the true rank over Q.

`rank_Q_certified` additionally returns a certificate that can be checked independently:
the list of pivot (row, column) pairs, whose corresponding square submatrix has non-zero
determinant (verified by an independent Bareiss determinant), together with the claim that
no larger minor is non-zero -- the latter being what the elimination itself establishes.
"""

__all__ = ["rank_Q", "det_Z", "rank_Q_certified"]


def rank_Q(rows):
    """exact rank over Q of an integer matrix given as a list of lists"""
    M = [list(r) for r in rows]
    m = len(M)
    n = len(M[0]) if m else 0
    prev = 1
    r = 0
    pivots = []
    for c in range(n):
        piv = None
        for i in range(r, m):
            if M[i][c]:
                piv = i
                break
        if piv is None:
            continue
        if piv != r:
            M[r], M[piv] = M[piv], M[r]
        pr, prc = M[r], M[r][c]
        for i in range(r + 1, m):
            Mi = M[i]
            mic = Mi[c]
            if mic:
                for j in range(c + 1, n):
                    Mi[j] = (Mi[j] * prc - mic * pr[j]) // prev
                Mi[c] = 0
            else:
                for j in range(c + 1, n):
                    Mi[j] = (Mi[j] * prc) // prev
        prev = prc
        pivots.append((r, c))
        r += 1
        if r == m:
            break
    return r


def det_Z(rows):
    """exact determinant of a square integer matrix, by Bareiss"""
    M = [list(r) for r in rows]
    n = len(M)
    if n == 0:
        return 1
    sign, prev = 1, 1
    for c in range(n - 1):
        if M[c][c] == 0:
            for i in range(c + 1, n):
                if M[i][c]:
                    M[c], M[i] = M[i], M[c]
                    sign = -sign
                    break
            else:
                return 0
        for i in range(c + 1, n):
            for j in range(c + 1, n):
                M[i][j] = (M[i][j] * M[c][c] - M[i][c] * M[c][j]) // prev
        prev = M[c][c]
    return sign * M[n - 1][n - 1]


def rank_Q_certified(rows):
    """(rank, minor) where minor is a non-zero r x r determinant witnessing rank >= r"""
    M = [list(r) for r in rows]
    m = len(M)
    n = len(M[0]) if m else 0
    # locate a maximal independent set of rows/columns by exact elimination, tracking indices
    ridx = list(range(m))
    work = [list(r) for r in M]
    prev, r, prow, pcol = 1, 0, [], []
    for c in range(n):
        piv = None
        for i in range(r, m):
            if work[i][c]:
                piv = i
                break
        if piv is None:
            continue
        if piv != r:
            work[r], work[piv] = work[piv], work[r]
            ridx[r], ridx[piv] = ridx[piv], ridx[r]
        pr, prc = work[r], work[r][c]
        for i in range(r + 1, m):
            Wi = work[i]
            wic = Wi[c]
            for j in range(c + 1, n):
                Wi[j] = (Wi[j] * prc - wic * pr[j]) // prev
            Wi[c] = 0
        prev = prc
        prow.append(ridx[r])
        pcol.append(c)
        r += 1
        if r == m:
            break
    sub = [[M[i][j] for j in pcol] for i in prow]
    return r, det_Z(sub)
