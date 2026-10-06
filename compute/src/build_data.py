"""Builds ../data.json, the precomputed table read by compute.html.

Sources (outside this repository):
  torus, punctured torus : uconf-torus/tables/homology-torus-k35.json (homology_optimized.cpp)
  plane, annulus         : colcx.py (Borel-Moore column complexes), one JSON line per k
Ordinary homology from Borel-Moore homology of the orientable 2k-manifold B_k(S):
  H_i = Z^{rank BM_{2k-i}} + tors BM_{2k-i-1}  (Poincare duality + universal coefficients).
Each group is stored as [rank, d_1, m_1, d_2, m_2, ...]: Z^rank + (Z/d_1)^m_1 + ..., with d_1 | d_2 | ...
"""
import json, sys
from collections import Counter

TORUS = '/home/fenu/idrissi-rocailucio/uconf-torus/tables/homology-torus-k35.json'

def runs(factors):
    out = []
    for d, m in sorted(Counter(int(x) for x in factors).items()):
        assert d < 2 ** 53
        out += [d, m]
    return out

def invariant_factors(primary):
    """{p: [valuations]} -> invariant factors d_1 | d_2 | ..."""
    n = max((len(v) for v in primary.values()), default=0)
    fs = [1] * n
    for p, vals in primary.items():
        for j, v in enumerate(sorted(vals, reverse=True)):
            fs[n - 1 - j] *= int(p) ** v
    return fs

def from_bm(k, bm):
    """bm: {n: (rank, invariant factors)} -> list over i of [rank, runs...]"""
    H = []
    for i in range(0, 2 * k + 1):
        r = bm.get(2 * k - i, (0, []))[0]
        t = bm.get(2 * k - i - 1, (0, []))[1]
        H.append([r] + runs(t))
    while len(H) > 1 and H[-1] == [0]:
        H.pop()
    return H

def main(plane_jsonl, annulus_jsonl, out):
    T = json.load(open(TORUS))
    data = {}
    for key, name in (('torus', 'closed'), ('punctured', 'punctured')):
        tab = {'0': [[1]]}
        for k, degs in T[name].items():
            n = max(int(i) for i in degs)
            tab[k] = [[degs[str(i)]['rank']] + runs(degs[str(i)]['invariant_factors']) if str(i) in degs else [0]
                      for i in range(n + 1)]
            while len(tab[k]) > 1 and tab[k][-1] == [0]: tab[k].pop()
        data[key] = {'kmax': T['kmax'], 'H': tab}
    for key, f in (('plane', plane_jsonl), ('annulus', annulus_jsonl)):
        tab = {'0': [[1]]}
        for line in open(f):
            d = json.loads(line); k = d['k']
            bm = {int(n): (b, invariant_factors(tor)) for n, (b, tor) in d['bm'].items()}
            tab[str(k)] = from_bm(k, bm)
        data[key] = {'kmax': max(int(k) for k in tab), 'H': tab}
    json.dump(data, open(out, 'w'), separators=(',', ':'))

if __name__ == '__main__':
    main(*sys.argv[1:4])
