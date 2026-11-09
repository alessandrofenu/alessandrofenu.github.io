/* Homotopy groups of spheres: everything known (here) about pi_{n+k}(S^n), with sources.
 * Pure logic, no DOM. Requires data.js (window.HOS_DATA). Exposes window.HOS. */
(function (root) {
    'use strict';
    const D = root.HOS_DATA;

    // ------------------------------------------------------------------ references
    // loc strings are the exact place in the source; every one was checked against the source itself.
    const REFS = {
        hatcher: { a: 'A. Hatcher', y: 2002, t: 'Algebraic Topology', v: 'Cambridge University Press, 2002', url: 'https://pi.math.cornell.edu/~hatcher/AT/AT.pdf', g: 'books' },
        ravenel: { a: 'D. C. Ravenel', y: 2004, t: 'Complex Cobordism and Stable Homotopy Groups of Spheres, 2nd ed.', v: 'AMS Chelsea Publishing 347, 2004', url: 'https://www.sas.rochester.edu/mth/sites/doug-ravenel/mybooks/ravenel.pdf', g: 'books' },
        toda62: { a: 'H. Toda', y: 1962, t: 'Composition Methods in Homotopy Groups of Spheres', v: 'Annals of Mathematics Studies 49, Princeton University Press, 1962', doi: '10.1515/9781400882625', g: 'books' },
        whitehead78: { a: 'G. W. Whitehead', y: 1978, t: 'Elements of Homotopy Theory', v: 'Graduate Texts in Mathematics, Springer, 1978', g: 'books' },
        hlmmoy26: { a: 'Y. Hirato, J. Lee, T. Miyauchi, J. Mukai, M. Ohara, J. Yang', y: 2026, t: 'The Art of Toda Brackets — A Journey through the Homotopy Groups of SO(n) in degrees 23–29', v: 'arXiv:2610.28447, 2026 (book)', url: 'https://arxiv.org/abs/2610.28447', g: 'unstable' },

        hopf31: { a: 'H. Hopf', y: 1931, t: 'Über die Abbildungen der dreidimensionalen Sphäre auf die Kugelfläche', v: 'Math. Ann. 104 (1931), 637–665', doi: '10.1007/BF01457962', g: 'classical' },
        freud38: { a: 'H. Freudenthal', y: 1938, t: 'Über die Klassen der Sphärenabbildungen I. Große Dimensionen', v: 'Compositio Math. 5 (1938), 299–314', url: 'http://www.numdam.org/item/CM_1938__5__299_0/', g: 'classical' },
        serre53: { a: 'J.-P. Serre', y: 1953, t: 'Groupes d’homotopie et classes de groupes abéliens', v: 'Ann. of Math. 58 (1953), 258–294', doi: '10.2307/1969789', g: 'classical' },
        james57: { a: 'I. M. James', y: 1957, t: 'On the suspension sequence', v: 'Ann. of Math. 65 (1957), 74–107', doi: '10.2307/1969666', g: 'classical' },
        adams60: { a: 'J. F. Adams', y: 1960, t: 'On the non-existence of elements of Hopf invariant one', v: 'Ann. of Math. 72 (1960), 20–104', doi: '10.2307/1970147', g: 'classical' },
        adams66: { a: 'J. F. Adams', y: 1966, t: 'On the groups J(X) — IV', v: 'Topology 5 (1966), 21–71', doi: '10.1016/0040-9383(66)90004-8', g: 'classical' },
        thomeier66: { a: 'S. Thomeier', y: 1966, t: 'Einige Ergebnisse über Homotopiegruppen von Sphären', v: 'Math. Ann. 164 (1966), 225–250', doi: '10.1007/BF01360248', g: 'classical' },
        mt67: { a: 'M. Mahowald, M. Tangora', y: 1967, t: 'Some differentials in the Adams spectral sequence', v: 'Topology 6 (1967), 349–369', doi: '10.1016/0040-9383(67)90023-7', g: 'classical' },
        browder69: { a: 'W. Browder', y: 1969, t: 'The Kervaire invariant of framed manifolds and its generalization', v: 'Ann. of Math. 90 (1969), 157–186', doi: '10.2307/1970686', g: 'classical' },
        gray69: { a: 'B. Gray', y: 1969, t: 'On the sphere of origin of infinite families in the homotopy groups of spheres', v: 'Topology 8 (1969), 219–232', doi: '10.1016/0040-9383(69)90012-3', g: 'classical' },
        smith70: { a: 'L. Smith', y: 1970, t: 'On realizing complex bordism modules. Applications to the stable homotopy of spheres', v: 'Amer. J. Math. 92 (1970), 793–856', doi: '10.2307/2373397', g: 'classical' },
        quillen71: { a: 'D. Quillen', y: 1971, t: 'The Adams conjecture', v: 'Topology 10 (1971), 67–80', doi: '10.1016/0040-9383(71)90018-8', g: 'classical' },
        toda71: { a: 'H. Toda', y: 1971, t: 'On spectra realizing exterior parts of the Steenrod algebra', v: 'Topology 10 (1971), 53–65', doi: '10.1016/0040-9383(71)90017-6', g: 'classical' },
        mahowald77: { a: 'M. Mahowald', y: 1977, t: 'A new infinite family in <sub>2</sub>π<sub>*</sub><sup>S</sup>', v: 'Topology 16 (1977), 249–256', doi: '10.1016/0040-9383(77)90005-2', g: 'classical' },
        selick78: { a: 'P. Selick', y: 1978, t: 'Odd primary torsion in π<sub>k</sub>(S<sup>3</sup>)', v: 'Topology 17 (1978), 407–412', doi: '10.1016/0040-9383(78)90007-1', g: 'classical' },
        cmn79: { a: 'F. R. Cohen, J. C. Moore, J. A. Neisendorfer', y: 1979, t: 'The double suspension and exponents of the homotopy groups of spheres', v: 'Ann. of Math. 110 (1979), 549–565', doi: '10.2307/1971238', g: 'classical' },
        neis81: { a: 'J. A. Neisendorfer', y: 1981, t: '3-primary exponents', v: 'Math. Proc. Cambridge Philos. Soc. 90 (1981), 63–83', doi: '10.1017/S0305004100058539', g: 'classical' },
        bjm84: { a: 'M. G. Barratt, J. D. S. Jones, M. E. Mahowald', y: 1984, t: 'Relations amongst Toda brackets and the Kervaire invariant in dimension 62', v: 'J. London Math. Soc. (2) 30 (1984), 533–550', doi: '10.1112/jlms/s2-30.3.533', g: 'classical' },

        mt63: { a: 'M. Mimura, H. Toda', y: 1963, t: 'The (n+20)-th homotopy groups of n-spheres', v: 'J. Math. Kyoto Univ. 3 (1963), 37–58', doi: '10.1215/kjm/1250524854', g: 'unstable' },
        mim65: { a: 'M. Mimura', y: 1965, t: 'On the generalized Hopf homomorphism and the higher composition. Part II. π<sub>n+i</sub>(S<sup>n</sup>) for i = 21 and 22', v: 'J. Math. Kyoto Univ. 4 (1965), 301–326', doi: '10.1215/kjm/1250524662', g: 'unstable' },
        mmo75: { a: 'M. Mimura, M. Mori, N. Oda', y: 1975, t: 'Determination of 2-components of the 23 and 24-stems in homotopy groups of spheres', v: 'Mem. Fac. Sci. Kyushu Univ. Ser. A 29 (1975), 1–42', doi: '10.2206/kyushumfs.29.1', g: 'unstable' },
        oda79: { a: 'N. Oda', y: 1979, t: 'Unstable homotopy groups of spheres', v: 'Bull. Inst. Adv. Res. Fukuoka Univ. 44 (1979), 49–152', g: 'unstable' },
        cm89: { a: 'E. B. Curtis, M. Mahowald', y: 1989, t: 'The unstable Adams spectral sequence for S<sup>3</sup>', v: 'Contemp. Math. 96 (1989), 125–162', doi: '10.1090/conm/096/1022678', g: 'unstable' },
        toda03: { a: 'H. Toda', y: 2003, t: 'Unstable 3-primary homotopy groups of spheres', v: 'Econoinformatics (Himeji Dokkyo Univ.) 29 (2003), 1–77', url: 'https://www.sas.rochester.edu/mth/sites/doug-ravenel/otherpapers/Toda03-ocr.pdf', g: 'unstable' },
        imm15: { a: 'T. Inoue, T. Miyauchi, J. Mukai', y: 2015, t: 'The 2-components of the 31-stem homotopy groups of the 9 and 10-spheres', v: 'J. Fac. Sci. Shinshu Univ. 46 (2015), 1–19; arXiv:1210.0995', url: 'https://arxiv.org/abs/1210.0995', g: 'unstable' },
        mm17: { a: 'T. Miyauchi, J. Mukai', y: 2017, t: 'Determination of the 2-primary components of the 32-stem homotopy groups of S<sup>n</sup>', v: 'Bol. Soc. Mat. Mex. (3) 23 (2017), 319–387', doi: '10.1007/s40590-016-0154-2', g: 'unstable' },
        yw24: { a: 'J. Yang, J. Wu', y: 2024, t: 'On the extension problems for three 33-stem homotopy groups', v: 'arXiv:2406.08621, 2024', url: 'https://arxiv.org/abs/2406.08621', g: 'unstable' },
        yw25: { a: 'J. Yang, J. Wu', y: 2025, t: 'Oda’s extension problem in the 33-stem homotopy group of the 6-sphere', v: 'Homology Homotopy Appl. 27(2) (2025), 53–60', doi: '10.4310/HHA.2025.v27.n2.a3', g: 'unstable' },
        baer26: { a: 'J. F. Baer', y: 2026, t: 'Automated proofs of unstable Adams differentials', v: 'arXiv:2609.23984, 2026; data: doi:10.5281/zenodo.22820890', url: 'https://arxiv.org/abs/2609.23984', g: 'unstable' },

        iwx23: { a: 'D. C. Isaksen, G. Wang, Z. Xu', y: 2023, t: 'Stable homotopy groups of spheres: from dimension 0 to 90', v: 'Publ. Math. IHÉS 137 (2023), 107–243', doi: '10.1007/s10240-023-00139-1', g: 'stable' },
        bix24: { a: 'R. Burklund, D. C. Isaksen, Z. Xu', y: 2024, t: 'Classical stable homotopy groups of spheres via 𝔽<sub>2</sub>-synthetic methods', v: 'arXiv:2408.00987, 2024', url: 'https://arxiv.org/abs/2408.00987', g: 'stable' },
        hhr16: { a: 'M. A. Hill, M. J. Hopkins, D. C. Ravenel', y: 2016, t: 'On the nonexistence of elements of Kervaire invariant one', v: 'Ann. of Math. 184 (2016), 1–262', doi: '10.4007/annals.2016.184.1.1', g: 'stable' },
        lwx24: { a: 'W. Lin, G. Wang, Z. Xu', y: 2024, t: 'On the last Kervaire invariant problem', v: 'arXiv:2412.10879, 2024 (revised 2025), preprint', url: 'https://arxiv.org/abs/2412.10879', g: 'stable' }
    };
    const short = id => {
        const r = REFS[id];
        const fam = r.a.split(',').map(s => s.trim().split(' ').pop());
        const names = fam.length > 3 ? fam[0] + ' et al.' : fam.length === 3 ? fam.join(', ').replace(/, ([^,]*)$/, ' and $1') : fam.join(' and ');
        return names + ' ' + r.y;
    };
    const C = (id, loc) => ({ id, loc: loc || '' });

    // ------------------------------------------------------------------ arithmetic
    const primeCache = [2];
    function primesUpTo(N) {
        for (let c = primeCache[primeCache.length - 1] + 1; primeCache[primeCache.length - 1] < N; c++) {
            if (primeCache.every(p => p * p > c || c % p)) primeCache.push(c);
            if (c > N) break;
        }
        return primeCache.filter(p => p <= N);
    }
    function vp(x, p) { let e = 0; while (x % p === 0) { x /= p; e++; } return e; }
    function ppart(orders, p) {         // p-power parts of a list of cyclic orders (0 = Z ignored)
        const out = [];
        for (const x of orders) { if (x > 1) { const e = vp(x, p); if (e) out.push(p ** e); } }
        return out.sort((a, b) => b - a);
    }
    function primesOf(orders) {
        const s = new Set();
        for (let x of orders) { for (let p = 2; p * p <= x; p++) while (x % p === 0) { s.add(p); x /= p; } if (x > 1) s.add(x); }
        return [...s].sort((a, b) => a - b);
    }
    const logp = (orders, p) => orders.reduce((s, x) => s + vp(x, p), 0);
    const parseStr = s => s === '.' ? [] : s.split('.').flatMap(t => { const [a, b] = t.split('^'); return Array(+(b || 1)).fill(+a); });

    // ------------------------------------------------------------------ p-primary knowledge objects
    // exact: orders known (p-powers, descending). Otherwise: summands (direct summands known),
    // elems (orders of known nonzero elements), log [lo, hi] (bounds on log_p of the order, hi may be Infinity),
    // alts (list of possible order lists), expo (p^expo kills the group).
    const exactP = (p, orders, src, how) => ({ p, exact: true, orders: orders.slice().sort((a, b) => b - a), src: src || [], how: how || '' });
    const partialP = (p, o) => Object.assign({ p, exact: false, summands: [], elems: [], log: [0, Infinity], alts: null, expo: Infinity, src: [], how: '' }, o);
    function logBounds(P) {
        if (P.exact) { const e = logp(P.orders, P.p); return [e, e]; }
        const lo = Math.max(P.log[0], logp(P.summands, P.p), P.elems.length ? Math.max(...P.elems.map(x => vp(x, P.p))) : 0);
        return [lo, P.log[1]];
    }
    function sumP(A, B) {               // knowledge of A ⊕ B at the same prime
        if (A.exact && B.exact) return exactP(A.p, A.orders.concat(B.orders), A.src.concat(B.src));
        const la = logBounds(A), lb = logBounds(B);
        const ea = A.exact ? (A.orders.length ? vp(A.orders[0], A.p) : 0) : A.expo;
        const eb = B.exact ? (B.orders.length ? vp(B.orders[0], B.p) : 0) : B.expo;
        let alts = null;
        if ((A.alts || A.exact) && (B.alts || B.exact) && (A.alts || B.alts)) {
            const aa = A.alts || [A.orders], bb = B.alts || [B.orders];
            alts = []; for (const x of aa) for (const y of bb) alts.push(x.concat(y).sort((u, v) => v - u));
        }
        return partialP(A.p, {
            summands: (A.exact ? A.orders : A.summands).concat(B.exact ? B.orders : B.summands),
            elems: A.elems ? A.elems.concat(B.elems || []) : (B.elems || []),
            log: [la[0] + lb[0], la[1] + lb[1]], alts, expo: Math.max(ea, eb), src: A.src.concat(B.src)
        });
    }
    const zeroP = (p, src, how) => exactP(p, [], src, how);

    // ------------------------------------------------------------------ stable stems
    const imJ2 = k => {                 // 2-primary v1-periodic summands (im J and the mu-family), Ravenel App. A3
        if (k === 1 || k === 2) return [2];
        if (k === 3) return [8];
        const r = k % 8;
        if (r === 7) { const t = (k + 1) / 8; return [2 ** (vp(t, 2) + 4)]; }
        if (r === 3) return [8];
        if (r === 0 || r === 2) return [2];
        if (r === 1) return [2, 2];
        return [];
    };
    const imJstrict2 = k => {           // 2-part of J(pi_k SO), Ravenel Thm 1.1.13
        const r = k % 8;
        if (r === 0 || r === 1) return k > 0 ? [2] : [];
        if (r === 3) return [8];
        if (r === 7) return [2 ** (vp((k + 1) / 8, 2) + 4)];
        return [];
    };
    const imJodd = (k, p) => { const q = 2 * p - 2; return (k + 1) % q === 0 ? [p ** (1 + vp((k + 1) / q, p))] : []; };
    function cokerSmall(k, p) {         // Ravenel Thm 1.1.14: p-component of coker J in stems <= 3pq-6, p odd
        const q = 2 * p - 2;
        const stems = [p * q - 2, (p + 1) * q - 3, 2 * p * q - 4, (2 * p + 1) * q - 5, (2 * p + 1) * q - 2, (2 * p + 2) * q - 3, 3 * p * q - 6];
        return stems.filter(s => s === k).map(() => p);
    }
    const STABLE_MAX = 90;
    function stableAlts(k) {            // 2-primary alternatives at uncertain stems (IWX Table 1)
        const s = D.stable[k];
        if (!s || s[0].indexOf(' or ') < 0) return null;
        return s[0].split(' or ').map(a => ppart(parseStr(a).concat(parseStr(s[2])), 2));
    }
    function stableP(k, p) {
        if (k === 0) return zeroP(p);
        if (2 * p - 3 > k) return zeroP(p, [C('serre53')], 'Serre');
        if (k <= STABLE_MAX) {
            const s = D.stable[k];
            const src = [C('iwx23', 'Table 1')];
            if (k <= 19) src.unshift(C('toda62'), C('ravenel', 'Table A3.6 (Toda’s table)'));
            if (k === 70 || k === 71) src.push(C('bix24', '§1.2 (correction of IWX Table 1)'));
            if (k === 82 || k === 83) src.push(C('bix24', 'Theorem 1.1'));
            if (p === 2) {
                const alts = stableAlts(k);
                if (alts) return partialP(2, { alts, log: [Math.min(...alts.map(a => logp(a, 2))), Math.max(...alts.map(a => logp(a, 2)))], summands: ppart(parseStr(s[2]), 2), src, expo: Math.max(...alts.map(a => vp(a[0] || 1, 2))) });
                return exactP(2, ppart(parseStr(s[0]).concat(parseStr(s[2])), 2), src);
            }
            if (p === 3 && k <= 55) src.push(C('toda03', '§5.6'));
            if (p === 3) src.push(C('ravenel', 'Table A3.4'));
            return exactP(p, ppart(parseStr(s[1]).concat(parseStr(s[2])), p), src);
        }
        if (p === 2) {                  // strict im J is a direct summand; Adams' mu-elements are recorded as elements
            const r = k % 8, strict = imJstrict2(k), elems = (r === 1 || r === 2) ? [2] : [];
            return partialP(2, { summands: strict, elems, elemNote: 'Adams’ element $' + (r === 1 ? '\\mu_{' + k + '}' : '\\eta\\mu_{' + (k - 1) + '}') + '$ (order 2, detected by real K-theory)', src: [C('adams66'), C('quillen71'), C('ravenel', 'Thm 1.1.13 and App. A3')], how: 'imJ' });
        }
        if (p === 3 && k <= 108) {
            return exactP(3, ((D.rav3[k] || []).concat(imJodd(k, 3))), [C('ravenel', 'Table A3.4 and Thm 1.1.13')]);
        }
        const q = 2 * p - 2;
        if (p >= 5 && k <= 3 * p * q - 6) {
            return exactP(p, cokerSmall(k, p).concat(imJodd(k, p)), [C('ravenel', k < p * q - 2 ? 'Thm 1.1.13 (im J; coker J starts in stem pq−2, Thm 1.1.14)' : 'Thms 1.1.13–1.1.14')]);
        }
        const elems = [];
        const src = [C('adams66'), C('quillen71'), C('ravenel', 'Thm 1.1.13')];
        if (p >= 5) {                   // beta_t (Smith) and gamma_t (Toda): elements of order p
            const t = (k + q + 2) / ((p + 1) * q);
            if (Number.isInteger(t) && t >= 1) { elems.push(p); src.push(C('smith70'), C('ravenel', 'Thm 1.3.15')); }
            if (p >= 7) {
                const u = (k + 3 + q * (p + 2)) / (q * (p * p + p + 1));
                if (Number.isInteger(u) && u >= 1) { elems.push(p); src.push(C('toda71'), C('ravenel', 'Thm 1.3.18(b)')); }
            }
        }
        return partialP(p, { summands: imJodd(k, p), elems, elemNote: (elems.length > 1 ? 'elements' : 'an element') + ' of order $' + p + '$ from Smith’s $\\beta_t$ / Toda’s $\\gamma_t$ families', src, how: 'imJ' });
    }
    function stableKnowledge(k) {       // {inf, P: {p: info}} over all primes p <= (k+3)/2
        const P = {};
        for (const p of primesUpTo(Math.floor((k + 3) / 2))) {
            const info = stableP(k, p);
            if (!info.exact || info.orders.length) P[p] = info;
        }
        return { inf: k === 0 ? 1 : 0, P };
    }

    // ------------------------------------------------------------------ unstable data access
    const toda = (k, n) => D.toda[k + ',' + n];
    const twoBook = (k, n) => D.two[k + ',' + n];
    const baer = (k, n) => D.baer[k + ',' + n];
    function three(k, n) {              // Toda 2003 table (odd n), returns orders or null
        if (k > 55 || n % 2 === 0) return null;
        const fams = D.three[k];
        if (!fams) return [];
        const out = [];
        for (const [d, st] of fams) {
            const ns = Object.keys(d).map(Number);
            if (d[n]) out.push(d[n]);
            else if (st && n > Math.max(...ns)) out.push(d[Math.max(...ns)]);
        }
        return out.sort((a, b) => b - a);
    }
    function twoSource(k, n) {          // original papers for the 2-primary group in stem k on S^n
        const b = C('hlmmoy26', k <= 19 ? 'App. B.1' : k === 20 ? 'App. B.2' : k <= 22 ? 'App. B.3' : k <= 24 ? 'App. B.4' : k <= 31 ? 'App. B.5' : k === 32 ? 'App. B.6' : 'App. B.7');
        if (k <= 19) return [C('toda62'), C('ravenel', 'Table A3.6 (Toda’s table)'), b];
        if (k === 20) return [C('mt63'), b];
        if (k <= 22) return [C('mim65'), b];
        if (k <= 24) return [C('mmo75'), b];
        if (k <= 31) return (k === 31 && (n === 9 || n === 10)) ? [C('oda79'), C('imm15', 'Theorem 1.1'), b] : [C('oda79'), b];
        if (k === 32) return [C('mm17'), b];
        if (n === 9) return [C('mm17'), b];
        if (n >= 6) return [C('oda79', 'Prop. 9.20'), C('yw25', 'Theorem 1.1, Corollary 1.2'), b];
        return [C('oda79', 'Prop. 9.20'), b];
    }
    // Thomeier, Satz 1.2 (case k = 8j+1), as quoted in Yang–Wu 2024, Lemma 2.6: integral groups just below the stable range
    function thomeier(k, n) {
        if (k < 9 || k % 8 !== 1 || k > STABLE_MAX) return null;
        const j = (k - 1) / 8, d = n - 8 * j;
        if (d < -4 || d > 2) return null;
        if (stableAlts(k)) return null;
        const extra = { '-4': [], '-3': [], '-2': [], '-1': [2], '0': [2, 2], '1': [2], '2': [] }[d];
        return { extra, inf: d === 2 ? 1 : 0 };
    }

    // ------------------------------------------------------------------ the core: knowledge of pi_{n+k}(S^n)
    // returns {inf: 0|1, P: {p: info}, notes: []}; primes absent from P have zero p-part (proved).
    const memo = new Map();
    function knowledge(n, k) {
        const key = n + ',' + k;
        if (memo.has(key)) return memo.get(key);
        const K = compute(n, k);
        memo.set(key, K);
        return K;
    }
    function compute(n, k) {
        if (k < 0) return { inf: 0, P: {}, kind: 'below' };
        if (k === 0) return { inf: 1, P: {}, kind: 'zero-stem' };
        if (n === 1) return { inf: 0, P: {}, kind: 'circle' };
        if (n === 2) { const K = knowledge(3, k - 1); return { inf: k === 1 ? 1 : K.inf, P: K.P, kind: 'hopf2' }; }
        if (n >= k + 2) { const S = stableKnowledge(k); return { inf: S.inf, P: S.P, kind: 'stable' }; }
        const inf = (n % 2 === 0 && k === n - 1) ? 1 : 0;
        const P = {};
        for (const p of primesUpTo(Math.floor((k + 3) / 2))) {
            const info = p === 2 ? twoInfo(n, k) : oddInfo(n, k, p);
            if (!info.exact || info.orders.length) P[p] = info;
        }
        return { inf, P, kind: 'unstable' };
    }
    function hopfSplitP(n, k, p) {      // n = 4, 8: pi_{n+k}(S^n) = pi_{n-1+k}(S^{n-1}) + pi_{n+k}(S^{2n-1})
        const A = getP(knowledge(n - 1, k), p);
        const s = k - n + 1;
        const B = s <= 0 ? zeroP(p) : getP(knowledge(2 * n - 1, s), p);
        const R = sumP(A, B);
        R.src = [C('hatcher', '§4.2, Exercise 31')].concat(R.src);
        R.how = 'hopf';
        const bb = p === 2 ? baer(k, n) : null;     // intersect with the direct machine bounds for S^n
        if (!R.exact && Array.isArray(bb)) {
            R.log = [Math.max(logBounds(R)[0], bb[0]), Math.min(R.log[1], bb[1])];
            R.src.push(C('baer26', 'Theorem 1.1, Tables 2–5'));
        }
        return R;
    }
    function getP(K, p) { return K.P[p] || zeroP(p); }

    function twoInfo(n, k) {
        const t = toda(k, n);
        if (t) return exactP(2, ppart(t, 2), twoSource(k, n));
        const b = twoBook(k, n);
        if (b) return exactP(2, b.filter(x => x > 0), twoSource(k, n));
        const th = thomeier(k, n);
        if (th) {
            const st = stableP(k, 2);
            if (st.exact) return exactP(2, st.orders.concat(th.extra), [C('thomeier66', 'Satz 1.2'), C('yw24', 'Lemma 2.6')].concat(st.src), 'thomeier');
        }
        if (n === 4 || n === 8) return hopfSplitP(n, k, 2);
        if (k === n - 1) {
            const st = stableP(k, 2);
            if (n % 2 === 0 && st.exact) return exactP(2, st.orders, [C('hatcher', '§4.J (before Cor. 4J.4)'), C('adams60')].concat(st.src), 'whitehead-even');
            if (n % 2 === 1 && st.exact) {
                if (n === 3 || n === 7) return exactP(2, st.orders, [C('hatcher', '§4.J and §4.B'), C('adams60')].concat(st.src), 'whitehead-odd-h');
                const e = logp(st.orders, 2) + 1;
                const R = partialP(2, { log: [e, e], src: [C('hatcher', '§4.J (before Cor. 4J.4)'), C('adams60'), C('whitehead78')].concat(st.src), how: 'whitehead-odd' });
                const bb = baer(k, n); if (bb && bb !== 'inf') R.src.push(C('baer26', 'Tables 2–5'));
                return R;
            }
        }
        const R = partialP(2, {});
        const bb = baer(k, n);
        if (bb && bb !== 'inf') { R.log = [bb[0], bb[1]]; R.src.push(C('baer26', 'Theorem 1.1, Tables 2–5')); R.how = 'baer'; }
        if (n % 2 === 1) { R.expo = n - 1; R.src.push(C('james57')); }
        if (n % 2 === 1 && R.how !== 'baer') R.how = 'james';
        return R;
    }

    function oddInfo(n, k, p) {
        if (2 * p - 3 > k) return zeroP(p, [C('serre53')], 'serre-first');
        if (n % 2 === 0) {               // Serre splitting at odd primes
            const A = getP(knowledge(n - 1, k), p);
            const s = k - n + 1;
            const B = s <= 0 ? zeroP(p) : getP(knowledge(2 * n - 1, s), p);
            const R = sumP(A, B);
            R.src = [C('ravenel', '§1.5, after (1.5.2)'), C('toda03', '(2.1.2)'), C('serre53')].concat(R.src);
            R.how = 'serre-split';
            return R;
        }
        const m = (n - 1) / 2;
        const t = toda(k, n);
        if (t) return exactP(p, ppart(t, p), [C('toda62'), C('ravenel', 'Table A3.6 (Toda’s table)')]);
        if (k < (p - 1) * (n + 1) - 2) {
            const st = stableP(k, p);
            const src = [C('toda03', 'Theorem 2.1 (Serre) and (2.1.7)'), C('ravenel', '§1.5, after Prop. 1.5.8')].concat(st.src);
            if (st.exact) return exactP(p, st.orders, src, 'p-stable');
            return Object.assign({}, st, { src, how: 'p-stable' });
        }
        if (p === 3 && k <= 55) return exactP(3, three(k, n), [C('toda03', '§5.6 (table for k ≤ 55)')], 'toda03');
        const bo = D.bookOdd[n + k];
        if (bo && n + k >= 23 && n + k <= 25) return exactP(p, ppart(bo[n] || [], p), [C('hlmmoy26', { 23: 'Prop. 8.27 (proof)', 24: 'Prop. 9.25 (proof)', 25: 'Prop. 10.28 (proof)' }[n + k])], 'book-odd');
        // partial: Toda's unstable alpha families + exponent theorems
        const q = 2 * p - 2, summands = [], src = [];
        let r = (k + 1) / q;
        if (Number.isInteger(r) && r % p !== 0) { summands.push(p); src.push(C('toda03', 'Theorem 4.2(1)')); }
        r = (k + 2) / q;
        if (Number.isInteger(r) && r % p !== 0 && m < r) { summands.push(p); src.push(C('toda03', 'Theorem 4.2(2)')); }
        let expo = m;
        if (n === 3) { expo = 1; src.push(C('selick78')); }
        else src.push(p === 3 ? C('neis81') : C('cmn79'));
        if (p === 3) src.unshift(C('toda03', '§6.5 has the table for 56 ≤ k < 80 (not transcribed here)'));
        return partialP(p, { summands, expo, src, how: 'alpha' });
    }

    // ------------------------------------------------------------------ formatting
    const TEX_Z = '\\mathbb{Z}';
    function cyclicTeX(o) { return o === 0 ? TEX_Z : TEX_Z + '/' + o; }
    function groupTeX(inf, orders) {    // orders: any list, printed grouped with exponents, descending
        const parts = [];
        if (inf) parts.push(inf > 1 ? TEX_Z + '^{' + inf + '}' : TEX_Z);
        const cnt = new Map();
        orders.slice().sort((a, b) => b - a).forEach(o => cnt.set(o, (cnt.get(o) || 0) + 1));
        for (const [o, c] of cnt) parts.push(c === 1 ? cyclicTeX(o) : '(' + cyclicTeX(o) + ')^{' + c + '}');
        return parts.length ? parts.join(' \\oplus ') : '0';
    }
    function invariantFactors(lists) {  // lists: array of p-power lists (desc) for distinct primes
        const L = Math.max(0, ...lists.map(l => l.length));
        const out = [];
        for (let i = 0; i < L; i++) out.push(lists.reduce((s, l) => s * (l[i] || 1), 1));
        return out;
    }
    function shortText(inf, orders) {  // Toda-style: "∞·120·2²"
        const parts = [];
        if (inf) parts.push('∞');
        const cnt = new Map();
        orders.forEach(o => cnt.set(o, (cnt.get(o) || 0) + 1));
        for (const [o, c] of cnt) parts.push(c === 1 ? String(o) : o + '<sup>' + c + '</sup>');
        return parts.length ? parts.join('·') : '0';
    }

    // summary of a knowledge object
    function summarize(K) {
        const ps = Object.keys(K.P).map(Number).sort((a, b) => a - b);
        const exactPs = ps.filter(p => K.P[p].exact), openPs = ps.filter(p => !K.P[p].exact);
        const inv = invariantFactors(exactPs.map(p => K.P[p].orders));
        return { exact: openPs.length === 0, inf: K.inf, inv, exactPs, openPs, primary: exactPs.flatMap(p => K.P[p].orders) };
    }

    // ------------------------------------------------------------------ cards
    // A card: {cls, icon, title, badge, body (TeX-in-text HTML), refs: [C]}
    function analyze(n, k) {
        const cards = [];
        const N = n + k;
        const pi = `\\pi_{${N}}(S^{${n}})`;
        const res = { n, k, N, cards };
        if (n < 1 || !Number.isInteger(n) || !Number.isInteger(k)) { res.headline = null; return res; }

        // --- trivial cases
        if (k < 0) {
            res.headline = `${pi} = 0`;
            cards.push({ cls: 'exact', icon: 'fa-check', title: 'Exact group', badge: 'Exact', body: `$${pi} = 0$, since $${N} < ${n}$: every map $S^{${N}} \\to S^{${n}}$ is homotopic to a constant one by cellular approximation.`, refs: [C('hatcher', 'Corollary 4.9'), C('ravenel', 'Thm 1.1.2')] });
            return res;
        }
        if (k === 0) {
            res.headline = `${pi} \\cong \\mathbb{Z}`;
            cards.push({ cls: 'exact', icon: 'fa-check', title: 'Exact group', badge: 'Exact', body: `$\\pi_{${n}}(S^{${n}}) \\cong \\mathbb{Z}$, generated by the identity; the isomorphism is the degree (Hurewicz).`, refs: [C('hatcher', 'Corollary 4.25'), C('ravenel', 'Thm 1.1.2 (Hurewicz)'), n === 1 ? C('hatcher', 'Theorem 1.7') : null].filter(Boolean) });
            return res;
        }
        if (n === 1) {
            res.headline = `${pi} = 0`;
            cards.push({ cls: 'exact', icon: 'fa-check', title: 'Exact group', badge: 'Exact', body: `$\\pi_{${N}}(S^1) = 0$: the universal cover $\\mathbb{R} \\to S^1$ induces isomorphisms on $\\pi_i$ for $i \\ge 2$, and $\\mathbb{R}$ is contractible.`, refs: [C('hatcher', 'Proposition 4.1'), C('ravenel', 'Thm 1.1.1')] });
            return res;
        }
        if (n === 2 && k === 1) {
            res.headline = `${pi} \\cong \\mathbb{Z}`;
            cards.push({ cls: 'exact', icon: 'fa-check', title: 'Exact group', badge: 'Exact', body: `$\\pi_3(S^2) \\cong \\mathbb{Z}$, generated by the Hopf map $\\eta\\colon S^3 \\to S^2$ (Hopf invariant one). From the Hopf fibration $S^1 \\to S^3 \\to S^2$: $\\pi_3(S^2) \\cong \\pi_3(S^3)$.`, refs: [C('hopf31'), C('hatcher', '§4.2, Example 4.45 and the discussion after Theorem 4.41'), C('ravenel', 'Example 1.1.7')] });
            return res;
        }

        const K = knowledge(n, k);
        const Sm = summarize(K);
        const stable = n >= k + 2;
        const nn = n === 2 ? 3 : n, kk = n === 2 ? k - 1 : k;   // n = 2 is read through S^3

        // --- headline
        if (Sm.exact) res.headline = `${pi} \\cong ${groupTeX(Sm.inf, Sm.inv)}`;
        else {
            const pieces = [];
            if (Sm.inf) pieces.push(TEX_Z);
            if (Sm.inv.length) pieces.push(groupTeX(0, Sm.inv));
            if (Sm.openPs.length <= 6) Sm.openPs.forEach(p => pieces.push(`G_{(${p})}`));
            else pieces.push(`\\textstyle\\bigoplus_{p \\in \\mathcal{P}} G_{(p)}`);
            res.headline = `${pi} \\cong ${pieces.join(' \\oplus ')}`;
            if (Sm.openPs.length > 6) res.openNote = `𝒫 = {${Sm.openPs.slice(0, 5).join(', ')}, …, ${Sm.openPs[Sm.openPs.length - 1]}} (${Sm.openPs.length} primes)`;
        }
        res.exact = Sm.exact;
        res.order = Sm.exact && !Sm.inf ? Sm.inv.reduce((a, b) => a * b, 1) : null;
        res.primaryTeX = Sm.exact && Sm.primary.length && Sm.inv.some(x => primesOf([x]).length > 1) ? groupTeX(Sm.inf, Sm.primary) : null;

        // --- reductions first (they explain where everything else comes from)
        if (n === 2) {
            cards.push({ cls: 'struct', icon: 'fa-link', title: 'Hopf fibration', badge: 'Reduction', body: `The fibration $S^1 \\to S^3 \\to S^2$ gives $\\pi_i(S^2) \\cong \\pi_i(S^3)$ for $i \\ge 3$, so $${pi} \\cong \\pi_{${N}}(S^3)$ (stem ${k - 1} on $S^3$). Everything below is about $\\pi_{${N}}(S^3)$.`, refs: [C('hopf31'), C('hatcher', '§4.2, discussion after Theorem 4.41')] });
        }
        if (n === 4 || n === 8) {
            const a = `\\pi_{${N - 1}}(S^{${n - 1}})`, b = `\\pi_{${N}}(S^{${2 * n - 1}})`;
            cards.push({ cls: 'struct', icon: 'fa-link', title: n === 4 ? 'Quaternionic Hopf fibration' : 'Octonionic Hopf fibration', badge: 'Reduction', body: `The Hopf bundle $S^{${n - 1}} \\to S^{${2 * n - 1}} \\to S^{${n}}$ has null-homotopic fiber inclusion, so its long exact sequence splits: $${pi} \\cong ${a} \\oplus ${b}$ (the second summand is stem ${k - n + 1} on $S^{${2 * n - 1}}$).`, refs: [C('hatcher', '§4.2, Exercise 31'), C('hopf31')] });
        }

        // --- where are we
        if (stable) {
            cards.push({ cls: 'struct', icon: 'fa-layer-group', title: 'Stable range', badge: 'Theorem', body: `Since $n = ${n} \\ge k + 2 = ${k + 2}$, the Freudenthal suspension theorem gives $${pi} \\cong \\pi_{${k}}^{s}$, the stable ${k}-stem: suspension $\\pi_i(S^n) \\to \\pi_{i+1}(S^{n+1})$ is an isomorphism for $i < 2n-1$.`, refs: [C('freud38'), C('hatcher', 'Corollary 4.24'), C('ravenel', 'Thm 1.1.4, Cor. 1.1.5')] });
        } else if (k === n - 1) {
            const evenTxt = n % 2 === 0
                ? (n === 2 || n === 4 || n === 8 ? `For $n = ${n}$ there is a map of Hopf invariant one, which splits off the $\\mathbb{Z}$ (see the Hopf fibration card).` : `For even $n \\ne 2, 4, 8$ the Whitehead square has Hopf invariant $\\pm 2$ and no element has Hopf invariant $\\pm 1$ (Adams), so $[\\iota_n,\\iota_n]$ generates a $\\mathbb{Z}$ summand and $${pi} \\cong \\mathbb{Z} \\oplus \\pi^s_{${k}}$.`)
                : (n === 3 || n === 7 ? `For $n = ${n}$, $S^{${n}}$ is an H-space, so $[\\iota_n,\\iota_n] = 0$ and $${pi} \\cong \\pi^s_{${k}}$.` : `For odd $n \\ne 1, 3, 7$, $[\\iota_n,\\iota_n] \\ne 0$ ($S^{${n}}$ is not an H-space, by Adams) and $2[\\iota_n,\\iota_n] = 0$ (graded symmetry of Whitehead products), so $|${pi}| = 2\\,|\\pi^s_{${k}}|$; the extension is not determined by this argument.`);
            cards.push({ cls: 'struct', icon: 'fa-angles-right', title: 'Edge of the stable range', badge: 'Theorem', body: `Here $k = n - 1$: suspension $${pi} \\to \\pi^s_{${k}}$ is onto (Freudenthal) and its kernel is generated by the Whitehead square $[\\iota_{${n}},\\iota_{${n}}]$. ${evenTxt}`, refs: [C('hatcher', 'Corollary 4.24 and §4.J (before Cor. 4J.4)'), C('adams60'), C('hatcher', '§4.B'), n % 2 ? C('whitehead78') : null].filter(Boolean) });
        } else {
            cards.push({ cls: 'struct', icon: 'fa-mountain', title: 'Unstable range', badge: 'Position', body: `Here $n = ${n} \\le k = ${k}$, so the group is unstable: it lies ${k - n + 2} step${k - n + 2 > 1 ? 's' : ''} before the stable range $n \\ge k+2$. Its image under iterated suspension lands in $\\pi^s_{${k}}$.`, refs: [C('freud38'), C('hatcher', 'Corollary 4.24')] });
        }

        // --- exact group card
        if (Sm.exact) {
            const src = uniqRefs(Object.values(K.P).flatMap(x => x.src));
            const fromToda = toda(k, n === 2 ? 3 : n) && !stable;
            let body = `$${res.headline}$`;
            if (res.primaryTeX) body += `<br><span class="muted">Primary decomposition: $${res.primaryTeX}$</span>`;
            if (!Object.keys(K.P).length) body += `<br><span class="muted">No torsion in any prime: see the vanishing card below.</span>`;
            body += assemblyNote(K, n, k, fromToda);
            cards.push({ cls: 'exact', icon: 'fa-check', title: 'Exact group', badge: 'Exact', body, refs: src });
        }

        // --- per-prime cards (only when something is open, or to show provenance)
        const ps = Object.keys(K.P).map(Number).sort((a, b) => a - b);
        const silent = [];                           // open primes about which nothing beyond an exponent is known
        for (const p of ps) {
            const I = K.P[p];
            if (I.exact && Sm.exact) continue;      // already in the exact card
            if (!I.exact && p > 7 && !I.alts && !I.summands.length && !I.elems.length && logBounds(I)[1] === Infinity) { silent.push(p); continue; }
            cards.push(primeCard(n, k, p, I));
        }
        if (silent.length) {
            const list = silent.length > 12 ? `${silent.slice(0, 6).join(', ')}, …, ${silent[silent.length - 1]} (${silent.length} primes)` : silent.join(', ');
            const ex = !stable && nn % 2 === 1 ? ` Each $G_{(p)}$ has exponent at most $p^{${nn === 3 ? 1 : (nn - 1) / 2}}$ (see the exponent card).` : '';
            cards.push({ cls: 'prime', icon: 'fa-cubes', title: 'Remaining odd primes', badge: 'Open', body: `For $p \\in \\{${list}\\}$ the $p$-primary part $G_{(p)}$ is not determined by the sources used here, and no nonzero element is known from them.${ex}`, refs: uniqRefs(silent.flatMap(p => K.P[p].src)).slice(0, 6) });
        }

        // --- vanishing primes / first torsion
        if (nn >= 3 && kk >= 1) {
            const firstP = (kk + 3) % 2 === 0 ? (kk + 3) / 2 : null;
            let body = `For a prime $p$, the $p$-primary part of $\\pi_{n+k}(S^n)$ ($n \\ge 3$) vanishes when $k < 2p-3$. ` + (n === 2 ? `Applied to $\\pi_{${N}}(S^3)$ (stem ${kk}), only` : 'Hence only') + ` primes $p \\le ${primesUpTo(Math.floor((kk + 3) / 2)).pop()}$ can appear here.`;
            if (firstP && primesUpTo(firstP).includes(firstP) && firstP > 2) body += ` Since $k = ${kk} = 2\\cdot${firstP}-3$, the ${firstP}-primary part is exactly $\\mathbb{Z}/${firstP}$ (its first appearance).`;
            cards.push({ cls: 'bound', icon: 'fa-filter', title: 'Which primes can occur', badge: 'Theorem', body, refs: [C('serre53')] });
        }

        // --- finiteness / rank
        {
            const infinite = K.inf;
            const body = infinite
                ? `$${pi}$ is infinite: it is $\\mathbb{Z} \\oplus F$ with $F$ finite (here $n$ is even and $k = n - 1$). All other $\\pi_{n+k}(S^n)$ with $k > 0$ are finite.`
                : `$${pi}$ is finite: for $k > 0$, $\\pi_{n+k}(S^n)$ is finite except when $n$ is even and $k = n-1$, where it is $\\mathbb{Z} \\oplus$ finite.`;
            cards.push({ cls: 'struct', icon: 'fa-infinity', title: 'Rank (Serre)', badge: 'Theorem', body, refs: [C('serre53'), C('ravenel', 'Thm 1.1.8')] });
        }

        // --- order card
        orderCard(cards, K, Sm, pi);

        // --- exponents (unstable only)
        if (!stable && nn >= 3) exponentCards(cards, nn, kk, n);

        // --- stable phenomena
        if (stable) stableSpecials(cards, k);
        if (!stable && k === n - 1 && (n & (n - 1)) === 0 && n >= 2) {
            const ok = n <= 8;
            cards.push({ cls: 'special', icon: 'fa-star', title: 'Hopf invariant one', badge: 'Theorem', body: ok ? `$${pi}$ contains an element of Hopf invariant one (a Hopf map), since $n = ${n} \\in \\{2, 4, 8\\}$.` : `$${pi}$ contains no element of Hopf invariant one: such elements exist only for $n = 2, 4, 8$. Equivalently, $h_{${Math.log2(n)}}$ does not survive the Adams spectral sequence ($d_2 h_{${Math.log2(n)}} = h_0 h_{${Math.log2(n) - 1}}^2$).`, refs: [C('adams60'), C('ravenel', 'Thms 1.2.11–1.2.12'), C('hatcher', '§4.B')] });
        }
        if (n === 3 && k <= 52 && !(K.P[2] && K.P[2].exact)) {
            cards.push({ cls: 'note', icon: 'fa-book', title: 'Further source for S³', badge: 'Note', body: `The 2-primary groups $\\pi_{3+k}(S^3)$ were determined through the 52-stem by Curtis and Mahowald with the unstable Adams spectral sequence (as recorded by Baer). Their values are not reproduced on this page.`, refs: [C('cm89'), C('baer26', '§1')] });
        }
        return res;
    }

    function assemblyNote(K, n, k, fromToda) {
        if (n >= k + 2) {
            if (k <= 19) return '<br><span class="muted">Stable stem: Toda’s tables, confirmed by the modern computation of Isaksen–Wang–Xu.</span>';
            const corr = (k === 70 || k === 71) ? ` The 2-primary part of this stem was <b>corrected</b> in 2024 by Burklund–Isaksen–Xu (the value in the published IWX Table 1 is ${shortText(0, ppart(parseStr(D.stableOrig[k]), 2))} for the 2-primary $v_1$-torsion).` : '';
            return `<br><span class="muted">Stable stem, assembled prime by prime (2-primary and odd parts from the sources listed).${corr}</span>`;
        }
        if (fromToda) return '<br><span class="muted">From Toda’s table of $\\pi_{n+k}(S^n)$, $k \\le 19$. The 2-primary part agrees with the 2026 summary of Hirato–Lee–Miyauchi–Mukai–Ohara–Yang, and the 3-primary part with Toda (2003).</span>';
        const hows = Object.values(K.P).map(I => I.how).filter(Boolean);
        let s = '<br><span class="muted">Assembled prime by prime';
        if (hows.includes('thomeier')) s = '<br><span class="muted">From Thomeier’s theorem (the eight groups below the stable range in stems $k \\equiv 1 \\bmod 8$ are determined by $\\pi^s_k$), with $\\pi^s_k$ from Isaksen–Wang–Xu; assembled prime by prime';
        return s + ': see the sources below for each prime.</span>';
    }

    function uniqRefs(list) {
        const seen = new Set(), out = [];
        for (const r of list) { const key = r.id + '|' + r.loc; if (!seen.has(key)) { seen.add(key); out.push(r); } }
        return out;
    }

    function primeCard(n, k, p, I) {
        const G = `G_{(${p})}`;
        let body = '', badge = 'Partial', cls = 'prime';
        if (I.exact) {
            badge = 'Exact';
            body = `$${G} \\cong ${groupTeX(0, I.orders)}$.`;
        } else {
            const lines = [];
            if (I.alts) {
                lines.push(`$${G}$ is one of: ` + I.alts.map(a => `$${groupTeX(0, a)}$`).join(', ') + '. The remaining uncertainty comes from unresolved Adams differentials or hidden extensions.');
            }
            const [lo, hi] = logBounds(I);
            if (!I.alts) {
                if (lo === hi && hi !== Infinity) lines.push(`Order: $|${G}| = ${pw(p, hi)}$` + (I.how === 'whitehead-odd' ? ` $= 2\\,|(\\pi^s_{${k}})_{(2)}|$ (Whitehead square, see the card above).` : '.'));
                else if (hi !== Infinity) lines.push(`Order: $${pw(p, lo)} \\le |${G}| \\le ${pw(p, hi)}$` + (I.how === 'baer' ? ' (machine bounds from the unstable Adams spectral sequence).' : '.'));
                else if (lo > 0) lines.push(`Order: $|${G}| \\ge ${pw(p, lo)}$.`);
            }
            if (I.summands.length && !I.alts) {
                const g = groupTeX(0, I.summands);
                lines.push(I.how === 'imJ' ? `Contains the image of $J$, $${g}$, as a direct summand.`
                    : I.how === 'alpha' ? `Contains $${g}$ as a direct summand, generated by elements of Toda’s unstable $\\alpha$-family.`
                    : `Contains $${g}$ as a direct summand.`);
            }
            if (I.elems.length) lines.push(`Contains ${I.elemNote || 'nonzero elements of order $' + I.elems.join(', ') + '$'} (nonzero; not claimed to be a direct summand).`);
            if (I.expo !== Infinity && !I.alts) lines.push(`Exponent: $${pw(p, I.expo)}\\cdot ${G} = 0$` + (I.expo === 1 && lo >= 1 ? `, so $${G} \\cong (\\mathbb{Z}/${p})^{a}$ for some $a \\ge ${lo}$.` : '.'));
            const nothing = !lines.length;
            if (nothing) lines.push(`Not determined by the sources used here, and no nonzero element is known from them.`);
            body = lines.join('<br>');
            if (nothing) return { cls, icon: 'fa-cube', title: (p === 2 ? '2' : p) + '-primary component', badge: 'Open', body, refs: [] };
        }
        const name = p === 2 ? '2-primary component' : p + '-primary component';
        return { cls, icon: 'fa-cube', title: name, badge, body, refs: uniqRefs(I.src) };
    }

    function orderCard(cards, K, Sm, pi) {
        if (K.inf && Sm.exact) {
            const t = Sm.inv.reduce((a, b) => a * b, 1);
            cards.push({ cls: 'order', icon: 'fa-hashtag', title: 'Torsion subgroup', badge: 'Exact', body: `The torsion subgroup has order $${fmtInt(t)}$` + (t > 1 ? ` $= ${factorTeX(t)}$.` : '.'), refs: [] });
            return;
        }
        if (Sm.exact) {
            const t = Sm.inv.reduce((a, b) => a * b, 1);
            cards.push({ cls: 'order', icon: 'fa-hashtag', title: 'Order', badge: 'Exact', body: `$|${pi}| = ${fmtInt(t)}$` + (t > 1 ? ` $= ${factorTeX(t)}$.` : '.'), refs: [] });
            return;
        }
        // partial: product of known exact parts times bounded parts
        let known = 1, loT = [], hiT = [], open = false;
        for (const p of Object.keys(K.P).map(Number)) {
            const I = K.P[p];
            if (I.exact) { known *= I.orders.reduce((a, b) => a * b, 1); continue; }
            const [lo, hi] = logBounds(I);
            if (hi === Infinity) open = true;
            loT.push(lo ? pw(p, lo) : null); hiT.push(hi !== Infinity ? pw(p, hi) : null);
        }
        const kn = known > 1 ? factorTeX(known) : '';
        let body;
        if (!open && loT.every((x, i) => x === hiT[i]))
            body = `$|${pi}| = ${[kn].concat(hiT).filter(Boolean).join(' \\cdot ')}$` + (K.inf ? ' (torsion part)' : '') + ', although the group structure is not fully determined.';
        else if (!open)
            body = `$${[kn].concat(loT).filter(Boolean).join(' \\cdot ') || '1'} \\le |${K.inf ? 'T' : pi}| \\le ${[kn].concat(hiT).filter(Boolean).join(' \\cdot ')}$` + (K.inf ? ' for the torsion subgroup $T$.' : '.');
        else {
            const lower = [kn].concat(loT).filter(Boolean).join(' \\cdot ');
            body = lower
                ? `$|${K.inf ? 'T' : pi}|$ is divisible by $${lower}$; no upper bound is available from these sources` + (K.inf ? ' ($T$ = torsion subgroup).' : '.')
                : 'The order is not determined by the sources used here.';
        }
        cards.push({ cls: 'order', icon: 'fa-hashtag', title: 'Order', badge: 'Bound', body, refs: [] });
    }

    function exponentCards(cards, n, k, nOrig) {
        if (n % 2 === 1) {
            const m = (n - 1) / 2;
            const lines = [];
            const refs = [];
            if (n === 3) { lines.push(`For every odd prime $p$, $p$ annihilates the $p$-torsion of $\\pi_*(S^3)$ (Selick).`); refs.push(C('selick78')); }
            lines.push(`For an odd prime $p$, $p^{${m}}$ annihilates the $p$-torsion of $\\pi_*(S^{${n}})$ (Cohen–Moore–Neisendorfer for $p \\ge 5$, Neisendorfer for $p = 3$). The bound is attained for suitable stems (Gray).`);
            refs.push(C('cmn79'), C('neis81'), C('ravenel', 'Thm 1.1.10'), C('gray69'));
            lines.push(`$4^{${m}} = 2^{${2 * m}}$ annihilates the 2-primary component of $\\pi_*(S^{${n}})$ (James).`);
            refs.push(C('james57'));
            cards.push({ cls: 'bound', icon: 'fa-compress', title: 'Exponent bounds' + (nOrig === 2 ? ' (for S³)' : ''), badge: 'Theorem', body: lines.join('<br>'), refs });
        } else if (n !== 4 && n !== 8) {
            const m = n / 2;
            cards.push({ cls: 'bound', icon: 'fa-compress', title: 'Exponent bounds (odd primes)', badge: 'Derived', body: `At an odd prime $p$, Serre’s splitting $\\pi_{i}(S^{${n}})_{(p)} \\cong \\pi_{i-1}(S^{${n - 1}})_{(p)} \\oplus \\pi_{i}(S^{${2 * n - 1}})_{(p)}$ and the Cohen–Moore–Neisendorfer theorem give: $p^{${2 * m - 1}}$ annihilates the $p$-torsion of $\\pi_*(S^{${n}})$${m > 1 ? `, and $p^{${m - 1}}$ suffices in stems $k < ${n - 1}$ (where the second summand vanishes)` : ''}.`, refs: [C('ravenel', '§1.5 after (1.5.2); Thm 1.1.10'), C('cmn79'), C('neis81')] });
        }
    }

    function stableSpecials(cards, k) {
        // image of J
        if (k >= 1) {
            const s2 = imJstrict2(k);
            const odd = primesUpTo(Math.floor((k + 3) / 2)).filter(p => p > 2).flatMap(p => imJodd(k, p));
            const all = s2.concat(odd);
            const r = k % 8, lines = [];
            if (all.length) {
                const ord = invariantFactors(primesOf(all).map(p => ppart(all, p)));
                lines.push(`$\\pi^s_{${k}}$ contains the image of $J\\colon \\pi_{${k}}(SO) \\to \\pi^s_{${k}}$ as a direct summand, here $${groupTeX(0, ord)}$` +
                    (k % 4 === 3 ? `, cyclic of order the denominator of $B_{${(k + 1) / 2}}/${k + 1}$ (Bernoulli number).` : ` ($J$ is injective on $\\pi_{${k}}(SO) \\cong \\mathbb{Z}/2$).`));
            }
            if (k > 2 && (r === 1 || r === 2)) lines.push(`It also contains Adams’ element $${r === 1 ? '\\mu_{' + k + '}' : '\\eta\\mu_{' + (k - 1) + '}'}$ of order 2, detected by real K-theory (not in im $J$).`);
            if (k === 2) lines.push('Here $\\pi_2(SO) = 0$; the group is generated by $\\eta^2$, which is $v_1$-periodic (detected by real K-theory).');
            if (lines.length) cards.push({ cls: 'special', icon: 'fa-wave-square', title: 'Image of J (v₁-periodic part)', badge: 'Theorem', body: lines.join('<br>'), refs: [C('adams66'), C('quillen71'), C('ravenel', 'Thm 1.1.13 and App. A3')] });
        }
        // Kervaire invariant one
        const j = Math.log2(k + 2) - 1;
        if (Number.isInteger(j) && j >= 1) {
            let body, refs = [C('browder69'), C('ravenel', 'Thm 3.4.4(a)')];
            const names = { 1: '\\eta^2', 2: '\\nu^2', 3: '\\sigma^2', 4: '\\theta_4', 5: '\\theta_5', 6: '\\theta_6' };
            if (j <= 5) {
                body = `$\\pi^s_{${k}}$ contains an element $${names[j]}$ of Kervaire invariant one (detected by $h_{${j}}^2$).`;
                if (j === 4) refs.push(C('mt67'));
                if (j === 5) refs.push(C('bjm84'));
                refs.push(C('lwx24', 'Abstract (summary of dims 2, 6, 14, 30, 62, 126)'));
            } else if (j === 6) {
                body = `$h_6^2$ is a permanent cycle, so $\\pi^s_{126}$ contains an element $\\theta_6$ of Kervaire invariant one: framed manifolds of Kervaire invariant one exist in dimension 126 (this settles the last open case; preprint).`;
                refs.push(C('lwx24', 'Theorem 1.1'), C('hhr16'));
            } else {
                body = `$\\pi^s_{${k}}$ contains no element of Kervaire invariant one: $h_{${j}}^2$ does not survive for $j \\ge 7$ (Hill–Hopkins–Ravenel).`;
                refs.push(C('hhr16'));
            }
            cards.push({ cls: 'special', icon: 'fa-star', title: 'Kervaire invariant one', badge: 'Theorem', body, refs });
        }
        // Mahowald's eta_j
        const jj = Math.log2(k);
        if (Number.isInteger(jj) && jj >= 3) {
            cards.push({ cls: 'special', icon: 'fa-star', title: 'Mahowald’s family', badge: 'Theorem', body: `$\\pi^s_{${k}}$ contains a nonzero element $\\eta_{${jj}}$, detected by $h_1 h_{${jj}}$ in the Adams spectral sequence.`, refs: [C('mahowald77'), C('ravenel', 'Thm 3.4.4(b)')] });
        }
        // Hopf invariant one (stable form)
        const hj = Math.log2(k + 1);
        if (Number.isInteger(hj) && hj >= 1) {
            cards.push({ cls: 'special', icon: 'fa-star', title: 'Hopf invariant one', badge: 'Theorem', body: hj <= 3 ? `The generator $${['', '\\eta', '\\nu', '\\sigma'][hj]}$ of im $J$ here is a Hopf map: it is detected by $h_{${hj}}$, an element of Hopf invariant one.` : `No element of $\\pi^s_{${k}}$ is detected by $h_{${hj}}$: elements of Hopf invariant one exist only in stems 1, 3, 7.`, refs: [C('adams60'), C('ravenel', 'Thms 1.2.11–1.2.12')] });
        }
        // corrections / uncertainty notes
        if (k === 70 || k === 71) {
            cards.push({ cls: 'warn', icon: 'fa-triangle-exclamation', title: 'Recent correction', badge: '2024', body: `The value of $\\pi^s_{${k}}$ in IWX’s Table 1, and in tables copied from it, is wrong at 2. Machine computations of Lin–Wang–Xu found a nonzero differential $d_5(h_1 p_1)$, and Burklund–Isaksen–Xu give the corrected 2-primary $v_1$-torsion: ${k === 70 ? '$(\\mathbb{Z}/2)^6 \\oplus \\mathbb{Z}/4$' : '$(\\mathbb{Z}/2)^5 \\oplus \\mathbb{Z}/4 \\oplus \\mathbb{Z}/8$'}. This page uses the corrected value.`, refs: [C('bix24', '§1.2'), C('iwx23', 'Table 1')] });
        }
        if (stableAlts(k)) {
            cards.push({ cls: 'warn', icon: 'fa-circle-question', title: 'Open problem', badge: 'Open', body: `The 2-primary part of $\\pi^s_{${k}}$ is not completely known: IWX list the alternatives shown in the 2-primary card. They come from unresolved Adams differentials (their Table 9) and possible hidden extensions (their Tables 17, 20, 23).`, refs: [C('iwx23', '§1.3, Table 1')] });
        }
        if (k === 82 || k === 83) {
            cards.push({ cls: 'note', icon: 'fa-circle-check', title: 'Recently settled', badge: '2024', body: `$\\pi^s_{${k}}$ was known only up to a factor of 2. Burklund–Isaksen–Xu proved that $h_6 g + h_2 e_2$ survives, which completes the computation: the 2-primary part is ${k === 82 ? '$(\\mathbb{Z}/2)^6 \\oplus \\mathbb{Z}/8$' : '$(\\mathbb{Z}/2)^3 \\oplus (\\mathbb{Z}/8)^2$'}.`, refs: [C('bix24', 'Theorem 1.1')] });
        }
        if (k > STABLE_MAX) {
            cards.push({ cls: 'note', icon: 'fa-circle-question', title: 'Beyond the computed range', badge: 'Note', body: `The 2-primary stable stems are computed only through stem 90 (with a few uncertainties above 83), so the 2-primary part of $\\pi^s_{${k}}$ is unknown beyond the summands listed. Odd primes are covered as far as the cited tables and theorems go (3-primary to stem 108).`, refs: [C('iwx23', 'Corollary 1.3'), C('ravenel', 'App. A3')] });
        }
    }

    const pw = (p, e) => e === 1 ? String(p) : `${p}^{${e}}`;
    const fmtInt = x => x.toLocaleString('en-US').replace(/,/g, '\\,');
    function factorTeX(x) {
        const f = [];
        for (let p = 2; p * p <= x; p++) { let e = 0; while (x % p === 0) { x /= p; e++; } if (e) f.push(e > 1 ? `${p}^{${e}}` : `${p}`); }
        if (x > 1) f.push(String(x));
        return f.join(' \\cdot ');
    }

    // ------------------------------------------------------------------ table cells
    // view: 'full' | '2' | '3'
    function cell(n, k, view) {
        if (k === 0) return { cls: 'inf', text: '∞', exact: true };
        if (n === 1) return { cls: 'zero', text: '·', exact: true };
        const stable = n >= k + 2;
        if (view === 'full') {
            if (n === 2 && k === 1) return { cls: 'inf', text: '∞', exact: true };
            const K = knowledge(n, k), S = summarize(K);
            const base = stable ? 'stable' : 'unstable';
            if (S.exact) {
                if (!S.inf && !S.inv.length) return { cls: 'zero', text: '·', exact: true };
                return { cls: S.inf ? 'inf' : base, text: shortText(S.inf, S.inv), exact: true };
            }
            const known = shortText(S.inf, S.inv);
            return { cls: 'partial ' + base, text: (S.inf || S.inv.length ? known + '·' : '') + '<i>?</i>', exact: false, open: S.openPs };
        }
        const p = +view;
        if (n === 2 && k === 1) return { cls: 'inf', text: '∞', exact: true };
        const K = knowledge(n, k);
        const I = K.P[p] || zeroP(p);
        const inf = K.inf;
        const base = stable ? 'stable' : 'unstable';
        if (I.exact) {
            if (!I.orders.length && !inf) return { cls: 'zero', text: '·', exact: true };
            return { cls: inf ? 'inf' : base, text: shortText(inf, I.orders), exact: true };
        }
        const [lo, hi] = logBounds(I);
        const infs = inf ? '∞·' : '';
        if (I.alts) return { cls: 'partial ' + base, text: infs + `${p}<sup>${lo}–${hi}</sup>`, exact: false, order: true };
        if (lo === hi && hi !== Infinity) return { cls: 'order ' + base, text: infs + `|${p}<sup>${hi}</sup>|`, exact: false, order: true };
        if (hi !== Infinity) return { cls: 'order ' + base, text: infs + `${p}<sup>${lo}–${hi}</sup>`, exact: false, order: true };
        return { cls: 'unknown ' + base, text: infs + '?', exact: false };
    }
    function cellTitle(n, k, view) {   // plain-text tooltip
        const K = knowledge(n, k), S = summarize(K);
        const N = n + k;
        if (S.exact) return `π${sub(N)}(S${sup(n)}) ≅ ` + plainGroup(S.inf, S.inv);
        return `π${sub(N)}(S${sup(n)}): partially known (open primes: ${S.openPs.join(', ')}) — click for details`;
    }
    const SUBS = '₀₁₂₃₄₅₆₇₈₉', SUPS = '⁰¹²³⁴⁵⁶⁷⁸⁹';
    const sub = x => String(x).split('').map(d => SUBS[+d]).join('');
    const sup = x => String(x).split('').map(d => SUPS[+d]).join('');
    function plainGroup(inf, orders) {
        const parts = [];
        if (inf) parts.push('ℤ');
        const cnt = new Map(); orders.forEach(o => cnt.set(o, (cnt.get(o) || 0) + 1));
        for (const [o, c] of cnt) parts.push(c === 1 ? `ℤ/${o}` : `(ℤ/${o})${sup(c)}`);
        return parts.length ? parts.join(' ⊕ ') : '0';
    }

    function stableRow(k) {             // for the stable stems list
        const S = summarize({ inf: k === 0 ? 1 : 0, P: stableKnowledge(k).P });
        return S;
    }

    root.HOS = { REFS, short, analyze, knowledge, summarize, cell, cellTitle, groupTeX, stableRow, stableAlts, STABLE_MAX, parseStr, ppart, invariantFactors, shortText };
})(typeof window !== 'undefined' ? window : globalThis);
