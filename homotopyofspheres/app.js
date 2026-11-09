/* UI for the homotopy groups of spheres page. Logic lives in engine.js (window.HOS). */
(() => {
    'use strict';
    const H = window.HOS;
    const $ = id => document.getElementById(id);
    const tex = (s, display = false) => katex.renderToString(s, { throwOnError: false, displayMode: display });
    // Render "text with $tex$" (our own strings only).
    const rich = s => s.split(/(\$[^$]+\$)/g).map(part => part.length > 1 && part[0] === '$' && part[part.length - 1] === '$' ? tex(part.slice(1, -1)) : part).join('');
    document.querySelectorAll('.tex').forEach(el => { el.innerHTML = tex(el.dataset.tex); });

    document.querySelector('.menu-toggle').addEventListener('click', () => document.querySelector('.site-nav').classList.toggle('open'));

    // ------------------------------------------------------------------ references
    const linkOf = r => r.doi ? 'https://doi.org/' + r.doi : r.url || null;
    function refItem(r, id, locs) {
        const link = linkOf(r);
        return `<li><a href="#ref-${id}" class="rl" data-ref="${id}"><b>${H.short(id)}</b></a>` +
            (locs.length ? `, <span class="loc">${locs.join('; ')}</span>` : '') +
            ` — <i class="t">${r.t}</i>` +
            (link ? `<a class="ext" href="${link}" target="_blank" rel="noopener" aria-label="Open source"><i class="fas fa-arrow-up-right-from-square"></i></a>` : '') + '</li>';
    }
    function renderRefs(refs) {
        const order = [], locs = new Map();
        for (const { id, loc } of refs) {
            if (!H.REFS[id]) continue;
            if (!locs.has(id)) { locs.set(id, []); order.push(id); }
            if (loc && !locs.get(id).includes(loc)) locs.get(id).push(loc);
        }
        if (!order.length) return '';
        return '<ul class="fact-refs">' + order.map(id => refItem(H.REFS[id], id, locs.get(id))).join('') + '</ul>';
    }
    (function renderBibliography() {
        const groups = [['books', 'Books and surveys'], ['classical', 'Classical theorems'], ['unstable', 'Computations of unstable groups'], ['stable', 'Stable stems']];
        let html = '';
        for (const [g, name] of groups) {
            const ids = Object.keys(H.REFS).filter(id => H.REFS[id].g === g).sort((a, b) => H.REFS[a].y - H.REFS[b].y);
            html += `<h3>${name}</h3><ol>` + ids.map(id => {
                const r = H.REFS[id], link = linkOf(r);
                const linkTxt = r.doi ? 'doi:' + r.doi : (r.url || '').replace(/^https?:\/\//, '');
                return `<li id="ref-${id}"><span class="key">${H.short(id)}</span>${r.a}, <span class="t">${r.t}</span>, ${r.v}.` +
                    (link ? ` <a href="${link}" target="_blank" rel="noopener">${linkTxt}</a>` : '') + '</li>';
            }).join('') + '</ol>';
        }
        $('ref-list').innerHTML = html;
    })();
    document.addEventListener('click', e => {           // flash the target reference
        const a = e.target.closest('a.rl');
        if (!a) return;
        const li = document.getElementById('ref-' + a.dataset.ref);
        if (!li) return;
        e.preventDefault();                               // keep #n=..&k=.. in the address bar
        li.scrollIntoView({ behavior: 'smooth', block: 'center' });
        li.classList.add('flash'); setTimeout(() => li.classList.remove('flash'), 1800);
    });

    // ------------------------------------------------------------------ calculator
    const inN = $('in-n'), inK = $('in-k'), supN = $('sup-n');
    const fit = el => { el.style.width = Math.max(1.3, el.value.length + 0.35) + 'ch'; };
    function sanitize(el, allowMinus) {
        let v = el.value.replace(allowMinus ? /[^\d-]/g : /\D/g, '');
        if (allowMinus) v = v.replace(/(?!^)-/g, '');
        if (v !== el.value) el.value = v;
        fit(el);
    }
    function mirror() {
        supN.textContent = inN.value || 'n';
        supN.classList.toggle('empty', !inN.value);
        supN.classList.toggle('live', !!inN.value);
    }
    function setInputs(n, k) {
        inN.value = String(n); inK.value = String(k);
        fit(inN); fit(inK); mirror(); run();
    }
    inN.addEventListener('input', () => { sanitize(inN, false); mirror(); run(); });
    inK.addEventListener('input', () => { sanitize(inK, true); run(); });
    inN.addEventListener('keydown', e => { if (e.key === '+' || (e.key === 'ArrowRight' && inN.selectionStart === inN.value.length)) { e.preventDefault(); inK.focus(); } });
    inK.addEventListener('keydown', e => { if ((e.key === 'Backspace' || e.key === 'ArrowLeft') && inK.selectionStart === 0 && inK.selectionEnd === 0) { e.preventDefault(); inN.focus(); } });

    let lastKey = '';
    function run() {
        const n = parseInt(inN.value, 10), k = parseInt(inK.value, 10);
        const result = $('result'), details = $('details'), facts = $('facts');
        if (!(n >= 1) || isNaN(k) || inK.value === '-') {
            result.className = 'result hint';
            result.textContent = !inN.value && !inK.value ? 'Type n and k, or click a cell of the table below.' : (n === 0 ? 'The sphere dimension must be at least 1.' : 'Type both n and k.');
            details.textContent = ''; facts.innerHTML = ''; markCell(null); return;
        }
        if (n + k < 1) {
            result.className = 'result hint'; result.textContent = 'The total degree n + k must be at least 1.';
            details.textContent = ''; facts.innerHTML = ''; markCell(null); return;
        }
        const key = n + ',' + k;
        if (key === lastKey) return;
        lastKey = key;
        const r = H.analyze(n, k);
        result.className = 'result';
        result.innerHTML = tex(r.headline, true);
        const bits = [];
        if (r.order) bits.push('order ' + r.order.toLocaleString('en-US'));
        if (r.primaryTeX) bits.push(tex(r.primaryTeX));
        if (r.exact === false) bits.push('partially known: each $G_{(p)}$ below is a $p$-primary part that the sources do not fully determine');
        if (r.openNote) bits.push(r.openNote);
        if (k > 0 && n > 1) bits.push(n >= k + 2 ? 'stable' : (k === n - 1 ? 'edge of the stable range' : 'unstable'));
        details.innerHTML = rich(bits.join(' · '));
        facts.innerHTML = r.cards.map(c =>
            `<div class="fact ${c.cls}"><div class="fact-head"><i class="fas ${c.icon}"></i><span class="title">${c.title}</span><span class="badge">${c.badge}</span></div>` +
            `<div class="fact-body">${rich(c.body)}</div>${renderRefs(c.refs)}</div>`).join('');
        markCell(n, k);
        const h = `#n=${n}&k=${k}`;
        if (location.hash !== h) history.replaceState(null, '', h);
    }

    // examples
    const EX = [[6, 33, '33-stem, settled 2025'], [72, 70, 'stable, corrected 2024'], [128, 126, 'Kervaire θ₆'], [12, 11, 'ℤ ⊕ im J'], [24, 40, 'only the order is bounded'], [3, 52, 'deep in the unstable range']];
    $('examples').innerHTML = EX.map(([n, k, lab]) => `<button class="chip" data-n="${n}" data-k="${k}">${tex(`\\pi_{${n + k}}(S^{${n}})`)} <small>${lab}</small></button>`).join('');
    $('examples').addEventListener('click', e => { const b = e.target.closest('.chip'); if (b) setInputs(+b.dataset.n, +b.dataset.k); });

    // ------------------------------------------------------------------ table
    const VIEWS = {
        full: { kmax: 33, note: 'Whole groups. Striped cells show the known part, and ? marks a p-primary part that the sources do not fully determine (5- and 7-primary parts on small spheres in stems 22–33, and the 2-primary 33-stem for 10 ≤ n ≤ 27).' },
        2: { kmax: 50, note: '2-primary components. Exact through stem 32, and in stem 33 for n ≤ 9 and n ≥ 28. Beyond that, gold striped cells give only the order 2^e, or bounds 2^a–2^b, from Baer’s 2026 machine computation.' },
        3: { kmax: 55, note: '3-primary components: Toda (2003) for odd n, Serre’s splitting for even n, and stable stems from Isaksen–Wang–Xu and Ravenel.' }
    };
    let view = 'full';
    function renderTable() {
        const { kmax, note } = VIEWS[view];
        $('table-note').textContent = note;
        const nmax = kmax + 2;
        let h = '<thead><tr><th class="corner">k \\ n</th>';
        for (let n = 1; n <= nmax; n++) h += `<th>${n}</th>`;
        h += '</tr></thead><tbody>';
        for (let k = 0; k <= kmax; k++) {
            h += `<tr><th>${k}</th>`;
            for (let n = 1; n <= nmax; n++) {
                if (n > k + 2) { h += '<td class="rep"></td>'; continue; }
                const c = H.cell(n, k, view);
                const len = c.text.replace(/<[^>]+>/g, '').length;
                const txt = len > 9 ? c.text.replace(/·/g, '·<wbr>') : c.text;
                h += `<td class="${c.cls}${len > 9 ? ' long' : ''}" data-n="${n}" data-k="${k}" title="${H.cellTitle(n, k, view).replace(/"/g, '&quot;')}">${txt}</td>`;
            }
            h += '</tr>';
        }
        $('hos-table').innerHTML = h + '</tbody>';
        drawLines();
        const n = parseInt(inN.value, 10), k = parseInt(inK.value, 10);
        if (n >= 1 && !isNaN(k)) markCell(n, k);
    }
    function drawLines() {
        const t = $('hos-table'), svg = $('stab-svg');
        const { kmax } = VIEWS[view];
        const td = (n, k) => t.querySelector(`td[data-n="${n}"][data-k="${k}"]`);
        const a = td(2, 0), b = td(kmax + 2, kmax);
        if (!a || !b) { svg.innerHTML = ''; return; }
        svg.setAttribute('width', t.offsetWidth); svg.setAttribute('height', t.offsetHeight);
        // staircase along the left edges of the first stable cells n = k + 2: one column per row, i.e. slope one
        let d = '';
        for (let k = 0; k <= kmax; k++) {
            const c = td(k + 2, k);
            if (!c) continue;
            const x = c.offsetLeft, y0 = c.offsetTop, y1 = y0 + c.offsetHeight;
            d += (d ? ` L${x},${y0}` : `M${x},${y0}`) + ` L${x},${y1}`;
        }
        const x2 = b.offsetLeft, y2 = b.offsetTop + b.offsetHeight;
        d += ` L${x2 + b.offsetWidth},${y2}`;
        svg.innerHTML = `<path d="${d}"/><text x="${Math.min(x2 + b.offsetWidth + 8, t.offsetWidth - 170)}" y="${y2 - 6}">n = k + 2: stable from here →</text>`;
    }
    function markCell(n, k) {
        document.querySelectorAll('#hos-table td.sel').forEach(c => c.classList.remove('sel'));
        if (n == null) return;
        const c = $('hos-table').querySelector(`td[data-n="${n}"][data-k="${k}"]`);
        if (c) c.classList.add('sel');
    }
    $('hos-table').addEventListener('click', e => {
        const c = e.target.closest('td[data-n]');
        if (!c) return;
        setInputs(+c.dataset.n, +c.dataset.k);
        $('calculator').scrollIntoView({ behavior: 'smooth', block: 'start' });
    });
    document.querySelectorAll('.seg').forEach(b => b.addEventListener('click', () => {
        view = b.dataset.view;
        document.querySelectorAll('.seg').forEach(x => x.setAttribute('aria-pressed', String(x === b)));
        renderTable();
    }));
    window.addEventListener('resize', () => { clearTimeout(drawLines.t); drawLines.t = setTimeout(drawLines, 120); });
    if (document.fonts && document.fonts.ready) document.fonts.ready.then(drawLines);

    // ------------------------------------------------------------------ stable stems list
    (function renderStable() {
        let h = '<thead><tr><th>k</th><th>π<sup>s</sup><sub>k</sub></th><th>Note</th></tr></thead><tbody>';
        for (let k = 0; k <= H.STABLE_MAX; k++) {
            let g, note = '';
            if (k === 0) g = tex('\\mathbb{Z}');
            else {
                const K = H.knowledge(k + 2, k), S = H.summarize(K);
                if (S.exact) g = tex(H.groupTeX(0, S.inv));
                else {
                    const alts = K.P[2].alts;
                    const odd = Object.keys(K.P).map(Number).filter(p => p > 2).map(p => K.P[p].orders);
                    g = alts.map(a => tex(H.groupTeX(0, H.invariantFactors([a].concat(odd))))).join(' <span class="muted">or</span> ');
                    note = '<span class="flag">open: 2-primary part not settled</span>';
                }
            }
            if (k === 70 || k === 71) note = '<span class="flag">corrected 2024 (Burklund–Isaksen–Xu)</span>';
            if (k === 82 || k === 83) note = 'settled 2024 (Burklund–Isaksen–Xu)';
            const kerv = { 2: 'θ₁ = η²', 6: 'θ₂ = ν²', 14: 'θ₃ = σ²', 30: 'θ₄', 62: 'θ₅' }[k];
            if (kerv) note = (note ? note + ' · ' : '') + 'Kervaire ' + kerv;
            h += `<tr data-k="${k}"><td>${k}</td><td>${g}</td><td>${note}</td></tr>`;
        }
        $('stable-list').innerHTML = h + '</tbody>';
        $('stable-list').addEventListener('click', e => {
            const tr = e.target.closest('tr[data-k]');
            if (!tr) return;
            const k = +tr.dataset.k;
            setInputs(Math.max(k + 2, 2), k);
            $('calculator').scrollIntoView({ behavior: 'smooth', block: 'start' });
        });
    })();

    // ------------------------------------------------------------------ start
    renderTable();
    function fromHash() {
        const m = location.hash.match(/n=(\d+)&k=(-?\d+)/);
        if (m) setInputs(+m[1], +m[2]); else { fit(inN); fit(inK); mirror(); }
    }
    window.addEventListener('hashchange', () => { lastKey = ''; fromHash(); });
    fromHash();
})();
