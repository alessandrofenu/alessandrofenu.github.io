<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Alessandro Fenu - New address</title>
    <!-- Fallback redirect if JavaScript is disabled (the script below redirects first, keeping any #section) -->
    <meta http-equiv="refresh" content="12;url=https://alefenu.com/">
    <link rel="canonical" href="https://alefenu.com/">
    <link href="https://fonts.googleapis.com/css2?family=Inter:wght@400;700&family=Oswald:wght@700&display=swap" rel="stylesheet">
    <style>
        /* Same palette as the main site */
        :root {
            --bg: #232223;
            --card: #363C3D;
            --text: #cccccc;
            --teal: #94A3A1;
            --accent: #847C6D;
            --track: #4A5152;
        }

        * { box-sizing: border-box; }

        body {
            margin: 0;
            min-height: 100vh;
            display: flex;
            align-items: center;
            justify-content: center;
            padding: 24px 16px;
            background-color: var(--bg);
            color: var(--text);
            font-family: 'Inter', system-ui, sans-serif;
            line-height: 1.7;
        }

        .card {
            width: 100%;
            max-width: 640px;
            background-color: var(--card);
            border-radius: 8px;
            box-shadow: 0 1px 5px rgba(0, 0, 0, 0.3);
            padding: 40px;
        }

        h1 {
            font-family: 'Oswald', sans-serif;
            text-transform: uppercase;
            letter-spacing: 0.5px;
            color: var(--accent);
            font-size: 2.2em;
            margin: 0 0 20px;
            padding-bottom: 5px;
            border-bottom: 2px solid var(--teal);
        }

        p { margin: 0 0 6px; }
        .lead { font-size: 1.1em; margin-bottom: 14px; }
        .translation { color: var(--teal); font-size: 0.95em; }

        a { color: var(--accent); text-decoration: none; transition: color 0.3s ease; }
        a:hover { color: var(--teal); }

        /* The path from the old address to the new one */
        svg { display: block; width: 100%; height: auto; margin: 32px 0 6px; overflow: visible; }
        .track { fill: none; stroke: var(--track); stroke-width: 3; stroke-dasharray: 2 8; stroke-linecap: round; }
        .trail {
            fill: none; stroke: var(--accent); stroke-width: 3; stroke-linecap: round;
            stroke-dasharray: 1; stroke-dashoffset: 1;
            animation: draw 10s linear forwards;
        }
        @keyframes draw { to { stroke-dashoffset: 0; } }
        .start { fill: var(--card); stroke: var(--teal); stroke-width: 3; }
        .end { fill: var(--teal); }
        .dot { fill: var(--accent); }

        .labels { display: flex; flex-wrap: wrap; justify-content: space-between; gap: 2px 12px; font-size: 0.85em; }
        .old { color: var(--teal); text-decoration: line-through; opacity: 0.8; }
        .new { font-weight: 700; margin-left: auto; }

        .actions {
            display: flex;
            flex-wrap: wrap;
            align-items: center;
            justify-content: space-between;
            gap: 16px;
            margin-top: 32px;
        }
        .countdown { color: var(--teal); }
        .button {
            font-family: 'Oswald', sans-serif;
            text-transform: uppercase;
            letter-spacing: 1px;
            background-color: var(--accent);
            color: var(--bg);
            padding: 10px 22px;
            border-radius: 6px;
            transition: background-color 0.3s ease;
        }
        .button:hover { background-color: var(--teal); color: var(--bg); }

        @media (prefers-reduced-motion: reduce) {
            .trail { animation: none; stroke-dashoffset: 0; }
        }

        @media (max-width: 480px) {
            .card { padding: 28px 20px; }
            h1 { font-size: 1.8em; }
        }
    </style>
</head>
<body>
    <main class="card">
        <h1>I've moved</h1>
        <p class="lead">This website has a new address: <a id="new-link" href="https://alefenu.com/">alefenu.com</a>. Up to homotopy, nothing else has changed.</p>
        <p class="translation" lang="it">Il sito ha cambiato indirizzo. A meno di omotopia, non è cambiato nient'altro.</p>
        <p class="translation" lang="fr">Le site a changé d'adresse. À homotopie près, rien d'autre n'a changé.</p>

        <svg viewBox="0 0 600 110" role="img" aria-label="A path from the old address to the new one">
            <path id="route" class="track" d="M 12 80 C 140 -10, 250 140, 370 55 S 530 20, 588 42"/>
            <path class="trail" pathLength="1" d="M 12 80 C 140 -10, 250 140, 370 55 S 530 20, 588 42"/>
            <circle class="start" cx="12" cy="80" r="7"/>
            <circle class="end" cx="588" cy="42" r="7"/>
            <circle id="dot" class="dot" r="6">
                <animateMotion dur="10s" fill="freeze" calcMode="paced">
                    <mpath href="#route"/>
                </animateMotion>
            </circle>
        </svg>
        <div class="labels">
            <span class="old">poisson.phc.dm.unipi.it/~afenu/</span>
            <span class="new">alefenu.com</span>
        </div>

        <div class="actions">
            <span class="countdown">Redirecting in <span id="seconds">10</span> s</span>
            <a class="button" id="go-now" href="https://alefenu.com/">Go now &rarr;</a>
        </div>
    </main>

    <script>
        // Keep the section of old links (e.g. .../tutto/#tutorato) when redirecting
        const target = 'https://alefenu.com/' + window.location.hash;
        document.getElementById('new-link').href = target;
        document.getElementById('go-now').href = target;

        // Without animations, just show the dot at the arrival point
        if (window.matchMedia('(prefers-reduced-motion: reduce)').matches) {
            const dot = document.getElementById('dot');
            dot.querySelector('animateMotion').remove();
            dot.setAttribute('cx', 588);
            dot.setAttribute('cy', 42);
        }

        let seconds = 10;
        const secondsEl = document.getElementById('seconds');
        const timer = setInterval(() => {
            seconds -= 1;
            secondsEl.textContent = seconds;
            if (seconds <= 0) {
                clearInterval(timer);
                window.location.replace(target);
            }
        }, 1000);
    </script>
</body>
</html>
