---
layout: default
title: Alessandro Fenu | Projects
---

<style>
    :root {
        --bg-color: #0f172a;
        --card-bg: #1e293b;
        --text-main: #f8fafc;
        --text-muted: #94a3b8;
        --accent: #38bdf8;
        --accent-gradient: linear-gradient(135deg, #38bdf8 0%, #818cf8 100%);
    }

    /* Target the GitHub Pages container if a theme is applied, or apply globally */
    body, .site-content, .wrapper {
        background-color: var(--bg-color) !important;
        color: var(--text-main) !important;
        font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
    }

    .homepage-container {
        display: flex;
        flex-direction: column;
        justify-content: center;
        align-items: center;
        padding: 2rem 1rem;
        min-height: 80vh;
    }

    header.custom-header {
        text-align: center;
        margin-bottom: 3rem;
        max-width: 600px;
    }

    header.custom-header h1 {
        font-size: 2.5rem;
        font-weight: 800;
        background: var(--accent-gradient);
        -webkit-background-clip: text;
        -webkit-text-fill-color: transparent;
        margin-bottom: 0.5rem;
        border: none !important;
        padding: 0 !important;
    }

    header.custom-header p {
        color: var(--text-muted);
        font-size: 1.1rem;
    }

    .grid {
        display: grid;
        grid-template-columns: repeat(auto-fit, minmax(280px, 1fr));
        gap: 1.5rem;
        width: 100%;
        max-width: 1000px;
    }

    .card {
        background-color: var(--card-bg);
        border: 1px solid #334155;
        border-radius: 12px;
        padding: 1.75rem;
        text-decoration: none;
        color: var(--text-main) !important;
        display: flex;
        flex-direction: column;
        justify-content: space-between;
        transition: transform 0.2s ease, border-color 0.2s ease, box-shadow 0.2s ease;
    }

    .card:hover {
        transform: translateY(-4px);
        border-color: var(--accent);
        box-shadow: 0 10px 25px -5px rgba(56, 189, 248, 0.15);
        text-decoration: none !important;
    }

    .card h2 {
        font-size: 1.35rem;
        font-weight: 700;
        margin-top: 0 !important;
        margin-bottom: 0.75rem;
        text-transform: capitalize;
        display: flex;
        align-items: center;
        justify-content: space-between;
        color: var(--text-main) !important;
        border: none !important;
    }

    .card h2::after {
        content: "→";
        font-size: 1.1rem;
        opacity: 0;
        transform: translateX(-5px);
        transition: all 0.2s ease;
        color: var(--accent);
    }

    .card:hover h2::after {
        opacity: 1;
        transform: translateX(0);
    }

    .card p {
        color: var(--text-muted);
        font-size: 0.95rem;
        line-height: 1.5;
        margin: 0;
    }

    footer.custom-footer {
        margin-top: 4rem;
        color: var(--text-muted);
        font-size: 0.85rem;
        text-align: center;
    }

    footer.custom-footer a {
        color: var(--accent);
        text-decoration: none;
    }
    
    footer.custom-footer a:hover {
        text-decoration: underline;
    }
</style>

<div class="homepage-container">

    <header class="custom-header">
        <h1>Alessandro Fenu</h1>
        <p>A collection of my web projects and deployments.</p>
    </header>

    <main class="grid">
        
        <a href="/mudare/" class="card">
            <div>
                <h2>Mudare</h2>
                <p>Explore the Mudare project dashboard and interactive interface.</p>
            </div>
        </a>

        <a href="/lotteria/" class="card">
            <div>
                <h2>Lotteria</h2>
                <p>A web-based lottery or serialization application engine.</p>
            </div>
        </a>

        <a href="/amsterdam/" class="card">
            <div>
                <h2>Amsterdam</h2>
                <p>Project hub focused on insights, mapping, or data related to Amsterdam.</p>
            </div>
        </a>

        <a href="/valencia/" class="card">
            <div>
                <h2>Valencia</h2>
                <p>Interactive web space or presentation dedicated to Valencia.</p>
            </div>
        </a>

        <a href="/visual-genealogy/" class="card">
            <div>
                <h2>Visual Genealogy</h2>
                <p>An interactive, visual data tool mapping out ancestral lineages and history.</p>
            </div>
        </a>

    </main>

    <footer class="custom-footer">
        <p>Built for speed. View main profile on <a href="https://github.com/alessandrofenu" target="_blank" rel="noopener noreferrer">GitHub</a>.</p>
    </footer>

</div>
