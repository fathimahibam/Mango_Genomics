/**
 * Mango Bioinformatics Lab — Core Application Logic
 *
 * Handles tab navigation, dashboard initialization (including 12-cultivar
 * resistance heatmap), gene explorer, pagination, and shared utility
 * functions for the web application.
 *
 * Data source: REAL NCBI GCF_011075055.1 (CATAS_Mindica_2.1)
 */

/* ───────────────────────────────────────────
 * Global State
 * ─────────────────────────────────────────── */

window.MangoApp = {
    diseases: [],
    cultivars: [],
    currentDisease: null,
    gwasData: null
};

let currentQuery = '';

/* ───────────────────────────────────────────
 * Utility Functions
 * ─────────────────────────────────────────── */

function formatNumber(n) {
    if (n == null || isNaN(n)) return '0';
    return Number(n).toLocaleString('en-US');
}

function showNotification(message, type = 'info') {
    const colorMap = {
        info:    { bg: 'rgba(59,130,246,0.92)',  icon: 'ℹ️' },
        success: { bg: 'rgba(34,197,94,0.92)',   icon: '✅' },
        error:   { bg: 'rgba(239,68,68,0.92)',   icon: '❌' },
        warning: { bg: 'rgba(245,158,11,0.92)',  icon: '⚠️' }
    };
    const cfg = colorMap[type] || colorMap.info;

    const toast = document.createElement('div');
    toast.className = 'toast-notification';
    toast.style.cssText = `
        position: fixed; top: 24px; right: 24px; z-index: 10000;
        padding: 14px 24px; border-radius: 12px;
        background: ${cfg.bg}; color: #fff;
        font-family: 'Inter', sans-serif; font-size: 14px;
        box-shadow: 0 8px 32px rgba(0,0,0,0.25);
        backdrop-filter: blur(12px);
        display: flex; align-items: center; gap: 10px;
        transform: translateX(120%); transition: transform 0.4s cubic-bezier(.23,1,.32,1);
        max-width: 400px;
    `;
    toast.innerHTML = `<span>${cfg.icon}</span><span>${message}</span>`;
    document.body.appendChild(toast);

    requestAnimationFrame(() => {
        toast.style.transform = 'translateX(0)';
    });

    setTimeout(() => {
        toast.style.transform = 'translateX(120%)';
        setTimeout(() => toast.remove(), 400);
    }, 3000);
}

function showLoading(element) {
    if (!element) return;
    if (element.querySelector('.mango-spinner')) return;

    const spinner = document.createElement('div');
    spinner.className = 'mango-spinner';
    spinner.style.cssText = `
        display: flex; align-items: center; justify-content: center;
        padding: 32px; width: 100%;
    `;
    spinner.innerHTML = `
        <div style="
            width: 40px; height: 40px;
            border: 3px solid rgba(255,255,255,0.15);
            border-top-color: #14b8a6;
            border-radius: 50%;
            animation: mangoSpin 0.8s linear infinite;
        "></div>
    `;
    element.appendChild(spinner);

    if (!document.getElementById('mango-spin-style')) {
        const style = document.createElement('style');
        style.id = 'mango-spin-style';
        style.textContent = `@keyframes mangoSpin { to { transform: rotate(360deg); } }`;
        document.head.appendChild(style);
    }
}

function hideLoading(element) {
    if (!element) return;
    const spinner = element.querySelector('.mango-spinner');
    if (spinner) spinner.remove();
}

function escapeHtml(str) {
    if (!str) return '';
    const div = document.createElement('div');
    div.textContent = str;
    return div.innerHTML;
}

/* ───────────────────────────────────────────
 * Animated Counter
 * ─────────────────────────────────────────── */

function animateCounter(el, target, duration = 1500) {
    if (!el || target == null) return;
    const start = performance.now();
    target = Number(target);

    function easeOutExpo(t) {
        return t === 1 ? 1 : 1 - Math.pow(2, -10 * t);
    }

    function tick(now) {
        const elapsed = now - start;
        const progress = Math.min(elapsed / duration, 1);
        const value = Math.floor(easeOutExpo(progress) * target);
        el.textContent = formatNumber(value);

        if (progress < 1) {
            requestAnimationFrame(tick);
        } else {
            el.textContent = formatNumber(target);
        }
    }

    requestAnimationFrame(tick);
}

/* ───────────────────────────────────────────
 * Tab Navigation
 * ─────────────────────────────────────────── */

function initTabNavigation() {
    const tabs = document.querySelectorAll('.tab-btn');
    const sections = document.querySelectorAll('.content-section');

    tabs.forEach(tab => {
        tab.addEventListener('click', () => {
            const targetId = tab.getAttribute('data-section');

            tabs.forEach(t => {
                t.classList.remove('active');
                t.setAttribute('aria-selected', 'false');
            });
            tab.classList.add('active');
            tab.setAttribute('aria-selected', 'true');

            sections.forEach(section => {
                section.classList.remove('active');
                section.style.display = 'none';
            });

            const target = document.getElementById(targetId);
            if (target) {
                target.classList.add('active');
                target.style.display = 'block';
            }
        });
    });
}

/* ───────────────────────────────────────────
 * Dashboard Initialization
 * ─────────────────────────────────────────── */

async function initDashboard() {
    // ── Real NCBI Genome Statistics ──
    try {
        const statsRes = await fetch('/api/genome/stats');
        if (!statsRes.ok) throw new Error(`Stats API error: ${statsRes.status}`);
        const stats = await statsRes.json();

        const statMapping = [
            { id: 'stat-total-genes',   key: 'total_genes' },
            { id: 'stat-total-cds',     key: 'total_cds' },
            { id: 'stat-chromosomes',   key: 'chromosome_count' },
            { id: 'stat-pseudogenes',   key: 'total_pseudogenes' },
            { id: 'stat-mrna',          key: 'total_mrna' },
        ];

        statMapping.forEach(({ id, key }) => {
            const el = document.getElementById(id);
            if (el && stats[key] != null) {
                animateCounter(el, stats[key], 1500);
            }
        });

        // Update assembly info if present
        showNotification(`Loaded real genome: ${stats.total_genes?.toLocaleString()} genes`, 'success');
    } catch (err) {
        console.error('[Dashboard] Failed to load genome stats:', err);
        showNotification('Failed to load genome statistics.', 'error');
    }

    // ── Load 12 Cultivar Data → build heatmap ──
    try {
        const cultRes = await fetch('/api/cultivars');
        if (!cultRes.ok) throw new Error(`Cultivars API error: ${cultRes.status}`);
        const cultivars = await cultRes.json();
        window.MangoApp.cultivars = cultivars;
        renderCultivarHeatmap(cultivars);
    } catch (err) {
        console.error('[Dashboard] Failed to load cultivars:', err);
    }

    // ── Disease Data ──
    try {
        const diseasesRes = await fetch('/api/diseases');
        if (!diseasesRes.ok) throw new Error(`Diseases API error: ${diseasesRes.status}`);
        const data = await diseasesRes.json();
        const diseases = Array.isArray(data) ? data : (data.diseases || []);

        window.MangoApp.diseases = diseases;

        // Wire up "View Gene" buttons on dashboard
        document.querySelectorAll('.view-gene-btn').forEach(btn => {
            btn.addEventListener('click', () => {
                const diseaseAttr = btn.getAttribute('data-disease');
                const idMap = {
                    'anthracnose':     'anthracnose',
                    'powdery-mildew':  'powdery_mildew',
                    'bacterial-canker': 'bacterial_canker',
                    'malformation':    'mango_malformation',
                    'fruit-rot':       'fruit_rot',
                    'stem-end-rot':    'stem_end_rot',
                };
                const diseaseId = idMap[diseaseAttr] || diseaseAttr;
                if (typeof window.selectDiseaseForCrispr === 'function') {
                    window.selectDiseaseForCrispr(diseaseId);
                }
            });
        });
    } catch (err) {
        console.error('[Dashboard] Failed to load diseases:', err);
        showNotification('Failed to load disease data.', 'error');
    }
}

/* ───────────────────────────────────────────
 * 12-Cultivar Resistance Heatmap
 * ─────────────────────────────────────────── */

/**
 * Render the 12-cultivar disease resistance heatmap table.
 * @param {Array<Object>} cultivars - Array of cultivar profiles from /api/cultivars
 */
function renderCultivarHeatmap(cultivars) {
    const tbody = document.getElementById('cultivar-heatmap-body');
    if (!tbody) return;

    const diseaseKeys = [
        'anthracnose',
        'powdery_mildew',
        'bacterial_canker',
        'mango_malformation',
        'fruit_rot',
        'stem_end_rot',
    ];

    tbody.innerHTML = '';

    cultivars.forEach(cult => {
        const row = document.createElement('tr');
        row.className = 'cultivar-heatmap-row';

        // Name cell
        const nameTd = document.createElement('td');
        nameTd.className = 'cultivar-name-cell';
        nameTd.innerHTML = `
            <div class="cultivar-heatmap-name">${escapeHtml(cult.name)}</div>
            <div class="cultivar-heatmap-type">${escapeHtml(cult.type || '')}</div>
        `;
        row.appendChild(nameTd);

        // Origin cell
        const originTd = document.createElement('td');
        originTd.className = 'cultivar-origin-cell';
        originTd.textContent = cult.origin || '—';
        row.appendChild(originTd);

        // Disease susceptibility cells
        diseaseKeys.forEach(key => {
            const level = (cult.disease_profile || {})[key] || '—';
            const td = document.createElement('td');
            td.className = 'heatmap-cell';

            let dotClass = 'heatmap-cell-dot';
            let levelClass = '';
            if (level === 'High') levelClass = 'heatmap-high';
            else if (level === 'Moderate') levelClass = 'heatmap-mod';
            else if (level === 'Low') levelClass = 'heatmap-low';

            td.innerHTML = `<span class="${dotClass} ${levelClass}">${escapeHtml(level)}</span>`;
            row.appendChild(td);
        });

        // Resistance score bar
        const scoreTd = document.createElement('td');
        scoreTd.className = 'score-cell';
        const score = cult.resistance_score || 0;
        const scoreColor = score >= 70 ? '#22c55e' : score >= 40 ? '#f59e0b' : '#ef4444';
        scoreTd.innerHTML = `
            <div class="score-bar-wrap">
                <div class="score-bar" style="width:${score}%; background:${scoreColor};"></div>
            </div>
            <span class="score-label" style="color:${scoreColor};">${score}</span>
        `;
        row.appendChild(scoreTd);

        tbody.appendChild(row);
    });
}

/* ───────────────────────────────────────────
 * Gene Explorer
 * ─────────────────────────────────────────── */

function initExplorer() {
    const searchBtn = document.getElementById('gene-search-btn');
    const searchInput = document.getElementById('gene-search-input');

    if (searchBtn) {
        searchBtn.addEventListener('click', () => {
            const query = searchInput?.value?.trim();
            if (query) {
                currentQuery = query;
                searchGenes(query, 1);
            }
        });
    }

    if (searchInput) {
        searchInput.addEventListener('keydown', (e) => {
            if (e.key === 'Enter') {
                const query = searchInput.value.trim();
                if (query) {
                    currentQuery = query;
                    searchGenes(query, 1);
                }
            }
        });
    }

    loadDiseaseGenes();
}

async function searchGenes(query, page = 1) {
    const tbody = document.getElementById('gene-table-body');
    const resultsCount = document.getElementById('results-count');

    if (tbody) {
        tbody.innerHTML = '';
        showLoading(tbody.parentElement);
    }

    try {
        const url = `/api/genes/search?query=${encodeURIComponent(query)}&page=${page}&per_page=20`;
        const res = await fetch(url);
        if (!res.ok) throw new Error(`Search API error: ${res.status}`);
        const data = await res.json();

        if (tbody) {
            hideLoading(tbody.parentElement);
            const genes = data.results || data.genes || [];

            if (genes.length === 0) {
                tbody.innerHTML = `
                    <tr><td colspan="8" style="text-align:center; padding:24px; color:#94a3b8;">
                        No genes found for "<strong>${escapeHtml(query)}</strong>"
                        <br><small style="color:#64748b;">Try: kinase, chitinase, LOC123, MYB, WRKY, invertase, NBS</small>
                    </td></tr>
                `;
            } else {
                genes.forEach(gene => {
                    const row = document.createElement('tr');
                    const productFull = gene.product || gene.gene_type || '—';
                    const productTrunc = productFull.length > 55
                        ? productFull.substring(0, 55) + '…'
                        : productFull;

                    // Badge for protein_coding vs pseudogene etc.
                    const typeColor = gene.gene_type === 'protein_coding' ? '#22c55e'
                        : gene.gene_type === 'pseudogene' ? '#f59e0b' : '#64748b';

                    row.innerHTML = `
                        <td style="font-family: 'Fira Code', 'Courier New', monospace; font-size: 0.82rem; color:#fbbf24;">
                            ${escapeHtml(gene.gene_id)}
                        </td>
                        <td style="${gene.name ? 'font-weight:700;' : 'color:#64748b;'}">
                            ${escapeHtml(gene.name || '—')}
                        </td>
                        <td><span style="color:#93c5fd;">${escapeHtml(gene.chromosome || '')}</span></td>
                        <td>${formatNumber(gene.start)}</td>
                        <td>${formatNumber(gene.end)}</td>
                        <td style="text-align:center; font-weight:600;
                            color:${gene.strand === '+' ? '#6ee7b7' : '#fca5a5'};">
                            ${gene.strand || ''}
                        </td>
                        <td>${formatNumber(gene.length_bp)}</td>
                        <td title="${escapeHtml(productFull)}" style="max-width:220px; overflow:hidden; text-overflow:ellipsis; white-space:nowrap;">
                            <span style="font-size:0.75rem; padding:2px 8px; border-radius:8px;
                                background:${typeColor}22; color:${typeColor}; margin-right:6px;">
                                ${escapeHtml(gene.gene_type || '—')}
                            </span>
                            ${escapeHtml(productTrunc)}
                        </td>
                    `;
                    tbody.appendChild(row);
                });
            }
        }

        if (resultsCount) {
            const shownCount = (data.results || data.genes || []).length;
            resultsCount.textContent = `Showing ${shownCount} of ${formatNumber(data.total || 0)} results from NCBI real genome`;
        }

        renderPagination(data.page || 1, data.pages || data.total_pages || 1);

    } catch (err) {
        console.error('[Explorer] Gene search failed:', err);
        if (tbody) hideLoading(tbody.parentElement);
        showNotification('Gene search failed. Please try again.', 'error');
    }
}

function renderPagination(currentPage, totalPages) {
    const prevBtn = document.getElementById('prev-page');
    const nextBtn = document.getElementById('next-page');
    const pageInfo = document.getElementById('page-info');

    if (pageInfo) {
        pageInfo.textContent = `Page ${currentPage} of ${totalPages}`;
    }

    if (prevBtn) {
        prevBtn.disabled = currentPage <= 1;
        prevBtn.onclick = () => {
            if (currentPage > 1) searchGenes(currentQuery, currentPage - 1);
        };
    }

    if (nextBtn) {
        nextBtn.disabled = currentPage >= totalPages;
        nextBtn.onclick = () => {
            if (currentPage < totalPages) searchGenes(currentQuery, currentPage + 1);
        };
    }
}

function loadDiseaseGenes() {
    const tbody = document.getElementById('candidate-gene-body');
    if (!tbody) return;

    const diseases = window.MangoApp.diseases;

    if (!diseases || diseases.length === 0) {
        setTimeout(loadDiseaseGenes, 800);
        return;
    }

    tbody.innerHTML = '';

    diseases.forEach(disease => {
        const row = document.createElement('tr');
        row.style.cssText = 'cursor: pointer; transition: background 0.2s;';
        row.className = 'disease-gene-row';

        row.innerHTML = `
            <td style="font-family: 'Fira Code', monospace; font-size: 0.82rem; color: #fbbf24;">
                ${escapeHtml(disease.susceptibility_gene || '')}
            </td>
            <td style="font-weight: 700;">${escapeHtml(disease.gene_full_name || '')}</td>
            <td><span style="color:#93c5fd;">${escapeHtml(disease.chromosome || '')}</span></td>
            <td>${disease.start ? formatNumber(disease.start) : '—'}</td>
            <td>${disease.end ? formatNumber(disease.end) : '—'}</td>
            <td style="text-align:center; color: #6ee7b7;">+</td>
            <td>${disease.start && disease.end ? formatNumber(disease.end - disease.start) : '—'}</td>
            <td>
                <span style="
                    padding: 4px 10px; border-radius: 12px; font-size: 0.78rem;
                    background: rgba(245,158,11,0.15); color: #fbbf24;
                ">${escapeHtml(disease.name)}</span>
            </td>
        `;

        row.addEventListener('click', () => {
            if (typeof window.selectDiseaseForCrispr === 'function') {
                window.selectDiseaseForCrispr(disease.id);
            }
        });

        row.addEventListener('mouseenter', () => { row.style.background = 'rgba(245,158,11,0.08)'; });
        row.addEventListener('mouseleave', () => { row.style.background = ''; });

        tbody.appendChild(row);
    });
}

/* ───────────────────────────────────────────
 * CSS for Cultivar Heatmap (injected)
 * ─────────────────────────────────────────── */

function injectHeatmapStyles() {
    if (document.getElementById('heatmap-inline-styles')) return;
    const style = document.createElement('style');
    style.id = 'heatmap-inline-styles';
    style.textContent = `
        .ncbi-source-banner {
            display: flex; align-items: flex-start; gap: 10px;
            background: rgba(99,102,241,0.1); border: 1px solid rgba(99,102,241,0.25);
            border-radius: 10px; padding: 12px 16px; margin-bottom: 20px;
            font-size: 0.88rem; color: #a5b4fc; line-height: 1.5;
        }
        .ncbi-icon { font-size: 1.2rem; flex-shrink: 0; }
        .ncbi-badge {
            display: flex; align-items: center; gap: 8px;
            background: rgba(99,102,241,0.12); border: 1px solid rgba(99,102,241,0.3);
            border-radius: 20px; padding: 4px 14px; font-size: 0.78rem;
            color: #a5b4fc; white-space: nowrap; margin-right: 12px;
        }
        .ncbi-dot {
            width: 8px; height: 8px; border-radius: 50%;
            background: #6366f1; animation: pulse 2s infinite;
        }
        @keyframes pulse { 0%,100%{opacity:1} 50%{opacity:0.4} }
        .ncbi-info-card { margin-bottom: 20px; }
        .cultivar-heatmap-card { overflow: hidden; }
        .heatmap-header {
            display: flex; gap: 20px; padding: 12px 0; margin-bottom: 8px;
            font-size: 0.82rem; flex-wrap: wrap;
        }
        .heatmap-legend-item { display: flex; align-items: center; gap: 6px; }
        .heatmap-dot { width: 12px; height: 12px; border-radius: 3px; }
        .heatmap-dot-high { background: rgba(239,68,68,0.8); }
        .heatmap-dot-mod  { background: rgba(245,158,11,0.8); }
        .heatmap-dot-low  { background: rgba(34,197,94,0.8); }
        .cultivar-heatmap-table-wrap { overflow-x: auto; }
        .cultivar-heatmap-table {
            width: 100%; border-collapse: collapse; font-size: 0.82rem;
        }
        .cultivar-heatmap-table th {
            padding: 10px 12px; text-align: center;
            background: rgba(255,255,255,0.04);
            border-bottom: 1px solid rgba(255,255,255,0.08);
            color: #94a3b8; font-weight: 600; white-space: nowrap;
        }
        .cultivar-heatmap-table th.cultivar-col { text-align: left; min-width: 160px; }
        .cultivar-heatmap-table th.origin-col { text-align: left; min-width: 120px; }
        .cultivar-heatmap-table th.score-col { min-width: 140px; }
        .cultivar-heatmap-row { transition: background 0.2s; }
        .cultivar-heatmap-row:hover { background: rgba(255,255,255,0.03); }
        .cultivar-heatmap-row td {
            padding: 10px 12px; border-bottom: 1px solid rgba(255,255,255,0.04);
            text-align: center; vertical-align: middle;
        }
        .cultivar-name-cell { text-align: left !important; }
        .cultivar-heatmap-name { font-weight: 700; color: #f1f5f9; font-size: 0.88rem; }
        .cultivar-heatmap-type { font-size: 0.72rem; color: #64748b; margin-top: 2px; }
        .cultivar-origin-cell { text-align: left !important; font-size: 0.76rem; color: #94a3b8; }
        .heatmap-cell { text-align: center; }
        .heatmap-cell-dot {
            display: inline-block; padding: 3px 10px;
            border-radius: 20px; font-size: 0.75rem; font-weight: 600;
        }
        .heatmap-high  { background: rgba(239,68,68,0.18); color: #f87171; }
        .heatmap-mod   { background: rgba(245,158,11,0.18); color: #fbbf24; }
        .heatmap-low   { background: rgba(34,197,94,0.18); color: #4ade80; }
        .score-cell { min-width: 120px; }
        .score-bar-wrap {
            width: 80px; height: 6px; background: rgba(255,255,255,0.1);
            border-radius: 3px; overflow: hidden; display: inline-block;
            margin-right: 8px; vertical-align: middle;
        }
        .score-bar { height: 100%; border-radius: 3px; transition: width 1s ease; }
        .score-label { font-size: 0.8rem; font-weight: 700; vertical-align: middle; }
        .gwas-stats-row {
            display: flex; gap: 24px; flex-wrap: wrap;
            padding: 12px 0; margin: 8px 0;
        }
        .gwas-stat { text-align: center; }
        .gwas-stat-num { display: block; font-size: 1.4rem; font-weight: 800; color: #fbbf24; }
        .gwas-stat-label { font-size: 0.72rem; color: #94a3b8; text-transform: uppercase; letter-spacing: 0.05em; }
        .cultivar-af-card .card-desc {
            color: #94a3b8; font-size: 0.88rem; margin-bottom: 16px;
        }
        .af-chart-container { width: 100%; max-height: 480px; overflow-y: auto; overflow-x: auto; border: 1px solid rgba(255,255,255,0.06); border-radius: 8px; }
        .af-chart-container canvas { max-width: 100%; }
        .af-chart-legend { display: flex; gap: 20px; flex-wrap: wrap; margin-top: 10px; font-size: 0.82rem; }
        .bar-fruitrot { background: linear-gradient(135deg, #f97316, #ea580c); }
        .bar-stemend  { background: linear-gradient(135deg, #8b5cf6, #7c3aed); }
    `;
    document.head.appendChild(style);
}

/* ───────────────────────────────────────────
 * Initialization
 * ─────────────────────────────────────────── */

document.addEventListener('DOMContentLoaded', () => {
    injectHeatmapStyles();
    initTabNavigation();
    initDashboard();
    initExplorer();
});
