/**
 * Mango Bioinformatics Lab — CRISPR-Cas9 Simulation Interface
 *
 * Provides gene selection, PAM site discovery, gRNA validation,
 * CRISPR simulation execution, and animated result visualization
 * including DNA cutting/repair animations.
 *
 * @author Mango Bioinformatics Lab
 * @version 1.0.0
 */

/**
 * @class CrisprSimulator
 * Manages the CRISPR-Cas9 simulation workflow: gene selection,
 * PAM site finding, gRNA input, simulation execution, and result display.
 */
class CrisprSimulator {
    constructor() {
        /** @type {string|null} Currently selected disease ID */
        this.currentDiseaseId = null;
    }

    /* ─────────────────────────────────
     * Initialization
     * ───────────────────────────────── */

    /**
     * Attach all event listeners for the CRISPR interface.
     * Uses card-based gene selection instead of dropdown.
     */
    init() {
        // Gene select cards (buttons with class "gene-select-card" and data-gene attr)
        const geneCards = document.querySelectorAll('.gene-select-card');
        geneCards.forEach(card => {
            card.addEventListener('click', () => {
                // Remove active from all cards
                geneCards.forEach(c => c.classList.remove('active'));
                card.classList.add('active');
                
                const geneAttr = card.getAttribute('data-gene') || '';
                // Map HTML data-gene values to backend disease IDs
                const idMap = {
                    'anthracnose':       'anthracnose',
                    'powdery_mildew':    'powdery_mildew',
                    'bacterial_canker':  'bacterial_canker',
                    'mango_malformation': 'mango_malformation',
                    'fruit_rot':         'fruit_rot',
                    'stem_end_rot':      'stem_end_rot',
                };
                const diseaseId = idMap[geneAttr] || geneAttr;
                this.selectGene(diseaseId);
            });
        });

        // Find PAM button - HTML ID: "find-pam-btn"
        const btnPam = document.getElementById('find-pam-btn');
        if (btnPam) {
            btnPam.addEventListener('click', () => this.findPamSites());
        }

        // Simulate button - HTML ID: "simulate-crispr-btn"
        const btnSim = document.getElementById('simulate-crispr-btn');
        if (btnSim) {
            btnSim.addEventListener('click', () => this.runSimulation());
        }

        // gRNA input character counter
        const grnaInput = document.getElementById('grna-input');
        const grnaLength = document.getElementById('grna-length');
        if (grnaInput && grnaLength) {
            grnaInput.addEventListener('input', () => {
                grnaLength.textContent = `${grnaInput.value.length}/20`;
            });
        }
    }

    /* ─────────────────────────────────
     * Gene Selection
     * ───────────────────────────────── */

    /**
     * Select a disease gene and display its information and sequences.
     * @param {string|number} diseaseId - The disease ID to select.
     */
    selectGene(diseaseId) {
        const disease = (window.MangoApp.diseases || []).find(
            d => String(d.id) === String(diseaseId)
        );

        if (!disease) {
            showNotification('Disease not found.', 'error');
            return;
        }

        this.currentDiseaseId = diseaseId;
        window.MangoApp.currentDisease = diseaseId;

        // ── Gene Info Panel ──
        // HTML uses gene-info-panel with inner fields crispr-gene-name, crispr-gene-chr, etc.
        const infoPanel = document.getElementById('gene-info-panel');
        if (infoPanel) {
            infoPanel.style.display = 'block';
            const setEl = (id, text) => {
                const el = document.getElementById(id);
                if (el) el.textContent = text;
            };
            setEl('crispr-gene-name', disease.gene_full_name || disease.susceptibility_gene || '—');
            setEl('crispr-gene-chr', disease.chromosome || '—');
            setEl('crispr-gene-coords', `${disease.start ? Number(disease.start).toLocaleString() : '—'} – ${disease.end ? Number(disease.end).toLocaleString() : '—'}`);
            setEl('crispr-gene-mutation', disease.mutation_description || '—');
        }

        // ── Sequences ──
        if (disease.susceptible_sequence) {
            this.colorizeSequence(disease.susceptible_sequence, 'susceptible-seq');
        }
        if (disease.resistant_sequence) {
            this.colorizeSequence(disease.resistant_sequence, 'resistant-seq');
        }
        if (disease.susceptible_sequence && disease.resistant_sequence) {
            this.showSequenceDiff(
                disease.susceptible_sequence,
                disease.resistant_sequence,
                'sequence-diff'
            );
        }

        // ── Clear previous results ──
        const pamList = document.getElementById('pam-sites-list');
        if (pamList) {
            pamList.style.display = 'none';
            const pamContainer = document.getElementById('pam-sites-container');
            if (pamContainer) pamContainer.innerHTML = '';
        }
        const results = document.getElementById('crispr-results');
        if (results) results.style.display = 'none';
        const grnaInput = document.getElementById('grna-input');
        if (grnaInput) grnaInput.value = '';
        const grnaLength = document.getElementById('grna-length');
        if (grnaLength) grnaLength.textContent = '0/20';

        // Enable Find PAM button
        const btnPam = document.getElementById('find-pam-btn');
        if (btnPam) btnPam.disabled = false;

        showNotification(`Selected: ${disease.name}`, 'info');
    }

    /* ─────────────────────────────────
     * Sequence Visualization
     * ───────────────────────────────── */

    /**
     * Render a DNA sequence with individually colored nucleotides.
     * Adds position markers every 10 characters.
     *
     * Nucleotide color scheme:
     * - A: #ff6b6b (red)
     * - T: #ffd93d (yellow)
     * - G: #6bcb77 (green)
     * - C: #4d96ff (blue)
     *
     * @param {string} sequence - The DNA sequence string (e.g. "ATCGATCG...").
     * @param {string} containerId - The id of the container element.
     */
    colorizeSequence(sequence, containerId) {
        const container = document.getElementById(containerId);
        if (!container || !sequence) return;

        const colorMap = {
            A: { color: '#ff6b6b', cls: 'nt-a' },
            T: { color: '#ffd93d', cls: 'nt-t' },
            G: { color: '#6bcb77', cls: 'nt-g' },
            C: { color: '#4d96ff', cls: 'nt-c' }
        };

        container.innerHTML = '';
        container.style.fontFamily = "'Fira Code', 'Courier New', monospace";
        container.style.fontSize = '13px';
        container.style.lineHeight = '1.8';
        container.style.letterSpacing = '1px';
        container.style.wordBreak = 'break-all';

        // Position number line
        const posLine = document.createElement('div');
        posLine.style.cssText = 'color:#64748b; font-size:10px; letter-spacing:1px; user-select:none;';
        let posText = '';
        for (let i = 0; i < sequence.length; i++) {
            if (i > 0 && i % 10 === 0) {
                const numStr = String(i);
                posText += numStr;
                // Pad to align with nucleotides
                i += numStr.length - 1;
            } else {
                posText += ' ';
            }
        }
        posLine.textContent = posText;
        container.appendChild(posLine);

        // Sequence line with colored spans
        const seqLine = document.createElement('div');
        for (let i = 0; i < sequence.length; i++) {
            const nt = sequence[i].toUpperCase();
            const span = document.createElement('span');
            span.textContent = nt;
            const mapping = colorMap[nt];
            if (mapping) {
                span.className = mapping.cls;
                span.style.color = mapping.color;
            } else {
                span.style.color = '#94a3b8';
            }

            // Add subtle separator every 10 chars
            if (i > 0 && i % 10 === 0) {
                span.style.marginLeft = '4px';
            }

            seqLine.appendChild(span);
        }
        container.appendChild(seqLine);
    }

    /**
     * Display a diff view of two sequences, highlighting mismatched positions.
     * Shows seq1 on top, a match/mismatch marker line, then seq2 below.
     *
     * @param {string} seq1 - First sequence (e.g. susceptible).
     * @param {string} seq2 - Second sequence (e.g. resistant).
     * @param {string} containerId - Container element id.
     */
    showSequenceDiff(seq1, seq2, containerId) {
        const container = document.getElementById(containerId);
        if (!container) return;

        container.innerHTML = '';
        container.style.fontFamily = "'Fira Code', 'Courier New', monospace";
        container.style.fontSize = '13px';
        container.style.lineHeight = '1.6';
        container.style.letterSpacing = '1px';

        const colorMap = {
            A: '#ff6b6b', T: '#ffd93d', G: '#6bcb77', C: '#4d96ff'
        };

        const maxLen = Math.max(seq1.length, seq2.length);

        // Seq1 line (top)
        const line1 = document.createElement('div');
        line1.style.marginBottom = '0';

        // Marker line
        const markerLine = document.createElement('div');
        markerLine.style.cssText = 'color:#64748b; margin:0;';

        // Seq2 line (bottom)
        const line2 = document.createElement('div');

        for (let i = 0; i < maxLen; i++) {
            const c1 = (seq1[i] || '-').toUpperCase();
            const c2 = (seq2[i] || '-').toUpperCase();
            const isMismatch = c1 !== c2;

            // Seq1 character
            const span1 = document.createElement('span');
            span1.textContent = c1;
            span1.style.color = colorMap[c1] || '#94a3b8';
            if (isMismatch) {
                span1.className = 'mutation-highlight';
                span1.style.background = 'rgba(245,158,11,0.3)';
                span1.style.borderRadius = '2px';
                span1.style.padding = '0 1px';
            }
            line1.appendChild(span1);

            // Marker
            const markerSpan = document.createElement('span');
            markerSpan.textContent = isMismatch ? '✱' : '│';
            markerSpan.style.color = isMismatch ? '#f59e0b' : '#334155';
            markerLine.appendChild(markerSpan);

            // Seq2 character
            const span2 = document.createElement('span');
            span2.textContent = c2;
            span2.style.color = colorMap[c2] || '#94a3b8';
            if (isMismatch) {
                span2.className = 'mutation-highlight';
                span2.style.background = 'rgba(245,158,11,0.3)';
                span2.style.borderRadius = '2px';
                span2.style.padding = '0 1px';
            }
            line2.appendChild(span2);

            // Spacer every 10 chars
            if (i > 0 && (i + 1) % 10 === 0) {
                [span1, markerSpan, span2].forEach(s => s.style.marginRight = '4px');
            }
        }

        // Labels
        const label1 = document.createElement('span');
        label1.textContent = ' Susceptible';
        label1.style.cssText = 'color:#fca5a5; font-size:10px; margin-left:8px;';
        line1.appendChild(label1);

        const label2 = document.createElement('span');
        label2.textContent = ' Resistant';
        label2.style.cssText = 'color:#6ee7b7; font-size:10px; margin-left:8px;';
        line2.appendChild(label2);

        container.appendChild(line1);
        container.appendChild(markerLine);
        container.appendChild(line2);
    }

    /* ─────────────────────────────────
     * PAM Site Discovery
     * ───────────────────────────────── */

    /**
     * Fetch PAM sites for the current disease gene and display them as selectable cards.
     * Clicking a card fills the gRNA input field.
     */
    async findPamSites() {
        if (!this.currentDiseaseId) {
            showNotification('Please select a gene first.', 'warning');
            return;
        }

        // HTML structure: pam-sites-list (outer, initially hidden) > h4 + pam-sites-container (inner)
        const pamList = document.getElementById('pam-sites-list');
        const pamContainer = document.getElementById('pam-sites-container');
        
        if (pamContainer) {
            pamContainer.innerHTML = '';
        }
        if (pamList) {
            pamList.style.display = 'block';
            showLoading(pamContainer || pamList);
        }

        try {
            const res = await fetch(`/api/crispr/find_pam?disease_id=${encodeURIComponent(this.currentDiseaseId)}`);
            if (!res.ok) throw new Error(`PAM API error: ${res.status}`);
            const data = await res.json();

            const targetContainer = pamContainer || pamList;
            if (targetContainer) {
                hideLoading(targetContainer);

                const sites = data.pam_sites || [];

                // Count badge
                const countBadge = document.createElement('div');
                countBadge.style.cssText = `
                    margin-bottom:12px; color:#94a3b8; font-size:0.85rem;
                `;
                countBadge.innerHTML = `Found <strong style="color:#fbbf24;">${sites.length}</strong> PAM sites in <em>${escapeHtml(data.gene || '')}</em> (${data.sequence_len || 0} bp)`;
                targetContainer.appendChild(countBadge);

                if (sites.length === 0) {
                    targetContainer.innerHTML += '<p style="color:#64748b; text-align:center; padding:16px;">No PAM sites found.</p>';
                    return;
                }

                sites.forEach((site, idx) => {
                    const card = document.createElement('div');
                    card.className = 'pam-site-card';
                    card.style.cssText = `
                        padding: 14px 18px; border-radius: 10px;
                        background: rgba(255,255,255,0.03);
                        border: 1px solid rgba(255,255,255,0.08);
                        margin-bottom: 8px; cursor: pointer;
                        transition: border-color 0.2s, background 0.2s, transform 0.15s;
                    `;

                    // Colorize gRNA inline
                    const grnaColored = this._colorizeInline(site.grna_target || '');

                    card.innerHTML = `
                        <div style="display:flex; align-items:center; justify-content:space-between; flex-wrap:wrap; gap:8px;">
                            <div style="display:flex; align-items:center; gap:12px;">
                                <span style="
                                    display:inline-flex; align-items:center; justify-content:center;
                                    width:28px; height:28px; border-radius:6px;
                                    background:rgba(245,158,11,0.15); color:#fbbf24;
                                    font-size:0.8rem; font-weight:700;
                                ">${idx + 1}</span>
                                <div>
                                    <div style="font-size:0.82rem; color:#94a3b8;">
                                        Position <strong style="color:#e2e8f0;">${site.position}</strong>
                                        &nbsp;·&nbsp; PAM: <strong style="color:#f59e0b; font-family:'Fira Code',monospace;">${escapeHtml(site.pam_sequence || '')}</strong>
                                        &nbsp;·&nbsp; Strand: <span style="color:${site.strand === '+' ? '#6ee7b7' : '#fca5a5'}; font-weight:600;">${site.strand || '?'}</span>
                                    </div>
                                    <div style="margin-top:6px; font-family:'Fira Code',monospace; font-size:12px; letter-spacing:1px;">
                                        gRNA: ${grnaColored}
                                    </div>
                                </div>
                            </div>
                            <span style="font-size:0.78rem; color:#64748b;">Click to select</span>
                        </div>
                    `;

                    // Click handler: fill gRNA input, highlight card
                    card.addEventListener('click', () => {
                        // Remove selected from all
                        targetContainer.querySelectorAll('.pam-site-card').forEach(c => {
                            c.classList.remove('selected');
                            c.style.borderColor = 'rgba(255,255,255,0.08)';
                            c.style.background = 'rgba(255,255,255,0.03)';
                        });

                        card.classList.add('selected');
                        card.style.borderColor = 'rgba(245,158,11,0.5)';
                        card.style.background = 'rgba(245,158,11,0.08)';

                        const grnaInput = document.getElementById('grna-input');
                        if (grnaInput) {
                            grnaInput.value = site.grna_target || '';
                            const grnaLength = document.getElementById('grna-length');
                            if (grnaLength) grnaLength.textContent = `${grnaInput.value.length}/20`;
                        }

                        const btnSim = document.getElementById('simulate-crispr-btn');
                        if (btnSim) btnSim.disabled = false;
                    });

                    // Hover
                    card.addEventListener('mouseenter', () => {
                        if (!card.classList.contains('selected')) {
                            card.style.borderColor = 'rgba(255,255,255,0.2)';
                            card.style.transform = 'translateX(4px)';
                        }
                    });
                    card.addEventListener('mouseleave', () => {
                        if (!card.classList.contains('selected')) {
                            card.style.borderColor = 'rgba(255,255,255,0.08)';
                        }
                        card.style.transform = 'translateX(0)';
                    });

                    targetContainer.appendChild(card);
                });
            }

            showNotification(`Found ${data.pam_sites?.length || 0} PAM sites`, 'success');

        } catch (err) {
            console.error('[CRISPR] PAM search failed:', err);
            if (targetContainer) hideLoading(targetContainer);
            showNotification('Failed to find PAM sites.', 'error');
        }
    }

    /**
     * Create inline colored HTML for a short DNA sequence.
     * @param {string} seq - DNA sequence string.
     * @returns {string} HTML string with colored spans.
     * @private
     */
    _colorizeInline(seq) {
        const colorMap = {
            A: '#ff6b6b', T: '#ffd93d', G: '#6bcb77', C: '#4d96ff'
        };
        return seq.split('').map(c => {
            const upper = c.toUpperCase();
            const color = colorMap[upper] || '#94a3b8';
            return `<span style="color:${color}">${upper}</span>`;
        }).join('');
    }

    /* ─────────────────────────────────
     * CRISPR Simulation
     * ───────────────────────────────── */

    /**
     * Validate the gRNA input and run the CRISPR simulation.
     * Displays animated results including DNA cutting animation,
     * step-by-step logs, status banner, alignment, and before/after comparison.
     */
    async runSimulation() {
        const grnaInput = document.getElementById('grna-input');
        const grnaValue = (grnaInput?.value || '').trim().toUpperCase();

        // ── Validation ──
        if (!grnaValue) {
            showNotification('Please enter a gRNA sequence.', 'warning');
            return;
        }
        if (grnaValue.length !== 20) {
            showNotification(`gRNA must be exactly 20 characters (got ${grnaValue.length}).`, 'error');
            return;
        }
        if (!/^[ATCG]+$/.test(grnaValue)) {
            showNotification('gRNA must contain only A, T, C, G nucleotides.', 'error');
            return;
        }
        if (!this.currentDiseaseId) {
            showNotification('Please select a gene first.', 'warning');
            return;
        }

        const resultsContainer = document.getElementById('crispr-results');
        if (resultsContainer) {
            resultsContainer.style.display = 'block';
            // Clear inner display areas
            const simLog = document.getElementById('simulation-log');
            if (simLog) simLog.innerHTML = '';
            const alignDisplay = document.getElementById('alignment-display');
            if (alignDisplay) alignDisplay.innerHTML = '';
            const repairedSeq = document.getElementById('repaired-sequence');
            if (repairedSeq) repairedSeq.innerHTML = '';
            const statusIndicator = document.getElementById('status-text');
            if (statusIndicator) statusIndicator.textContent = 'Running simulation...';
        }

        try {
            const res = await fetch('/api/crispr/simulate', {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify({
                    disease_id: this.currentDiseaseId,
                    grna_sequence: grnaValue
                })
            });
            if (!res.ok) throw new Error(`Simulation API error: ${res.status}`);
            const data = await res.json();

            if (resultsContainer) {
                // Populate simulation log
                const simLog = document.getElementById('simulation-log');
                if (simLog && data.logs) {
                    simLog.innerHTML = '';
                    data.logs.forEach((msg, i) => {
                        const lower = msg.toLowerCase();
                        let color = '#fbbf24';
                        let icon = 'ℹ';
                        if (lower.includes('error') || lower.includes('fail') || lower.includes('not found')) {
                            color = '#f87171'; icon = '✗';
                        } else if (lower.includes('found') || lower.includes('success') || lower.includes('complete') || lower.includes('valid')) {
                            color = '#4ade80'; icon = '✓';
                        }
                        simLog.innerHTML += `<div style="padding:4px 0; font-size:0.88rem;"><span style="color:${color}; font-weight:700; margin-right:8px;">${i+1}. ${icon}</span><span style="color:#e2e8f0;">${escapeHtml(msg)}</span></div>`;
                    });
                }

                // Alignment display
                const alignDisplay = document.getElementById('alignment-display');
                if (alignDisplay && data.alignment_display) {
                    alignDisplay.textContent = data.alignment_display;
                }

                // Status
                const statusIndicator = document.getElementById('status-indicator');
                const statusText = document.getElementById('status-text');
                if (statusIndicator && statusText) {
                    const dot = statusIndicator.querySelector('.status-dot');
                    if (data.success) {
                        statusText.textContent = '✅ CRISPR Edit Successful!';
                        statusText.style.color = '#4ade80';
                        if (dot) dot.style.background = '#4ade80';
                    } else {
                        let errorMsg = 'Simulation did not succeed.';
                        if (!data.target_found) errorMsg = 'gRNA target sequence not found in gene.';
                        else if (!data.pam_found) errorMsg = 'No valid PAM site adjacent to target.';
                        statusText.textContent = '❌ ' + errorMsg;
                        statusText.style.color = '#f87171';
                        if (dot) dot.style.background = '#f87171';
                    }
                }

                // Repaired sequence
                const repairedSeq = document.getElementById('repaired-sequence');
                if (repairedSeq && data.success && data.repaired_sequence) {
                    this.colorizeSequence(data.repaired_sequence, 'repaired-sequence');
                }
            }

        } catch (err) {
            console.error('[CRISPR] Simulation failed:', err);
            if (resultsContainer) hideLoading(resultsContainer);
            showNotification('CRISPR simulation failed.', 'error');
        }
    }

    /* ─────────────────────────────────
     * Result Display Helpers
     * ───────────────────────────────── */

    /**
     * Show a brief DNA cutting and repair animation before results.
     * @param {HTMLElement} container - Results container.
     * @param {Object} data - Simulation result data.
     * @returns {Promise<void>}
     * @private
     */
    async _showDnaAnimation(container, data) {
        // Inject animation keyframes
        if (!document.getElementById('crispr-anim-styles')) {
            const style = document.createElement('style');
            style.id = 'crispr-anim-styles';
            style.textContent = `
                @keyframes slideScissors {
                    0% { left: 0%; opacity: 0; }
                    20% { opacity: 1; }
                    60% { left: ${data.cut_position ? Math.min((data.cut_position / 60) * 100, 85) : 50}%; }
                    70% { left: ${data.cut_position ? Math.min((data.cut_position / 60) * 100, 85) : 50}%; transform: rotate(0deg) scale(1.3); }
                    80% { transform: rotate(20deg) scale(1.3); }
                    90% { transform: rotate(-10deg) scale(1.1); }
                    100% { left: ${data.cut_position ? Math.min((data.cut_position / 60) * 100, 85) : 50}%; opacity: 0.3; transform: rotate(0deg) scale(1); }
                }
                @keyframes dnaSplit {
                    0% { gap: 0px; }
                    50% { gap: 20px; }
                    100% { gap: 4px; }
                }
                @keyframes repairSlide {
                    0% { opacity: 0; transform: translateY(-8px); }
                    100% { opacity: 1; transform: translateY(0); }
                }
                @keyframes pulseGlow {
                    0%, 100% { box-shadow: 0 0 20px rgba(34,197,94,0.3); }
                    50% { box-shadow: 0 0 40px rgba(34,197,94,0.6); }
                }
                @keyframes shakeError {
                    0%, 100% { transform: translateX(0); }
                    20% { transform: translateX(-8px); }
                    40% { transform: translateX(8px); }
                    60% { transform: translateX(-6px); }
                    80% { transform: translateX(6px); }
                }
            `;
            document.head.appendChild(style);
        }

        const animDiv = document.createElement('div');
        animDiv.style.cssText = `
            position: relative; height: 80px; margin-bottom: 20px;
            background: rgba(255,255,255,0.02); border-radius: 12px;
            overflow: hidden; display: flex; align-items: center;
            justify-content: center; padding: 0 20px;
        `;

        // DNA backbone representation
        const dnaStrand = document.createElement('div');
        dnaStrand.style.cssText = `
            display: flex; align-items: center; gap: 0px;
            font-family: 'Fira Code', monospace; font-size: 14px;
            letter-spacing: 2px; animation: dnaSplit 1.5s ease-in-out forwards;
            animation-delay: 0.8s;
        `;

        // Build a visual DNA snippet
        const seqSnippet = (data.original_sequence || 'ATCGATCGATCGATCGATCG').substring(0, 30);
        const cutPos = Math.min(data.cut_position || 15, seqSnippet.length);

        const leftPart = document.createElement('span');
        leftPart.innerHTML = this._colorizeInline(seqSnippet.substring(0, cutPos));
        const rightPart = document.createElement('span');
        rightPart.innerHTML = this._colorizeInline(seqSnippet.substring(cutPos));

        dnaStrand.appendChild(leftPart);
        dnaStrand.appendChild(rightPart);

        // Scissors emoji sliding along
        const scissors = document.createElement('span');
        scissors.textContent = '✂️';
        scissors.style.cssText = `
            position: absolute; top: 50%; transform: translateY(-50%);
            font-size: 24px; animation: slideScissors 1.5s ease-in-out forwards;
            z-index: 2;
        `;

        animDiv.appendChild(dnaStrand);
        animDiv.appendChild(scissors);
        container.appendChild(animDiv);

        // Wait for animation to complete
        await new Promise(resolve => setTimeout(resolve, 2000));

        // Fade out animation
        animDiv.style.transition = 'opacity 0.5s';
        animDiv.style.opacity = '0.3';
    }

    /**
     * Display simulation log messages one by one with a typing-like effect.
     * @param {HTMLElement} container - Results container.
     * @param {string[]} logs - Array of log message strings.
     * @returns {Promise<void>}
     * @private
     */
    async _showSimulationLog(container, logs) {
        if (!logs || logs.length === 0) return;

        const logSection = document.createElement('div');
        logSection.style.cssText = `
            background: rgba(0,0,0,0.3); border-radius: 10px;
            padding: 16px 20px; margin-bottom: 20px;
            border: 1px solid rgba(255,255,255,0.06);
        `;

        const logTitle = document.createElement('div');
        logTitle.style.cssText = 'color:#94a3b8; font-size:0.8rem; text-transform:uppercase; letter-spacing:1px; margin-bottom:12px;';
        logTitle.textContent = 'Simulation Log';
        logSection.appendChild(logTitle);

        container.appendChild(logSection);

        for (let i = 0; i < logs.length; i++) {
            const msg = logs[i];
            const logLine = document.createElement('div');
            logLine.style.cssText = `
                opacity: 0; transform: translateX(-10px);
                transition: opacity 0.3s, transform 0.3s;
                padding: 6px 0; font-size: 0.88rem;
                font-family: 'Fira Code', monospace;
            `;

            // Categorize message
            const lower = msg.toLowerCase();
            let icon, color;
            if (lower.includes('error') || lower.includes('fail') || lower.includes('not found')) {
                icon = '✗';
                color = '#f87171';
            } else if (lower.includes('found') || lower.includes('success') || lower.includes('complete') || lower.includes('valid')) {
                icon = '✓';
                color = '#4ade80';
            } else {
                icon = 'ℹ';
                color = '#fbbf24';
            }

            logLine.innerHTML = `
                <span style="color:${color}; font-weight:700; margin-right:8px;">${i + 1}. ${icon}</span>
                <span style="color:#e2e8f0;">${escapeHtml(msg)}</span>
            `;

            logSection.appendChild(logLine);

            // Animate in with delay
            await new Promise(resolve => setTimeout(resolve, 200));
            logLine.style.opacity = '1';
            logLine.style.transform = 'translateX(0)';
        }
    }

    /**
     * Show a large status banner (success or failure) with animation.
     * @param {HTMLElement} container - Results container.
     * @param {Object} data - Simulation result data.
     * @private
     */
    _showStatusBanner(container, data) {
        const banner = document.createElement('div');

        if (data.success) {
            banner.style.cssText = `
                background: linear-gradient(135deg, rgba(34,197,94,0.15), rgba(16,185,129,0.1));
                border: 1px solid rgba(34,197,94,0.3);
                border-radius: 12px; padding: 24px; text-align: center;
                margin-bottom: 20px;
                animation: pulseGlow 2s ease-in-out infinite;
            `;
            banner.innerHTML = `
                <div style="font-size: 48px; margin-bottom: 8px;">✅</div>
                <div style="font-size: 1.4rem; font-weight: 700; color: #4ade80;">CRISPR Edit Successful!</div>
                <div style="color: #94a3b8; margin-top: 8px; font-size: 0.88rem;">
                    Cut at position ${data.cut_position || '—'} · PAM: ${escapeHtml(data.pam_sequence || '—')}
                </div>
            `;
        } else {
            banner.style.cssText = `
                background: linear-gradient(135deg, rgba(239,68,68,0.15), rgba(220,38,38,0.1));
                border: 1px solid rgba(239,68,68,0.3);
                border-radius: 12px; padding: 24px; text-align: center;
                margin-bottom: 20px;
                animation: shakeError 0.5s ease-in-out;
            `;

            let errorMsg = 'Simulation did not succeed.';
            if (!data.target_found) errorMsg = 'gRNA target sequence not found in gene.';
            else if (!data.pam_found) errorMsg = 'No valid PAM site adjacent to target.';

            banner.innerHTML = `
                <div style="font-size: 48px; margin-bottom: 8px;">❌</div>
                <div style="font-size: 1.4rem; font-weight: 700; color: #f87171;">Simulation Failed</div>
                <div style="color: #fca5a5; margin-top: 8px; font-size: 0.88rem;">${escapeHtml(errorMsg)}</div>
            `;
        }

        container.appendChild(banner);
    }

    /**
     * Display the alignment visualization in a pre-formatted block.
     * @param {HTMLElement} container - Results container.
     * @param {string} alignmentText - Pre-formatted alignment string.
     * @private
     */
    _showAlignment(container, alignmentText) {
        const section = document.createElement('div');
        section.style.cssText = 'margin-bottom:20px;';
        section.innerHTML = `
            <div style="color:#94a3b8; font-size:0.8rem; text-transform:uppercase; letter-spacing:1px; margin-bottom:8px;">
                Alignment Display
            </div>
            <pre style="
                background: rgba(0,0,0,0.4); border-radius: 10px;
                padding: 16px 20px; overflow-x: auto; color: #e2e8f0;
                font-family: 'Fira Code', 'Courier New', monospace;
                font-size: 12px; line-height: 1.6;
                border: 1px solid rgba(255,255,255,0.06);
            ">${escapeHtml(alignmentText)}</pre>
        `;
        container.appendChild(section);
    }

    /**
     * Show the repaired sequence with colorized nucleotides and a success badge.
     * @param {HTMLElement} container - Results container.
     * @param {string} sequence - The repaired DNA sequence.
     * @private
     */
    _showRepairedSequence(container, sequence) {
        const section = document.createElement('div');
        section.style.cssText = `
            margin-bottom: 20px; padding: 16px 20px;
            background: rgba(34,197,94,0.05);
            border: 1px solid rgba(34,197,94,0.2);
            border-radius: 10px;
        `;
        section.innerHTML = `
            <div style="display:flex; align-items:center; gap:10px; margin-bottom:10px;">
                <span style="color:#94a3b8; font-size:0.8rem; text-transform:uppercase; letter-spacing:1px;">
                    Repaired Sequence
                </span>
                <span style="
                    background:rgba(34,197,94,0.2); color:#4ade80;
                    padding:2px 10px; border-radius:12px; font-size:0.75rem; font-weight:600;
                ">✅ Repaired</span>
            </div>
            <div id="crispr-repaired-seq-display"></div>
        `;
        container.appendChild(section);

        this.colorizeSequence(sequence, 'crispr-repaired-seq-display');
    }

    /**
     * Show a before/after comparison of the original and repaired sequences.
     * @param {HTMLElement} container - Results container.
     * @param {string} original - Original (susceptible) sequence.
     * @param {string} repaired - Repaired (resistant) sequence.
     * @private
     */
    _showBeforeAfter(container, original, repaired) {
        const section = document.createElement('div');
        section.style.cssText = `
            margin-bottom: 20px; display: grid;
            grid-template-columns: 1fr 1fr; gap: 16px;
        `;

        // Before
        const beforeBox = document.createElement('div');
        beforeBox.style.cssText = `
            padding: 16px; border-radius: 10px;
            background: rgba(239,68,68,0.05);
            border: 1px solid rgba(239,68,68,0.2);
        `;
        beforeBox.innerHTML = `
            <div style="color:#fca5a5; font-size:0.8rem; text-transform:uppercase; letter-spacing:1px; margin-bottom:10px;">
                Before (Susceptible)
            </div>
            <div id="crispr-before-display" style="overflow-x:auto;"></div>
        `;

        // After
        const afterBox = document.createElement('div');
        afterBox.style.cssText = `
            padding: 16px; border-radius: 10px;
            background: rgba(34,197,94,0.05);
            border: 1px solid rgba(34,197,94,0.2);
        `;
        afterBox.innerHTML = `
            <div style="color:#6ee7b7; font-size:0.8rem; text-transform:uppercase; letter-spacing:1px; margin-bottom:10px;">
                After (Resistant)
            </div>
            <div id="crispr-after-display" style="overflow-x:auto;"></div>
        `;

        section.appendChild(beforeBox);
        section.appendChild(afterBox);
        container.appendChild(section);

        this.colorizeSequence(original, 'crispr-before-display');
        this.colorizeSequence(repaired, 'crispr-after-display');
    }
}

/* ───────────────────────────────────────────
 * Global Function: selectDiseaseForCrispr
 * ─────────────────────────────────────────── */

/**
 * Navigate to the CRISPR tab and auto-select a disease gene.
 * Called from other modules (Dashboard, GWAS) for cross-tab navigation.
 * @param {string|number} diseaseId - The disease ID to select.
 */
function selectDiseaseForCrispr(diseaseId) {
    // Switch to CRISPR tab — HTML uses class "tab-btn" with data-section="section-crispr"
    const crisprTab = document.querySelector('.tab-btn[data-section="section-crispr"]');
    if (crisprTab) {
        crisprTab.click();
    }

    // After a short delay for DOM update, select the gene
    setTimeout(() => {
        // Highlight the matching gene card
        const geneCards = document.querySelectorAll('.gene-select-card');
        geneCards.forEach(card => {
            card.classList.remove('active');
            const geneAttr = card.getAttribute('data-gene') || '';
            // Check if this card maps to the requested disease ID
            const idMap = {
                'anthracnose': 'anthracnose',
                'powdery_mildew': 'powdery_mildew',
                'bacterial_canker': 'bacterial_canker',
                'malformation': 'mango_malformation'
            };
            if ((idMap[geneAttr] || geneAttr) === diseaseId) {
                card.classList.add('active');
            }
        });

        // Call selectGene on the simulator
        if (window.crisprSim) {
            window.crisprSim.selectGene(diseaseId);
        }
    }, 200);
}

/* ───────────────────────────────────────────
 * Initialization
 * ─────────────────────────────────────────── */

document.addEventListener('DOMContentLoaded', () => {
    /** @type {CrisprSimulator} Global CRISPR simulator instance */
    window.crisprSim = new CrisprSimulator();
    window.crisprSim.init();

    /** @type {Function} Global function for cross-module CRISPR navigation */
    window.selectDiseaseForCrispr = selectDiseaseForCrispr;
});
