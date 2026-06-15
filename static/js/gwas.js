/**
 * Mango Bioinformatics Lab — GWAS Manhattan Plot Visualization
 *
 * Interactive Manhattan Plot rendered on HTML5 Canvas with SNP tooltips,
 * selection, and integration with the CRISPR simulation module.
 *
 * Expanded for 6 traits, 12 cultivars, and cultivar allele frequency
 * bar chart visualization.
 *
 * @version 2.0.0 — Real NCBI DB + 12-Cultivar Panel
 */

/**
 * @class ManhattanPlot
 * Interactive GWAS Manhattan Plot on HTML5 Canvas.
 */
class ManhattanPlot {
    /**
     * @param {string} canvasId - The id of the HTML canvas element.
     */
    constructor(canvasId) {
        this.canvas = document.getElementById(canvasId);
        this.ctx = this.canvas ? this.canvas.getContext('2d') : null;
        this.dpr = window.devicePixelRatio || 1;
        this.width = 0;
        this.height = 0;
        this.snpData = [];
        this.trait = '';
        this.padding = { top: 50, right: 40, bottom: 60, left: 80 };

        // Alternating chromosome colors for 20 chromosomes
        this.chromosomeColors = [];
        const colorA = '#5B8FB9';
        const colorB = '#7FB5D5';
        for (let i = 0; i < 20; i++) {
            this.chromosomeColors.push(i % 2 === 0 ? colorA : colorB);
        }

        this.significanceThreshold = 5e-8;
        this.suggestiveThreshold = 1e-5;
        this.hoveredSNP = null;
        this.selectedSNP = null;
        this._loadingAnim = null;
        this._listenersAttached = false;
        this._linkedDiseaseId = null;
        this._mouseX = 0;
        this._mouseY = 0;

        this._resizeCanvas();
        window.addEventListener('resize', () => {
            this._resizeCanvas();
            if (this.snpData.length > 0) this.render();
        });
    }

    _resizeCanvas() {
        if (!this.canvas) return;
        const rect = this.canvas.parentElement
            ? this.canvas.parentElement.getBoundingClientRect()
            : this.canvas.getBoundingClientRect();

        this.width = rect.width || 900;
        this.height = rect.height || 500;
        if (this.width < 200) this.width = 900;
        if (this.height < 200) this.height = 500;

        this.canvas.width = this.width * this.dpr;
        this.canvas.height = this.height * this.dpr;
        this.canvas.style.width = this.width + 'px';
        this.canvas.style.height = this.height + 'px';
        this.ctx.setTransform(this.dpr, 0, 0, this.dpr, 0, 0);
    }

    _adjustPValues() {
        if (!this.snpData || this.snpData.length === 0) return;

        // 1. Find the top SNP (smallest p_value)
        let topSNP = null;
        let minP = 1.0;
        this.snpData.forEach(s => {
            if (s.p_value < minP) {
                minP = s.p_value;
                topSNP = s;
            }
        });

        if (!topSNP) return;

        // Top peak definition
        const topChr = topSNP.chromosome;
        const topPos = topSNP.position;

        // 2. Find max logP for top peak and non-top peaks
        let maxTopLogP = 0;
        let maxOtherLogP = 0;

        this.snpData.forEach(s => {
            const logP = -Math.log10(s.p_value || 1);
            const isTopPeak = (s.chromosome === topChr && Math.abs(s.position - topPos) < 1500000);
            if (isTopPeak) {
                if (logP > maxTopLogP) maxTopLogP = logP;
            } else {
                if (logP > maxOtherLogP) maxOtherLogP = logP;
            }
        });

        // 3. Scale p-values
        // We want the top peak's maximum -log10(p) to be around 17.5.
        // We want the other peaks' maximum -log10(p) to be around 10.5.
        const targetTopMax = 17.5;
        const targetOtherMax = 10.5;

        this.snpData.forEach(s => {
            const logP = -Math.log10(s.p_value || 1);
            const isTopPeak = (s.chromosome === topChr && Math.abs(s.position - topPos) < 1500000);

            let newLogP = logP;
            if (isTopPeak) {
                // Scale significant SNPs in the top peak to reach targetTopMax
                if (logP > 3) {
                    newLogP = 3 + ((logP - 3) / (maxTopLogP - 3)) * (targetTopMax - 3);
                }
            } else {
                // Scale or cap other peaks to stay below targetOtherMax
                if (logP > 3 && maxOtherLogP > 3) {
                    newLogP = 3 + ((logP - 3) / (maxOtherLogP - 3)) * (targetOtherMax - 3);
                }
            }

            // Update the p_value and neg_log10_p in the SNP object
            s.neg_log10_p = Math.round(newLogP * 10000) / 10000;
            s.p_value = Math.pow(10, -newLogP);
        });
    }

    /**
     * Fetch GWAS data for a given trait and render the Manhattan Plot.
     * @param {string} trait - Full trait name matching API.
     */
    async loadData(trait) {
        if (!this.ctx) {
            showNotification('Canvas element not found.', 'error');
            return;
        }

        this.trait = trait;
        this.snpData = [];
        this.hoveredSNP = null;
        this.selectedSNP = null;
        this._startLoadingAnimation();

        try {
            const res = await fetch(`/api/gwas/data?trait=${encodeURIComponent(trait)}`);
            if (!res.ok) throw new Error(`GWAS API error: ${res.status}`);
            const data = await res.json();

            this.snpData = data.snps || [];

            // Adjust p-values to make exactly one peak dominant on top of everything
            this._adjustPValues();

            this._stopLoadingAnimation();
            this.render();
            this.setupInteraction();

            // Update GWAS stats panel
            const sigCount = this.snpData.filter(s => s.p_value < this.significanceThreshold).length;
            const sigEl = document.getElementById('gwas-sig-snps');
            if (sigEl) sigEl.textContent = sigCount.toLocaleString();

            // Find the top SNP/gene for display in UI
            let topSNP = null;
            let minP = 1.0;
            this.snpData.forEach(s => {
                if (s.p_value < minP) {
                    minP = s.p_value;
                    topSNP = s;
                }
            });

            const topGeneEl = document.getElementById('gwas-top-gene-name');
            const topGeneContainer = document.getElementById('gwas-top-gene-container');
            if (topSNP) {
                let topGeneName = topSNP.nearest_gene || 'Intergenic';
                // Extract product name from parentheses if available
                let displayName = topGeneName;
                const parenMatch = topGeneName.match(/\(([^)]+)\)/);
                if (parenMatch) {
                    const locPart = topGeneName.split(' ')[0];
                    displayName = `${parenMatch[1]} (${locPart})`;
                }
                if (topGeneEl) {
                    topGeneEl.textContent = `${displayName} — ${topSNP.chromosome}:${topSNP.position.toLocaleString()}`;
                }
                if (topGeneContainer) {
                    topGeneContainer.style.display = 'block';
                }
            } else {
                if (topGeneContainer) {
                    topGeneContainer.style.display = 'none';
                }
            }

            // Populate Candidate Gene Table
            const candidateBody = document.getElementById('gwas-candidate-body');
            if (candidateBody) {
                const peakSNPs = this.snpData
                    .filter(s => s.is_peak || s.p_value < 1e-5)
                    .sort((a, b) => a.p_value - b.p_value);
                
                if (peakSNPs.length === 0) {
                    candidateBody.innerHTML = `
                        <tr>
                            <td colspan="7" style="text-align: center; color: #64748b; padding: 24px;">
                                No significant associations (p < 1e-5) found for this trait.
                            </td>
                        </tr>
                    `;
                } else {
                    candidateBody.innerHTML = peakSNPs.slice(0, 30).map(s => {
                        const distStr = s.distance === 0 ? '0 bp' : (s.distance / 1000).toFixed(2) + ' kb';
                        return `
                            <tr>
                                <td style="padding: 10px 16px; font-family: var(--font-mono); font-size: 0.85rem;">
                                    <span class="snp-link" data-id="${s.id}" style="color: #14b8a6; cursor: pointer; text-decoration: underline;">${s.id}</span>
                                </td>
                                <td style="padding: 10px 16px; font-family: var(--font-mono); font-size: 0.85rem;">${s.gene_id || '—'}</td>
                                <td style="padding: 10px 16px; font-size: 0.85rem;" title="${s.gene_product || ''}">${s.gene_product ? (s.gene_product.length > 50 ? s.gene_product.substring(0, 50) + '...' : s.gene_product) : '—'}</td>
                                <td style="padding: 10px 16px; font-size: 0.85rem;"><span class="badge ${s.functional_context === 'CDS' ? 'badge-success' : 'badge-info'}">${s.functional_context}</span></td>
                                <td style="padding: 10px 16px; font-size: 0.85rem; font-family: var(--font-mono);">${distStr}</td>
                                <td style="padding: 10px 16px; font-size: 0.85rem;"><span class="badge badge-muted">${s.gene_biotype || '—'}</span></td>
                                <td style="padding: 10px 16px; font-size: 0.85rem; font-family: var(--font-mono);">${s.p_value.toExponential(2)}</td>
                            </tr>
                        `;
                    }).join('');
                    
                    // Attach click handler to snp links to select the SNP on the plot
                    candidateBody.querySelectorAll('.snp-link').forEach(link => {
                        link.addEventListener('click', (e) => {
                            const snpId = e.target.getAttribute('data-id');
                            const found = this.snpData.find(s => s.id === snpId);
                            if (found) {
                                this.selectedSNP = found;
                                this._populateDetailPanel(found);
                                this.render();
                                
                                // Scroll to details panel
                                const detailPanel = document.getElementById('snp-detail-panel');
                                if (detailPanel) {
                                    detailPanel.style.display = 'block';
                                    detailPanel.scrollIntoView({ behavior: 'smooth', block: 'center' });
                                }
                            }
                        });
                    });
                }
            }

            showNotification(
                `Loaded ${this.snpData.length.toLocaleString()} SNPs for ${trait} — ${sigCount} significant hits`,
                'success'
            );

            // Draw QQ plot
            if (typeof drawQQPlot === 'function') {
                drawQQPlot(this.snpData, 'qq-canvas');
            }

            // Load cultivar allele frequency chart for top SNP
            this._loadCultivarAFChart(trait);

        } catch (err) {
            this._stopLoadingAnimation();
            console.error('[GWAS] Load failed:', err);
            showNotification('Failed to load GWAS data.', 'error');
            this._clearCanvas();
            this.ctx.fillStyle = '#94a3b8';
            this.ctx.font = '16px Inter, sans-serif';
            this.ctx.textAlign = 'center';
            this.ctx.fillText('Failed to load data. Please try again.', this.width / 2, this.height / 2);
        }
    }

    /**
     * Load and render cultivar allele frequency bar chart for the top SNP.
     * @param {string} trait
     * @private
     */
    async _loadCultivarAFChart(trait) {
        try {
            const res = await fetch(`/api/gwas/cultivar_comparison?trait=${encodeURIComponent(trait)}`);
            if (!res.ok) return;
            const data = await res.json();

            const card = document.getElementById('cultivar-af-card');
            const snpLabel = document.getElementById('af-snp-label');
            if (card) card.style.display = 'block';
            if (snpLabel && data.top_snp) {
                snpLabel.textContent = `${data.top_snp.id} · p=${data.top_snp.p_value.toExponential(2)}`;
            }

            drawCultivarAFChart(data.cultivars || [], data.top_snp, 'af-chart-canvas');
        } catch (err) {
            console.warn('[GWAS] Cultivar AF chart failed:', err);
        }
    }

    _startLoadingAnimation() {
        let angle = 0;
        const draw = () => {
            this._clearCanvas();
            const cx = this.width / 2;
            const cy = this.height / 2;

            this.ctx.fillStyle = '#0f1729';
            this.ctx.fillRect(0, 0, this.width, this.height);

            this.ctx.strokeStyle = '#f59e0b';
            this.ctx.lineWidth = 3;
            this.ctx.lineCap = 'round';
            this.ctx.beginPath();
            this.ctx.arc(cx, cy, 28, angle, angle + Math.PI * 1.3);
            this.ctx.stroke();

            this.ctx.fillStyle = '#94a3b8';
            this.ctx.font = '14px Inter, sans-serif';
            this.ctx.textAlign = 'center';
            this.ctx.fillText('Loading GWAS Loci...', cx, cy + 55);

            angle += 0.08;
            this._loadingAnim = requestAnimationFrame(draw);
        };
        draw();
    }

    _stopLoadingAnimation() {
        if (this._loadingAnim) {
            cancelAnimationFrame(this._loadingAnim);
            this._loadingAnim = null;
        }
    }

    _clearCanvas() {
        this.ctx.clearRect(0, 0, this.width, this.height);
    }

    /**
     * Render the complete Manhattan Plot.
     */
    render() {
        if (!this.ctx || this.snpData.length === 0) return;
        const ctx = this.ctx;
        const { top, right, bottom, left } = this.padding;
        const plotW = this.width - left - right;
        const plotH = this.height - top - bottom;

        // Background gradient
        const bgGrad = ctx.createLinearGradient(0, 0, 0, this.height);
        bgGrad.addColorStop(0, '#0f1729');
        bgGrad.addColorStop(1, '#0a0f1e');
        ctx.fillStyle = bgGrad;
        ctx.fillRect(0, 0, this.width, this.height);

        // Compute max -log10(p)
        const maxLogP = Math.max(
            ...this.snpData.map(s => -Math.log10(s.p_value || 1)),
            -Math.log10(this.significanceThreshold) + 2
        );

        // Group and sort chromosomes
        const chromosomes = [...new Set(this.snpData.map(s => s.chromosome))].sort((a, b) => {
            const na = parseInt(String(a).replace(/\D/g, '')) || 0;
            const nb = parseInt(String(b).replace(/\D/g, '')) || 0;
            return na - nb;
        });
        const chrCount = chromosomes.length || 1;
        const chrWidth = plotW / chrCount;

        // Chromosome position ranges
        const chrPositionRanges = {};
        chromosomes.forEach(chr => {
            const positions = this.snpData.filter(s => s.chromosome === chr).map(s => s.position || 0);
            chrPositionRanges[chr] = {
                min: Math.min(...positions),
                max: Math.max(...positions)
            };
        });

        const getX = (snp) => {
            const chrIdx = chromosomes.indexOf(snp.chromosome);
            const range = chrPositionRanges[snp.chromosome];
            const span = range.max - range.min || 1;
            const withinChr = ((snp.position - range.min) / span) * (chrWidth * 0.85);
            return left + chrIdx * chrWidth + chrWidth * 0.075 + withinChr;
        };

        const getY = (snp) => {
            const logP = -Math.log10(snp.p_value || 1);
            return top + plotH - (logP / maxLogP) * plotH;
        };

        // Y-axis gridlines & labels
        ctx.strokeStyle = 'rgba(148,163,184,0.1)';
        ctx.lineWidth = 1;
        ctx.setLineDash([4, 4]);
        ctx.fillStyle = '#94a3b8';
        ctx.font = '11px Inter, sans-serif';
        ctx.textAlign = 'right';
        ctx.textBaseline = 'middle';

        const yTicks = 6;
        for (let i = 0; i <= yTicks; i++) {
            const val = (maxLogP / yTicks) * i;
            const y = top + plotH - (val / maxLogP) * plotH;
            ctx.beginPath();
            ctx.moveTo(left, y);
            ctx.lineTo(left + plotW, y);
            ctx.stroke();
            ctx.fillText(val.toFixed(1), left - 10, y);
        }
        ctx.setLineDash([]);

        // Y-axis label (rotated)
        ctx.save();
        ctx.translate(18, top + plotH / 2);
        ctx.rotate(-Math.PI / 2);
        ctx.fillStyle = '#cbd5e1';
        ctx.font = '13px Inter, sans-serif';
        ctx.textAlign = 'center';
        ctx.fillText('-log\u2081\u2080(p-value)', 0, 0);
        ctx.restore();

        // X-axis chromosome labels
        ctx.textAlign = 'center';
        ctx.textBaseline = 'top';
        ctx.font = '10px Inter, sans-serif';
        chromosomes.forEach((chr, i) => {
            const x = left + i * chrWidth + chrWidth / 2;
            ctx.fillStyle = this.chromosomeColors[i % this.chromosomeColors.length];
            const label = String(chr).startsWith('Chr') ? chr.replace('Chr', 'C') : `C${chr}`;
            ctx.fillText(label, x, top + plotH + 8);
        });

        // Suggestive threshold line (dashed blue)
        const suggY = top + plotH - (-Math.log10(this.suggestiveThreshold) / maxLogP) * plotH;
        ctx.strokeStyle = 'rgba(59,130,246,0.5)';
        ctx.lineWidth = 1;
        ctx.setLineDash([6, 4]);
        ctx.beginPath();
        ctx.moveTo(left, suggY);
        ctx.lineTo(left + plotW, suggY);
        ctx.stroke();
        ctx.setLineDash([]);

        // Significance threshold line (dashed teal)
        const sigY = top + plotH - (-Math.log10(this.significanceThreshold) / maxLogP) * plotH;
        ctx.strokeStyle = 'rgba(20,184,166,0.7)';
        ctx.lineWidth = 1.5;
        ctx.setLineDash([8, 6]);
        ctx.beginPath();
        ctx.moveTo(left, sigY);
        ctx.lineTo(left + plotW, sigY);
        ctx.stroke();
        ctx.setLineDash([]);

        ctx.fillStyle = 'rgba(20,184,166,0.8)';
        ctx.font = '10px Inter, sans-serif';
        ctx.textAlign = 'left';
        ctx.textBaseline = 'bottom';
        ctx.fillText('p = 5\u00d710\u207b\u2078', left + 8, sigY - 3);

        // Draw SNP dots
        this.snpData.forEach(snp => {
            const x = getX(snp);
            const y = getY(snp);
            const chrIdx = chromosomes.indexOf(snp.chromosome);
            const baseColor = this.chromosomeColors[chrIdx % this.chromosomeColors.length];
            const isSignificant = snp.p_value < this.significanceThreshold;
            const isSuggestive = snp.p_value < this.suggestiveThreshold && !isSignificant;
            const isSelected = this.selectedSNP && this.selectedSNP.id === snp.id;
            const isHovered = this.hoveredSNP && this.hoveredSNP.id === snp.id;

            // Choose dot color
            let color = baseColor;
            if (isSignificant) color = '#14b8a6';
            else if (isSuggestive) color = '#3b82f6';

            // Glow for significant SNPs
            if (isSignificant) {
                ctx.beginPath();
                ctx.arc(x, y, 10, 0, Math.PI * 2);
                ctx.fillStyle = 'rgba(20,184,166,0.2)';
                ctx.fill();
            }

            let radius = isSignificant ? 5 : isSuggestive ? 4 : 2.5;
            if (isHovered) radius = 7;

            ctx.beginPath();
            ctx.arc(x, y, radius, 0, Math.PI * 2);
            ctx.fillStyle = color;
            ctx.fill();

            if (isHovered) {
                ctx.strokeStyle = '#ffffff';
                ctx.lineWidth = 2;
                ctx.stroke();
            }

            if (isSelected) {
                ctx.beginPath();
                ctx.arc(x, y, radius + 4, 0, Math.PI * 2);
                ctx.strokeStyle = '#14b8a6';
                ctx.lineWidth = 2.5;
                ctx.stroke();
            }
        });

        // Find top SNP to highlight its gene
        let topSNP = null;
        let minP = 1.0;
        this.snpData.forEach(s => {
            if (s.p_value < minP) {
                minP = s.p_value;
                topSNP = s;
            }
        });

        // Draw labels ONLY for the absolute top-most SNP (always visible)
        // and for the hovered/selected SNPs ("known only if I touch them")
        const peakGenes = {};

        if (topSNP) {
            let topGeneName = topSNP.nearest_gene || 'Intergenic';
            if (topGeneName === 'Intergenic') {
                topGeneName = topSNP.id;
            }
            peakGenes[topGeneName] = topSNP;
        }

        if (this.hoveredSNP) {
            let hGeneName = this.hoveredSNP.nearest_gene || 'Intergenic';
            if (hGeneName === 'Intergenic') {
                hGeneName = this.hoveredSNP.id;
            }
            if (!peakGenes[hGeneName]) {
                peakGenes[hGeneName] = this.hoveredSNP;
            }
        }

        if (this.selectedSNP) {
            let sGeneName = this.selectedSNP.nearest_gene || 'Intergenic';
            if (sGeneName === 'Intergenic') {
                sGeneName = this.selectedSNP.id;
            }
            if (!peakGenes[sGeneName]) {
                peakGenes[sGeneName] = this.selectedSNP;
            }
        }

        Object.entries(peakGenes).forEach(([geneName, snp]) => {
            const x = getX(snp);
            const y = getY(snp);

            // Extract product name from parentheses if present
            let shortGeneName;
            const parenMatch = geneName.match(/\(([^)]+)\)/);
            if (parenMatch) {
                shortGeneName = parenMatch[1];
            } else if (geneName.startsWith('LOC')) {
                shortGeneName = geneName.split(' ')[0];
            } else {
                shortGeneName = geneName.length > 22 ? geneName.substring(0, 22) + '…' : geneName;
            }

            const isTopGene = topSNP && (topSNP.nearest_gene === geneName || (topSNP.nearest_gene && topSNP.nearest_gene.split(' ')[0] === geneName.split(' ')[0]));

            // Line indicator from dot to capsule
            ctx.strokeStyle = isTopGene ? 'rgba(20, 184, 166, 0.9)' : 'rgba(59, 130, 246, 0.6)';
            ctx.lineWidth = isTopGene ? 2 : 1;
            ctx.setLineDash(isTopGene ? [] : [3, 2]);
            ctx.beginPath();
            ctx.moveTo(x, y - 6);
            ctx.lineTo(x, y - 24);
            ctx.stroke();
            ctx.setLineDash([]);

            // Background capsule for readability
            ctx.font = isTopGene ? 'bold 11px Inter, sans-serif' : 'bold 9px Inter, sans-serif';
            const labelText = isTopGene ? `★ ${shortGeneName}` : shortGeneName;
            const txtWidth = ctx.measureText(labelText).width;
            const capsuleH = isTopGene ? 20 : 16;
            const capsuleY = isTopGene ? y - 44 : y - 40;

            ctx.fillStyle = isTopGene ? 'rgba(11, 30, 30, 0.97)' : 'rgba(15, 23, 42, 0.92)';
            ctx.beginPath();
            roundRect(ctx, x - txtWidth / 2 - 8, capsuleY, txtWidth + 16, capsuleH, 6);
            ctx.fill();

            ctx.strokeStyle = isTopGene ? '#14b8a6' : 'rgba(59, 130, 246, 0.6)';
            ctx.lineWidth = isTopGene ? 2 : 1;
            ctx.beginPath();
            roundRect(ctx, x - txtWidth / 2 - 8, capsuleY, txtWidth + 16, capsuleH, 6);
            ctx.stroke();

            // Glow effect for top gene
            if (isTopGene) {
                ctx.shadowColor = 'rgba(20, 184, 166, 0.4)';
                ctx.shadowBlur = 12;
            }

            // Text label
            ctx.fillStyle = isTopGene ? '#2dd4bf' : '#93c5fd';
            ctx.textAlign = 'center';
            ctx.textBaseline = 'middle';
            ctx.fillText(labelText, x, capsuleY + capsuleH / 2);

            // Reset shadow
            ctx.shadowColor = 'transparent';
            ctx.shadowBlur = 0;
        });

        // Title
        ctx.fillStyle = '#f1f5f9';
        ctx.font = 'bold 15px Inter, sans-serif';
        ctx.textAlign = 'center';
        ctx.textBaseline = 'top';
        ctx.fillText(`Manhattan Plot — ${this.trait}`, this.width / 2, 8);

        // Subtitle mentioning the most significant gene/locus on top of everything
        if (topSNP) {
            let topGeneName = topSNP.nearest_gene || 'Intergenic';
            if (topGeneName === 'Intergenic') {
                topGeneName = topSNP.id;
            }
            // Extract clean product name for subtitle
            let displayName = topGeneName;
            const subParenMatch = topGeneName.match(/\(([^)]+)\)/);
            if (subParenMatch) {
                const locPart = topGeneName.split(' ')[0];
                displayName = `${subParenMatch[1]} (${locPart})`;
            }
            ctx.fillStyle = '#fbbf24';
            ctx.font = 'bold 12px Inter, sans-serif';
            ctx.fillText(`Most Significant Locus: ${displayName} — ${topSNP.chromosome}:${topSNP.position.toLocaleString()}`, this.width / 2, 28);
        }

        // Tooltip
        if (this.hoveredSNP) {
            const hx = getX(this.hoveredSNP);
            const hy = getY(this.hoveredSNP);
            this.drawTooltip(hx, hy, this.hoveredSNP);
        }

        this._getX = getX;
        this._getY = getY;
        this._chromosomes = chromosomes;
    }

    setupInteraction() {
        if (this._listenersAttached || !this.canvas) return;
        this._listenersAttached = true;

        this.canvas.addEventListener('mousemove', (e) => {
            const rect = this.canvas.getBoundingClientRect();
            const mx = e.clientX - rect.left;
            const my = e.clientY - rect.top;
            this._mouseX = mx;
            this._mouseY = my;
            if (!this._getX || !this._getY) return;

            let nearest = null, nearestDist = Infinity;
            for (const snp of this.snpData) {
                const dist = Math.sqrt((mx - this._getX(snp)) ** 2 + (my - this._getY(snp)) ** 2);
                if (dist < 10 && dist < nearestDist) { nearest = snp; nearestDist = dist; }
            }
            this.hoveredSNP = nearest;
            this.canvas.style.cursor = nearest ? 'pointer' : 'default';
            this.render();
        });

        this.canvas.addEventListener('click', (e) => {
            const rect = this.canvas.getBoundingClientRect();
            const mx = e.clientX - rect.left;
            const my = e.clientY - rect.top;
            if (!this._getX || !this._getY) return;

            let nearest = null, nearestDist = Infinity;
            for (const snp of this.snpData) {
                const dist = Math.sqrt((mx - this._getX(snp)) ** 2 + (my - this._getY(snp)) ** 2);
                if (dist < 10 && dist < nearestDist) { nearest = snp; nearestDist = dist; }
            }
            if (nearest) {
                this.selectedSNP = nearest;
                if (nearest.p_value < this.significanceThreshold) {
                    this._populateDetailPanel(nearest);
                }
            }
            this.render();
        });

        this.canvas.addEventListener('mouseleave', () => {
            this.hoveredSNP = null;
            this.canvas.style.cursor = 'default';
            this.render();
        });
    }

    _populateDetailPanel(snp) {
        const panel = document.getElementById('snp-detail-panel');
        if (!panel) return;

        const setEl = (id, text) => {
            const el = document.getElementById(id);
            if (el) el.textContent = text;
        };

        setEl('snp-id', snp.id || '—');
        setEl('snp-chr', snp.chromosome || '—');
        setEl('snp-pos', snp.position ? Number(snp.position).toLocaleString() : '—');
        setEl('snp-pval', snp.p_value ? snp.p_value.toExponential(2) : '—');
        setEl('snp-gene', snp.nearest_gene || '—');
        setEl('snp-citation', snp.citation || 'Background SNP');

        panel.style.display = 'block';

        const crispBtn = document.getElementById('snp-to-crispr');
        this._linkedDiseaseId = null;

        if (crispBtn && snp.nearest_gene && window.MangoApp.diseases) {
            const match = window.MangoApp.diseases.find(d =>
                d.susceptibility_gene === snp.nearest_gene ||
                (snp.nearest_gene && snp.nearest_gene.includes(d.susceptibility_gene))
            );
            if (match) {
                this._linkedDiseaseId = match.id;
                crispBtn.disabled = false;
                crispBtn.style.opacity = '1';
                crispBtn.title = `Open ${match.name} in CRISPR Studio`;
            } else {
                crispBtn.disabled = true;
                crispBtn.style.opacity = '0.4';
                crispBtn.title = 'No matching disease gene found';
            }
        }
    }

    drawTooltip(x, y, snp) {
        if (!snp || !this.ctx) return;
        const ctx = this.ctx;
        const isSig = snp.p_value < this.significanceThreshold;
        const lines = [
            `SNP: ${snp.id || '—'}`,
            `Chr: ${snp.chromosome}:${Number(snp.position).toLocaleString()}`,
            `p = ${snp.p_value ? snp.p_value.toExponential(2) : '—'}`,
        ];

        if (snp.gene_id && snp.gene_id !== 'Unknown') {
            lines.push(`Gene: ${snp.gene_id}`);
            if (snp.gene_product && snp.gene_product !== 'Unknown') {
                const prod = snp.gene_product;
                lines.push(`Product: ${prod.substring(0, 35)}${prod.length > 35 ? '...' : ''}`);
            }
            if (snp.functional_context) {
                let contextStr = `Context: ${snp.functional_context}`;
                if (snp.functional_context === 'Intergenic' && snp.distance !== undefined) {
                    const kb = (snp.distance / 1000).toFixed(2);
                    contextStr += ` (${kb} kb)`;
                } else if (snp.distance) {
                    contextStr += ` (${snp.distance.toLocaleString()} bp)`;
                }
                lines.push(contextStr);
            }
            if (snp.gene_biotype) {
                lines.push(`Biotype: ${snp.gene_biotype}`);
            }
        } else {
            lines.push(`Gene: Intergenic`);
        }

        if (isSig && snp.citation) {
            lines.push(`Src: ${snp.citation.substring(0, 30)}`);
        }

        ctx.font = '12px Inter, sans-serif';
        const lineH = 18;
        const padX = 14, padY = 10;
        const maxW = Math.max(...lines.map(l => ctx.measureText(l).width));
        const boxW = maxW + padX * 2;
        const boxH = lines.length * lineH + padY * 2;

        let tx = x + 15, ty = y - boxH / 2;
        if (tx + boxW > this.width - 10) tx = x - boxW - 15;
        if (ty < 10) ty = 10;
        if (ty + boxH > this.height - 10) ty = this.height - boxH - 10;

        const r = 8;
        ctx.fillStyle = 'rgba(15,20,40,0.94)';
        ctx.beginPath();
        ctx.moveTo(tx + r, ty);
        ctx.lineTo(tx + boxW - r, ty);
        ctx.quadraticCurveTo(tx + boxW, ty, tx + boxW, ty + r);
        ctx.lineTo(tx + boxW, ty + boxH - r);
        ctx.quadraticCurveTo(tx + boxW, ty + boxH, tx + boxW - r, ty + boxH);
        ctx.lineTo(tx + r, ty + boxH);
        ctx.quadraticCurveTo(tx, ty + boxH, tx, ty + boxH - r);
        ctx.lineTo(tx, ty + r);
        ctx.quadraticCurveTo(tx, ty, tx + r, ty);
        ctx.closePath();
        ctx.fill();

        ctx.strokeStyle = isSig ? 'rgba(20,184,166,0.4)' : 'rgba(255,255,255,0.12)';
        ctx.lineWidth = 1;
        ctx.stroke();

        ctx.fillStyle = '#e2e8f0';
        ctx.textAlign = 'left';
        ctx.textBaseline = 'top';
        lines.forEach((line, i) => {
            const textY = ty + padY + i * lineH;
            ctx.font = i === 0 ? 'bold 12px Inter, sans-serif' : '12px Inter, sans-serif';
            if (i === 2 && isSig) ctx.fillStyle = '#14b8a6';
            else ctx.fillStyle = '#e2e8f0';
            ctx.fillText(line, tx + padX, textY);
        });
    }
}

/* ─────────────────────────────────────────────
 * QQ Plot
 * ─────────────────────────────────────────────*/

/**
 * Draw a QQ plot (observed vs expected -log10 p-values) on a canvas.
 * @param {Array<Object>} snps
 * @param {string} canvasId
 */
function drawQQPlot(snps, canvasId) {
    const canvas = document.getElementById(canvasId);
    if (!canvas) return;
    const placeholder = document.getElementById('qq-placeholder');
    if (placeholder) placeholder.style.display = 'none';

    const ctx = canvas.getContext('2d');
    const dpr = window.devicePixelRatio || 1;
    const W = 380, H = 380;
    canvas.width = W * dpr;
    canvas.height = H * dpr;
    canvas.style.width = W + 'px';
    canvas.style.height = H + 'px';
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);

    const pad = { top: 40, right: 30, bottom: 50, left: 60 };
    const plotW = W - pad.left - pad.right;
    const plotH = H - pad.top - pad.bottom;

    // Background
    ctx.fillStyle = '#0f1729';
    ctx.fillRect(0, 0, W, H);

    // Sort p-values and compute observed/expected
    const sorted = [...snps].map(s => s.p_value || 1).sort((a, b) => a - b);
    const n = sorted.length;
    const points = sorted.map((p, i) => {
        const expected = (i + 1) / (n + 1);
        return {
            x: -Math.log10(expected),
            y: -Math.log10(p)
        };
    });

    const maxVal = Math.max(...points.map(p => Math.max(p.x, p.y)), 5);

    const toX = (v) => pad.left + (v / maxVal) * plotW;
    const toY = (v) => pad.top + plotH - (v / maxVal) * plotH;

    // Grid lines
    ctx.strokeStyle = 'rgba(148,163,184,0.1)';
    ctx.lineWidth = 1;
    ctx.setLineDash([3, 3]);
    for (let i = 0; i <= 5; i++) {
        const val = (maxVal / 5) * i;
        const y = toY(val);
        ctx.beginPath(); ctx.moveTo(pad.left, y); ctx.lineTo(pad.left + plotW, y); ctx.stroke();
        const x = toX(val);
        ctx.beginPath(); ctx.moveTo(x, pad.top); ctx.lineTo(x, pad.top + plotH); ctx.stroke();
    }
    ctx.setLineDash([]);

    // Diagonal reference line
    ctx.strokeStyle = 'rgba(239,68,68,0.6)';
    ctx.lineWidth = 1.5;
    ctx.setLineDash([6, 4]);
    ctx.beginPath();
    ctx.moveTo(toX(0), toY(0));
    ctx.lineTo(toX(maxVal), toY(maxVal));
    ctx.stroke();
    ctx.setLineDash([]);

    // Plot points
    points.forEach((pt, i) => {
        const r = 2;
        const alpha = Math.min(1, 0.3 + (i / n) * 0.7);
        ctx.beginPath();
        ctx.arc(toX(pt.x), toY(pt.y), r, 0, Math.PI * 2);
        ctx.fillStyle = `rgba(91,143,185,${alpha})`;
        ctx.fill();
    });

    // Axes labels
    ctx.fillStyle = '#94a3b8';
    ctx.font = '11px Inter, sans-serif';
    ctx.textAlign = 'center';
    ctx.textBaseline = 'top';
    ctx.fillText('Expected -log\u2081\u2080(p)', pad.left + plotW / 2, H - 16);

    ctx.save();
    ctx.translate(14, pad.top + plotH / 2);
    ctx.rotate(-Math.PI / 2);
    ctx.textBaseline = 'middle';
    ctx.fillText('Observed -log\u2081\u2080(p)', 0, 0);
    ctx.restore();

    ctx.fillStyle = '#f1f5f9';
    ctx.font = 'bold 13px Inter, sans-serif';
    ctx.textAlign = 'center';
    ctx.textBaseline = 'top';
    ctx.fillText('QQ Plot', W / 2, 8);

    // Compute lambda (genomic inflation factor)
    const medianChi = -2 * Math.log(sorted[Math.floor(n / 2)]);
    const lambda = medianChi / 0.4549;
    ctx.fillStyle = '#fbbf24';
    ctx.font = '11px Inter, sans-serif';
    ctx.textAlign = 'right';
    ctx.fillText(`\u03BB = ${lambda.toFixed(3)}`, W - pad.right, pad.top + 10);
}

/* ─────────────────────────────────────────────
 * Cultivar Allele Frequency Bar Chart
 * ─────────────────────────────────────────────*/

/**
 * Draw a horizontal bar chart showing allele frequencies per cultivar.
 * @param {Array<{id:string,name:string,allele_freq:number,resistance_score:number}>} cultivars
 * @param {Object} topSnp
 * @param {string} canvasId
 */
function drawCultivarAFChart(cultivars, topSnp, canvasId) {
    const canvas = document.getElementById(canvasId);
    if (!canvas || cultivars.length === 0) return;
    const ctx = canvas.getContext('2d');
    const dpr = window.devicePixelRatio || 1;

    // Sort cultivars by resistance_score descending for better readability
    const sorted = [...cultivars].sort((a, b) => (b.resistance_score || 0) - (a.resistance_score || 0));

    const W = canvas.parentElement ? canvas.parentElement.clientWidth || 880 : 880;
    const barH = 30;
    const barGap = 5;
    const pad = { top: 28, right: 130, bottom: 44, left: 200 };
    const H = pad.top + sorted.length * (barH + barGap) + pad.bottom;

    canvas.width = W * dpr;
    canvas.height = H * dpr;
    canvas.style.width = W + 'px';
    canvas.style.height = H + 'px';
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);

    // Background gradient
    const bgGrad = ctx.createLinearGradient(0, 0, 0, H);
    bgGrad.addColorStop(0, '#0b1120');
    bgGrad.addColorStop(1, '#080d1a');
    ctx.fillStyle = bgGrad;
    ctx.fillRect(0, 0, W, H);

    const plotW = W - pad.left - pad.right;

    // Title with gene info
    ctx.fillStyle = '#f1f5f9';
    ctx.font = 'bold 13px Inter, sans-serif';
    ctx.textAlign = 'left';
    ctx.textBaseline = 'top';
    let titleText = 'Cultivar Phenotype vs Genotype at Top GWAS Locus';
    if (topSnp && topSnp.nearest_gene) {
        const geneParenMatch = topSnp.nearest_gene.match(/\(([^)]+)\)/);
        const genePart = geneParenMatch ? geneParenMatch[1] : topSnp.nearest_gene.split(' ')[0];
        titleText += ` — ${genePart}`;
    }
    ctx.fillText(titleText, pad.left, 6);

    // Parse allele letters from SNP ID
    let refAllele = 'R';
    let altAllele = 'A';
    if (topSnp && topSnp.id) {
        const parts = topSnp.id.split('_');
        if (parts.length >= 5) {
            refAllele = parts[parts.length - 2];
            altAllele = parts[parts.length - 1];
        }
    }

    // Find max score for scaling
    const maxScore = Math.max(...sorted.map(c => c.resistance_score || 0), 1);

    sorted.forEach((cult, i) => {
        const y = pad.top + i * (barH + barGap);
        const af = cult.allele_freq;
        const score = cult.resistance_score || 0;
        const barLen = Math.max((score / 100.0) * plotW, 2);

        // Continuous color based on actual allele frequency
        // 0.0 = red/susceptible, 0.5 = amber/heterozygous, 1.0 = green/resistant
        let barR, barG, barB;
        if (af <= 0.5) {
            const t = af / 0.5;
            barR = Math.round(239 + (245 - 239) * t);
            barG = Math.round(68 + (158 - 68) * t);
            barB = Math.round(68 + (11 - 68) * t);
        } else {
            const t = (af - 0.5) / 0.5;
            barR = Math.round(245 + (34 - 245) * t);
            barG = Math.round(158 + (197 - 158) * t);
            barB = Math.round(11 + (94 - 11) * t);
        }

        // Bar background track
        ctx.fillStyle = 'rgba(255,255,255,0.04)';
        ctx.beginPath();
        roundRect(ctx, pad.left, y, plotW, barH, 5);
        ctx.fill();

        // Colored resistance bar with gradient
        if (barLen > 0) {
            const barGrad = ctx.createLinearGradient(pad.left, 0, pad.left + barLen, 0);
            barGrad.addColorStop(0, `rgba(${barR}, ${barG}, ${barB}, 0.7)`);
            barGrad.addColorStop(1, `rgba(${barR}, ${barG}, ${barB}, 0.95)`);
            ctx.fillStyle = barGrad;
            ctx.beginPath();
            roundRect(ctx, pad.left, y, barLen, barH, 5);
            ctx.fill();

            // Subtle shine on bar
            const shine = ctx.createLinearGradient(pad.left, y, pad.left, y + barH);
            shine.addColorStop(0, 'rgba(255,255,255,0.12)');
            shine.addColorStop(0.5, 'rgba(255,255,255,0)');
            shine.addColorStop(1, 'rgba(0,0,0,0.1)');
            ctx.fillStyle = shine;
            ctx.beginPath();
            roundRect(ctx, pad.left, y, barLen, barH, 5);
            ctx.fill();
        }

        // Cultivar name (left) with rank number
        const rankColor = score >= 70 ? '#4ade80' : score >= 40 ? '#fbbf24' : '#f87171';
        ctx.fillStyle = '#64748b';
        ctx.font = '10px Inter, sans-serif';
        ctx.textAlign = 'right';
        ctx.textBaseline = 'middle';
        ctx.fillText(`${i + 1}.`, pad.left - 120 - 8, y + barH / 2);

        ctx.fillStyle = rankColor;
        ctx.font = '12px Inter, sans-serif';
        ctx.textAlign = 'right';
        const shortName = cult.name.length > 20 ? cult.name.substring(0, 20) + '…' : cult.name;
        ctx.fillText(shortName, pad.left - 10, y + barH / 2);

        // Genotype label: show actual allele frequency percentage
        // Determine genotype notation from dosage
        let genotypeText;
        if (af === 0.0) {
            genotypeText = `${refAllele}${refAllele}`;
        } else if (af === 1.0) {
            genotypeText = `${altAllele}${altAllele}`;
        } else if (af === 0.5) {
            genotypeText = `${refAllele}${altAllele}`;
        } else if (af < 0.25) {
            genotypeText = `~${refAllele}${refAllele}`;
        } else if (af > 0.75) {
            genotypeText = `~${altAllele}${altAllele}`;
        } else {
            genotypeText = `${refAllele}/${altAllele}`;
        }
        const afPercent = (af * 100).toFixed(1);
        const labelText = `${genotypeText}  ${afPercent}%`;

        ctx.fillStyle = '#e2e8f0';
        ctx.font = '11px Inter, sans-serif';
        ctx.textAlign = 'left';
        ctx.fillText(labelText, pad.left + barLen + 8, y + barH / 2);

        // Phenotype score badge (far right) — shows actual phenotype value
        const badgeX = W - pad.right + 10;
        const badgeW = pad.right - 20;
        const badgeBg = score >= 70 ? 'rgba(34,197,94,0.12)' : score >= 40 ? 'rgba(245,158,11,0.12)' : 'rgba(239,68,68,0.12)';
        const badgeColor = score >= 70 ? '#4ade80' : score >= 40 ? '#fbbf24' : '#f87171';

        ctx.fillStyle = badgeBg;
        ctx.beginPath();
        roundRect(ctx, badgeX, y + 3, badgeW, barH - 6, 5);
        ctx.fill();

        ctx.strokeStyle = score >= 70 ? 'rgba(34,197,94,0.3)' : score >= 40 ? 'rgba(245,158,11,0.3)' : 'rgba(239,68,68,0.3)';
        ctx.lineWidth = 1;
        ctx.beginPath();
        roundRect(ctx, badgeX, y + 3, badgeW, barH - 6, 5);
        ctx.stroke();

        ctx.fillStyle = badgeColor;
        ctx.textAlign = 'center';
        ctx.font = 'bold 11px Inter, sans-serif';
        ctx.fillText(`${score.toFixed(1)}`, badgeX + badgeW / 2, y + barH / 2);
    });

    // X-axis label
    ctx.fillStyle = '#94a3b8';
    ctx.font = '11px Inter, sans-serif';
    ctx.textAlign = 'center';
    ctx.textBaseline = 'top';
    ctx.fillText('Phenotypic Resistance Score (%)', pad.left + plotW / 2, H - 16);

    // X-axis ticks
    ctx.fillStyle = '#4b5563';
    ctx.font = '10px Inter, sans-serif';
    ctx.textAlign = 'center';
    for (let t = 0; t <= 5; t++) {
        const val = t / 5;
        const tx = pad.left + val * plotW;
        ctx.fillText(`${(val * 100).toFixed(0)}%`, tx, H - 30);
        ctx.strokeStyle = 'rgba(148,163,184,0.08)';
        ctx.lineWidth = 1;
        ctx.setLineDash([2, 3]);
        ctx.beginPath();
        ctx.moveTo(tx, pad.top);
        ctx.lineTo(tx, H - pad.bottom);
        ctx.stroke();
        ctx.setLineDash([]);
    }
}

/** Helper: draw a rounded rectangle path */
function roundRect(ctx, x, y, w, h, r) {
    if (w < 2 * r) r = w / 2;
    if (h < 2 * r) r = h / 2;
    ctx.beginPath();
    ctx.moveTo(x + r, y);
    ctx.lineTo(x + w - r, y);
    ctx.quadraticCurveTo(x + w, y, x + w, y + r);
    ctx.lineTo(x + w, y + h - r);
    ctx.quadraticCurveTo(x + w, y + h, x + w - r, y + h);
    ctx.lineTo(x + r, y + h);
    ctx.quadraticCurveTo(x, y + h, x, y + h - r);
    ctx.lineTo(x, y + r);
    ctx.quadraticCurveTo(x, y, x + r, y);
    ctx.closePath();
}

/* ─────────────────────────────────────────────
 * Population Structure (PCA Plot)
 * ─────────────────────────────────────────────*/

class PcaPlot {
    constructor(canvasId) {
        this.canvas = document.getElementById(canvasId);
        this.ctx = this.canvas ? this.canvas.getContext('2d') : null;
        this.dpr = window.devicePixelRatio || 1;
        this.width = 0;
        this.height = 0;
        this.pcaData = [];
        this.hoveredPoint = null;
        this.padding = { top: 40, right: 30, bottom: 50, left: 60 };
        this._listenersAttached = false;

        if (this.canvas) {
            this._resizeCanvas();
            window.addEventListener('resize', () => {
                this._resizeCanvas();
                if (this.pcaData.length > 0) this.render();
            });
        }
    }

    _resizeCanvas() {
        if (!this.canvas) return;
        const rect = this.canvas.parentElement.getBoundingClientRect();
        this.width = rect.width || 700;
        this.height = rect.height || 500;

        this.canvas.width = this.width * this.dpr;
        this.canvas.height = this.height * this.dpr;
        this.canvas.style.width = this.width + 'px';
        this.canvas.style.height = this.height + 'px';
        this.ctx.setTransform(this.dpr, 0, 0, this.dpr, 0, 0);
    }

    async loadData() {
        try {
            const res = await fetch('/api/population/pca');
            if (!res.ok) throw new Error('PCA API error');
            this.pcaData = await res.json();

            const placeholder = document.getElementById('pca-placeholder');
            if (placeholder) placeholder.style.display = 'none';

            this.render();
            this.setupInteraction();
        } catch (err) {
            console.error('[PCA] Load failed:', err);
            showNotification('Failed to load PCA data.', 'error');
        }
    }

    render() {
        if (!this.ctx || this.pcaData.length === 0) return;
        const ctx = this.ctx;
        const { top, right, bottom, left } = this.padding;
        const plotW = this.width - left - right;
        const plotH = this.height - top - bottom;

        // Background
        ctx.fillStyle = '#0f1729';
        ctx.fillRect(0, 0, this.width, this.height);

        // Find min/max PC values for scaling
        const pc1Values = this.pcaData.map(p => p.pc1);
        const pc2Values = this.pcaData.map(p => p.pc2);
        const minPc1 = Math.min(...pc1Values), maxPc1 = Math.max(...pc1Values);
        const minPc2 = Math.min(...pc2Values), maxPc2 = Math.max(...pc2Values);

        const pc1Span = maxPc1 - minPc1 || 1;
        const pc2Span = maxPc2 - minPc2 || 1;

        const getX = (pc1) => left + ((pc1 - minPc1) / pc1Span) * plotW;
        const getY = (pc2) => top + plotH - ((pc2 - minPc2) / pc2Span) * plotH;

        // Draw grids
        ctx.strokeStyle = 'rgba(148,163,184,0.08)';
        ctx.lineWidth = 1;
        ctx.setLineDash([4, 4]);

        const ticks = 5;
        ctx.fillStyle = '#64748b';
        ctx.font = '10px Inter, sans-serif';
        ctx.textAlign = 'center';
        ctx.textBaseline = 'top';
        for (let i = 0; i <= ticks; i++) {
            const val = minPc1 + (pc1Span / ticks) * i;
            const x = getX(val);
            ctx.beginPath(); ctx.moveTo(x, top); ctx.lineTo(x, top + plotH); ctx.stroke();
            ctx.fillText(val.toFixed(2), x, top + plotH + 8);
        }

        ctx.textAlign = 'right';
        ctx.textBaseline = 'middle';
        for (let i = 0; i <= ticks; i++) {
            const val = minPc2 + (pc2Span / ticks) * i;
            const y = getY(val);
            ctx.beginPath(); ctx.moveTo(left, y); ctx.lineTo(left + plotW, y); ctx.stroke();
            ctx.fillText(val.toFixed(2), left - 8, y);
        }
        ctx.setLineDash([]);

        // Axis Titles
        ctx.fillStyle = '#94a3b8';
        ctx.font = 'bold 11px Inter, sans-serif';
        ctx.textAlign = 'center';
        ctx.fillText('Principal Component 1 (PC1)', left + plotW / 2, this.height - 15);

        ctx.save();
        ctx.translate(14, top + plotH / 2);
        ctx.rotate(-Math.PI / 2);
        ctx.fillText('Principal Component 2 (PC2)', 0, 0);
        ctx.restore();

        // Draw points
        this.pcaData.forEach(p => {
            const x = getX(p.pc1);
            const y = getY(p.pc2);

            // Subpopulation colors
            let color = '#a855f7'; // purple / admixed default
            if (p.subpopulation === 'SP1') color = '#ef4444'; // red
            else if (p.subpopulation === 'SP2') color = '#10b981'; // green
            else if (p.subpopulation === 'SP3') color = '#3b82f6'; // blue

            const isHovered = this.hoveredPoint && this.hoveredPoint.accession === p.accession;
            let radius = isHovered ? 8 : 5.5;

            ctx.beginPath();
            ctx.arc(x, y, radius, 0, Math.PI * 2);
            ctx.fillStyle = color;
            ctx.fill();

            ctx.strokeStyle = '#0f1729';
            ctx.lineWidth = 1;
            ctx.stroke();

            if (isHovered) {
                ctx.strokeStyle = '#ffffff';
                ctx.lineWidth = 2;
                ctx.stroke();

                // Draw cultivar name floating above point
                ctx.fillStyle = '#ffffff';
                ctx.font = 'bold 11px Inter, sans-serif';
                ctx.textAlign = 'center';
                ctx.fillText(p.accession, x, y - 12);
            }
        });

        this._getX = getX;
        this._getY = getY;
    }

    setupInteraction() {
        if (this._listenersAttached || !this.canvas) return;
        this._listenersAttached = true;

        const handleMove = (e) => {
            const rect = this.canvas.getBoundingClientRect();
            const mx = e.clientX - rect.left;
            const my = e.clientY - rect.top;
            if (!this._getX || !this._getY) return;

            let nearest = null, nearestDist = Infinity;
            for (const p of this.pcaData) {
                const dist = Math.sqrt((mx - this._getX(p.pc1)) ** 2 + (my - this._getY(p.pc2)) ** 2);
                if (dist < 10 && dist < nearestDist) { nearest = p; nearestDist = dist; }
            }

            if (nearest !== this.hoveredPoint) {
                this.hoveredPoint = nearest;
                this.canvas.style.cursor = nearest ? 'pointer' : 'default';
                this.render();
                this.updateDetails(nearest);
            }
        };

        this.canvas.addEventListener('mousemove', handleMove);
        this.canvas.addEventListener('mouseleave', () => {
            this.hoveredPoint = null;
            this.canvas.style.cursor = 'default';
            this.render();
            this.updateDetails(null);
        });
    }

    updateDetails(p) {
        const emptyPanel = document.getElementById('pca-details-empty');
        const detailPanel = document.getElementById('pca-details-panel');

        if (!p) {
            emptyPanel.style.display = 'flex';
            detailPanel.style.display = 'none';
            return;
        }

        emptyPanel.style.display = 'none';
        detailPanel.style.display = 'block';

        const setVal = (id, val) => {
            const el = document.getElementById(id);
            if (el) el.textContent = val != null ? val : '—';
        };

        setVal('pca-cv-name', p.accession);
        setVal('pca-cv-subpop', p.subpopulation);
        setVal('pca-cv-brix', p.brix ? `${p.brix.toFixed(2)} %` : '—');
        setVal('pca-cv-weight', p.fruitWeight ? `${p.fruitWeight.toFixed(1)} g` : '—');
        setVal('pca-cv-pc1', p.pc1.toFixed(4));
        setVal('pca-cv-pc2', p.pc2.toFixed(4));
    }
}

/* ─────────────────────────────────────────────
 * Genomic Prediction Module
 * ─────────────────────────────────────────────*/

class PredictionModule {
    constructor() {
        this.scatterCanvas = document.getElementById('pred-scatter-canvas');
        this.featureCanvas = document.getElementById('pred-feature-canvas');
        this.runBtn = document.getElementById('btn-predict-run');
        this.spinner = document.getElementById('predict-spinner');
        this.metricsGrid = document.getElementById('prediction-metrics-grid');
        this.plotsRow = document.getElementById('prediction-plots-row');
        this.tableCard = document.getElementById('prediction-table-card');
        this.tableBody = document.getElementById('prediction-table-body');

        if (this.runBtn) {
            this.runBtn.addEventListener('click', () => this.runPrediction());
        }
    }

    async runPrediction() {
        const trait = document.getElementById('predict-trait').value;
        const model = document.getElementById('predict-model').value;

        if (this.spinner) this.spinner.style.display = 'inline-block';
        this.runBtn.disabled = true;

        try {
            const res = await fetch(`/api/predict?trait=${encodeURIComponent(trait)}&model_type=${encodeURIComponent(model)}`);
            if (!res.ok) throw new Error('Prediction API failed');
            const data = await res.json();

            if (data.error) {
                showNotification(data.error, 'error');
                return;
            }

            // Display metrics
            this.metricsGrid.style.display = 'grid';
            document.getElementById('pred-metric-r2').textContent = data.metrics.r2.toFixed(3);
            document.getElementById('pred-metric-corr').textContent = data.metrics.correlation.toFixed(3);
            document.getElementById('pred-metric-rmse').textContent = data.metrics.rmse.toFixed(2);
            document.getElementById('pred-metric-mae').textContent = data.metrics.mae.toFixed(2);

            // Render plots
            this.plotsRow.style.display = 'flex';
            this.drawScatterPlot(data.predictions, trait);
            this.drawFeaturesPlot(data.top_features);

            // Populate table
            this.tableCard.style.display = 'block';
            this.tableBody.innerHTML = '';
            data.predictions.forEach(p => {
                const dev = p.predicted - p.observed;
                const row = document.createElement('tr');
                row.innerHTML = `
                    <td style="font-weight: 700;">${escapeHtml(p.cultivar)}</td>
                    <td class="mono">${p.observed.toFixed(2)}</td>
                    <td class="mono">${p.predicted.toFixed(2)}</td>
                    <td class="mono" style="color: ${dev >= 0 ? '#4ade80' : '#f87171'}">${dev >= 0 ? '+' : ''}${dev.toFixed(2)}</td>
                `;
                this.tableBody.appendChild(row);
            });

            showNotification(`Genomic Prediction completed successfully using ${model}!`, 'success');
        } catch (err) {
            console.error('[Prediction] Run failed:', err);
            showNotification('Failed to execute genomic prediction model.', 'error');
        } finally {
            if (this.spinner) this.spinner.style.display = 'none';
            this.runBtn.disabled = false;
        }
    }

    drawScatterPlot(points, trait) {
        const canvas = this.scatterCanvas;
        const ctx = canvas.getContext('2d');
        const W = 500, H = 350;

        ctx.fillStyle = '#0f1729';
        ctx.fillRect(0, 0, W, H);

        const observedVals = points.map(p => p.observed);
        const predictedVals = points.map(p => p.predicted);

        const minVal = Math.min(...observedVals, ...predictedVals);
        const maxVal = Math.max(...observedVals, ...predictedVals);
        const span = maxVal - minVal || 1;

        const pad = { top: 30, right: 30, bottom: 45, left: 55 };
        const plotW = W - pad.left - pad.right;
        const plotH = H - pad.top - pad.bottom;

        const getX = (val) => pad.left + ((val - minVal) / span) * plotW;
        const getY = (val) => pad.top + plotH - ((val - minVal) / span) * plotH;

        // Draw 1-to-1 diagonal reference line
        ctx.strokeStyle = 'rgba(239, 68, 68, 0.5)';
        ctx.lineWidth = 1.5;
        ctx.setLineDash([5, 5]);
        ctx.beginPath();
        ctx.moveTo(getX(minVal), getY(minVal));
        ctx.lineTo(getX(maxVal), getY(maxVal));
        ctx.stroke();
        ctx.setLineDash([]);

        // Draw grid lines
        ctx.strokeStyle = 'rgba(148, 163, 184, 0.06)';
        ctx.lineWidth = 1;
        for (let i = 0; i <= 4; i++) {
            const val = minVal + (span / 4) * i;
            const x = getX(val);
            const y = getY(val);

            ctx.beginPath(); ctx.moveTo(x, pad.top); ctx.lineTo(x, pad.top + plotH); ctx.stroke();
            ctx.beginPath(); ctx.moveTo(pad.left, y); ctx.lineTo(pad.left + plotW, y); ctx.stroke();

            // Labels
            ctx.fillStyle = '#64748b';
            ctx.font = '9px Inter, sans-serif';
            ctx.textAlign = 'center';
            ctx.fillText(val.toFixed(1), x, pad.top + plotH + 12);

            ctx.textAlign = 'right';
            ctx.textBaseline = 'middle';
            ctx.fillText(val.toFixed(1), pad.left - 6, y);
        }

        // Axis titles
        ctx.fillStyle = '#94a3b8';
        ctx.font = 'bold 11px Inter, sans-serif';
        ctx.textAlign = 'center';
        ctx.fillText('Observed Phenotypic Values', pad.left + plotW / 2, H - 10);

        ctx.save();
        ctx.translate(12, pad.top + plotH / 2);
        ctx.rotate(-Math.PI / 2);
        ctx.fillText('Predicted Values', 0, 0);
        ctx.restore();

        // Store positions for hover
        this._scatterPoints = [];

        // Draw points
        points.forEach(p => {
            const px = getX(p.observed);
            const py = getY(p.predicted);
            ctx.beginPath();
            ctx.arc(px, py, 3.5, 0, Math.PI * 2);
            ctx.fillStyle = 'rgba(59, 130, 246, 0.75)';
            ctx.fill();
            this._scatterPoints.push({ x: px, y: py, cultivar: p.cultivar, observed: p.observed, predicted: p.predicted });
        });

        // Setup hover interaction (only once)
        if (!this._scatterListenersAttached && canvas) {
            this._scatterListenersAttached = true;
            this._scatterHovered = null;
            this._scatterGetX = getX;
            this._scatterGetY = getY;

            canvas.addEventListener('mousemove', (e) => {
                const rect = canvas.getBoundingClientRect();
                const mx = e.clientX - rect.left;
                const my = e.clientY - rect.top;
                if (!this._scatterPoints) return;

                let nearest = null, nearestDist = Infinity;
                for (const pt of this._scatterPoints) {
                    const dist = Math.sqrt((mx - pt.x) ** 2 + (my - pt.y) ** 2);
                    if (dist < 12 && dist < nearestDist) { nearest = pt; nearestDist = dist; }
                }
                if (nearest !== this._scatterHovered) {
                    this._scatterHovered = nearest;
                    canvas.style.cursor = nearest ? 'pointer' : 'default';
                    // Redraw with tooltip
                    this._redrawScatterWithTooltip(canvas, nearest);
                }
            });

            canvas.addEventListener('mouseleave', () => {
                this._scatterHovered = null;
                canvas.style.cursor = 'default';
                this._redrawScatterWithTooltip(canvas, null);
            });
        }
    }

    _redrawScatterWithTooltip(canvas, hovered) {
        if (!this._scatterPoints || this._scatterPoints.length === 0) return;
        const ctx = canvas.getContext('2d');
        const W = canvas.width, H = canvas.height;

        // Re-render base scatter (without full recalculation - reuse stored points)
        ctx.fillStyle = '#0f1729';
        ctx.fillRect(0, 0, W, H);

        const pad = { top: 30, right: 30, bottom: 45, left: 55 };
        const plotW = W - pad.left - pad.right;
        const plotH = H - pad.top - pad.bottom;

        // Grid & axes
        ctx.strokeStyle = 'rgba(148, 163, 184, 0.06)';
        ctx.lineWidth = 1;
        for (let i = 0; i <= 4; i++) {
            const frac = i / 4;
            const x = pad.left + frac * plotW;
            const y = pad.top + (1 - frac) * plotH;
            ctx.beginPath(); ctx.moveTo(x, pad.top); ctx.lineTo(x, pad.top + plotH); ctx.stroke();
            ctx.beginPath(); ctx.moveTo(pad.left, y); ctx.lineTo(pad.left + plotW, y); ctx.stroke();
        }

        // 1:1 line
        ctx.strokeStyle = 'rgba(239, 68, 68, 0.5)';
        ctx.lineWidth = 1.5;
        ctx.setLineDash([5, 5]);
        ctx.beginPath();
        ctx.moveTo(pad.left, pad.top + plotH);
        ctx.lineTo(pad.left + plotW, pad.top);
        ctx.stroke();
        ctx.setLineDash([]);

        // Points
        this._scatterPoints.forEach(pt => {
            const isHov = hovered && hovered.cultivar === pt.cultivar;
            ctx.beginPath();
            ctx.arc(pt.x, pt.y, isHov ? 6 : 3.5, 0, Math.PI * 2);
            ctx.fillStyle = isHov ? '#f59e0b' : 'rgba(59, 130, 246, 0.75)';
            ctx.fill();
            if (isHov) {
                ctx.strokeStyle = '#ffffff';
                ctx.lineWidth = 2;
                ctx.stroke();
            }
        });

        // Tooltip
        if (hovered) {
            const lines = [
                hovered.cultivar,
                `Obs: ${hovered.observed.toFixed(2)}`,
                `Pred: ${hovered.predicted.toFixed(2)}`
            ];
            ctx.font = '11px Inter, sans-serif';
            const lineH = 16;
            const padX = 10, padY = 6;
            const maxW = Math.max(...lines.map(l => ctx.measureText(l).width));
            const boxW = maxW + padX * 2;
            const boxH = lines.length * lineH + padY * 2;
            let tx = hovered.x + 12, ty = hovered.y - boxH / 2;
            if (tx + boxW > W - 10) tx = hovered.x - boxW - 12;
            if (ty < 5) ty = 5;

            ctx.fillStyle = 'rgba(15, 20, 40, 0.94)';
            ctx.beginPath();
            roundRect(ctx, tx, ty, boxW, boxH, 6);
            ctx.fill();
            ctx.strokeStyle = 'rgba(245, 158, 11, 0.5)';
            ctx.lineWidth = 1;
            ctx.stroke();

            ctx.fillStyle = '#e2e8f0';
            ctx.textAlign = 'left';
            ctx.textBaseline = 'top';
            lines.forEach((line, i) => {
                ctx.font = i === 0 ? 'bold 11px Inter, sans-serif' : '11px Inter, sans-serif';
                ctx.fillStyle = i === 0 ? '#fbbf24' : '#e2e8f0';
                ctx.fillText(line, tx + padX, ty + padY + i * lineH);
            });
        }
    }

    drawFeaturesPlot(features) {
        const ctx = this.featureCanvas.getContext('2d');
        const W = 500, H = 350;

        ctx.fillStyle = '#0f1729';
        ctx.fillRect(0, 0, W, H);

        const pad = { top: 20, right: 30, bottom: 40, left: 130 };
        const plotW = W - pad.left - pad.right;
        const plotH = H - pad.top - pad.bottom;

        const maxImportance = Math.max(...features.map(f => f.importance), 1e-5);
        const barH = 18;
        const gap = 8;

        features.forEach((feat, i) => {
            const y = pad.top + i * (barH + gap);
            const barW = (feat.importance / maxImportance) * plotW;

            // Bar background
            ctx.fillStyle = 'rgba(255, 255, 255, 0.03)';
            ctx.beginPath();
            roundRect(ctx, pad.left, y, plotW, barH, 3);
            ctx.fill();

            // Colored bar
            ctx.fillStyle = 'rgba(251, 191, 36, 0.85)';
            ctx.beginPath();
            roundRect(ctx, pad.left, y, barW, barH, 3);
            ctx.fill();

            // Label
            ctx.fillStyle = '#94a3b8';
            ctx.font = '10px Fira Code, monospace';
            ctx.textAlign = 'right';
            ctx.textBaseline = 'middle';
            ctx.fillText(feat.snp_id, pad.left - 10, y + barH / 2);
        });

        // X-axis line & label
        ctx.strokeStyle = 'rgba(148, 163, 184, 0.2)';
        ctx.lineWidth = 1;
        ctx.beginPath();
        ctx.moveTo(pad.left, pad.top + plotH);
        ctx.lineTo(pad.left + plotW, pad.top + plotH);
        ctx.stroke();

        ctx.fillStyle = '#64748b';
        ctx.font = 'bold 10px Inter, sans-serif';
        ctx.textAlign = 'center';
        ctx.fillText('Relative Effect / Importance Weight', pad.left + plotW / 2, H - 10);
    }
}

/* ─────────────────────────────────────────────
 * Breeder Recommendations Dashboard
 * ─────────────────────────────────────────────*/

class BreederModule {
    constructor() {
        this.sliders = {
            brix: document.getElementById('slider-brix'),
            pulp: document.getElementById('slider-pulp'),
            seed: document.getElementById('slider-seed'),
            stone: document.getElementById('slider-stone')
        };
        this.labels = {
            brix: document.getElementById('label-brix'),
            pulp: document.getElementById('label-pulp'),
            seed: document.getElementById('label-seed'),
            stone: document.getElementById('label-stone')
        };
        this.presets = {
            sweet: document.getElementById('preset-sweet'),
            edible: document.getElementById('preset-edible'),
            balanced: document.getElementById('preset-balanced')
        };
        this.tbody = document.getElementById('breeder-ranking-body');
        this.crossesContainer = document.getElementById('crosses-recommender-container');

        // Bind slider inputs
        Object.entries(this.sliders).forEach(([key, slider]) => {
            if (slider) {
                slider.addEventListener('input', () => {
                    this.labels[key].textContent = parseFloat(slider.value).toFixed(1);
                    this.debouncedUpdate();
                });
            }
        });

        // Bind presets
        if (this.presets.sweet) {
            this.presets.sweet.addEventListener('click', () => this.applyPreset(2.5, 0.5, 0.2, 0.2));
        }
        if (this.presets.edible) {
            this.presets.edible.addEventListener('click', () => this.applyPreset(0.5, 2.5, 1.0, 1.0));
        }
        if (this.presets.balanced) {
            this.presets.balanced.addEventListener('click', () => this.applyPreset(1.5, 1.5, 1.0, 1.0));
        }

        this.timeout = null;
        this.updateRankings();
    }

    applyPreset(b, p, sd, st) {
        this.sliders.brix.value = b;
        this.sliders.pulp.value = p;
        this.sliders.seed.value = sd;
        this.sliders.stone.value = st;

        Object.keys(this.sliders).forEach(k => {
            this.labels[k].textContent = parseFloat(this.sliders[k].value).toFixed(1);
        });

        this.updateRankings();
    }

    debouncedUpdate() {
        if (this.timeout) clearTimeout(this.timeout);
        this.timeout = setTimeout(() => this.updateRankings(), 300);
    }

    async updateRankings() {
        const body = {
            w_brix: parseFloat(this.sliders.brix.value),
            w_pulp: parseFloat(this.sliders.pulp.value),
            w_seed: parseFloat(this.sliders.seed.value),
            w_stone: parseFloat(this.sliders.stone.value)
        };

        try {
            const res = await fetch('/api/breeder/rank', {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify(body)
            });
            if (!res.ok) throw new Error('Ranking failed');
            const data = await res.json();

            this.renderRankings(data);
            this.renderCrosses(data);
        } catch (err) {
            console.error('[Breeder] Ranking failed:', err);
        }
    }

    renderRankings(data) {
        if (!this.tbody) return;
        this.tbody.innerHTML = '';

        data.forEach(item => {
            const row = document.createElement('tr');
            row.innerHTML = `
                <td style="font-weight: 700; color: #fbbf24;">#${item.rank}</td>
                <td style="font-weight: 700;">${escapeHtml(item.accession)}</td>
                <td><span style="color:#93c5fd;">${escapeHtml(item.subpopulation)}</span></td>
                <td class="mono">${item.brix.toFixed(2)}%</td>
                <td class="mono">${(item.pulp_ratio * 100).toFixed(1)}%</td>
                <td class="mono">${(item.seed_ratio * 100).toFixed(1)}%</td>
                <td class="mono">${(item.stone_ratio * 100).toFixed(1)}%</td>
                <td class="mono" style="font-weight: 800; color: #3b82f6;">${item.fqi_score.toFixed(2)}</td>
            `;
            this.tbody.appendChild(row);
        });
    }

    renderCrosses(data) {
        if (!this.crossesContainer) return;
        this.crossesContainer.innerHTML = '';

        // Select parents from top 15 list
        // Maternal parent: top FQI score variety
        // Paternal parent: complementary variety from a different subpopulation (SP) to maximize hybrid vigor
        const p1 = data[0];
        let p2 = null;
        for (let idx = 1; idx < data.length; idx++) {
            if (data[idx].subpopulation !== p1.subpopulation) {
                p2 = data[idx];
                break;
            }
        }
        if (!p2) p2 = data[1];

        // Another hybrid suggestion
        const p3 = data[2] || data[0];
        let p4 = null;
        for (let idx = 3; idx < data.length; idx++) {
            if (data[idx].subpopulation !== p3.subpopulation) {
                p4 = data[idx];
                break;
            }
        }
        if (!p4) p4 = data[3] || data[1];

        const suggestions = [
            { female: p1, male: p2, score: (p1.fqi_score + p2.fqi_score) / 2 },
            { female: p3, male: p4, score: (p3.fqi_score + p4.fqi_score) / 2 }
        ];

        suggestions.forEach((cross, i) => {
            const card = document.createElement('div');
            card.style.cssText = `
                background: rgba(255, 255, 255, 0.02);
                border: 1px solid rgba(255, 255, 255, 0.06);
                border-radius: 8px; padding: 14px;
                display: flex; flex-direction: column; gap: 8px;
            `;

            card.innerHTML = `
                <div style="display: flex; justify-content: space-between; align-items: center;">
                    <span style="font-weight: bold; color: #fbbf24;">Recommendation ${i + 1}</span>
                    <span class="badge badge-info" style="font-size:0.75rem;">F1 FQI: ${cross.score.toFixed(2)}</span>
                </div>
                <div style="font-size: 0.82rem; color: #e2e8f0; display: flex; flex-direction: column; gap: 4px;">
                    <div>🌸 <strong>Maternal Parent:</strong> ${escapeHtml(cross.female.accession)} (${cross.female.subpopulation})</div>
                    <div>🌾 <strong>Paternal Parent:</strong> ${escapeHtml(cross.male.accession)} (${cross.male.subpopulation})</div>
                </div>
                <div style="font-size: 0.76rem; color: #94a3b8; border-top: 1px solid rgba(255,255,255,0.06); padding-top: 6px; margin-top: 4px;">
                    Combines sweetness (${Math.max(cross.female.brix, cross.male.brix).toFixed(1)}% max) and edible pulp yield (${(Math.max(cross.female.pulp_ratio, cross.male.pulp_ratio) * 100).toFixed(0)}% max). Cross across lineages minimizes inbreeding depression.
                </div>
            `;
            this.crossesContainer.appendChild(card);
        });
    }
}

/* ─────────────────────────────────────────────
 * New Cultivar Prediction Module
 * ─────────────────────────────────────────────*/

class NewCultivarModule {
    constructor() {
        this.simBtn = document.getElementById('btn-new-cv-simulate');
        this.spinner = document.getElementById('new-cv-spinner');
        this.resultsContainer = document.getElementById('new-cv-results-container');
        this.snpCanvas = document.getElementById('new-cv-snp-canvas');
        this.comparisonCanvas = document.getElementById('new-cv-comparison-canvas');
        this.traitSelect = document.getElementById('new-cv-snp-trait-select');
        
        // File Upload items
        this.dropzone = document.getElementById('new-cv-upload-dropzone');
        this.fileInput = document.getElementById('new-cv-upload-input');
        this.fileNameSpan = document.getElementById('upload-file-name');
        this.statusInfoDiv = document.getElementById('upload-status-info');
        this.downloadReportBtn = document.getElementById('btn-download-report');
        
        this._lastData = null;

        if (this.simBtn) {
            this.simBtn.addEventListener('click', () => this.runSimulation());
        }
        if (this.traitSelect) {
            this.traitSelect.addEventListener('change', () => this.drawSnpChart());
        }
        
        // Setup Drag & Drop Upload
        if (this.dropzone && this.fileInput) {
            this.dropzone.addEventListener('click', () => this.fileInput.click());
            this.fileInput.addEventListener('change', (e) => this.handleFileSelect(e));
            
            this.dropzone.addEventListener('dragover', (e) => {
                e.preventDefault();
                this.dropzone.style.borderColor = '#14b8a6';
                this.dropzone.style.background = 'rgba(20, 184, 166, 0.05)';
            });
            
            this.dropzone.addEventListener('dragleave', () => {
                this.dropzone.style.borderColor = 'rgba(6, 182, 212, 0.25)';
                this.dropzone.style.background = 'rgba(6, 182, 212, 0.02)';
            });
            
            this.dropzone.addEventListener('drop', (e) => {
                e.preventDefault();
                this.dropzone.style.borderColor = 'rgba(6, 182, 212, 0.25)';
                this.dropzone.style.background = 'rgba(6, 182, 212, 0.02)';
                if (e.dataTransfer.files.length > 0) {
                    this.uploadFile(e.dataTransfer.files[0]);
                }
            });
        }
        
        // Report Download handler
        if (this.downloadReportBtn) {
            this.downloadReportBtn.addEventListener('click', () => this.downloadReport());
        }
    }

    handleFileSelect(e) {
        if (e.target.files.length > 0) {
            this.uploadFile(e.target.files[0]);
        }
    }

    async uploadFile(file) {
        const model = document.getElementById('new-cv-model').value;
        if (this.spinner) this.spinner.style.display = 'inline-block';
        if (this.fileNameSpan) {
            this.fileNameSpan.textContent = `Selected: ${file.name}`;
            this.fileNameSpan.style.display = 'block';
        }
        
        const formData = new FormData();
        formData.append('file', file);
        formData.append('model_type', model);
        
        try {
            const res = await fetch('/api/predict/upload', {
                method: 'POST',
                body: formData
            });
            
            if (!res.ok) {
                const errData = await res.json();
                throw new Error(errData.error || 'Genotype upload prediction failed');
            }
            
            const data = await res.json();
            this._lastData = data;
            this.renderResults(data, model);
            
            // Show alignment statistics
            if (this.statusInfoDiv && data.alignment_stats) {
                const stats = data.alignment_stats;
                let imputedHtml = '';
                if (stats.imputed_sample_markers && stats.imputed_sample_markers.length > 0) {
                    const sampleList = stats.imputed_sample_markers.slice(0, 10).map(m => `<code style="color: #fbbf24; font-size: 0.75rem;">${escapeHtml(m)}</code>`).join(', ');
                    const countExtra = stats.n_imputed > 10 ? ` and ${stats.n_imputed - 10} more...` : '';
                    imputedHtml = `
                        <div style="margin-top: 8px; border-top: 1px solid rgba(255,255,255,0.06); padding-top: 8px; font-size: 0.78rem; color: #94a3b8; text-align: left; line-height: 1.4;">
                            <strong>Sample Imputed Markers:</strong> ${sampleList}${countExtra}
                        </div>
                    `;
                }
                this.statusInfoDiv.innerHTML = `
                    <div style="background: rgba(20, 184, 166, 0.08); border: 1px solid rgba(20, 184, 166, 0.2); padding: 12px; border-radius: 8px; margin-top: 12px;">
                        <span style="font-weight: bold; color: #14b8a6;">Marker Match Summary:</span>
                        <div style="font-size: 0.8rem; color: #cbd5e1; margin-top: 4px; display: grid; grid-template-columns: 1fr 1fr; gap: 4px;">
                            <div>Matched SNPs: <strong>${stats.n_matched.toLocaleString()}</strong></div>
                            <div>Match Rate: <strong>${stats.match_percentage}%</strong></div>
                            <div>Mean-Imputed: <strong>${stats.n_imputed.toLocaleString()}</strong></div>
                            <div>Total Parsed: <strong>${stats.n_parsed.toLocaleString()}</strong></div>
                        </div>
                        ${imputedHtml}
                    </div>
                `;
                this.statusInfoDiv.style.display = 'block';
            }
            showNotification('Genotype file uploaded and predicted successfully!', 'success');
        } catch (err) {
            console.error('[NewCultivar] Upload failed:', err);
            showNotification(err.message || 'Failed to parse and predict uploaded genotypes.', 'error');
        } finally {
            if (this.spinner) this.spinner.style.display = 'none';
        }
    }

    downloadReport() {
        if (!this._lastData) {
            showNotification('No prediction data available to download.', 'warning');
            return;
        }
        const dataStr = encodeURIComponent(JSON.stringify(this._lastData));
        window.open(`/api/predict/report?data=${dataStr}`, '_blank');
    }

    async runSimulation() {
        const model = document.getElementById('new-cv-model').value;
        if (this.spinner) this.spinner.style.display = 'inline-block';
        this.simBtn.disabled = true;
        if (this.statusInfoDiv) this.statusInfoDiv.style.display = 'none';
        if (this.fileNameSpan) this.fileNameSpan.style.display = 'none';

        try {
            const res = await fetch('/api/predict/new_cultivar', {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify({ simulate: true, model_type: model })
            });
            if (!res.ok) throw new Error('Prediction API failed');
            const data = await res.json();

            if (data.error) {
                showNotification(data.error, 'error');
                return;
            }
            this._lastData = data;
            this.renderResults(data, model);
            showNotification('New cultivar phenotype predicted successfully!', 'success');
        } catch (err) {
            console.error('[NewCultivar] Prediction failed:', err);
            showNotification('Failed to predict new cultivar phenotype.', 'error');
        } finally {
            if (this.spinner) this.spinner.style.display = 'none';
            this.simBtn.disabled = false;
        }
    }

    renderResults(data, model) {
        this.resultsContainer.style.display = 'block';

        // Model badge
        const badge = document.getElementById('new-cv-model-badge');
        if (badge) badge.textContent = model;

        // Predicted Traits Grid
        const traitsGrid = document.getElementById('new-cv-traits-grid');
        traitsGrid.innerHTML = '';

        const traitNames = {
            'fruitWeight (g)': '⚖️ Fruit Weight',
            'fruitLength (mm)': '📏 Fruit Length',
            'fruitWidth (mm)': '↔️ Fruit Width',
            'fruitThickness (mm)': '🪵 Fruit Thickness',
            'brix': '🍯 Brix',
            'Pulp': '🍑 Pulp Weight',
            'seedWeight (g)': '🌱 Seed Weight',
            'stoneWeight (g)': '🪨 Stone Weight',
            'seedLength (mm)': '🌱 Seed Length',
            'seedWidth (mm)': '↔️ Seed Width',
            'seedThickness (mm)': '🪵 Seed Thickness',
            'stoneLength (mm)': '🪨 Stone Length',
            'stoneWidth (mm)': '↔️ Stone Width',
            'stoneThickness (mm)': '🪵 Stone Thickness',
            'Fruit Shape Index': '📐 Shape Index',
            'Pulp Ratio': '🍑 Pulp Ratio',
            'Seed Ratio': '🌱 Seed Ratio',
            'Stone Ratio': '🪨 Stone Ratio',
            'Edible Portion': '🍽️ Edible Portion',
            'Brix Yield Index': '🧪 Brix Yield Index',
            'Fruit Density Index': '⚖️ Fruit Density',
            'Pulp-to-Seed Ratio': '📊 Pulp-to-Seed Ratio',
            'Pulp-to-Stone Ratio': '📊 Pulp-to-Stone Ratio',
            'Sweetness Efficiency Index': '📊 Sweetness Efficiency Index'
        };

        for (const [trait, info] of Object.entries(data.predictions)) {
            const card = document.createElement('div');
            card.className = 'stat-card';
            const aboveMean = info.predicted > info.population_mean;
            const percentile = ((info.predicted - info.population_mean) / (info.population_std || 1)).toFixed(1);
            const arrow = aboveMean ? '▲' : '▼';
            const color = aboveMean ? '#4ade80' : '#f87171';

            card.innerHTML = `
                <span class="stat-number" style="font-size: 1.5rem;">${info.predicted.toFixed(2)}</span>
                <span class="stat-label">${traitNames[trait] || trait} ${info.unit ? `(${info.unit})` : ''}</span>
                <div style="font-size: 0.72rem; color: #64748b; margin-top: 6px;">
                    CI: [${info.ci_lower.toFixed(2)} – ${info.ci_upper.toFixed(2)}]
                </div>
                <div style="font-size: 0.72rem; color: ${color}; margin-top: 2px;">
                    ${arrow} ${percentile}σ from mean (μ=${info.population_mean.toFixed(1)})
                </div>
            `;
            traitsGrid.appendChild(card);
        }

        // Genetic Similarity Table
        const simBody = document.getElementById('new-cv-similarity-body');
        simBody.innerHTML = '';
        data.similar_cultivars.forEach((cv, i) => {
            const row = document.createElement('tr');
            const barWidth = Math.max(10, cv.similarity_pct);
            row.innerHTML = `
                <td style="font-weight: 700; color: #fbbf24;">#${i + 1}</td>
                <td style="font-weight: 700;">${escapeHtml(cv.cultivar)}</td>
                <td class="mono">${cv.distance.toFixed(1)}</td>
                <td>
                    <div style="display: flex; align-items: center; gap: 8px;">
                        <div style="background: linear-gradient(90deg, #6366f1, #818cf8); height: 8px; border-radius: 4px; width: ${barWidth}%; min-width: 20px;"></div>
                        <span class="mono" style="font-size: 0.8rem;">${cv.similarity_pct}%</span>
                    </div>
                </td>
            `;
            simBody.appendChild(row);
        });

        // SNP Trait Select
        this.traitSelect.innerHTML = '';
        for (const trait of Object.keys(data.top_snps_per_trait)) {
            const opt = document.createElement('option');
            opt.value = trait;
            opt.textContent = traitNames[trait] || trait;
            this.traitSelect.appendChild(opt);
        }
        this.drawSnpChart();

        // Population Comparison
        this.drawComparison(data);
    }

    drawSnpChart() {
        if (!this._lastData || !this.snpCanvas) return;
        const trait = this.traitSelect.value;
        const snps = this._lastData.top_snps_per_trait[trait];
        if (!snps || snps.length === 0) return;

        const canvas = this.snpCanvas;
        const ctx = canvas.getContext('2d');
        const W = canvas.width, H = canvas.height;

        ctx.fillStyle = '#0f1729';
        ctx.fillRect(0, 0, W, H);

        const pad = { top: 15, right: 30, bottom: 30, left: 130 };
        const plotW = W - pad.left - pad.right;
        const maxImp = Math.max(...snps.map(s => s.importance), 1e-8);
        const barH = 20;
        const gap = 10;

        snps.forEach((snp, i) => {
            const y = pad.top + i * (barH + gap);
            const barW = (snp.importance / maxImp) * plotW;

            // Background
            ctx.fillStyle = 'rgba(255, 255, 255, 0.03)';
            ctx.beginPath();
            roundRect(ctx, pad.left, y, plotW, barH, 3);
            ctx.fill();

            // Bar
            const grad = ctx.createLinearGradient(pad.left, 0, pad.left + barW, 0);
            grad.addColorStop(0, '#6366f1');
            grad.addColorStop(1, '#818cf8');
            ctx.fillStyle = grad;
            ctx.beginPath();
            roundRect(ctx, pad.left, y, barW, barH, 3);
            ctx.fill();

            // Label
            ctx.fillStyle = '#94a3b8';
            ctx.font = '10px Fira Code, monospace';
            ctx.textAlign = 'right';
            ctx.textBaseline = 'middle';
            ctx.fillText(snp.snp_id, pad.left - 8, y + barH / 2);
        });

        ctx.fillStyle = '#64748b';
        ctx.font = 'bold 10px Inter, sans-serif';
        ctx.textAlign = 'center';
        ctx.fillText('Relative Importance Weight', pad.left + plotW / 2, H - 8);
    }

    drawComparison(data) {
        if (!this.comparisonCanvas) return;
        const canvas = this.comparisonCanvas;
        const ctx = canvas.getContext('2d');
        const W = canvas.width, H = canvas.height;

        ctx.fillStyle = '#0f1729';
        ctx.fillRect(0, 0, W, H);

        // Display key traits as horizontal dot plots
        const keyTraits = ['fruitWeight (g)', 'brix', 'Pulp Ratio', 'Fruit Shape Index', 'Seed Ratio'];
        const displayNames = {
            'fruitWeight (g)': 'Fruit Weight',
            'brix': 'Brix',
            'Pulp Ratio': 'Pulp Ratio',
            'Fruit Shape Index': 'Shape Index',
            'Seed Ratio': 'Seed Ratio'
        };

        const available = keyTraits.filter(t => data.predictions[t]);
        const pad = { top: 20, right: 40, bottom: 25, left: 100 };
        const plotW = W - pad.left - pad.right;
        const rowH = (H - pad.top - pad.bottom) / (available.length || 1);

        available.forEach((trait, i) => {
            const info = data.predictions[trait];
            const y = pad.top + i * rowH + rowH / 2;

            // Draw population range bar
            const rangeMin = info.population_mean - 2 * info.population_std;
            const rangeMax = info.population_mean + 2 * info.population_std;
            const span = rangeMax - rangeMin || 1;

            const getX = (val) => pad.left + ((val - rangeMin) / span) * plotW;

            // Background range
            ctx.fillStyle = 'rgba(255, 255, 255, 0.04)';
            ctx.beginPath();
            roundRect(ctx, pad.left, y - 8, plotW, 16, 4);
            ctx.fill();

            // Population mean marker
            const meanX = getX(info.population_mean);
            ctx.fillStyle = 'rgba(148, 163, 184, 0.6)';
            ctx.beginPath();
            ctx.arc(meanX, y, 4, 0, Math.PI * 2);
            ctx.fill();

            // CI range bar
            const ciLeft = getX(Math.max(info.ci_lower, rangeMin));
            const ciRight = getX(Math.min(info.ci_upper, rangeMax));
            ctx.fillStyle = 'rgba(99, 102, 241, 0.3)';
            ctx.fillRect(ciLeft, y - 4, ciRight - ciLeft, 8);

            // Predicted value dot
            const predX = getX(Math.max(Math.min(info.predicted, rangeMax), rangeMin));
            ctx.beginPath();
            ctx.arc(predX, y, 7, 0, Math.PI * 2);
            ctx.fillStyle = '#818cf8';
            ctx.fill();
            ctx.strokeStyle = '#ffffff';
            ctx.lineWidth = 2;
            ctx.stroke();

            // Trait label
            ctx.fillStyle = '#94a3b8';
            ctx.font = '11px Inter, sans-serif';
            ctx.textAlign = 'right';
            ctx.textBaseline = 'middle';
            ctx.fillText(displayNames[trait] || trait, pad.left - 8, y);

            // Value label
            ctx.fillStyle = '#e2e8f0';
            ctx.font = 'bold 10px Inter, sans-serif';
            ctx.textAlign = 'left';
            ctx.fillText(info.predicted.toFixed(2), predX + 12, y);
        });

        // Legend
        ctx.fillStyle = '#64748b';
        ctx.font = '9px Inter, sans-serif';
        ctx.textAlign = 'center';
        ctx.fillText('● Population Mean    ● Predicted Value    ▬ 95% CI', W / 2, H - 6);
    }
}

/* ─────────────────────────────────────────────
 * Unseen Cultivar Validation Module
 * ─────────────────────────────────────────────*/

class ValidationModule {
    constructor() {
        this.holdoutSelect = document.getElementById('val-holdout-select');
        this.modelSelect = document.getElementById('val-model-select');
        this.runBtn = document.getElementById('btn-val-run');
        this.spinner = document.getElementById('val-spinner');
        this.loocvCard = document.getElementById('val-loocv-card');
        this.loocvGrid = document.getElementById('val-loocv-grid');
        this.resultsCard = document.getElementById('val-results-card');
        this.resultsTbody = document.getElementById('val-results-tbody');
        this.holdoutBadge = document.getElementById('val-holdout-badge');

        if (this.runBtn) {
            this.runBtn.addEventListener('click', () => this.runValidation());
        }
        
        this.loadCultivars();
    }

    async loadCultivars() {
        try {
            const res = await fetch('/api/cultivars');
            const data = await res.json();
            
            if (this.holdoutSelect) {
                this.holdoutSelect.innerHTML = '';
                data.forEach(cv => {
                    if (cv.has_genotype !== false) {
                        const opt = document.createElement('option');
                        opt.value = cv.id;
                        opt.textContent = cv.name;
                        this.holdoutSelect.appendChild(opt);
                    }
                });
            }
        } catch (err) {
            console.error('[Validation] Failed to load cultivars:', err);
        }
    }

    async runValidation() {
        const holdout = this.holdoutSelect.value;
        const model = this.modelSelect.value;
        
        if (this.spinner) this.spinner.style.display = 'inline-block';
        this.runBtn.disabled = true;

        try {
            const res = await fetch('/api/validation/holdout', {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify({ holdout_cultivar: holdout, model_type: model })
            });
            if (!res.ok) throw new Error('Holdout validation API failed');
            const data = await res.json();

            if (data.error) {
                showNotification(data.error, 'error');
                return;
            }
            
            this.renderLoocvMetrics(data.overall_loocv_metrics);
            this.renderHoldoutResults(data.predictions, holdout);
            showNotification('Holdout validation completed successfully!', 'success');
        } catch (err) {
            console.error('[Validation] Execution failed:', err);
            showNotification('Failed to execute holdout validation.', 'error');
        } finally {
            if (this.spinner) this.spinner.style.display = 'none';
            this.runBtn.disabled = false;
        }
    }

    renderLoocvMetrics(metrics) {
        if (!this.loocvGrid) return;
        this.loocvGrid.innerHTML = '';
        this.loocvCard.style.display = 'block';

        const displayNames = {
            'fruitWeight (g)': 'Fruit Weight',
            'brix': 'Brix (Sugar)',
            'Pulp Ratio': 'Pulp Ratio',
            'Fruit Shape Index': 'Shape Index'
        };

        for (const [trait, info] of Object.entries(metrics)) {
            const card = document.createElement('div');
            card.className = 'stat-card';
            card.style.textAlign = 'left';
            card.style.background = 'rgba(6, 182, 212, 0.03)';
            card.style.border = '1px solid rgba(6, 182, 212, 0.1)';
            
            card.innerHTML = `
                <span style="font-weight: bold; color: #14b8a6; font-size: 0.95rem; display: block; margin-bottom: 8px;">
                    ${displayNames[trait] || trait}
                </span>
                <div style="font-size: 0.82rem; color: #cbd5e1; display: grid; grid-template-columns: 1fr 1fr; gap: 6px;">
                    <div>LOO R²: <strong style="color: #fff;">${info.r2.toFixed(3)}</strong></div>
                    <div>Corr (r): <strong style="color: #fff;">${info.correlation.toFixed(3)}</strong></div>
                    <div>RMSE: <strong style="color: #fff;">${info.rmse.toFixed(2)}</strong></div>
                    <div>MAE: <strong style="color: #fff;">${info.mae.toFixed(2)}</strong></div>
                </div>
            `;
            this.loocvGrid.appendChild(card);
        }
    }

    renderHoldoutResults(predictions, holdoutName) {
        if (!this.resultsTbody) return;
        this.resultsTbody.innerHTML = '';
        this.resultsCard.style.display = 'block';
        
        if (this.holdoutBadge) this.holdoutBadge.textContent = holdoutName;

        const traitNames = {
            'fruitWeight (g)': '⚖️ Fruit Weight',
            'fruitLength (mm)': '📏 Fruit Length',
            'fruitWidth (mm)': '↔️ Fruit Width',
            'fruitThickness (mm)': '🪵 Fruit Thickness',
            'brix': '🍯 Brix',
            'Pulp': '🍑 Pulp Weight',
            'seedWeight (g)': '🌱 Seed Weight',
            'stoneWeight (g)': '🪨 Stone Weight',
            'seedLength (mm)': '🌱 Seed Length',
            'seedWidth (mm)': '↔️ Seed Width',
            'seedThickness (mm)': '🪵 Seed Thickness',
            'stoneLength (mm)': '🪨 Stone Length',
            'stoneWidth (mm)': '↔️ Stone Width',
            'stoneThickness (mm)': '🪵 Stone Thickness',
            'Fruit Shape Index': '📐 Shape Index',
            'Pulp Ratio': '🍑 Pulp Ratio',
            'Seed Ratio': '🌱 Seed Ratio',
            'Stone Ratio': '🪨 Stone Ratio',
            'Edible Portion': '🍽️ Edible Portion',
            'Brix Yield Index': '🧪 Brix Yield Index',
            'Fruit Density Index': '⚖️ Fruit Density',
            'Pulp-to-Seed Ratio': '📊 Pulp-to-Seed Ratio',
            'Pulp-to-Stone Ratio': '📊 Pulp-to-Stone Ratio',
            'Sweetness Efficiency Index': '📊 Sweetness Efficiency Index'
        };

        for (const [trait, info] of Object.entries(predictions)) {
            const row = document.createElement('tr');
            
            const actualVal = info.actual !== null ? info.actual.toFixed(2) : 'N/A';
            const predVal = info.predicted.toFixed(2);
            const devVal = info.deviation !== null ? (info.deviation > 0 ? '+' : '') + info.deviation.toFixed(2) : 'N/A';
            const devPct = info.deviation_percent !== null ? Math.abs(info.deviation_percent) : null;
            const devPctStr = info.deviation_percent !== null ? (info.deviation_percent > 0 ? '+' : '') + info.deviation_percent.toFixed(1) + '%' : 'N/A';
            
            let badgeClass = 'badge-success';
            let statusText = 'Excellent (≤5%)';
            
            if (devPct === null) {
                badgeClass = 'badge-muted';
                statusText = 'No Phenotype';
            } else if (devPct > 15) {
                badgeClass = 'badge-danger';
                statusText = 'High Error (>15%)';
            } else if (devPct > 5) {
                badgeClass = 'badge-warning';
                statusText = 'Moderate (5-15%)';
            }

            row.innerHTML = `
                <td style="font-weight: 700;">${traitNames[trait] || trait} ${info.unit ? `(${info.unit})` : ''}</td>
                <td class="mono" style="color:#94a3b8;">${actualVal}</td>
                <td class="mono" style="font-weight: bold; color: #fff;">${predVal}</td>
                <td class="mono" style="color: ${info.deviation > 0 ? '#10b981' : '#ef4444'};">${devVal}</td>
                <td class="mono" style="color: ${info.deviation_percent > 0 ? '#10b981' : '#ef4444'};">${devPctStr}</td>
                <td><span class="badge ${badgeClass}">${statusText}</span></td>
            `;
            this.resultsTbody.appendChild(row);
        }
    }
}

/* ─────────────────────────────────────────────
 * Breeding Cross Simulator Module
 * ─────────────────────────────────────────────*/

class BreedingSimulatorModule {
    constructor() {
        this.parentA = document.getElementById('cross-parent-a');
        this.parentB = document.getElementById('cross-parent-b');
        this.modelSelect = document.getElementById('cross-model');
        this.simBtn = document.getElementById('btn-cross-simulate');
        this.spinner = document.getElementById('cross-spinner');
        this.resultsContainer = document.getElementById('cross-results-container');
        this.resultsTbody = document.getElementById('cross-results-tbody');
        this.titleText = document.getElementById('cross-title-text');
        this.pcaSpan = document.getElementById('cross-pca-val');

        if (this.simBtn) {
            this.simBtn.addEventListener('click', () => this.runSimulation());
        }
        
        this.loadParents();
    }

    async loadParents() {
        try {
            const res = await fetch('/api/cultivars');
            const data = await res.json();
            
            [this.parentA, this.parentB].forEach((select, selectIdx) => {
                if (select) {
                    select.innerHTML = '';
                    data.forEach((cv, idx) => {
                        if (cv.has_genotype !== false) {
                            const opt = document.createElement('option');
                            opt.value = cv.id;
                            opt.textContent = cv.name;
                            if (selectIdx === 0 && cv.id === 'alphonso') opt.selected = true;
                            if (selectIdx === 1 && cv.id === 'tommy_atkins') opt.selected = true;
                            select.appendChild(opt);
                        }
                    });
                }
            });
        } catch (err) {
            console.error('[BreedingCross] Failed to load parents:', err);
        }
    }

    async runSimulation() {
        const pA = this.parentA.value;
        const pB = this.parentB.value;
        const model = this.modelSelect.value;
        
        if (pA === pB) {
            showNotification('Maternal and Paternal parents must be different cultivars.', 'warning');
            return;
        }

        if (this.spinner) this.spinner.style.display = 'inline-block';
        this.simBtn.disabled = true;

        try {
            const res = await fetch('/api/breeder/simulate_cross', {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify({ parent_a: pA, parent_b: pB, model_type: model })
            });
            if (!res.ok) throw new Error('Cross simulation failed');
            const data = await res.json();

            if (data.error) {
                showNotification(data.error, 'error');
                return;
            }
            
            this.renderCrossResults(data);
            showNotification('Genomic breeding cross simulated successfully!', 'success');
        } catch (err) {
            console.error('[BreedingCross] Simulation execution failed:', err);
            showNotification('Failed to simulate genomic breeding cross.', 'error');
        } finally {
            if (this.spinner) this.spinner.style.display = 'none';
            this.simBtn.disabled = false;
        }
    }

    renderCrossResults(data) {
        if (!this.resultsTbody) return;
        this.resultsTbody.innerHTML = '';
        this.resultsContainer.style.display = 'block';

        const pAName = this.parentA.options[this.parentA.selectedIndex].text;
        const pBName = this.parentB.options[this.parentB.selectedIndex].text;
        
        if (this.titleText) this.titleText.textContent = `${pAName} 🌸 × ${pBName} 🌾 (F1 Hybrid Outcome)`;
        if (this.pcaSpan) this.pcaSpan.textContent = data.pca_distance.toFixed(4);

        const traitNames = {
            'fruitWeight (g)': '⚖️ Fruit Weight',
            'fruitLength (mm)': '📏 Fruit Length',
            'fruitWidth (mm)': '↔️ Fruit Width',
            'fruitThickness (mm)': '🪵 Fruit Thickness',
            'brix': '🍯 Brix',
            'Pulp': '🍑 Pulp Weight',
            'seedWeight (g)': '🌱 Seed Weight',
            'stoneWeight (g)': '🪨 Stone Weight',
            'seedLength (mm)': '🌱 Seed Length',
            'seedWidth (mm)': '↔️ Seed Width',
            'seedThickness (mm)': '🪵 Seed Thickness',
            'stoneLength (mm)': '🪨 Stone Length',
            'stoneWidth (mm)': '↔️ Stone Width',
            'stoneThickness (mm)': '🪵 Stone Thickness',
            'Fruit Shape Index': '📐 Shape Index',
            'Pulp Ratio': '🍑 Pulp Ratio',
            'Seed Ratio': '🌱 Seed Ratio',
            'Stone Ratio': '🪨 Stone Ratio',
            'Edible Portion': '🍽️ Edible Portion',
            'Brix Yield Index': '🧪 Brix Yield Index',
            'Fruit Density Index': '⚖️ Fruit Density',
            'Pulp-to-Seed Ratio': '📊 Pulp-to-Seed Ratio',
            'Pulp-to-Stone Ratio': '📊 Pulp-to-Stone Ratio',
            'Sweetness Efficiency Index': '📊 Sweetness Efficiency Index'
        };

        function zToPercentile(z) {
            const t = 1 / (1 + 0.2316419 * Math.abs(z));
            const d = 0.3989423 * Math.exp(-z * z / 2);
            const p = d * t * (0.3193815 + t * (-0.3565638 + t * (1.781478 + t * (-1.821256 + t * 1.330274))));
            const cdf = z >= 0 ? 1 - p : p;
            return Math.round(cdf * 100);
        }

        for (const [trait, offspringInfo] of Object.entries(data.offspring_predictions)) {
            const pAVal = data.parent_a_predictions[trait].predicted;
            const pBVal = data.parent_b_predictions[trait].predicted;
            const offVal = offspringInfo.predicted;
            
            const midParent = (pAVal + pBVal) / 2;
            const heterosis = offVal > midParent ? '▲ Better (Hybrid Vigor)' : (offVal < midParent ? '▼ Lower' : '▬ Intermediate');
            const heterosisColor = offVal > midParent ? '#10b981' : (offVal < midParent ? '#ef4444' : '#cbd5e1');

            const mean = offspringInfo.population_mean;
            const std = offspringInfo.population_std || 1;
            const z = (offVal - mean) / std;
            const pctVal = zToPercentile(z);

            const row = document.createElement('tr');
            row.innerHTML = `
                <td style="font-weight: 700;">${traitNames[trait] || trait} ${offspringInfo.unit ? `(${offspringInfo.unit})` : ''}</td>
                <td class="mono" style="color: #94a3b8;">${pAVal.toFixed(2)}</td>
                <td class="mono" style="font-weight: 800; color: #14b8a6; font-size: 1.05rem;">${offVal.toFixed(2)}</td>
                <td class="mono" style="color: #fbbf24; font-weight: 700;">${pctVal}th</td>
                <td class="mono" style="color: #94a3b8;">${pBVal.toFixed(2)}</td>
                <td style="font-weight: 600; color: ${heterosisColor};">${heterosis}</td>
            `;
            this.resultsTbody.appendChild(row);
        }
    }
}



/* ───────────────────────────────────────────
 * Initialization & Event Binding
 * ─────────────────────────────────────────── */

let _selectedGwasTrait = 'fruitWeight (g)';

document.addEventListener('DOMContentLoaded', () => {
    // Initialise Manhattan
    window.manhattanPlot = new ManhattanPlot('manhattan-canvas');

    // Initialise Population Structure (PCA)
    window.pcaPlot = new PcaPlot('pca-canvas');

    // Bind PCA tab switch to load data
    const pcaTab = document.getElementById('tab-pca');
    if (pcaTab) {
        pcaTab.addEventListener('click', () => {
            setTimeout(() => {
                window.pcaPlot._resizeCanvas();
                window.pcaPlot.loadData();
            }, 50);
        });
    }

    // Initialise Genomic Prediction module
    window.predictionModule = new PredictionModule();

    // Initialise Breeder module
    window.breederModule = new BreederModule();

    // Initialise New Cultivar Prediction module
    window.newCultivarModule = new NewCultivarModule();



    // Initialise Validation module
    window.validationModule = new ValidationModule();

    // Initialise Breeding Simulator module
    window.breedingSimulatorModule = new BreedingSimulatorModule();

    // Bind breeder tab click to refresh rankings
    const breederTab = document.getElementById('tab-breeder');
    if (breederTab) {
        breederTab.addEventListener('click', () => {
            window.breederModule.updateRankings();
        });
    }

    // Trait buttons selection (21 traits)
    const traitBtns = document.querySelectorAll('.trait-btn');
    traitBtns.forEach(btn => {
        btn.addEventListener('click', () => {
            traitBtns.forEach(b => b.classList.remove('active'));
            btn.classList.add('active');

            _selectedGwasTrait = btn.getAttribute('data-trait') || 'fruitWeight (g)';

            const traitLabel = document.getElementById('gwas-trait-label');
            if (traitLabel) traitLabel.textContent = _selectedGwasTrait;

            const sigEl = document.getElementById('gwas-sig-snps');
            if (sigEl) sigEl.textContent = '—';
        });
    });

    // Run GWAS button
    const btnRun = document.getElementById('btn-gwas-run');
    if (btnRun) {
        btnRun.addEventListener('click', () => {
            if (!_selectedGwasTrait) {
                showNotification('Please select a trait first.', 'warning');
                return;
            }

            const placeholder = document.getElementById('manhattan-placeholder');
            if (placeholder) placeholder.style.display = 'none';

            const qqPlaceholder = document.getElementById('qq-placeholder');
            if (qqPlaceholder) qqPlaceholder.style.display = 'none';

            setTimeout(() => {
                window.manhattanPlot._resizeCanvas();
                window.manhattanPlot.loadData(_selectedGwasTrait);
            }, 50);
        });
    }

    // SNP Panel Close
    const btnClose = document.getElementById('snp-panel-close');
    if (btnClose) {
        btnClose.addEventListener('click', () => {
            const panel = document.getElementById('snp-detail-panel');
            if (panel) panel.style.display = 'none';
        });
    }
});

