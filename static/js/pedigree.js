/**
 * Mango Bioinformatics Lab — Genetic Pedigree & Ancestry Module
 *
 * Implements historical pedigree/ancestry tree explorer and
 * multi-generational hybrid cross (F1 -> F2) genomic prediction simulation.
 */

class PedigreeModule {
    constructor() {
        // Mode Selectors
        this.btnExplorer = document.getElementById('btn-pedigree-mode-explorer');
        this.btnDesigner = document.getElementById('btn-pedigree-mode-designer');
        this.viewExplorer = document.getElementById('pedigree-view-explorer');
        this.viewDesigner = document.getElementById('pedigree-view-designer');

        // Explorer Elements
        this.selectCultivar = document.getElementById('pedigree-select-cultivar');
        this.explorerNodes = document.getElementById('pedigree-explorer-nodes');
        this.explorerSvg = document.getElementById('pedigree-explorer-svg');

        // Designer Elements
        this.parentA = document.getElementById('pedigree-parent-a');
        this.parentB = document.getElementById('pedigree-parent-b');
        this.parentC = document.getElementById('pedigree-parent-c');
        this.modelSelect = document.getElementById('pedigree-model');
        this.btnSimulate = document.getElementById('btn-pedigree-simulate');
        this.spinner = document.getElementById('pedigree-spinner');
        this.designerNodes = document.getElementById('pedigree-designer-nodes');
        this.designerSvg = document.getElementById('pedigree-designer-svg');

        // State
        this.lineageData = null;
        this.currentMode = 'explorer'; // 'explorer' or 'designer'
        this.simulatedData = null;
        this.selectedExplorerCultivar = '';

        this.init();
    }

    async init() {
        // Toggle view events
        if (this.btnExplorer && this.btnDesigner) {
            this.btnExplorer.addEventListener('click', () => this.switchMode('explorer'));
            this.btnDesigner.addEventListener('click', () => this.switchMode('designer'));
        }

        // Run simulation event
        if (this.btnSimulate) {
            this.btnSimulate.addEventListener('click', () => this.runDesignerCross());
        }

        // Select cultivar change event
        if (this.selectCultivar) {
            this.selectCultivar.addEventListener('change', (e) => {
                this.selectedExplorerCultivar = e.target.value;
                this.renderExplorerTree();
            });
        }

        // Window resize event to redraw SVG connection lines
        window.addEventListener('resize', () => {
            if (this.currentMode === 'explorer') {
                this.renderExplorerTree();
            } else {
                this.renderDesignerTree();
            }
        });

        // Load baseline data
        await this.loadLineageAndCultivars();
    }

    switchMode(mode) {
        this.currentMode = mode;
        if (mode === 'explorer') {
            this.btnExplorer.className = 'btn btn-primary';
            this.btnDesigner.className = 'btn btn-outline';
            this.viewExplorer.style.display = 'block';
            this.viewDesigner.style.display = 'none';
            this.renderExplorerTree();
        } else {
            this.btnExplorer.className = 'btn btn-outline';
            this.btnDesigner.className = 'btn btn-primary';
            this.viewExplorer.style.display = 'none';
            this.viewDesigner.style.display = 'block';
            this.renderDesignerTree();
        }
    }

    async loadLineageAndCultivars() {
        try {
            // Fetch pedigree historical facts
            const res = await fetch('/api/pedigree/lineage');
            if (res.ok) {
                this.lineageData = await res.json();
            }

            // Populate explorer select
            if (this.selectCultivar) {
                this.selectCultivar.innerHTML = '';
                
                // Get cultivars from window.MangoApp
                let cultivars = window.MangoApp.cultivars;
                if (!cultivars || cultivars.length === 0) {
                    // Fallback fetch
                    const cRes = await fetch('/api/cultivars');
                    if (cRes.ok) cultivars = await cRes.json();
                }

                if (cultivars && cultivars.length > 0) {
                    cultivars.forEach(cv => {
                        const opt = document.createElement('option');
                        opt.value = cv.id;
                        opt.textContent = cv.name;
                        this.selectCultivar.appendChild(opt);
                    });

                    // Add options to designer dropdowns (filtering for has_genotype !== false)
                    const parentOptions = cultivars.filter(c => c.has_genotype !== false);
                    [this.parentA, this.parentB, this.parentC].forEach(select => {
                        if (select) {
                            select.innerHTML = '';
                            parentOptions.forEach(cv => {
                                const opt = document.createElement('option');
                                opt.value = cv.id;
                                opt.textContent = cv.name;
                                select.appendChild(opt);
                            });
                        }
                    });

                    // Set default choices
                    if (this.parentA && parentOptions.length > 0) this.parentA.value = 'sindhri';
                    if (this.parentB && parentOptions.length > 1) this.parentB.value = 'earlygold';
                    if (this.parentC && parentOptions.length > 2) this.parentC.value = 'alphonso';

                    this.selectedExplorerCultivar = cultivars[0].id;
                    this.renderExplorerTree();
                }
            }
        } catch (err) {
            console.error('[Pedigree] Load error:', err);
        }
    }

    renderExplorerTree() {
        if (!this.lineageData || !this.selectedExplorerCultivar) return;

        const cvId = this.selectedExplorerCultivar;
        const cv = this.lineageData[cvId] || { name: cvId, parents: ['unknown', 'unknown'], history: 'N/A', details: 'N/A' };
        
        // Clear previous SVGs and Nodes
        this.explorerSvg.innerHTML = '';
        this.explorerNodes.innerHTML = '';

        const width = this.explorerNodes.clientWidth;
        const height = this.explorerNodes.clientHeight;

        // Position nodes:
        // Offspring is center bottom
        const childNode = { id: 'child', name: cv.name, x: width / 2, y: height * 0.72, type: 'Selected Cultivar', details: cv.details };
        
        // Parents
        const parentAName = cv.parents[0];
        const parentBName = cv.parents[1];
        
        const parentANode = { id: 'parentA', name: this.formatName(parentAName), x: width * 0.28, y: height * 0.32, type: 'Maternal Ancestor' };
        const parentBNode = { id: 'parentB', name: this.formatName(parentBName), x: width * 0.72, y: height * 0.32, type: 'Paternal Ancestor' };

        const nodes = [childNode, parentANode, parentBNode];

        // If parents are known, draw grandparent links conceptually if available in lineageData
        if (this.lineageData[parentAName]) {
            const gp = this.lineageData[parentAName];
            nodes.push({ id: 'gpA1', name: this.formatName(gp.parents[0]), x: width * 0.12, y: height * 0.08, type: 'Grandparent' });
            nodes.push({ id: 'gpA2', name: this.formatName(gp.parents[1]), x: width * 0.38, y: height * 0.08, type: 'Grandparent' });
            this.drawConnection(width * 0.12, height * 0.08, parentANode.x, parentANode.y, this.explorerSvg);
            this.drawConnection(width * 0.38, height * 0.08, parentANode.x, parentANode.y, this.explorerSvg);
        }
        if (this.lineageData[parentBName]) {
            const gp = this.lineageData[parentBName];
            nodes.push({ id: 'gpB1', name: this.formatName(gp.parents[0]), x: width * 0.62, y: height * 0.08, type: 'Grandparent' });
            nodes.push({ id: 'gpB2', name: this.formatName(gp.parents[1]), x: width * 0.88, y: height * 0.08, type: 'Grandparent' });
            this.drawConnection(width * 0.62, height * 0.08, parentBNode.x, parentBNode.y, this.explorerSvg);
            this.drawConnection(width * 0.88, height * 0.08, parentBNode.x, parentBNode.y, this.explorerSvg);
        }

        // Draw connections for parents to child
        this.drawConnection(parentANode.x, parentANode.y, childNode.x, childNode.y, this.explorerSvg);
        this.drawConnection(parentBNode.x, parentBNode.y, childNode.x, childNode.y, this.explorerSvg);

        // Render nodes
        nodes.forEach(node => {
            const card = document.createElement('div');
            card.className = `pedigree-node-card ${node.id === 'child' ? 'node-highlight' : ''}`;
            card.style.left = `${node.x}px`;
            card.style.top = `${node.y}px`;
            
            // Generate some mock/real traits to display
            let traitsHtml = '';
            if (node.id === 'child') {
                // Try to find matching profile details
                const appCv = window.MangoApp.cultivars?.find(c => c.id === cvId);
                if (appCv) {
                    traitsHtml = `
                        <div class="pedigree-node-metrics">
                            <div class="pedigree-node-metric-row"><span>⚖️ Weight:</span><span class="pedigree-node-metric-val">${appCv.fruit_weight_g || 250}g</span></div>
                            <div class="pedigree-node-metric-row"><span>🍯 Brix:</span><span class="pedigree-node-metric-val">${appCv.brix?.toFixed(1) || 18.0}</span></div>
                            <div class="pedigree-node-metric-row"><span>🛡️ Score:</span><span class="pedigree-node-metric-val" style="color:#22c55e;">${appCv.resistance_score || 50}</span></div>
                        </div>
                    `;
                }
            }

            card.innerHTML = `
                <div class="pedigree-node-header" title="${node.name}">${node.name}</div>
                <div class="pedigree-node-type">${node.type}</div>
                ${traitsHtml}
            `;
            
            // Click to trace ancestors if they exist in lineageData
            if (node.id !== 'child' && this.lineageData[node.name.toLowerCase().replace(/ \(.+\)/, '')]) {
                card.style.borderStyle = 'dashed';
                card.addEventListener('click', () => {
                    this.selectCultivar.value = node.name.toLowerCase().replace(/ \(.+\)/, '');
                    this.selectedExplorerCultivar = this.selectCultivar.value;
                    this.renderExplorerTree();
                });
            }

            this.explorerNodes.appendChild(card);
        });
    }

    async runDesignerCross() {
        const pA = this.parentA.value;
        const pB = this.parentB.value;
        const pC = this.parentC.value;
        const model = this.modelSelect.value;

        if (pA === pB || pB === pC || pA === pC) {
            showNotification('Please select distinct parents for each node.', 'warning');
            return;
        }

        if (this.spinner) this.spinner.style.display = 'inline-block';
        this.btnSimulate.disabled = true;

        try {
            const res = await fetch('/api/pedigree/cross_f2', {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify({ parent_a: pA, parent_b: pB, parent_c: pC, model_type: model })
            });

            if (!res.ok) throw new Error('Simulation failed');
            const data = await res.json();
            
            if (data.error) {
                showNotification(data.error, 'error');
                return;
            }

            this.simulatedData = data;
            this.renderDesignerTree();
            showNotification('F2 generation cross tree simulated successfully!', 'success');
        } catch (err) {
            console.error('[Pedigree] Designer simulation failed:', err);
            showNotification('Genomic simulation failed.', 'error');
        } finally {
            if (this.spinner) this.spinner.style.display = 'none';
            this.btnSimulate.disabled = false;
        }
    }

    renderDesignerTree() {
        if (!this.simulatedData) {
            // Draw placeholder tree
            this.designerSvg.innerHTML = '';
            this.designerNodes.innerHTML = `
                <div style="text-align: center; color: #64748b; padding-top: 150px; font-size: 0.95rem;">
                    👈 Select three parents and click "Simulate F2 Tree" to generate the multi-generation genomic pedigree.
                </div>
            `;
            return;
        }

        const data = this.simulatedData;
        this.designerSvg.innerHTML = '';
        this.designerNodes.innerHTML = '';

        const width = this.designerNodes.clientWidth;
        const height = this.designerNodes.clientHeight;

        // Position coordinates:
        // Tier 1 (Grandparents)
        const nodeA = { id: 'parent_a', name: this.formatName(data.parent_a), x: width * 0.18, y: height * 0.16, type: 'Grandmother', metrics: data.parent_a_predictions };
        const nodeB = { id: 'parent_b', name: this.formatName(data.parent_b), x: width * 0.46, y: height * 0.16, type: 'Grandfather', metrics: data.parent_b_predictions };
        const nodeC = { id: 'parent_c', name: this.formatName(data.parent_c), x: width * 0.82, y: height * 0.16, type: 'Donor Parent', metrics: data.parent_c_predictions };

        // Tier 2 (F1 Hybrid)
        const nodeF1 = { id: 'f1', name: 'F1 Hybrid (A × B)', x: width * 0.32, y: height * 0.50, type: 'F1 Generation', metrics: data.f1_predictions };

        // Tier 3 (F2 Hybrid)
        const nodeF2 = { id: 'f2', name: 'F2 Hybrid (F1 × C)', x: width * 0.57, y: height * 0.82, type: 'F2 Generation', metrics: data.f2_predictions };

        const nodes = [nodeA, nodeB, nodeC, nodeF1, nodeF2];

        // Draw connections
        this.drawConnection(nodeA.x, nodeA.y, nodeF1.x, nodeF1.y, this.designerSvg);
        this.drawConnection(nodeB.x, nodeB.y, nodeF1.x, nodeF1.y, this.designerSvg);
        this.drawConnection(nodeF1.x, nodeF1.y, nodeF2.x, nodeF2.y, this.designerSvg, true);
        this.drawConnection(nodeC.x, nodeC.y, nodeF2.x, nodeF2.y, this.designerSvg, true);

        // Render nodes
        nodes.forEach(node => {
            const card = document.createElement('div');
            const isHybrid = (node.id === 'f1' || node.id === 'f2');
            card.className = `pedigree-node-card ${isHybrid ? 'node-hybrid' : ''} ${node.id === 'f2' ? 'node-highlight' : ''}`;
            card.style.left = `${node.x}px`;
            card.style.top = `${node.y}px`;

            // Read predicted metrics
            const weightVal = node.metrics['fruitWeight (g)']?.predicted || 0;
            const brixVal = node.metrics['brix']?.predicted || 0;
            const edibleVal = node.metrics['Edible Portion']?.predicted || 0;

            card.innerHTML = `
                <div class="pedigree-node-header" title="${node.name}">${node.name}</div>
                <div class="pedigree-node-type">${node.type}</div>
                <div class="pedigree-node-metrics">
                    <div class="pedigree-node-metric-row"><span>⚖️ Fruit Weight:</span><span class="pedigree-node-metric-val">${weightVal.toFixed(0)}g</span></div>
                    <div class="pedigree-node-metric-row"><span>🍯 Brix sugar:</span><span class="pedigree-node-metric-val">${brixVal.toFixed(1)}</span></div>
                    <div class="pedigree-node-metric-row"><span>🍽️ Edible %:</span><span class="pedigree-node-metric-val">${edibleVal.toFixed(1)}%</span></div>
                </div>
            `;

            this.designerNodes.appendChild(card);
        });
    }

    drawConnection(x1, y1, x2, y2, svg, isHybrid = false) {
        // Create an SVG line or path
        const path = document.createElementNS('http://www.w3.org/2000/svg', 'path');
        
        // Draw an elegant S-curve path connecting nodes
        const midY = (y1 + y2) / 2;
        const d = `M ${x1} ${y1} C ${x1} ${midY}, ${x2} ${midY}, ${x2} ${y2}`;
        
        path.setAttribute('d', d);
        path.setAttribute('class', `pedigree-line ${isHybrid ? 'line-hybrid' : 'pedigree-tree-line'}`);
        path.style.cssText = `
            stroke: ${isHybrid ? 'rgba(168, 85, 247, 0.45)' : 'rgba(20, 184, 166, 0.45)'};
            stroke-width: 2.5px;
            fill: none;
            stroke-dasharray: 6;
            animation: pedigreeDash 30s linear infinite;
        `;

        // Add arrowhead marker (conceptual, using SVG line marker or a small circle indicator)
        svg.appendChild(path);
    }

    formatName(str) {
        if (!str) return 'Unknown';
        // Capitalize words
        return str
            .replace(/_/g, ' ')
            .split(' ')
            .map(w => w.charAt(0).toUpperCase() + w.slice(1))
            .join(' ');
    }
}

// Bind to DOMContentLoaded to initialize
document.addEventListener('DOMContentLoaded', () => {
    window.pedigreeModule = new PedigreeModule();
});
