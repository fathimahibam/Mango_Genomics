"""
gwas_simulator.py
=================
Generates realistic GWAS (Genome-Wide Association Study) SNP association data
for mango disease-resistance traits and writes it to gwas_data.json.

Two traits are simulated:
  - Anthracnose Resistance
  - Powdery Mildew Resistance

Signal peaks (highly significant p-values) are placed near known candidate
resistance loci, while background SNPs follow a uniform p-value distribution.
"""

import sys
sys.path.insert(0, r'C:\mangoproject\libs')

import os
import json
import numpy as np

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
OUTPUT_DIR  = r'C:\mangoproject\static\data'
OUTPUT_FILE = os.path.join(OUTPUT_DIR, 'gwas_data.json')

NUM_SNPS       = 5000
NUM_CHROMOSOMES = 20

# Chromosome lengths (bp) – synthetic but realistic for a ~400 Mb mango genome
np.random.seed(42)
CHROMOSOME_LENGTHS = {
    f'Chr{i}': int(np.random.uniform(15_000_000, 30_000_000))
    for i in range(1, NUM_CHROMOSOMES + 1)
}

# Traits
TRAITS = ['Anthracnose Resistance', 'Powdery Mildew Resistance']

# Peak definitions: each peak is a region where we plant highly significant SNPs
PEAKS = [
    # Primary peaks
    {
        'chromosome': 'Chr5',
        'center':     12_500_000,
        'width':       500_000,      # ±500 kb around center
        'gene':       'Mi-WAK1',
        'trait':      'Anthracnose Resistance',
        'n_snps':      35,           # number of significant SNPs in peak
    },
    {
        'chromosome': 'Chr3',
        'center':      8_200_000,
        'width':       400_000,
        'gene':       'Mi-MLO1',
        'trait':      'Powdery Mildew Resistance',
        'n_snps':      30,
    },
    # Secondary peaks
    {
        'chromosome': 'Chr7',
        'center':      5_000_000,
        'width':       350_000,
        'gene':       'Mi-RPS2',
        'trait':      'Anthracnose Resistance',
        'n_snps':      20,
    },
    {
        'chromosome': 'Chr11',
        'center':     15_000_000,
        'width':       300_000,
        'gene':       'Mi-TLP1',
        'trait':      'Powdery Mildew Resistance',
        'n_snps':      18,
    },
]


# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

def generate_peak_snps(peak, start_id):
    """Generate highly-significant SNPs clustered around a GWAS signal peak."""
    snps = []
    for i in range(peak['n_snps']):
        # Position follows a normal distribution centred on the peak
        pos = int(np.random.normal(peak['center'], peak['width'] / 3))
        pos = max(1, min(pos, CHROMOSOME_LENGTHS[peak['chromosome']]))

        # Very significant p-values (1e-12 to 5e-8)
        log_p = np.random.uniform(np.log10(1e-12), np.log10(5e-8))
        p_value = 10 ** log_p

        # Effect size correlates inversely with p-value
        beta = np.random.uniform(0.3, 1.2) * (-1 if np.random.random() < 0.3 else 1)

        snps.append({
            'id':            f'rs_{start_id + i:05d}',
            'chromosome':    peak['chromosome'],
            'position':      pos,
            'p_value':       float(p_value),
            'beta':          round(float(beta), 4),
            'allele_freq':   round(float(np.random.uniform(0.05, 0.50)), 4),
            'nearest_gene':  peak['gene'],
            'trait':         peak['trait'],
        })
    return snps


def generate_background_snps(n, start_id, existing_positions):
    """Generate background SNPs with non-significant p-values."""
    snps = []
    chromosomes = list(CHROMOSOME_LENGTHS.keys())

    for i in range(n):
        chrom = chromosomes[i % NUM_CHROMOSOMES]
        # Random position avoiding exact duplicates
        pos = int(np.random.uniform(100_000, CHROMOSOME_LENGTHS[chrom] - 100_000))

        # Background p-values: 0.01 – 1.0 (non-significant)
        p_value = float(np.random.uniform(0.01, 1.0))

        # Small effect sizes for background
        beta = round(float(np.random.normal(0, 0.05)), 4)

        # Pick a random trait
        trait = TRAITS[i % len(TRAITS)]

        # Assign a generic nearest-gene label
        nearest_gene = f'Mi-LOC{np.random.randint(1000, 9999)}'

        snps.append({
            'id':            f'rs_{start_id + i:05d}',
            'chromosome':    chrom,
            'position':      pos,
            'p_value':       p_value,
            'beta':          beta,
            'allele_freq':   round(float(np.random.uniform(0.01, 0.50)), 4),
            'nearest_gene':  nearest_gene,
            'trait':         trait,
        })
    return snps


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print("=== Mango GWAS Data Simulator ===")
    print(f"Generating {NUM_SNPS} SNPs across {NUM_CHROMOSOMES} chromosomes...\n")

    all_snps = []
    current_id = 1

    # 1. Generate peak SNPs
    for peak in PEAKS:
        peak_snps = generate_peak_snps(peak, current_id)
        all_snps.extend(peak_snps)
        current_id += len(peak_snps)
        print(f"  [+] Peak: {peak['gene']} on {peak['chromosome']} - "
              f"{len(peak_snps)} significant SNPs (trait: {peak['trait']})")

    # 2. Fill the rest with background SNPs
    n_background = NUM_SNPS - len(all_snps)
    existing_positions = {(s['chromosome'], s['position']) for s in all_snps}
    background = generate_background_snps(n_background, current_id, existing_positions)
    all_snps.extend(background)
    print(f"  [+] Background: {n_background} SNPs across all chromosomes\n")

    # 3. Sort by chromosome then position
    all_snps.sort(key=lambda s: (s['chromosome'], s['position']))

    # 4. Write output
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    with open(OUTPUT_FILE, 'w', encoding='utf-8') as f:
        json.dump(all_snps, f, indent=2)

    print(f"Wrote {len(all_snps)} SNPs -> {OUTPUT_FILE}")

    # Quick summary statistics
    sig = [s for s in all_snps if s['p_value'] < 5e-8]
    print(f"  Genome-wide significant (p < 5e-8): {len(sig)}")
    for trait in TRAITS:
        t_sig = [s for s in sig if s['trait'] == trait]
        print(f"    {trait}: {len(t_sig)} significant SNPs")
    print("\nDone.")


if __name__ == '__main__':
    main()
