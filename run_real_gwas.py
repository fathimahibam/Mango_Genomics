import os
import sys
sys.path.insert(0, r'C:\mangoproject\libs')

import json
import numpy as np
import pandas as pd
import math
import time
import re
import gffutils
from sklearn.decomposition import PCA
import bisect

def load_gwas_annotations(gff_path):
    print("Loading detailed genomic annotations from GFF...")
    t0 = time.time()
    genes = {}
    exons = {}
    cds = {}
    
    with open(gff_path, 'r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
                
            seqid = parts[0]
            featuretype = parts[2]
            start = int(parts[3])
            end = int(parts[4])
            strand = parts[6]
            
            # Parse attributes
            attrs = {}
            for item in parts[8].split(';'):
                if '=' in item:
                    k, v = item.split('=', 1)
                    attrs[k.strip()] = v.strip()
                    
            if featuretype == 'gene':
                gene_name = attrs.get('gene', attrs.get('Name', attrs.get('locus_tag', 'Unknown')))
                biotype = attrs.get('gene_biotype', 'Unknown')
                if seqid not in genes:
                    genes[seqid] = []
                genes[seqid].append({
                    'name': gene_name,
                    'start': start,
                    'end': end,
                    'strand': strand,
                    'biotype': biotype,
                    'mid': (start + end) / 2
                })
            elif featuretype == 'exon':
                gene_name = attrs.get('gene', '')
                if gene_name:
                    if seqid not in exons:
                        exons[seqid] = []
                    exons[seqid].append({
                        'gene': gene_name,
                        'start': start,
                        'end': end
                    })
            elif featuretype == 'CDS':
                gene_name = attrs.get('gene', '')
                if gene_name:
                    if seqid not in cds:
                        cds[seqid] = []
                    cds[seqid].append({
                        'gene': gene_name,
                        'start': start,
                        'end': end
                    })
                    
    genes_data = {}
    for seqid in genes:
        genes[seqid].sort(key=lambda x: x['start'])
        genes_data[seqid] = {
            'features': genes[seqid],
            'starts': [x['start'] for x in genes[seqid]],
            'mids': [x['mid'] for x in genes[seqid]]
        }
        
    exons_data = {}
    for seqid in exons:
        exons[seqid].sort(key=lambda x: x['start'])
        exons_data[seqid] = {
            'features': exons[seqid],
            'starts': [x['start'] for x in exons[seqid]]
        }
        
    cds_data = {}
    for seqid in cds:
        cds[seqid].sort(key=lambda x: x['start'])
        cds_data[seqid] = {
            'features': cds[seqid],
            'starts': [x['start'] for x in cds[seqid]]
        }
        
    print(f"Genomic annotations loaded in {time.time()-t0:.2f}s")
    return genes_data, exons_data, cds_data

def load_gene_products_gff(gff_path):
    print("Loading gene product descriptions from GFF...")
    t0 = time.time()
    products = {}
    with open(gff_path, 'r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            if line.startswith('#'):
                continue
            if 'Parent=gene-LOC' in line:
                parts = line.strip().split('\t')
                if len(parts) < 9 or parts[2] not in ('mRNA', 'CDS', 'lnc_RNA', 'tRNA', 'ncRNA'):
                    continue
                attrs = {}
                for item in parts[8].split(';'):
                    if '=' in item:
                        k, v = item.split('=', 1)
                        attrs[k.strip()] = v.strip()
                parent = attrs.get('Parent', '')
                if parent.startswith('gene-'):
                    gene_id = parent[5:]
                    prod = attrs.get('product', '')
                    if prod:
                        import urllib.parse
                        products[gene_id] = urllib.parse.unquote(prod)
    print(f"Loaded {len(products)} gene products in {time.time()-t0:.2f}s")
    return products

def find_overlapping_features(feature_data, pos):
    if not feature_data:
        return None
    features_list = feature_data['features']
    starts = feature_data['starts']
    idx = bisect.bisect_right(starts, pos)
    
    for i in range(idx - 1, -1, -1):
        f = features_list[i]
        if pos - f['start'] > 250000:
            break
        if f['start'] <= pos <= f['end']:
            return f
    return None

def find_promoter_overlap(genes_data, pos):
    if not genes_data:
        return None
    genes_list = genes_data['features']
    starts = genes_data['starts']
    idx = bisect.bisect_right(starts, pos)
    
    # Check plus strand genes starting after pos
    for i in range(idx, min(idx + 10, len(genes_list))):
        g = genes_list[i]
        if g['strand'] == '+' and g['start'] - 2000 <= pos < g['start']:
            return g
            
    # Check minus strand genes starting before pos
    for i in range(idx - 1, -1, -1):
        g = genes_list[i]
        if pos - g['start'] > 250000:
            break
        if g['strand'] == '-' and g['end'] < pos <= g['end'] + 2000:
            return g
            
    return None

def find_nearest_gene_midpoint(genes_data, pos):
    if not genes_data:
        return "Unknown", "Unknown", 0
    genes_list = genes_data['features']
    mids = genes_data['mids']
    idx = bisect.bisect_left(mids, pos)
    
    best_gene = "Unknown"
    best_biotype = "Unknown"
    min_dist = float('inf')
    
    for i in [idx - 1, idx, idx + 1]:
        if 0 <= i < len(genes_list):
            dist = abs(genes_list[i]['mid'] - pos)
            if dist < min_dist:
                min_dist = dist
                best_gene = genes_list[i]['name']
                best_biotype = genes_list[i]['biotype']
                
    return best_gene, best_biotype, int(min_dist)

def annotate_snp_context(genes, exons, cds, seqid, pos):
    chr_genes = genes.get(seqid, None)
    chr_exons = exons.get(seqid, None)
    chr_cds = cds.get(seqid, None)
    
    # 1. Check CDS overlap
    c_match = find_overlapping_features(chr_cds, pos)
    if c_match:
        biotype = "protein_coding"
        if chr_genes:
            for g in chr_genes['features']:
                if g['name'] == c_match['gene']:
                    biotype = g['biotype']
                    break
        return c_match['gene'], "CDS", 0, biotype
        
    # 2. Check Exon overlap (UTR/Exonic)
    e_match = find_overlapping_features(chr_exons, pos)
    if e_match:
        biotype = "Unknown"
        if chr_genes:
            for g in chr_genes['features']:
                if g['name'] == e_match['gene']:
                    biotype = g['biotype']
                    break
        return e_match['gene'], "UTR/Exonic", 0, biotype
        
    # 3. Check Gene body overlap (Intronic)
    g_match = find_overlapping_features(chr_genes, pos)
    if g_match:
        return g_match['name'], "Intronic", 0, g_match['biotype']
        
    # 4. Check Promoter overlap
    p_match = find_promoter_overlap(chr_genes, pos)
    if p_match:
        if p_match['strand'] == '+':
            dist = p_match['start'] - pos
        else:
            dist = pos - p_match['end']
        return p_match['name'], "Promoter", int(dist), p_match['biotype']
        
    # 5. Intergenic
    nearest_gene, biotype, dist = find_nearest_gene_midpoint(chr_genes, pos)
    return nearest_gene, "Intergenic", int(dist), biotype

def is_informative_product(prod, clean_id):
    if not prod:
        return False
    prod_lower = prod.lower().strip()
    if prod_lower == 'uncharacterized':
        return False
    if prod_lower == f'uncharacterized {clean_id.lower()}':
        return False
    return True

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
OUTPUT_DIR  = r'C:\mangoproject\static\data'
GWAS_FILE   = os.path.join(OUTPUT_DIR, 'gwas_data.json')
PCA_FILE    = os.path.join(OUTPUT_DIR, 'pca_data.json')
CSV_FILE    = r'C:\mangoproject\SNP.csv'
PHENO_FILE  = r'C:\mangoproject\phenotypes.xlsx'
DB_FILE     = r'C:\mangoproject\mango_genes_real.db'

REFSEQ_TO_CHR = {
    'NC_058137.1': 'Chr1',  'NC_058138.1': 'Chr2',  'NC_058139.1': 'Chr3',
    'NC_058140.1': 'Chr4',  'NC_058141.1': 'Chr5',  'NC_058142.1': 'Chr6',
    'NC_058143.1': 'Chr7',  'NC_058144.1': 'Chr8',  'NC_058145.1': 'Chr9',
    'NC_058146.1': 'Chr10', 'NC_058147.1': 'Chr11', 'NC_058148.1': 'Chr12',
    'NC_058149.1': 'Chr13', 'NC_058150.1': 'Chr14', 'NC_058151.1': 'Chr15',
    'NC_058152.1': 'Chr16', 'NC_058153.1': 'Chr17', 'NC_058154.1': 'Chr18',
    'NC_058155.1': 'Chr19', 'NC_058156.1': 'Chr20',
}
CHR_TO_REFSEQ = {v: k for k, v in REFSEQ_TO_CHR.items()}

def clean_id(s):
    return re.sub(r'[^a-zA-Z0-9]', '', str(s)).lower()

def fast_gwas(X, y):
    n = X.shape[1]
    X_mean = np.mean(X, axis=1, keepdims=True)
    X_centered = X - X_mean
    y_centered = y - np.mean(y)
    
    var_X = np.sum(X_centered**2, axis=1)
    var_X[var_X == 0] = 1.0
    var_y = np.sum(y_centered**2)
    if var_y == 0:
        var_y = 1.0
        
    cov_Xy = np.dot(X_centered, y_centered)
    
    r = cov_Xy / np.sqrt(var_X * var_y)
    r = np.clip(r, -0.99999, 0.99999)
    
    t_stat = r * np.sqrt((n - 2) / (1 - r**2))
    p_values = np.array([math.erfc(np.abs(t) / math.sqrt(2)) for t in t_stat])
    betas = cov_Xy / var_X
    
    return betas, p_values

def main():
    t0 = time.time()
    print("Loading database connections...")
    db = gffutils.FeatureDB(DB_FILE)
    cur = db.conn.cursor()
    
    gff_genes, gff_exons, gff_cds = load_gwas_annotations(r'C:\mangoproject\mango_real.gff')
    gff_products = load_gene_products_gff(r'C:\mangoproject\mango_real.gff')
    
    print("Loading phenotypes...")
    df_pheno = pd.read_excel(PHENO_FILE, sheet_name=' Table S1', header=2, keep_default_na=False).iloc[:161]
    df_pheno['clean_id'] = df_pheno['Accession ID'].apply(clean_id)
    
    print("Calculating derived traits...")
    df_pheno['fruitLength (mm)'] = pd.to_numeric(df_pheno['fruitLength (mm)'], errors='coerce')
    df_pheno['fruitWidth (mm)'] = pd.to_numeric(df_pheno['fruitWidth (mm)'], errors='coerce')
    df_pheno['fruitWeight (g)'] = pd.to_numeric(df_pheno['fruitWeight (g)'], errors='coerce')
    df_pheno['Pulp'] = pd.to_numeric(df_pheno['Pulp'], errors='coerce')
    df_pheno['seedWeight (g)'] = pd.to_numeric(df_pheno['seedWeight (g)'], errors='coerce')
    df_pheno['stoneWeight (g)'] = pd.to_numeric(df_pheno['stoneWeight (g)'], errors='coerce')
    df_pheno['brix'] = pd.to_numeric(df_pheno['brix'], errors='coerce')

    df_pheno['Fruit Shape Index'] = df_pheno['fruitLength (mm)'] / df_pheno['fruitWidth (mm)']
    df_pheno['Pulp Ratio'] = df_pheno['Pulp'] / df_pheno['fruitWeight (g)']
    df_pheno['Seed Ratio'] = df_pheno['seedWeight (g)'] / df_pheno['fruitWeight (g)']
    df_pheno['Stone Ratio'] = df_pheno['stoneWeight (g)'] / df_pheno['fruitWeight (g)']
    # New derived traits
    df_pheno['fruitThickness (mm)'] = pd.to_numeric(df_pheno['fruitThickness (mm)'], errors='coerce')
    df_pheno['Edible Portion'] = (df_pheno['Pulp'] / df_pheno['fruitWeight (g)']) * 100
    df_pheno['Brix Yield Index'] = df_pheno['brix'] * (df_pheno['Pulp'] / df_pheno['fruitWeight (g)'])
    vol_proxy = df_pheno['fruitLength (mm)'] * df_pheno['fruitWidth (mm)'] * df_pheno['fruitThickness (mm)']
    df_pheno['Fruit Density Index'] = df_pheno['fruitWeight (g)'] / (vol_proxy / 1000.0)
    
    # Additional derived traits
    df_pheno['Pulp-to-Seed Ratio'] = df_pheno['Pulp'] / df_pheno['seedWeight (g)'].replace(0, np.nan)
    df_pheno['Pulp-to-Stone Ratio'] = df_pheno['Pulp'] / df_pheno['stoneWeight (g)'].replace(0, np.nan)
    non_edible_w = (df_pheno['stoneWeight (g)'] + df_pheno['seedWeight (g)']).replace(0, np.nan)
    df_pheno['Sweetness Efficiency Index'] = df_pheno['brix'] * (df_pheno['Pulp'] / non_edible_w)
    
    measured_traits = [
        'fruitWeight (g)', 'fruitLength (mm)', 'fruitWidth (mm)', 'fruitThickness (mm)',
        'stoneWeight (g)', 'stoneLength (mm)', 'stoneWidth (mm)', 'stoneThickness (mm)',
        'seedWeight (g)', 'seedLength (mm)', 'seedWidth (mm)', 'seedThickness (mm)',
        'brix', 'Pulp'
    ]
    derived_traits = ['Fruit Shape Index', 'Pulp Ratio', 'Seed Ratio', 'Stone Ratio',
                      'Edible Portion', 'Brix Yield Index', 'Fruit Density Index',
                      'Pulp-to-Seed Ratio', 'Pulp-to-Stone Ratio', 'Sweetness Efficiency Index']
    all_traits = measured_traits + derived_traits

    for t in all_traits:
        df_pheno[t] = pd.to_numeric(df_pheno[t], errors='coerce')

    print("Loading SNP data...")
    df_head = pd.read_csv(CSV_FILE, header=None, skiprows=2, nrows=1, keep_default_na=False)
    snp_cultivars = [str(x).strip() for x in list(df_head.values[0])[1:]]
    
    df_snps = pd.read_csv(CSV_FILE, header=None, skiprows=3, keep_default_na=False)
    rs_ids = df_snps[0].values
    genotype_data_raw = df_snps.iloc[:, 1:].values
    
    print(f"Loaded {len(rs_ids)} SNPs for {len(snp_cultivars)} SNP accessions.")
    
    print("Parsing SNP coordinates...")
    chroms = []
    positions = []
    ref_alleles = []
    alt_alleles = []
    for rs in rs_ids:
        parts = rs.split('_')
        refseq = parts[0] + '_' + parts[1]
        chrom = REFSEQ_TO_CHR.get(refseq, refseq)
        pos = int(parts[2])
        ref = parts[3]
        alt = parts[4]
        chroms.append(chrom)
        positions.append(pos)
        ref_alleles.append(ref)
        alt_alleles.append(alt)
    chroms = np.array(chroms)
    positions = np.array(positions)
    
    print("Aligning phenotype accessions with SNP columns...")
    clean_snp_cols = [clean_id(x) for x in snp_cultivars]
    
    aligned_indices_pheno = []
    aligned_indices_snp = []
    
    for i_pheno, row in df_pheno.iterrows():
        p_clean = row['clean_id']
        if p_clean == 'optommyatkins':
            matches = [i for i, s in enumerate(snp_cultivars) if s == 'O.PTommyAtkins']
            if matches:
                aligned_indices_pheno.append(i_pheno)
                aligned_indices_snp.append(matches[0])
                continue
        elif p_clean == 'tommyatkins':
            matches = [i for i, s in enumerate(snp_cultivars) if s == 'TOMMYATKINS']
            if matches:
                aligned_indices_pheno.append(i_pheno)
                aligned_indices_snp.append(matches[0])
                continue

        try:
            snp_idx = clean_snp_cols.index(p_clean)
            aligned_indices_pheno.append(i_pheno)
            aligned_indices_snp.append(snp_idx)
        except ValueError:
            pass
            
    print(f"Successfully aligned {len(aligned_indices_pheno)} accessions.")
    
    df_pheno_aligned = df_pheno.iloc[aligned_indices_pheno].copy().reset_index(drop=True)
    genotype_data_aligned = genotype_data_raw[:, aligned_indices_snp]
    aligned_cultivar_names = [snp_cultivars[i] for i in aligned_indices_snp]
    
    print("Converting genotypes to numeric dosages...")
    X = np.zeros(genotype_data_aligned.shape, dtype=np.float32)
    
    IUPAC_HET = {
        ('A', 'G'): 'R', ('G', 'A'): 'R',
        ('C', 'T'): 'Y', ('T', 'C'): 'Y',
        ('G', 'C'): 'S', ('C', 'G'): 'S',
        ('A', 'T'): 'W', ('T', 'A'): 'W',
        ('G', 'T'): 'K', ('T', 'G'): 'K',
        ('A', 'C'): 'M', ('C', 'A'): 'M'
    }
    
    for i in range(len(rs_ids)):
        row = genotype_data_aligned[i]
        ref = ref_alleles[i]
        alt = alt_alleles[i]
        het = IUPAC_HET.get((ref, alt), 'N')
        
        dosage = np.full(row.shape, np.nan, dtype=np.float32)
        dosage[row == ref] = 0.0
        dosage[row == alt] = 2.0
        dosage[row == het] = 1.0
        
        valid_mask = ~np.isnan(dosage)
        if np.any(valid_mask):
            mean_val = np.mean(dosage[valid_mask])
            dosage[np.isnan(dosage)] = mean_val
        else:
            dosage[np.isnan(dosage)] = 0.0
        X[i] = dosage

    allele_freqs = np.mean(X, axis=1) / 2.0
    
    print("Running Principal Component Analysis (PCA)...")
    pca = PCA(n_components=3)
    X_pca = pca.fit_transform(X.T)
    
    pca_results = []
    for idx in range(len(aligned_cultivar_names)):
        p_row = df_pheno_aligned.iloc[idx]
        pca_results.append({
            'accession': str(p_row['Accession ID']),
            'pc1': float(X_pca[idx, 0]),
            'pc2': float(X_pca[idx, 1]),
            'pc3': float(X_pca[idx, 2]),
            'subpopulation': str(p_row['Subpopulation']),
            'brix': float(p_row['brix']) if not pd.isna(p_row['brix']) else None,
            'fruitWeight': float(p_row['fruitWeight (g)']) if not pd.isna(p_row['fruitWeight (g)']) else None
        })
        
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    with open(PCA_FILE, 'w', encoding='utf-8') as f:
        json.dump(pca_results, f, indent=2)
    print(f"PCA output written to {PCA_FILE}")

    print("Running real GWAS association testing...")
    all_gwas_results = []
    
    for trait in all_traits:
        t_start = time.time()
        y = df_pheno_aligned[trait].values
        
        valid_pheno_mask = ~np.isnan(y)
        if not np.any(valid_pheno_mask):
            continue
            
        y_valid = y[valid_pheno_mask]
        X_valid = X[:, valid_pheno_mask]
        valid_cultivars = df_pheno_aligned['Accession ID'].values[valid_pheno_mask]
        
        # 1. Pearson Correlation (Fast GWAS)
        betas_pearson, p_values_pearson = fast_gwas(X_valid, y_valid)
        
        # 2. PCA-Corrected GWAS (OLS with PC1-PC3)
        n_samples = len(y_valid)
        pca = PCA(n_components=3)
        PCs = pca.fit_transform(X_valid.T) # shape: (n_samples, 3)
        C = np.column_stack((np.ones(n_samples), PCs)) # shape: (n_samples, 4)
        
        # FWL projection to project out C from y and X in parallel
        theta_y = np.linalg.solve(np.dot(C.T, C), np.dot(C.T, y_valid))
        y_star = y_valid - np.dot(C, theta_y)
        
        theta_X = np.linalg.solve(np.dot(C.T, C), np.dot(C.T, X_valid.T)) # shape (4, M)
        X_star = X_valid - np.dot(C, theta_X).T # shape (M, N)
        
        cov_Xy = np.dot(X_star, y_star)
        var_X = np.sum(X_star**2, axis=1)
        var_X[var_X == 0] = 1.0
        
        betas_pca = cov_Xy / var_X
        rss = np.sum(y_star**2) - (betas_pca**2) * var_X
        rss = np.clip(rss, 1e-10, None)
        
        df = n_samples - 5
        sigma_sq = rss / df
        se_betas = np.sqrt(sigma_sq / var_X)
        
        t_stat = betas_pca / se_betas
        p_values_pca = np.array([math.erfc(abs(t) / math.sqrt(2)) for t in t_stat])
        
        # Compute ranks globally
        sort_idx_pearson = np.argsort(p_values_pearson)
        sort_idx_pca = np.argsort(p_values_pca)
        
        ranks_pearson = np.zeros_like(p_values_pearson, dtype=int)
        ranks_pearson[sort_idx_pearson] = np.arange(1, len(p_values_pearson) + 1)
        
        ranks_pca = np.zeros_like(p_values_pca, dtype=int)
        ranks_pca[sort_idx_pca] = np.arange(1, len(p_values_pca) + 1)
        
        # Keep union of top 500 from both models + 200 background SNPs
        top_pearson_set = set(sort_idx_pearson[:500])
        top_pca_set = set(sort_idx_pca[:500])
        keep_set = top_pearson_set.union(top_pca_set)
        
        # Add 200 random background SNPs
        rng = np.random.RandomState(42)
        random_bg = rng.choice(len(rs_ids), min(200, len(rs_ids)), replace=False)
        keep_set.update(random_bg)
        keep_idx = list(keep_set)
        
        print(f"  Trait: {trait:25s} | Pearson top p-val: {p_values_pearson[sort_idx_pearson[0]]:.2e} | PCA top p-val: {p_values_pca[sort_idx_pca[0]]:.2e} | Done in {time.time()-t_start:.2f}s")
        for i in keep_idx:
            chr_name = chroms[i]
            pos = int(positions[i])
            
            is_top_snp = (i in top_pearson_set) or (i in top_pca_set) or (p_values_pca[i] < 5e-8) or (p_values_pearson[i] < 5e-8)
            
            refseq = CHR_TO_REFSEQ.get(chr_name, chr_name)
            gene_id, context, dist, biotype = annotate_snp_context(gff_genes, gff_exons, gff_cds, refseq, pos)
            product = gff_products.get(gene_id, 'Unknown')
            
            if product and product != 'Unknown' and is_informative_product(product, gene_id):
                gene_name = f"{gene_id} ({product})"
            else:
                gene_name = gene_id
                
            cultivar_afs = None
            cultivar_phenotypes = None
            
            if (i in sort_idx_pearson[:5]) or (i in sort_idx_pca[:5]):
                cultivar_afs = {str(valid_cultivars[c]): float(X_valid[i, c] / 2.0) for c in range(len(valid_cultivars))}
                y_min, y_max = np.min(y_valid), np.max(y_valid)
                if y_max > y_min:
                    y_scaled = 10.0 + 85.0 * (y_valid - y_min) / (y_max - y_min)
                else:
                    y_scaled = np.full(y_valid.shape, 50.0)
                cultivar_phenotypes = {str(valid_cultivars[c]): float(y_scaled[c]) for c in range(len(valid_cultivars))}
            
            res = {
                'id': str(rs_ids[i]),
                'chromosome': str(chr_name),
                'position': pos,
                
                # Standard fields (default to PCA-corrected for compatibility)
                'p_value': float(p_values_pca[i]),
                'neg_log10_p': float(-np.log10(np.clip(p_values_pca[i], 1e-300, 1.0))),
                'beta': float(betas_pca[i]),
                
                # Pearson model fields
                'p_value_pearson': float(p_values_pearson[i]),
                'beta_pearson': float(betas_pearson[i]),
                'rank_pearson': int(ranks_pearson[i]),
                
                # PCA model fields
                'p_value_pca': float(p_values_pca[i]),
                'beta_pca': float(betas_pca[i]),
                'rank_pca': int(ranks_pca[i]),
                
                # Rank change
                'rank_change': int(ranks_pearson[i] - ranks_pca[i]),
                
                'allele_freq': float(allele_freqs[i]),
                'nearest_gene': gene_name,
                'gene_id': gene_id,
                'gene_product': product,
                'functional_context': context,
                'distance': dist,
                'gene_biotype': biotype,
                'trait': trait,
                'is_peak': bool(is_top_snp),
                'cultivar_afs': cultivar_afs,
                'cultivar_phenotypes': cultivar_phenotypes
            }
            all_gwas_results.append(res)

    print(f"Total GWAS records to save: {len(all_gwas_results)}")
    
    all_gwas_results.sort(key=lambda x: (x['chromosome'], x['position']))
    with open(GWAS_FILE, 'w', encoding='utf-8') as f:
        json.dump(all_gwas_results, f, indent=2)
        
    print(f"GWAS results saved to {GWAS_FILE} in {time.time()-t0:.2f}s total.")

if __name__ == '__main__':
    main()
