import os
import sys
import json
import sqlite3
import numpy as np
import pandas as pd
import time
import random

sys.path.insert(0, r'C:\mangoproject')
from prediction_engine import engine

BASE_DIR = r'C:\mangoproject'
DB_FILE = os.path.join(BASE_DIR, 'synthetic_cultivars.db')
EXCEL_FILE = os.path.join(BASE_DIR, 'synthetic_cultivars.xlsx')

def init_db():
    print("Initializing synthetic_cultivars.db...")
    conn = sqlite3.connect(DB_FILE)
    c = conn.cursor()
    
    # Table for cultivars metadata and genotype
    c.execute('''
        CREATE TABLE IF NOT EXISTS cultivars (
            id TEXT PRIMARY KEY,
            name TEXT,
            parent_a TEXT,
            parent_b TEXT,
            type TEXT, -- 'hybrid', 'virtual', 'uploaded'
            genotype TEXT, -- compact string of 0, 1, 2
            mutation_rate REAL,
            created_at TEXT
        )
    ''')
    
    # Table for cached predictions
    c.execute('''
        CREATE TABLE IF NOT EXISTS predictions (
            cultivar_id TEXT,
            trait TEXT,
            model_type TEXT, -- 'rrBLUP', 'Random Forest', 'XGBoost'
            predicted REAL,
            ci_lower REAL,
            ci_upper REAL,
            PRIMARY KEY (cultivar_id, trait, model_type)
        )
    ''')
    
    conn.commit()
    conn.close()

def simulate_gamete(genotype_vector, rng):
    # For each locus, a parent transmits one of its alleles
    # 0 (Ref/Ref) transmits 0 (Ref allele)
    # 2 (Alt/Alt) transmits 1 (Alt allele)
    # 1 (Ref/Alt) transmits 0 or 1 with 50% probability each
    gamete = np.zeros_like(genotype_vector)
    
    # Alt/Alt (2.0) transmits 1.0
    gamete[genotype_vector == 2.0] = 1.0
    
    # Heterozygous (1.0) transmits 0.0 or 1.0 with 50% probability
    het_mask = (genotype_vector == 1.0)
    n_hets = np.sum(het_mask)
    if n_hets > 0:
        gamete[het_mask] = rng.choice([0.0, 1.0], size=n_hets)
        
    return gamete

def make_hybrid_cross(genotype_a, genotype_b, rng):
    gamete_a = simulate_gamete(genotype_a, rng)
    gamete_b = simulate_gamete(genotype_b, rng)
    offspring = gamete_a + gamete_b
    return offspring

def add_mutation(genotype_vector, mutation_rate, rng):
    mutated = genotype_vector.copy()
    n_markers = len(mutated)
    n_mutations = int(mutation_rate * n_markers)
    if n_mutations > 0:
        mut_idx = rng.choice(n_markers, size=n_mutations, replace=False)
        for idx in mut_idx:
            mutated[idx] = rng.choice([0.0, 1.0, 2.0])
    return mutated

def main():
    t_start = time.time()
    init_db()
    
    # Ensure prediction engine data is loaded
    print("Loading BreedingEngine training data...")
    engine.load_data()
    n_snps = engine.X.shape[1]
    n_ref_cultivars = len(engine.cultivars)
    print(f"Genomic panel loaded: {n_ref_cultivars} accessions, {n_snps} SNPs.")
    
    # Map from cultivars names to index in panel
    cv_indices = {c.lower(): idx for idx, c in enumerate(engine.cultivars)}
    
    # 1. Define Benchmark Hybrids to generate
    benchmarks = [
        ('Alphonso', 'Keitt', 'TEST_HYBRID_Alphonso_Keitt_F1', 0.0),
        ('Tommy Atkins', 'Kent', 'TEST_HYBRID_TommyAtkins_Kent_F1', 0.0),
        ('Langra', 'Haden', 'TEST_HYBRID_Langra_Haden_F1', 0.0),
        ('Chaunsa', 'Keitt', 'TEST_HYBRID_Chaunsa_Keitt_F1', 0.0),
    ]
    
    generated_cultivars = []
    rng = np.random.RandomState(42)
    
    print("\nGenerating F1 hybrid benchmarks...")
    for parent_a_name, parent_b_name, child_id, mut_rate in benchmarks:
        idx_a = cv_indices.get(parent_a_name.lower().replace(' ', '').replace('_', ''))
        idx_b = cv_indices.get(parent_b_name.lower().replace(' ', '').replace('_', ''))
        
        # Fallback if names are slightly different (e.g. tommy_atkins vs tommyatkins)
        if idx_a is None:
            # try fuzzy matching
            for cv_name, idx in cv_indices.items():
                if parent_a_name.lower().replace(' ', '') in cv_name or cv_name in parent_a_name.lower():
                    idx_a = idx
                    break
        if idx_b is None:
            for cv_name, idx in cv_indices.items():
                if parent_b_name.lower().replace(' ', '') in cv_name or cv_name in parent_b_name.lower():
                    idx_b = idx
                    break
                    
        if idx_a is None or idx_b is None:
            print(f"Error: parents {parent_a_name} ({idx_a}) or {parent_b_name} ({idx_b}) not found.")
            continue
            
        genotype_a = engine.X[idx_a]
        genotype_b = engine.X[idx_b]
        
        offspring_gt = make_hybrid_cross(genotype_a, genotype_b, rng)
        if mut_rate > 0:
            offspring_gt = add_mutation(offspring_gt, mut_rate, rng)
            
        generated_cultivars.append({
            'id': child_id,
            'name': child_id.replace('TEST_HYBRID_', '').replace('_F1', '').replace('_', ' × '),
            'parent_a': engine.cultivars[idx_a],
            'parent_b': engine.cultivars[idx_b],
            'type': 'hybrid',
            'genotype': offspring_gt,
            'mutation_rate': mut_rate
        })
        print(f"  Generated hybrid {child_id} (Parent A: {engine.cultivars[idx_a]}, Parent B: {engine.cultivars[idx_b]})")

    # 2. Generate Virtual Cultivars Library (50 virtual cultivars for speed and size, expandable)
    # The prompt says "Generate hundreds of synthetic cultivars...". Let's do 100 virtual cultivars.
    n_virtual = 100
    print(f"\nGenerating {n_virtual} virtual cultivars from panel...")
    for i in range(1, n_virtual + 1):
        idx_a = rng.randint(0, n_ref_cultivars)
        idx_b = rng.randint(0, n_ref_cultivars)
        while idx_a == idx_b:
            idx_b = rng.randint(0, n_ref_cultivars)
            
        genotype_a = engine.X[idx_a]
        genotype_b = engine.X[idx_b]
        
        offspring_gt = make_hybrid_cross(genotype_a, genotype_b, rng)
        mut_rate = 0.005 # 0.5% mutation rate
        offspring_gt = add_mutation(offspring_gt, mut_rate, rng)
        
        cv_id = f"TEST_VIRTUAL_cultivar_{i:03d}"
        generated_cultivars.append({
            'id': cv_id,
            'name': f"Virtual Cultivar {i:03d}",
            'parent_a': engine.cultivars[idx_a],
            'parent_b': engine.cultivars[idx_b],
            'type': 'virtual',
            'genotype': offspring_gt,
            'mutation_rate': mut_rate
        })
        if i % 20 == 0:
            print(f"  Generated {i}/{n_virtual} virtual cultivars...")

    # Save to SQLite Database and Predict
    print("\nRunning multi-model phenotypic predictions on all generated cultivars...")
    conn = sqlite3.connect(DB_FILE)
    c = conn.cursor()
    
    metadata_list = []
    predictions_excel_list = []
    
    all_traits = [
        'fruitWeight (g)', 'fruitLength (mm)', 'fruitWidth (mm)', 'fruitThickness (mm)',
        'brix', 'Pulp', 'seedWeight (g)', 'stoneWeight (g)',
        'seedLength (mm)', 'seedWidth (mm)', 'seedThickness (mm)',
        'stoneLength (mm)', 'stoneWidth (mm)', 'stoneThickness (mm)',
        'Fruit Shape Index', 'Pulp Ratio', 'Seed Ratio', 'Stone Ratio',
        'Edible Portion', 'Brix Yield Index', 'Fruit Density Index',
        'Pulp-to-Seed Ratio', 'Pulp-to-Stone Ratio', 'Sweetness Efficiency Index'
    ]
    
    gen_date = time.strftime("%Y-%m-%d %H:%M:%S")
    
    for idx, cv in enumerate(generated_cultivars):
        cv_id = cv['id']
        name = cv['name']
        p_a = cv['parent_a']
        p_b = cv['parent_b']
        cv_type = cv['type']
        gt = cv['genotype']
        mut_rate = cv['mutation_rate']
        
        # 1. Encode genotype to compact string
        gt_str = "".join([str(int(round(x))) for x in gt])
        
        # Write cultivar to SQLite
        c.execute('''
            INSERT OR REPLACE INTO cultivars (id, name, parent_a, parent_b, type, genotype, mutation_rate, created_at)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?)
        ''', (cv_id, name, p_a, p_b, cv_type, gt_str, mut_rate, gen_date))
        
        # 2. Run predictions under rrBLUP, RF, and XGBoost
        # We can optimize: train each model once per trait and predict all cultivars!
        # But prediction_engine predict_new_cultivar is already optimized and trains internally.
        # Let's call predict_new_cultivar for the cultivar
        pred_rrblup = engine.predict_new_cultivar(gt.tolist(), 'rrBLUP')
        pred_rf = engine.predict_new_cultivar(gt.tolist(), 'Random Forest')
        pred_xgboost = engine.predict_new_cultivar(gt.tolist(), 'XGBoost')
        
        metadata_list.append({
            'Synthetic Cultivar ID': cv_id,
            'Cultivar Name': name,
            'Parent A': p_a,
            'Parent B': p_b,
            'Hybrid Type': cv_type,
            'Mutation Rate': f"{mut_rate*100:.1f}%",
            'Generation Date': gen_date
        })
        
        for trait in all_traits:
            res_rr = pred_rrblup['predictions'].get(trait, None)
            res_rf = pred_rf['predictions'].get(trait, None)
            res_xgb = pred_xgboost['predictions'].get(trait, None)
            
            if res_rr and res_rf and res_xgb:
                pred_rr_val = res_rr['predicted']
                pred_rf_val = res_rf['predicted']
                pred_xgb_val = res_xgb['predicted']
                
                mean_pred = (pred_rr_val + pred_rf_val + pred_xgb_val) / 3.0
                sd_pred = np.std([pred_rr_val, pred_rf_val, pred_xgb_val])
                
                # Write predictions to SQLite (for fast dashboard lookups)
                c.execute('''
                    INSERT OR REPLACE INTO predictions (cultivar_id, trait, model_type, predicted, ci_lower, ci_upper)
                    VALUES (?, ?, ?, ?, ?, ?)
                ''', (cv_id, trait, 'rrBLUP', pred_rr_val, res_rr['ci_lower'], res_rr['ci_upper']))
                
                c.execute('''
                    INSERT OR REPLACE INTO predictions (cultivar_id, trait, model_type, predicted, ci_lower, ci_upper)
                    VALUES (?, ?, ?, ?, ?, ?)
                ''', (cv_id, trait, 'Random Forest', pred_rf_val, res_rf['ci_lower'], res_rf['ci_upper']))
                
                c.execute('''
                    INSERT OR REPLACE INTO predictions (cultivar_id, trait, model_type, predicted, ci_lower, ci_upper)
                    VALUES (?, ?, ?, ?, ?, ?)
                ''', (cv_id, trait, 'XGBoost', pred_xgb_val, res_xgb['ci_lower'], res_xgb['ci_upper']))
                
                # Append to Excel records
                predictions_excel_list.append({
                    'Synthetic Cultivar ID': cv_id,
                    'Trait': trait,
                    'rrBLUP_Prediction': pred_rr_val,
                    'rrBLUP_CI_Lower': res_rr['ci_lower'],
                    'rrBLUP_CI_Upper': res_rr['ci_upper'],
                    'RF_Prediction': pred_rf_val,
                    'XGBoost_Prediction': pred_xgb_val,
                    'Model_Mean': round(mean_pred, 3),
                    'Model_SD': round(sd_pred, 3)
                })
                
        if (idx + 1) % 20 == 0 or cv_id.startswith('TEST_HYBRID'):
            print(f"  Predicted {idx + 1}/{len(generated_cultivars)} cultivars...")

    conn.commit()
    conn.close()
    print("SQLite database successfully populated and closed.")
    
    # 3. Export to Excel (synthetic_cultivars.xlsx)
    print("\nExporting synthetic cultivars dataset to Excel...")
    df_meta = pd.DataFrame(metadata_list)
    df_pred = pd.DataFrame(predictions_excel_list)
    
    with pd.ExcelWriter(EXCEL_FILE, engine='openpyxl') as writer:
        df_meta.to_excel(writer, sheet_name='Metadata', index=False)
        df_pred.to_excel(writer, sheet_name='Predicted_Phenotypes', index=False)
        
    print(f"Excel dataset successfully written to {EXCEL_FILE}")
    print(f"Process completed in {time.time()-t_start:.2f}s total.")

if __name__ == '__main__':
    main()
