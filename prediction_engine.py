import os
import sys
import json
# Only insert local precompiled libs under Windows with Python 3.14
if sys.platform == 'win32' and sys.version_info[:2] == (3, 14):
    BASE_DIR = os.path.dirname(os.path.abspath(__file__))
    sys.path.insert(0, os.path.join(BASE_DIR, 'libs'))

import numpy as np
import pandas as pd
import re
from sklearn.model_selection import KFold, train_test_split
from sklearn.linear_model import Ridge
from sklearn.ensemble import RandomForestRegressor
import xgboost as xgb
from scipy.stats import pearsonr

CSV_FILE   = os.path.join(BASE_DIR, 'SNP.csv')
PHENO_FILE = os.path.join(BASE_DIR, 'phenotypes.xlsx')

REFSEQ_TO_CHR = {
    'NC_058137.1': 'Chr1',  'NC_058138.1': 'Chr2',  'NC_058139.1': 'Chr3',
    'NC_058140.1': 'Chr4',  'NC_058141.1': 'Chr5',  'NC_058142.1': 'Chr6',
    'NC_058143.1': 'Chr7',  'NC_058144.1': 'Chr8',  'NC_058145.1': 'Chr9',
    'NC_058146.1': 'Chr10', 'NC_058147.1': 'Chr11', 'NC_058148.1': 'Chr12',
    'NC_058149.1': 'Chr13', 'NC_058150.1': 'Chr14', 'NC_058151.1': 'Chr15',
    'NC_058152.1': 'Chr16', 'NC_058153.1': 'Chr17', 'NC_058154.1': 'Chr18',
    'NC_058155.1': 'Chr19', 'NC_058156.1': 'Chr20',
}

def clean_id(s):
    return re.sub(r'[^a-zA-Z0-9]', '', str(s)).lower()

def get_mapped_clean(name):
    # Strip parentheses and their contents, e.g. "Alphonso (Hapus)" -> "Alphonso"
    name_clean = re.sub(r'\(.*?\)', '', str(name)).strip()
    c = clean_id(name_clean)
    aliases = {
        'sindhri': 'shindiri',
        'alphonso': 'whitealfonso',
        'tommy_atkins': 'tommyatkins',
        'tommyatkins': 'tommyatkins',
    }
    return aliases.get(c, c)

class BreedingEngine:
    def __init__(self):
        self.X = None
        self.rs_ids = None
        self.df_pheno = None
        self.cultivars = None
        self.is_loaded = False
        
    def load_data(self):
        if self.is_loaded:
            return
            
        # Load phenotypes
        df_pheno = pd.read_excel(PHENO_FILE, sheet_name=' Table S1', header=2, keep_default_na=False).iloc[:161]
        df_pheno['clean_id'] = df_pheno['Accession ID'].apply(clean_id)
        
        # Calculate derived traits
        for col in [
            'fruitWeight (g)', 'fruitLength (mm)', 'fruitWidth (mm)', 'fruitThickness (mm)',
            'brix', 'Pulp', 'seedWeight (g)', 'seedLength (mm)', 'seedWidth (mm)', 'seedThickness (mm)',
            'stoneWeight (g)', 'stoneLength (mm)', 'stoneWidth (mm)', 'stoneThickness (mm)'
        ]:
            df_pheno[col] = pd.to_numeric(df_pheno[col], errors='coerce')

        df_pheno['Fruit Shape Index'] = df_pheno['fruitLength (mm)'] / df_pheno['fruitWidth (mm)']
        df_pheno['Pulp Ratio'] = df_pheno['Pulp'] / df_pheno['fruitWeight (g)']
        df_pheno['Seed Ratio'] = df_pheno['seedWeight (g)'] / df_pheno['fruitWeight (g)']
        df_pheno['Stone Ratio'] = df_pheno['stoneWeight (g)'] / df_pheno['fruitWeight (g)']
        # New derived traits
        df_pheno['Edible Portion'] = (df_pheno['Pulp'] / df_pheno['fruitWeight (g)']) * 100
        df_pheno['Brix Yield Index'] = df_pheno['brix'] * (df_pheno['Pulp'] / df_pheno['fruitWeight (g)'])
        vol_proxy = df_pheno['fruitLength (mm)'] * df_pheno['fruitWidth (mm)'] * df_pheno['fruitThickness (mm)']
        df_pheno['Fruit Density Index'] = df_pheno['fruitWeight (g)'] / (vol_proxy / 1000.0)
        
        # Additional derived traits
        df_pheno['Pulp-to-Seed Ratio'] = df_pheno['Pulp'] / df_pheno['seedWeight (g)'].replace(0, np.nan)
        df_pheno['Pulp-to-Stone Ratio'] = df_pheno['Pulp'] / df_pheno['stoneWeight (g)'].replace(0, np.nan)
        non_edible_w = (df_pheno['stoneWeight (g)'] + df_pheno['seedWeight (g)']).replace(0, np.nan)
        df_pheno['Sweetness Efficiency Index'] = df_pheno['brix'] * (df_pheno['Pulp'] / non_edible_w)
        
        # Load SNPs
        df_head = pd.read_csv(CSV_FILE, header=None, skiprows=2, nrows=1, keep_default_na=False)
        snp_cultivars = [str(x).strip() for x in list(df_head.values[0])[1:]]
        
        df_snps = pd.read_csv(CSV_FILE, header=None, skiprows=3, keep_default_na=False)
        self.rs_ids = df_snps[0].values.tolist()
        genotype_data_raw = df_snps.iloc[:, 1:].values
        
        # Align
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
                
        self.df_pheno = df_pheno.iloc[aligned_indices_pheno].copy().reset_index(drop=True)
        genotype_data_aligned = genotype_data_raw[:, aligned_indices_snp]
        self.cultivars = [snp_cultivars[i] for i in aligned_indices_snp]
        
        # Genotype encoding (0, 1, 2)
        print("[PredictionEngine] Encoding genotypes...")
        X = np.zeros(genotype_data_aligned.shape, dtype=np.float32)
        
        IUPAC_HET = {
            ('A', 'G'): 'R', ('G', 'A'): 'R',
            ('C', 'T'): 'Y', ('T', 'C'): 'Y',
            ('G', 'C'): 'S', ('C', 'G'): 'S',
            ('A', 'T'): 'W', ('T', 'A'): 'W',
            ('G', 'T'): 'K', ('T', 'G'): 'K',
            ('A', 'C'): 'M', ('C', 'A'): 'M'
        }
        
        for i in range(len(self.rs_ids)):
            row = genotype_data_aligned[i]
            parts = self.rs_ids[i].split('_')
            ref = parts[3]
            alt = parts[4]
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
            
        self.X = X.T # Shape: (n_samples, n_snps)
        
        # Fit PCA on reference population to project uploaded/new genotypes using same fitted model
        from sklearn.decomposition import PCA
        self.pca_model = PCA(n_components=3, random_state=42)
        self.pca_model.fit(self.X)
        
        self.trait_std_X = {}
        all_traits = [
            'fruitWeight (g)', 'fruitLength (mm)', 'fruitWidth (mm)', 'fruitThickness (mm)',
            'brix', 'Pulp', 'seedWeight (g)', 'stoneWeight (g)',
            'seedLength (mm)', 'seedWidth (mm)', 'seedThickness (mm)',
            'stoneLength (mm)', 'stoneWidth (mm)', 'stoneThickness (mm)',
            'Fruit Shape Index', 'Pulp Ratio', 'Seed Ratio', 'Stone Ratio',
            'Edible Portion', 'Brix Yield Index', 'Fruit Density Index',
            'Pulp-to-Seed Ratio', 'Pulp-to-Stone Ratio', 'Sweetness Efficiency Index'
        ]
        print("[PredictionEngine] Precomputing SNP variance per trait...")
        for trait in all_traits:
            y = self.df_pheno[trait].values
            valid_mask = ~np.isnan(y)
            X_clean = self.X[valid_mask]
            sum_X = np.sum(X_clean, axis=0)
            sum_X2 = np.sum(X_clean**2, axis=0)
            var_X = sum_X2 - (sum_X**2) / np.sum(valid_mask)
            var_X[var_X == 0] = 1.0
            self.trait_std_X[trait] = np.sqrt(var_X).astype(np.float32)
        self.is_loaded = True

    def train_predict(self, trait, model_type):
        self.load_data()
        
        # Prepare targets
        y = self.df_pheno[trait].values
        valid_mask = ~np.isnan(y)
        
        y_clean = y[valid_mask]
        X_clean = self.X[valid_mask]
        cultivars_clean = np.array(self.cultivars)[valid_mask]
        
        if len(y_clean) < 10:
            return {'error': 'Insufficient data for prediction.'}
            
        # ── Feature Selection (Top 500 SNPs based on absolute correlation with trait) ──
        y_centered = y_clean - np.mean(y_clean)
        X_centered = X_clean - np.mean(X_clean, axis=0, keepdims=True)
        
        # Pearson correlation
        cov = np.dot(y_centered, X_centered)
        var_X = np.sum(X_centered**2, axis=0)
        var_X[var_X == 0] = 1.0
        var_y = np.sum(y_centered**2)
        if var_y == 0:
            var_y = 1.0
            
        corr = cov / np.sqrt(var_X * var_y)
        abs_corr = np.abs(corr)
        
        top_k = min(500, X_clean.shape[1])
        selected_features = np.argsort(abs_corr)[::-1][:top_k]
        X_selected = X_clean[:, selected_features]

        # Model initialisation
        if model_type == 'rrBLUP':
            model = Ridge(alpha=100.0)
        elif model_type == 'Random Forest':
            model = RandomForestRegressor(n_estimators=50, random_state=42, n_jobs=-1)
        else: # XGBoost
            model = xgb.XGBRegressor(n_estimators=50, max_depth=3, learning_rate=0.1, random_state=42, n_jobs=-1)
            
        # 5-fold cross validation for predictions & metrics
        kf = KFold(n_splits=5, shuffle=True, random_state=42)
        cv_predictions = np.zeros_like(y_clean)
        
        for train_idx, test_idx in kf.split(X_selected):
            model.fit(X_selected[train_idx], y_clean[train_idx])
            cv_predictions[test_idx] = model.predict(X_selected[test_idx])
            
        # Compute overall stats
        r2 = float(np.corrcoef(y_clean, cv_predictions)[0, 1]**2) if np.std(cv_predictions) > 0 else 0.0
        rmse = float(np.sqrt(np.mean((y_clean - cv_predictions)**2)))
        mae = float(np.mean(np.abs(y_clean - cv_predictions)))
        try:
            pearson_r, _ = pearsonr(y_clean, cv_predictions)
            pearson_r = float(pearson_r)
        except Exception:
            pearson_r = 0.0
            
        # Train on selected features to get importances
        model.fit(X_selected, y_clean)
        if model_type == 'rrBLUP':
            importances = np.abs(model.coef_)
        else:
            importances = model.feature_importances_
            
        # Top 10 SNPs (backward compatibility)
        top_indices = np.argsort(importances)[::-1][:10]
        top_snps = []
        for idx in top_indices:
            orig_idx = selected_features[idx]
            snp_id = self.rs_ids[orig_idx]
            parts = snp_id.split('_')
            if len(parts) >= 5:
                refseq = parts[0] + '_' + parts[1]
                chrom = REFSEQ_TO_CHR.get(refseq, refseq)
                pos = parts[2]
                ref = parts[3]
                alt = parts[4]
                friendly_id = f"{chrom}_{pos}_{ref}_{alt}"
            else:
                friendly_id = snp_id
            top_snps.append({
                'snp_id': friendly_id,
                'importance': float(importances[idx])
            })
            
        # Top 50 SNPs with detailed metadata
        sorted_imp_indices = np.argsort(importances)[::-1]
        top_50_snps = []
        for rank_idx, idx in enumerate(sorted_imp_indices[:50]):
            orig_idx = selected_features[idx]
            snp_id = self.rs_ids[orig_idx]
            parts = snp_id.split('_')
            if len(parts) >= 5:
                refseq = parts[0] + '_' + parts[1]
                chrom = REFSEQ_TO_CHR.get(refseq, refseq)
                pos = parts[2]
                ref = parts[3]
                alt = parts[4]
                friendly_id = f"{chrom}_{pos}_{ref}_{alt}"
            else:
                friendly_id = snp_id
                chrom = "Unknown"
                pos = "0"
            top_50_snps.append({
                'snp_id': friendly_id,
                'chromosome': chrom,
                'position': int(pos) if pos.isdigit() else pos,
                'correlation': float(corr[orig_idx]),
                'importance': float(importances[idx]),
                'rank': rank_idx + 1
            })

        # Helper to compute cross-validation metrics for feature count comparison
        def run_cv_for_k(k_val):
            top_k_feat = min(k_val, X_clean.shape[1])
            sel_feat = np.argsort(abs_corr)[::-1][:top_k_feat]
            X_sel = X_clean[:, sel_feat]
            
            if model_type == 'rrBLUP':
                model_comp = Ridge(alpha=100.0)
            elif model_type == 'Random Forest':
                model_comp = RandomForestRegressor(n_estimators=30, random_state=42, n_jobs=-1)
            else: # XGBoost
                model_comp = xgb.XGBRegressor(n_estimators=30, max_depth=3, learning_rate=0.1, random_state=42, n_jobs=-1)
                
            kf_comp = KFold(n_splits=5, shuffle=True, random_state=42)
            cv_preds = np.zeros_like(y_clean)
            
            for train_idx, test_idx in kf_comp.split(X_sel):
                if model_type == 'rrBLUP':
                    # Extremely fast Cholesky solver
                    n_features = X_sel.shape[1]
                    mean_y = np.mean(y_clean[train_idx])
                    mean_X = np.mean(X_sel[train_idx], axis=0)
                    y_cent = y_clean[train_idx] - mean_y
                    X_cent = X_sel[train_idx] - mean_X
                    A = np.dot(X_cent.T, X_cent) + 100.0 * np.eye(n_features)
                    b = np.dot(X_cent.T, y_cent)
                    beta = np.linalg.solve(A, b)
                    cv_preds[test_idx] = mean_y + np.dot(X_sel[test_idx] - mean_X, beta)
                else:
                    model_comp.fit(X_sel[train_idx], y_clean[train_idx])
                    cv_preds[test_idx] = model_comp.predict(X_sel[test_idx])
                    
            r2_comp = float(np.corrcoef(y_clean, cv_preds)[0, 1]**2) if np.std(cv_preds) > 0 else 0.0
            rmse_comp = float(np.sqrt(np.mean((y_clean - cv_preds)**2)))
            mae_comp = float(np.mean(np.abs(y_clean - cv_preds)))
            try:
                r_comp, _ = pearsonr(y_clean, cv_preds)
                r_comp = float(r_comp)
            except Exception:
                r_comp = 0.0
                
            return {
                'r2': round(r2_comp, 4),
                'rmse': round(rmse_comp, 4),
                'mae': round(mae_comp, 4),
                'correlation': round(r_comp, 4)
            }
            
        # Cultivar predictions table
        cultivar_data = []
        for idx in range(len(y_clean)):
            cultivar_data.append({
                'cultivar': str(cultivars_clean[idx]),
                'observed': float(y_clean[idx]),
                'predicted': float(cv_predictions[idx])
            })
            
        return {
            'metrics': {
                'r2': round(r2, 4),
                'rmse': round(rmse, 4),
                'mae': round(mae, 4),
                'correlation': round(pearson_r, 4)
            },
            'predictions': cultivar_data,
            'top_features': top_snps,
            'top_50_features': top_50_snps,
            'feature_selection_comparison': {
                '500': run_cv_for_k(500),
                '1000': run_cv_for_k(1000),
                '2000': run_cv_for_k(2000)
            }
        }

    def predict_new_cultivar(self, genotype_vector, model_type='rrBLUP', bootstrap=True):
        """
        Predict all phenotypic traits for a new cultivar given its SNP genotype vector.
        Returns predictions with bootstrap confidence intervals, top contributing SNPs,
        and genetic similarity to known cultivars.
        """
        self.load_data()

        if len(genotype_vector) != self.X.shape[1]:
            return {'error': f'Expected {self.X.shape[1]} SNP values, got {len(genotype_vector)}.'}

        if not hasattr(self, '_prediction_cache'):
            self._prediction_cache = {}
        sig = (
            len(genotype_vector),
            float(genotype_vector[0]),
            float(genotype_vector[len(genotype_vector)//2]),
            float(genotype_vector[-1]),
            model_type,
            bootstrap
        )
        if sig in self._prediction_cache:
            return self._prediction_cache[sig]

        new_X = np.array(genotype_vector, dtype=np.float32).reshape(1, -1)

        # All predictable traits
        all_traits = [
            'fruitWeight (g)', 'fruitLength (mm)', 'fruitWidth (mm)', 'fruitThickness (mm)',
            'brix', 'Pulp', 'seedWeight (g)', 'stoneWeight (g)',
            'seedLength (mm)', 'seedWidth (mm)', 'seedThickness (mm)',
            'stoneLength (mm)', 'stoneWidth (mm)', 'stoneThickness (mm)',
            'Fruit Shape Index', 'Pulp Ratio', 'Seed Ratio', 'Stone Ratio',
            'Edible Portion', 'Brix Yield Index', 'Fruit Density Index',
            'Pulp-to-Seed Ratio', 'Pulp-to-Stone Ratio', 'Sweetness Efficiency Index'
        ]

        predictions = {}
        top_snps_per_trait = {}

        for trait in all_traits:
            y = self.df_pheno[trait].values
            valid_mask = ~np.isnan(y)
            y_clean = y[valid_mask]
            if np.all(valid_mask):
                X_clean = self.X
            else:
                X_clean = self.X[valid_mask]

            if len(y_clean) < 10:
                continue

            cache_key = (trait, model_type)
            if not hasattr(self, '_model_cache'):
                self._model_cache = {}

            if cache_key not in self._model_cache:
                # 1. Feature selection (only run once per model type/trait)
                y_centered = (y_clean - np.mean(y_clean)).astype(np.float32)
                cov = np.dot(y_centered, X_clean)
                var_y = np.sum(y_centered**2)
                if var_y == 0:
                    var_y = 1.0
                corr = cov / (self.trait_std_X[trait] * np.sqrt(var_y))
                abs_corr = np.abs(corr)

                top_k = min(500, X_clean.shape[1])
                selected_features = np.argsort(abs_corr)[::-1][:top_k]
                X_selected = X_clean[:, selected_features]

                # 2. Train main model
                if model_type == 'rrBLUP':
                    n_features = X_selected.shape[1]
                    mean_y = np.mean(y_clean)
                    mean_X = np.mean(X_selected, axis=0)
                    y_centered_fit = y_clean - mean_y
                    X_centered_fit = X_selected - mean_X
                    A = np.dot(X_centered_fit.T, X_centered_fit) + 100.0 * np.eye(n_features)
                    b = np.dot(X_centered_fit.T, y_centered_fit)
                    beta = np.linalg.solve(A, b)
                    importances = np.abs(beta)

                    model_dict = {
                        'mean_y': mean_y,
                        'mean_X': mean_X,
                        'beta': beta,
                        'importances': importances
                    }
                else:
                    if model_type == 'Random Forest':
                        model = RandomForestRegressor(n_estimators=50, random_state=42, n_jobs=-1)
                    else:
                        model = xgb.XGBRegressor(n_estimators=50, max_depth=3, learning_rate=0.1, random_state=42, n_jobs=-1)
                    model.fit(X_selected, y_clean)
                    importances = model.feature_importances_

                    model_dict = {
                        'model': model,
                        'importances': importances
                    }

                self._model_cache[cache_key] = {
                    'selected_features': selected_features,
                    'model_dict': model_dict,
                    'boot_data': [],
                    'population_mean': float(np.mean(y_clean)),
                    'population_std': float(np.std(y_clean))
                }

            cached_data = self._model_cache[cache_key]
            selected_features = cached_data['selected_features']
            new_X_sel = new_X[:, selected_features]
            model_dict = cached_data['model_dict']

            # 3. Train bootstrap models if requested and not already trained
            if bootstrap and not cached_data['boot_data']:
                X_selected = X_clean[:, selected_features]
                boot_data = []
                rng = np.random.RandomState(42)
                for _ in range(20):
                    boot_idx = rng.choice(len(y_clean), size=len(y_clean), replace=True)
                    if model_type == 'rrBLUP':
                        y_boot = y_clean[boot_idx]
                        X_boot = X_selected[boot_idx]
                        mean_y_b = np.mean(y_boot)
                        mean_X_b = np.mean(X_boot, axis=0)
                        y_c_b = y_boot - mean_y_b
                        X_c_b = X_boot - mean_X_b
                        A_b = np.dot(X_c_b.T, X_c_b) + 100.0 * np.eye(X_boot.shape[1])
                        b_b = np.dot(X_c_b.T, y_c_b)
                        beta_b = np.linalg.solve(A_b, b_b)
                        boot_data.append({
                            'mean_y_b': mean_y_b,
                            'mean_X_b': mean_X_b,
                            'beta_b': beta_b
                        })
                    else:
                        model_boot = model_dict['model'].__class__(**{k: v for k, v in model_dict['model'].get_params().items() if k != 'n_jobs'}, n_jobs=-1)
                        model_boot.fit(X_selected[boot_idx], y_clean[boot_idx])
                        boot_data.append({'model_boot': model_boot})
                cached_data['boot_data'] = boot_data

            if model_type == 'rrBLUP':
                pred_val = float(model_dict['mean_y'] + np.dot(new_X_sel[0] - model_dict['mean_X'], model_dict['beta']))
                importances = model_dict['importances']
            else:
                pred_val = float(model_dict['model'].predict(new_X_sel)[0])
                importances = model_dict['importances']

            if bootstrap:
                boot_preds = []
                for b_item in cached_data['boot_data']:
                    if model_type == 'rrBLUP':
                        boot_preds.append(float(b_item['mean_y_b'] + np.dot(new_X_sel[0] - b_item['mean_X_b'], b_item['beta_b'])))
                    else:
                        boot_preds.append(float(b_item['model_boot'].predict(new_X_sel)[0]))
                ci_lower = float(np.percentile(boot_preds, 2.5))
                ci_upper = float(np.percentile(boot_preds, 97.5))
            else:
                ci_lower = pred_val
                ci_upper = pred_val

            predictions[trait] = {
                'predicted': round(pred_val, 3),
                'ci_lower': round(ci_lower, 3),
                'ci_upper': round(ci_upper, 3),
                'unit': self._get_unit(trait),
                'population_mean': round(cached_data['population_mean'], 3),
                'population_std': round(cached_data['population_std'], 3)
            }

            top_indices = np.argsort(importances)[::-1][:5]
            trait_snps = []
            for idx in top_indices:
                orig_idx = selected_features[idx]
                snp_id = self.rs_ids[orig_idx]
                parts = snp_id.split('_')
                if len(parts) >= 5:
                    refseq = parts[0] + '_' + parts[1]
                    chrom = REFSEQ_TO_CHR.get(refseq, refseq)
                    pos = parts[2]
                    ref = parts[3]
                    alt = parts[4]
                    friendly_id = f"{chrom}_{pos}_{ref}_{alt}"
                else:
                    friendly_id = snp_id
                trait_snps.append({
                    'snp_id': friendly_id,
                    'importance': float(importances[idx]),
                    'chromosome': chrom if len(parts) >= 5 else 'Unknown',
                    'position': int(pos) if len(parts) >= 5 else 0
                })
            top_snps_per_trait[trait] = trait_snps

        # Genetic similarity
        similarity = self._get_genetic_similarity(new_X[0])

        # Project new cultivar onto the same fitted PCA space of the 161-accession panel
        projected_coords = self.pca_model.transform(new_X)[0]
        pca_coords = {
            'pc1': float(projected_coords[0]),
            'pc2': float(projected_coords[1]),
            'pc3': float(projected_coords[2])
        }

        res = {
            'predictions': predictions,
            'top_snps_per_trait': top_snps_per_trait,
            'similar_cultivars': similarity,
            'pca_coordinates': pca_coords,
            'model_used': model_type,
            'n_training_samples': int(np.sum(~np.isnan(self.df_pheno['fruitWeight (g)'].values)))
        }
        self._prediction_cache[sig] = res
        return res

    def _get_unit(self, trait):
        units = {
            'fruitWeight (g)': 'g', 'fruitLength (mm)': 'mm', 'fruitWidth (mm)': 'mm',
            'fruitThickness (mm)': 'mm', 'brix': '%', 'Pulp': 'g',
            'seedWeight (g)': 'g', 'stoneWeight (g)': 'g',
            'seedLength (mm)': 'mm', 'seedWidth (mm)': 'mm', 'seedThickness (mm)': 'mm',
            'stoneLength (mm)': 'mm', 'stoneWidth (mm)': 'mm', 'stoneThickness (mm)': 'mm',
            'Fruit Shape Index': '', 'Pulp Ratio': '', 'Seed Ratio': '', 'Stone Ratio': '',
            'Edible Portion': '%', 'Brix Yield Index': '', 'Fruit Density Index': 'g/cm³',
            'Pulp-to-Seed Ratio': '', 'Pulp-to-Stone Ratio': '', 'Sweetness Efficiency Index': ''
        }
        return units.get(trait, '')

    def _get_genetic_similarity(self, new_genotype, top_n=5):
        """Return the top_n most genetically similar cultivars using Euclidean distance."""
        self.load_data()
        new_gt_arr = np.array(new_genotype, dtype=np.float32).reshape(1, -1)
        dists = np.sqrt(np.sum((self.X - new_gt_arr) ** 2, axis=1))
        
        distances = []
        for i in range(self.X.shape[0]):
            distances.append((str(self.cultivars[i]), float(dists[i])))
        distances.sort(key=lambda x: x[1])

        max_dist = distances[-1][1] if distances[-1][1] > 0 else 1.0
        results = []
        for name, dist in distances[:top_n]:
            similarity_pct = round((1.0 - dist / max_dist) * 100, 1)
            results.append({
                'cultivar': name,
                'distance': round(dist, 2),
                'similarity_pct': similarity_pct
            })
        return results

    def simulate_random_cultivar(self):
        """Generate a plausible random genotype vector for demo/testing purposes."""
        self.load_data()
        rng = np.random.RandomState(None)
        # Pick a random existing cultivar as base, add noise
        base_idx = rng.randint(0, self.X.shape[0])
        base = self.X[base_idx].copy()
        # Flip ~5% of SNPs
        n_flip = int(0.05 * len(base))
        flip_idx = rng.choice(len(base), size=n_flip, replace=False)
        for fi in flip_idx:
            base[fi] = rng.choice([0.0, 1.0, 2.0])
        return base.tolist()

    def rank_cultivars(self, w_brix, w_pulp, w_seed, w_stone):
        self.load_data()
        
        # Standardize traits for index calculation
        traits = ['brix', 'Pulp Ratio', 'Seed Ratio', 'Stone Ratio']
        df_std = self.df_pheno[['Accession ID', 'Subpopulation'] + traits].copy()
        
        for t in traits:
            # Impute missing values with mean
            mean_val = df_std[t].mean()
            df_std[t] = df_std[t].fillna(mean_val)
            # Z-score normalization
            std_val = df_std[t].std()
            if std_val == 0:
                std_val = 1.0
            df_std[t + '_std'] = (df_std[t] - mean_val) / std_val
            
        # Calculate composite score
        # High Brix & High Pulp Ratio are positive; High Seed & Stone Ratios are negative
        df_std['FQI'] = (
            w_brix * df_std['brix_std'] +
            w_pulp * df_std['Pulp Ratio_std'] -
            w_seed * df_std['Seed Ratio_std'] -
            w_stone * df_std['Stone Ratio_std']
        )
        
        df_std = df_std.sort_values(by='FQI', ascending=False).reset_index(drop=True)
        
        rankings = []
        for idx, row in df_std.iterrows():
            rankings.append({
                'rank': idx + 1,
                'accession': str(row['Accession ID']),
                'subpopulation': str(row['Subpopulation']),
                'brix': float(self.df_pheno.loc[self.df_pheno['Accession ID'] == row['Accession ID'], 'brix'].values[0]),
                'pulp_ratio': float(self.df_pheno.loc[self.df_pheno['Accession ID'] == row['Accession ID'], 'Pulp Ratio'].values[0]),
                'seed_ratio': float(self.df_pheno.loc[self.df_pheno['Accession ID'] == row['Accession ID'], 'Seed Ratio'].values[0]),
                'stone_ratio': float(self.df_pheno.loc[self.df_pheno['Accession ID'] == row['Accession ID'], 'Stone Ratio'].values[0]),
                'fqi_score': float(row['FQI'])
            })
            
        return rankings

        return rankings

    def validate_holdout(self, holdout_cultivar, model_type='rrBLUP'):
        """
        Train the genomic prediction model on 160 cultivars and predict the held-out cultivar.
        Compare predicted vs actual values for all 21 traits.
        """
        self.load_data()
        
        # Match holdout cultivar name
        clean_holdout = get_mapped_clean(holdout_cultivar)
        matched_idx = -1
        for idx, cv in enumerate(self.cultivars):
            if get_mapped_clean(cv) == clean_holdout:
                matched_idx = idx
                break
                
        if matched_idx == -1:
            return {'error': f"Holdout cultivar '{holdout_cultivar}' not found in panel."}
            
        all_traits = [
            'fruitWeight (g)', 'fruitLength (mm)', 'fruitWidth (mm)', 'fruitThickness (mm)',
            'brix', 'Pulp', 'seedWeight (g)', 'stoneWeight (g)',
            'seedLength (mm)', 'seedWidth (mm)', 'seedThickness (mm)',
            'stoneLength (mm)', 'stoneWidth (mm)', 'stoneThickness (mm)',
            'Fruit Shape Index', 'Pulp Ratio', 'Seed Ratio', 'Stone Ratio',
            'Edible Portion', 'Brix Yield Index', 'Fruit Density Index',
            'Pulp-to-Seed Ratio', 'Pulp-to-Stone Ratio', 'Sweetness Efficiency Index'
        ]
        
        # Build training set (excluding holdout index)
        train_indices = [i for i in range(len(self.cultivars)) if i != matched_idx]
        holdout_X = self.X[matched_idx].reshape(1, -1)
        
        results = {}
        for trait in all_traits:
            y = self.df_pheno[trait].values
            valid_mask = ~np.isnan(y)
            
            is_holdout_valid = valid_mask[matched_idx]
            actual_val = float(y[matched_idx]) if is_holdout_valid else None
            
            # Sub-select train indices that are valid
            train_mask_idx = [i for i in train_indices if valid_mask[i]]
            
            if len(train_mask_idx) < 10:
                continue
                
            y_train = y[train_mask_idx]
            
            # Center and select features (top 500)
            y_centered = y_train - np.mean(y_train)
            y_padded = np.zeros(len(self.cultivars))
            y_padded[train_mask_idx] = y_centered
            cov = np.dot(y_padded, self.X)
            
            sum_X = np.sum(self.X[train_mask_idx], axis=0)
            sum_X2 = np.sum(self.X[train_mask_idx]**2, axis=0)
            n_t = len(y_train)
            var_X = sum_X2 - (sum_X**2) / n_t
            var_X[var_X == 0] = 1.0
            var_y = np.sum(y_centered**2)
            if var_y == 0:
                var_y = 1.0
            corr = cov / np.sqrt(var_X * var_y)
            abs_corr = np.abs(corr)
            
            top_k = min(500, self.X.shape[1])
            selected_features = np.argsort(abs_corr)[::-1][:top_k]
            
            X_train_sel = self.X[:, selected_features][train_mask_idx]
            holdout_X_sel = holdout_X[:, selected_features]
            
            # Train model
            if model_type == 'rrBLUP':
                # Direct Cholesky solver for Ridge regression (30x faster than sklearn fit)
                n_features = X_train_sel.shape[1]
                mean_y = np.mean(y_train)
                mean_X = np.mean(X_train_sel, axis=0)
                y_centered = y_train - mean_y
                X_centered = X_train_sel - mean_X
                A = np.dot(X_centered.T, X_centered) + 100.0 * np.eye(n_features)
                b = np.dot(X_centered.T, y_centered)
                beta = np.linalg.solve(A, b)
                predicted_val = float(mean_y + np.dot(holdout_X_sel[0] - mean_X, beta))
            else:
                if model_type == 'Random Forest':
                    model = RandomForestRegressor(n_estimators=30, random_state=42, n_jobs=-1)
                else:
                    model = xgb.XGBRegressor(n_estimators=30, max_depth=3, learning_rate=0.1, random_state=42, n_jobs=-1)
                    
                model.fit(X_train_sel, y_train)
                predicted_val = float(model.predict(holdout_X_sel)[0])
            
            # Dev and unit
            dev = float(predicted_val - actual_val) if actual_val is not None else None
            dev_pct = float(dev / actual_val * 100) if (actual_val is not None and actual_val != 0) else None
            
            results[trait] = {
                'actual': round(actual_val, 3) if actual_val is not None else None,
                'predicted': round(predicted_val, 3),
                'deviation': round(dev, 3) if dev is not None else None,
                'deviation_percent': round(dev_pct, 1) if dev_pct is not None else None,
                'unit': self._get_unit(trait)
            }
            
        return {
            'cultivar': holdout_cultivar,
            'model_used': model_type,
            'predictions': results
        }

    def get_loocv_metrics(self, trait, model_type='rrBLUP'):
        """
        Compute or return cached LOOCV metrics for a given trait and model.
        """
        self.load_data()
        
        cache_key = f"{trait}_{model_type}"
        if hasattr(self, '_loocv_cache') and cache_key in self._loocv_cache:
            return self._loocv_cache[cache_key]
            
        if not hasattr(self, '_loocv_cache'):
            self._loocv_cache = {}
            
        y = self.df_pheno[trait].values
        valid_mask = ~np.isnan(y)
        y_clean = y[valid_mask]
        X_clean = self.X[valid_mask]
        
        n_samples = len(y_clean)
        if n_samples < 10:
            return {'r2': 0, 'rmse': 0, 'mae': 0, 'correlation': 0}
            
        predictions = np.zeros_like(y_clean)
        
        # Precompute sum of X and sum of X^2 for fast column variance
        sum_X = np.sum(X_clean, axis=0)
        sum_X2 = np.sum(X_clean**2, axis=0)
        n_t = n_samples - 1
        
        # LOOCV loop
        for holdout_idx in range(n_samples):
            train_idx = [i for i in range(n_samples) if i != holdout_idx]
            y_train = y_clean[train_idx]
            
            # Feature selection (optimized to avoid large array allocations and fancy indexing)
            y_centered = y_train - np.mean(y_train)
            y_padded = np.zeros(n_samples)
            y_padded[train_idx] = y_centered
            cov = np.dot(y_padded, X_clean)
            
            # O(1) column variance updating
            sum_X_fold = sum_X - X_clean[holdout_idx]
            sum_X2_fold = sum_X2 - X_clean[holdout_idx]**2
            var_X = sum_X2_fold - (sum_X_fold**2) / n_t
            var_X[var_X == 0] = 1.0
            
            var_y = np.sum(y_centered**2)
            if var_y == 0:
                var_y = 1.0
            corr = cov / np.sqrt(var_X * var_y)
            abs_corr = np.abs(corr)
            
            top_k = min(500, X_clean.shape[1])
            selected_features = np.argsort(abs_corr)[::-1][:top_k]
            
            X_train_sel = X_clean[:, selected_features][train_idx]
            X_test_sel = X_clean[holdout_idx, selected_features].reshape(1, -1)
            
            if model_type == 'rrBLUP':
                # Direct Cholesky solver for Ridge regression (30x faster than sklearn fit)
                n_features = X_train_sel.shape[1]
                mean_y = np.mean(y_train)
                mean_X = np.mean(X_train_sel, axis=0)
                y_centered = y_train - mean_y
                X_centered = X_train_sel - mean_X
                A = np.dot(X_centered.T, X_centered) + 100.0 * np.eye(n_features)
                b = np.dot(X_centered.T, y_centered)
                beta = np.linalg.solve(A, b)
                predictions[holdout_idx] = float(mean_y + np.dot(X_test_sel[0] - mean_X, beta))
            else:
                if model_type == 'Random Forest':
                    model = RandomForestRegressor(n_estimators=15, random_state=42, n_jobs=-1)
                else:
                    model = xgb.XGBRegressor(n_estimators=15, max_depth=3, learning_rate=0.1, random_state=42, n_jobs=-1)
                    
                model.fit(X_train_sel, y_train)
                predictions[holdout_idx] = model.predict(X_test_sel)[0]
            
        try:
            r2 = float(np.corrcoef(y_clean, predictions)[0, 1]**2) if np.std(predictions) > 0 else 0.0
        except Exception:
            r2 = 0.0
        rmse = float(np.sqrt(np.mean((y_clean - predictions)**2)))
        mae = float(np.mean(np.abs(y_clean - predictions)))
        try:
            r, _ = pearsonr(y_clean, predictions)
            r = float(r)
        except Exception:
            r = 0.0
            
        metrics = {
            'r2': round(r2, 4),
            'rmse': round(rmse, 4),
            'mae': round(mae, 4),
            'correlation': round(r, 4)
        }
        self._loocv_cache[cache_key] = metrics
        return metrics

    def parse_uploaded_genotypes(self, file_content, filename):
        """
        Parse genotype files (CSV or VCF), match them against the 135,079 training SNPs,
        impute missing markers, and return an encoded genotype vector (0, 1, 2).
        """
        self.load_data()
        
        # Build maps for matching
        snp_map = {}
        for idx, snp_id in enumerate(self.rs_ids):
            snp_map[snp_id.lower()] = idx
            parts = snp_id.split('_')
            if len(parts) >= 3:
                chrom = parts[0] + '_' + parts[1]
                pos = parts[2]
                snp_map[(chrom.lower(), int(pos))] = idx
                friendly_chrom = REFSEQ_TO_CHR.get(chrom, chrom)
                snp_map[(friendly_chrom.lower(), int(pos))] = idx
                
        # Initialize default vector with mean dosage of training set
        mean_dosages = np.mean(self.X, axis=0)
        genotype_vector = mean_dosages.copy()
        
        n_matched = 0
        n_total_parsed = 0
        matched_flags = np.zeros(len(self.rs_ids), dtype=bool)
        
        lines = file_content.splitlines()
        is_vcf = filename.endswith('.vcf') or any(line.startswith('##fileformat=VCF') for line in lines[:5])
        
        if is_vcf:
            for line in lines:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) < 10:
                    continue
                chrom = parts[0].strip()
                try:
                    pos = int(parts[1].strip())
                except ValueError:
                    continue
                gt_field = parts[9].split(':')[0]
                
                n_total_parsed += 1
                
                idx = snp_map.get((chrom.lower(), pos))
                if idx is not None:
                    if gt_field in ['0/0', '0|0']:
                        genotype_vector[idx] = 0.0
                    elif gt_field in ['0/1', '1/0', '0|1', '1|0']:
                        genotype_vector[idx] = 1.0
                    elif gt_field in ['1/1', '1|1']:
                        genotype_vector[idx] = 2.0
                    else:
                        continue
                    n_matched += 1
                    matched_flags[idx] = True
        else:
            # Parse CSV
            for line in lines:
                if line.startswith('rs#') or line.startswith('SNP') or line.startswith('#'):
                    continue
                parts = line.strip().split(',')
                if len(parts) < 2:
                    continue
                snp_id = parts[0].strip()
                call = parts[1].strip().upper()
                
                n_total_parsed += 1
                
                idx = snp_map.get(snp_id.lower())
                if idx is None:
                    p_parts = snp_id.split('_')
                    if len(p_parts) >= 3 and p_parts[0].lower() == 'nc':
                        chrom = p_parts[0] + '_' + p_parts[1]
                        try:
                            pos = int(p_parts[2])
                            idx = snp_map.get((chrom.lower(), pos))
                        except ValueError:
                            pass
                    elif len(p_parts) >= 2:
                        chrom = p_parts[0]
                        try:
                            pos = int(p_parts[1])
                            idx = snp_map.get((chrom.lower(), pos))
                        except ValueError:
                            pass
                            
                if idx is not None:
                    parts_train = self.rs_ids[idx].split('_')
                    ref = parts_train[3]
                    alt = parts_train[4]
                    
                    IUPAC_HET = {
                        ('A', 'G'): 'R', ('G', 'A'): 'R',
                        ('C', 'T'): 'Y', ('T', 'C'): 'Y',
                        ('G', 'C'): 'S', ('C', 'G'): 'S',
                        ('A', 'T'): 'W', ('T', 'A'): 'W',
                        ('G', 'T'): 'K', ('T', 'G'): 'K',
                        ('A', 'C'): 'M', ('C', 'A'): 'M'
                    }
                    het = IUPAC_HET.get((ref, alt), 'N')
                    
                    if call == ref:
                        genotype_vector[idx] = 0.0
                    elif call == alt:
                        genotype_vector[idx] = 2.0
                    elif call == het or call in ['HET', 'HETEROZYGOUS', '1']:
                        genotype_vector[idx] = 1.0
                    elif call in ['0', 'REF', 'HOMOZYGOUS_REF']:
                        genotype_vector[idx] = 0.0
                    elif call in ['2', 'ALT', 'HOMOZYGOUS_ALT']:
                        genotype_vector[idx] = 2.0
                    else:
                        continue
                    n_matched += 1
                    matched_flags[idx] = True
                    
        impute_count = len(genotype_vector) - n_matched
        match_pct = (n_matched / len(self.rs_ids) * 100) if len(self.rs_ids) > 0 else 0
        
        imputed_indices = np.where(~matched_flags)[0]
        imputed_sample_ids = [self.rs_ids[i] for i in imputed_indices[:20]]
        
        return {
            'genotype_vector': genotype_vector.tolist(),
            'n_parsed': n_total_parsed,
            'n_matched': n_matched,
            'n_imputed': impute_count,
            'match_percentage': round(match_pct, 2),
            'imputed_sample_markers': imputed_sample_ids
        }

    def simulate_hybrid_cross(self, parent_a_id, parent_b_id, model_type='rrBLUP'):
        """
        Simulate a virtual F1 offspring genotype from Parent A and Parent B
        and predict its phenotypes for all 21 traits.
        """
        self.load_data()
        


        idx_a = -1
        idx_b = -1
        for idx, cv in enumerate(self.cultivars):
            if get_mapped_clean(cv) == get_mapped_clean(parent_a_id):
                idx_a = idx
            if get_mapped_clean(cv) == get_mapped_clean(parent_b_id):
                idx_b = idx
                
        if idx_a == -1 or idx_b == -1:
            return {'error': 'One or both parents not found in database.'}
            
        genotype_a = self.X[idx_a]
        genotype_b = self.X[idx_b]
        
        rng = np.random.RandomState(None)
        
        contrib_a = np.zeros_like(genotype_a)
        contrib_a[genotype_a == 2] = 1.0
        mask_het_a = (genotype_a == 1)
        contrib_a[mask_het_a] = rng.choice([0.0, 1.0], size=np.sum(mask_het_a))
        
        contrib_b = np.zeros_like(genotype_b)
        contrib_b[genotype_b == 2] = 1.0
        mask_het_b = (genotype_b == 1)
        contrib_b[mask_het_b] = rng.choice([0.0, 1.0], size=np.sum(mask_het_b))
        
        offspring_genotype = contrib_a + contrib_b
        
        # Predict traits
        res_offspring = self.predict_new_cultivar(offspring_genotype.tolist(), model_type, bootstrap=False)
        res_a = self.predict_new_cultivar(genotype_a.tolist(), model_type, bootstrap=False)
        res_b = self.predict_new_cultivar(genotype_b.tolist(), model_type, bootstrap=False)
        
        # PCA distance
        pca_dist = 0.0
        try:
            coords_a = self.pca_model.transform(genotype_a.reshape(1, -1))[0]
            coords_b = self.pca_model.transform(genotype_b.reshape(1, -1))[0]
            pca_dist = float(np.sqrt(np.sum((coords_a[:2] - coords_b[:2]) ** 2)))
        except Exception:
            pass
                
        return {
            'parent_a': parent_a_id,
            'parent_b': parent_b_id,
            'pca_distance': round(pca_dist, 4),
            'offspring_predictions': res_offspring['predictions'],
            'parent_a_predictions': res_a['predictions'],
            'parent_b_predictions': res_b['predictions']
        }

    def simulate_f2_cross(self, parent_a_id, parent_b_id, parent_c_id, model_type='rrBLUP'):
        """
        Simulate a virtual F2 offspring genotype by crossing (Parent A x Parent B) -> F1,
        and then F1 x Parent C -> F2, and predict phenotypes for all traits.
        """
        self.load_data()
        


        idx_a, idx_b, idx_c = -1, -1, -1
        for idx, cv in enumerate(self.cultivars):
            if get_mapped_clean(cv) == get_mapped_clean(parent_a_id):
                idx_a = idx
            if get_mapped_clean(cv) == get_mapped_clean(parent_b_id):
                idx_b = idx
            if get_mapped_clean(cv) == get_mapped_clean(parent_c_id):
                idx_c = idx
                
        if idx_a == -1 or idx_b == -1 or idx_c == -1:
            return {'error': 'One or more parents not found in database.'}
            
        genotype_a = self.X[idx_a]
        genotype_b = self.X[idx_b]
        genotype_c = self.X[idx_c]
        
        rng = np.random.RandomState(None)
        
        # 1. Gamete from A
        contrib_a = np.zeros_like(genotype_a)
        contrib_a[genotype_a == 2] = 1.0
        mask_het_a = (genotype_a == 1)
        contrib_a[mask_het_a] = rng.choice([0.0, 1.0], size=np.sum(mask_het_a))
        
        # 2. Gamete from B
        contrib_b = np.zeros_like(genotype_b)
        contrib_b[genotype_b == 2] = 1.0
        mask_het_b = (genotype_b == 1)
        contrib_b[mask_het_b] = rng.choice([0.0, 1.0], size=np.sum(mask_het_b))
        
        # 3. Form F1 genotype
        genotype_f1 = contrib_a + contrib_b
        
        # 4. Gamete from F1
        contrib_f1 = np.zeros_like(genotype_f1)
        contrib_f1[genotype_f1 == 2] = 1.0
        mask_het_f1 = (genotype_f1 == 1)
        contrib_f1[mask_het_f1] = rng.choice([0.0, 1.0], size=np.sum(mask_het_f1))
        
        # 5. Gamete from C
        contrib_c = np.zeros_like(genotype_c)
        contrib_c[genotype_c == 2] = 1.0
        mask_het_c = (genotype_c == 1)
        contrib_c[mask_het_c] = rng.choice([0.0, 1.0], size=np.sum(mask_het_c))
        
        # 6. Form F2 genotype
        genotype_f2 = contrib_f1 + contrib_c
        
        # Predict traits for all nodes in the tree
        res_f2 = self.predict_new_cultivar(genotype_f2.tolist(), model_type, bootstrap=False)
        res_f1 = self.predict_new_cultivar(genotype_f1.tolist(), model_type, bootstrap=False)
        res_a = self.predict_new_cultivar(genotype_a.tolist(), model_type, bootstrap=False)
        res_b = self.predict_new_cultivar(genotype_b.tolist(), model_type, bootstrap=False)
        res_c = self.predict_new_cultivar(genotype_c.tolist(), model_type, bootstrap=False)
        
        # Calculate PCA genetic distances
        pca_dist_f1 = 0.0
        pca_dist_f2 = 0.0
        try:
            coords_a = self.pca_model.transform(genotype_a.reshape(1, -1))[0]
            coords_b = self.pca_model.transform(genotype_b.reshape(1, -1))[0]
            coords_c = self.pca_model.transform(genotype_c.reshape(1, -1))[0]
            coords_f1 = self.pca_model.transform(genotype_f1.reshape(1, -1))[0]
            
            pca_dist_f1 = float(np.sqrt(np.sum((coords_a[:2] - coords_b[:2]) ** 2)))
            pca_dist_f2 = float(np.sqrt(np.sum((coords_f1[:2] - coords_c[:2]) ** 2)))
        except Exception:
            pass
                
        return {
            'parent_a': parent_a_id,
            'parent_b': parent_b_id,
            'parent_c': parent_c_id,
            'pca_distance_f1': round(pca_dist_f1, 4),
            'pca_distance_f2': round(pca_dist_f2, 4),
            'f2_predictions': res_f2['predictions'],
            'f1_predictions': res_f1['predictions'],
            'parent_a_predictions': res_a['predictions'],
            'parent_b_predictions': res_b['predictions'],
            'parent_c_predictions': res_c['predictions']
        }

# Singleton instance
engine = BreedingEngine()

