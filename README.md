# MangoGenomics Lab

A premium bioinformatics research platform for **Mangifera indica** (mango) genomics, featuring real NCBI-annotated genome data, GWAS analysis, genomic prediction, CRISPR simulation, and pedigree exploration.

## Features

- **Genome Dashboard** — Real NCBI GCF_011075055.1 (CATAS_Mindica_2.1) with 15,414+ genes
- **GWAS Lab** — Genome-wide association studies across 161 accessions with 135,079 SNPs
- **Genomic Prediction** — rrBLUP, Random Forest, and XGBoost models for 29 traits
- **Holdout Validation** — Leave-one-out cross-validation for prediction accuracy
- **CRISPR Studio** — PAM site finder and HDR gene-editing simulation
- **Genetic Pedigree** — Breeding cross simulator (F1, F2) with phenotype prediction
- **Breeding Report** — Printable PDF breeding value reports

## Tech Stack

- **Backend:** Python Flask + Waitress WSGI server
- **Data:** SQLite (gffutils), NumPy, Pandas, scikit-learn, XGBoost
- **Frontend:** Vanilla HTML/CSS/JS with interactive visualizations
- **Deployment:** Render (Docker-free, Git-based)

## Local Development

```bash
pip install -r requirements.txt
python run.py
# Visit http://localhost:5000
```

## Deployment (Render)

1. Fork/clone this repository
2. Create a new Web Service on [Render](https://render.com)
3. Connect your GitHub repo
4. Use these settings:
   - **Runtime:** Python 3.11
   - **Build Command:** `bash build.sh`
   - **Start Command:** `python run.py`
   - **Plan:** Standard (2GB+ RAM recommended)
5. Render auto-detects `render.yaml` for Blueprint deployments

## Data Sources

- **Genome:** NCBI GCF_011075055.1 (CATAS_Mindica_2.1)
- **SNPs:** 135,079 markers across 161 mango accessions
- **Phenotypes:** 29 traits (fruit weight, brix, pulp ratio, etc.)

## License

MIT License
