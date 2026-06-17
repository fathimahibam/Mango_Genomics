# MangoGenomics Lab - Production Deployment Guide

This document describes how to deploy the MangoGenomics Lab application to cloud platforms (Render, Railway, AWS, Azure) and run it in a production environment.

---

## Architecture Overview

1. **Backend**: Python Flask app served via the multi-threaded **Waitress WSGI** server (`run.py`).
2. **Databases**: SQLite (`mango_genes_real.db`, `synthetic_cultivars.db`).
3. **Frontend**: Vanilla HTML/CSS/JavaScript under `templates/` and `static/`.
4. **Environment Modes**:
   - **Local Dev (Windows)**: Run with Python 3.14.5, utilizing precompiled C-extensions in the [libs](file:///C:/mangoproject/libs) folder.
   - **Cloud Prod (Linux/Windows)**: Install standard dependencies from [requirements.txt](file:///C:/mangoproject/requirements.txt) via `pip` on standard Python environments (Python 3.10 to 3.12). The precompiled `libs` folder will be skipped automatically on non-Windows/non-3.14 platforms.

---

## 1. Preparing the Deployment Package

Ensure the following files are present in the root folder of your repository:
- `requirements.txt` (List of required packages)
- `run.py` (WSGI server entrypoint)
- `app.py` (Main Flask application logic)
- `prediction_engine.py` (Genomic prediction engine)
- `SNP.csv` (SNP dataset)
- `phenotypes.xlsx` (Phenotype data)
- `mango_genes_real.db` (Gene lookup SQLite database)
- `static/data/loocv_metrics.json` (Precomputed LOOCV metrics)

---

## 2. Cloud Deployment Options

### Option A: Render (Recommended)
Render is the easiest way to deploy python web services.
1. Sign in to [Render](https://render.com) and create a **New Web Service**.
2. Connect your Git repository.
3. Configure the following settings:
   - **Environment**: `Python`
   - **Build Command**: `pip install -r requirements.txt`
   - **Start Command**: `python run.py`
   - **Instance Type**: Select an instance with at least **2GB RAM** (due to the large SQLite database and SNP file memory requirements).
4. Click **Deploy**. Render will bind to the `$PORT` environment variable automatically.

### Option B: Railway
Railway offers fast deployments and automatic scaling.
1. Sign in to [Railway](https://railway.app) and create a **New Project**.
2. Click **Deploy from GitHub repo** and select your repository.
3. Railway will automatically detect the Python environment.
4. Go to **Variables** and add:
   - `PORT`: `5000` (or leave empty, Railway handles dynamic ports).
5. Go to **Settings** and ensure the start command is:
   - `python run.py`
6. Click **Deploy**.

### Option C: Microsoft Azure Web App
Deploy as a Python Linux App Service.
1. Create a **Web App** in Azure Portal:
   - **Publish**: Code
   - **Runtime Stack**: `Python 3.10` or `Python 3.11`
   - **Operating System**: Linux
   - **Sku**: Basic (B1) or higher (recommended >= 2GB memory).
2. Deploy using Azure CLI, Local Git, or GitHub Actions.
3. In App Service **Configuration** settings:
   - Add Application Setting: `SCM_DO_BUILD_DURING_DEPLOYMENT` = `true`
   - Set Startup Command: `python run.py`

### Option D: AWS (Elastic Beanstalk)
1. Initialize Elastic Beanstalk using the CLI: `eb init -p python-3.11 mango-genomics`
2. Create environment: `eb create mango-genomics-env`
3. Configure environment properties to set the command to: `python run.py`
4. Deploy updates via `eb deploy`.

---

## 3. Post-Deployment Verification

Once deployed, visit your cloud URL and confirm:
1. The **GWAS Lab** and **Genomic Prediction** tabs load instantly (under 1s).
2. The **Hybrid Cross Simulation** (F1) completes under **1 second** (Target: < 5s).
3. The **Pedigree Explorer** (F2) completes under **2 seconds** (Target: < 10s).
4. The uploaded genotype parser works using `dummy_genotype.csv` from the web interface.
