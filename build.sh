#!/usr/bin/env bash
set -e

echo "=== MangoGenomics Lab - Render Build ==="

# Install Python dependencies
echo "[1/3] Installing Python dependencies..."
pip install --upgrade pip
pip install -r requirements.txt

# Decompress GFF if needed (the .gz is in the repo, the uncompressed may not be)
GFF_DIR="/opt/render/project/src"
if [ ! -f "$GFF_DIR/mango_real.gff" ] && [ -f "$GFF_DIR/mango_real.gff.gz" ]; then
    echo "[2/3] Decompressing mango_real.gff.gz (this may take a few minutes)..."
    gzip -dk "$GFF_DIR/mango_real.gff.gz"
    echo "GFF decompressed successfully."
else
    echo "[2/3] GFF file already exists or no .gz found, skipping decompression."
fi

echo "[3/3] Build complete!"
echo "Files ready:"
ls -lh "$GFF_DIR/mango_genes_real.db" 2>/dev/null || echo "  WARNING: mango_genes_real.db not found"
ls -lh "$GFF_DIR/mango_real.gff" 2>/dev/null || echo "  WARNING: mango_real.gff not found"
ls -lh "$GFF_DIR/SNP.csv" 2>/dev/null || echo "  WARNING: SNP.csv not found"
