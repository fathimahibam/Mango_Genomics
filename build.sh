#!/usr/bin/env bash
set -e

echo "=== MangoGenomics Lab - Render Build ==="

# Install Python dependencies
echo "[1/3] Installing Python dependencies..."
pip install --upgrade pip
pip install -r requirements.txt

# Determine project directory (Render sets this, or fall back to current dir)
PROJECT_DIR="${RENDER_SRC_DIR:-$(pwd)}"
echo "[2/3] Project directory: $PROJECT_DIR"

# Decompress GFF if needed
if [ ! -f "$PROJECT_DIR/mango_real.gff" ] && [ -f "$PROJECT_DIR/mango_real.gff.gz" ]; then
    echo "Decompressing mango_real.gff.gz (this may take a few minutes)..."
    gzip -dk "$PROJECT_DIR/mango_real.gff.gz"
    echo "GFF decompressed successfully."
else
    echo "GFF file already exists or no .gz found, skipping."
fi

echo "[3/3] Build complete!"
ls -lh "$PROJECT_DIR/mango_genes_real.db" 2>/dev/null || echo "  NOTE: mango_genes_real.db not found (LFS may still be pulling)"
ls -lh "$PROJECT_DIR/mango_real.gff" 2>/dev/null || echo "  NOTE: mango_real.gff not found"
ls -lh "$PROJECT_DIR/SNP.csv" 2>/dev/null || echo "  NOTE: SNP.csv not found"
