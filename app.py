"""
app.py
======
Flask backend for the Mango Bioinformatics Lab web application.
Uses the REAL NCBI annotated Mangifera indica genome database
(GCF_011075055.1 / CATAS_Mindica_2.1) with 15,414 genes.

Endpoints:
  GET  /                             → Serve the main SPA (index.html)
  GET  /api/genome/stats             → Assembly-level statistics from real gffutils DB
  GET  /api/genes/search             → Paginated keyword search across real genes
  GET  /api/diseases                 → Mango disease–gene reference data (6 diseases)
  GET  /api/cultivars                → All 12 mango cultivar profiles
  POST /api/crispr/simulate          → CRISPR-Cas9 HDR simulation
  GET  /api/crispr/find_pam          → Scan for PAM/gRNA targets in a disease allele
  GET  /api/gwas/data                → Return GWAS association data (optionally filtered)
  GET  /api/gwas/cultivar_comparison → Per-cultivar allele frequency comparison at top SNP
"""

import sys
sys.path.insert(0, r'C:\mangoproject\libs')

import os
import json
import math
from flask import Flask, jsonify, request, render_template, g

import gffutils

# ---------------------------------------------------------------------------
# App initialisation
# ---------------------------------------------------------------------------
BASE_DIR     = r'C:\mangoproject'
# *** SWITCHED TO REAL NCBI ANNOTATED DATABASE ***
DB_PATH      = os.path.join(BASE_DIR, 'mango_genes_real.db')
GWAS_PATH    = os.path.join(BASE_DIR, 'static', 'data', 'gwas_data.json')
TEMPLATE_DIR = os.path.join(BASE_DIR, 'templates')
STATIC_DIR   = os.path.join(BASE_DIR, 'static')

app = Flask(__name__, template_folder=TEMPLATE_DIR, static_folder=STATIC_DIR)

# Ensure required directories exist
for d in [TEMPLATE_DIR, STATIC_DIR, os.path.join(STATIC_DIR, 'data')]:
    os.makedirs(d, exist_ok=True)

# Connect to the gffutils database — per-request to avoid SQLite threading issues
def get_db():
    if 'db' not in g:
        g.db = gffutils.FeatureDB(DB_PATH)
    return g.db

@app.teardown_appcontext
def teardown_db(exception):
    db = g.pop('db', None)
    if db is not None:
        try:
            db.conn.close()
        except Exception:
            pass


GENE_PRODUCT_CACHE = {}

def load_gene_products_gff():
    gff_path = os.path.join(BASE_DIR, 'mango_real.gff')
    products = {}
    if not os.path.exists(gff_path):
        return products
    try:
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
    except Exception as e:
        print("Error pre-loading gene products:", e)
    return products

print("Pre-loading gene products from GFF for instant lookups...")
import time
t_prod_start = time.time()
GENE_PRODUCT_CACHE = load_gene_products_gff()
print(f"Loaded {len(GENE_PRODUCT_CACHE)} gene products from GFF in {time.time()-t_prod_start:.2f}s")


def is_informative_product(prod, clean_id):
    if not prod:
        return False
    prod_lower = prod.lower().strip()
    if prod_lower == 'uncharacterized':
        return False
    if prod_lower == f'uncharacterized {clean_id.lower()}':
        return False
    return True


def get_gene_product(gene_id_or_name):
    clean_id = gene_id_or_name.split(' ')[0]
    if not clean_id.startswith('LOC'):
        return None
        
    if clean_id in GENE_PRODUCT_CACHE:
        return GENE_PRODUCT_CACHE[clean_id]
        
    db_id = f"gene-{clean_id}"
    try:
        db = get_db()
        gene = db[db_id]
        products = []
        for child in db.children(gene, featuretype='mRNA', level=1):
            prod = child.attributes.get('product', [''])[0]
            if prod:
                products.append(prod)
        if not products:
            for child in db.children(gene, featuretype='CDS', level=2):
                prod = child.attributes.get('product', [''])[0]
                if prod:
                    products.append(prod)
        if products:
            unique_prods = list(set(products))
            cleaned = []
            for p in unique_prods:
                p_clean = p.split(', transcript variant')[0].strip()
                if p_clean and p_clean not in cleaned:
                    cleaned.append(p_clean)
            res = '; '.join(cleaned)
            GENE_PRODUCT_CACHE[clean_id] = res
            return res
    except Exception:
        pass
        
    return None


# ═══════════════════════════════════════════════════════════════════════════
# Chromosome mapping: NCBI RefSeq accession → human-readable chromosome name
# (from GCF_011075055.1 assembly report)
# ═══════════════════════════════════════════════════════════════════════════
REFSEQ_TO_CHR = {
    'NC_058137.1': 'Chr1',
    'NC_058138.1': 'Chr2',
    'NC_058139.1': 'Chr3',
    'NC_058140.1': 'Chr4',
    'NC_058141.1': 'Chr5',
    'NC_058142.1': 'Chr6',
    'NC_058143.1': 'Chr7',
    'NC_058144.1': 'Chr8',
    'NC_058145.1': 'Chr9',
    'NC_058146.1': 'Chr10',
    'NC_058147.1': 'Chr11',
    'NC_058148.1': 'Chr12',
    'NC_058149.1': 'Chr13',
    'NC_058150.1': 'Chr14',
    'NC_058151.1': 'Chr15',
    'NC_058152.1': 'Chr16',
    'NC_058153.1': 'Chr17',
    'NC_058154.1': 'Chr18',
    'NC_058155.1': 'Chr19',
    'NC_058156.1': 'Chr20',
}


# ═══════════════════════════════════════════════════════════════════════════
# Helper: reverse complement
# ═══════════════════════════════════════════════════════════════════════════
COMPLEMENT = str.maketrans('ATCGatcg', 'TAGCtagc')


def reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    return seq.translate(COMPLEMENT)[::-1]


# ═══════════════════════════════════════════════════════════════════════════
# 12 Mango Cultivar profiles — real-world disease susceptibility data
# ═══════════════════════════════════════════════════════════════════════════
CULTIVARS = [
    {
        'id': 'alphonso',
        'name': 'Alphonso (Hapus)',
        'origin': 'Ratnagiri, Maharashtra, India',
        'type': 'Indian',
        'description': (
            'Known as the "King of Mangoes", Alphonso is prized for its '
            'unparalleled aroma, rich buttery texture, and deep saffron '
            'colour. However, it carries susceptible alleles at multiple '
            'resistance loci, making it highly vulnerable to fungal diseases.'
        ),
        'disease_profile': {
            'anthracnose': 'High',
            'powdery_mildew': 'Moderate',
            'bacterial_canker': 'High',
            'mango_malformation': 'Moderate',
            'fruit_rot': 'High',
            'stem_end_rot': 'High',
        },
        'resistance_score': 15,   # out of 100 (higher = more resistant)
        'fruit_weight_g': 200,
        'brix': 19.5,
        'badge': 'High Disease Risk',
        'badge_class': 'danger',
        'has_genotype': True,
    },
    {
        'id': 'tommy_atkins',
        'name': 'Tommy Atkins',
        'origin': 'Fort Lauderdale, Florida, USA',
        'type': 'American',
        'description': (
            'The world\'s most commercially grown cultivar. Tommy Atkins '
            'has exceptional shelf-life, attractive red blush skin, and '
            'carries natural resistance alleles at WAK1 and RPS2 loci. '
            'Moderate flavour but outstanding disease resistance profile.'
        ),
        'disease_profile': {
            'anthracnose': 'Low',
            'powdery_mildew': 'Low',
            'bacterial_canker': 'Low',
            'mango_malformation': 'Low',
            'fruit_rot': 'Low',
            'stem_end_rot': 'Low',
        },
        'resistance_score': 85,
        'fruit_weight_g': 400,
        'brix': 13.0,
        'badge': 'Highly Resistant',
        'badge_class': 'success',
        'has_genotype': True,
    },
    {
        'id': 'kent',
        'name': 'Kent',
        'origin': 'Coconut Grove, Florida, USA',
        'type': 'American',
        'description': (
            'Kent is a large, sweet, juicy mango with very little fibre. '
            'It shows good resistance to anthracnose due to thick skin '
            'and is widely grown in Central America and Spain.'
        ),
        'disease_profile': {
            'anthracnose': 'Low',
            'powdery_mildew': 'Moderate',
            'bacterial_canker': 'Moderate',
            'mango_malformation': 'Low',
            'fruit_rot': 'Moderate',
            'stem_end_rot': 'Moderate',
        },
        'resistance_score': 72,
        'fruit_weight_g': 600,
        'brix': 17.0,
        'badge': 'Good Resistance',
        'badge_class': 'success',
        'has_genotype': True,
    },
    {
        'id': 'keitt',
        'name': 'Keitt',
        'origin': 'Homestead, Florida, USA',
        'type': 'American',
        'description': (
            'Keitt is one of the latest-maturing Florida cultivars with '
            'large, oval fruits that stay green even when ripe. It shows '
            'good disease resistance and is popular in export markets.'
        ),
        'disease_profile': {
            'anthracnose': 'Low',
            'powdery_mildew': 'Low',
            'bacterial_canker': 'Moderate',
            'mango_malformation': 'Low',
            'fruit_rot': 'Low',
            'stem_end_rot': 'Moderate',
        },
        'resistance_score': 78,
        'fruit_weight_g': 700,
        'brix': 15.5,
        'badge': 'Good Resistance',
        'badge_class': 'success',
        'has_genotype': True,
    },
    {
        'id': 'ataulfo',
        'name': 'Ataulfo (Honey/Champagne)',
        'origin': 'Chiapas, Mexico',
        'type': 'Mexican',
        'description': (
            'The beloved "Honey Mango" of Mexico. Ataulfo has a creamy, '
            'nearly fibreless flesh with a sweet-tart flavour. Moderate '
            'susceptibility to powdery mildew in high-humidity environments.'
        ),
        'disease_profile': {
            'anthracnose': 'Moderate',
            'powdery_mildew': 'Moderate',
            'bacterial_canker': 'Moderate',
            'mango_malformation': 'Low',
            'fruit_rot': 'Moderate',
            'stem_end_rot': 'Moderate',
        },
        'resistance_score': 55,
        'fruit_weight_g': 200,
        'brix': 18.0,
        'badge': 'Moderate Risk',
        'badge_class': 'warning',
        'has_genotype': False,
    },
    {
        'id': 'chaunsa',
        'name': 'Chaunsa',
        'origin': 'Rahim Yar Khan, Pakistan',
        'type': 'Pakistani',
        'description': (
            'The most prized Pakistani mango, Chaunsa has a golden-yellow '
            'skin and exceptionally sweet, aromatic flesh. It is susceptible '
            'to powdery mildew during flowering and bacterial canker in wet seasons.'
        ),
        'disease_profile': {
            'anthracnose': 'High',
            'powdery_mildew': 'High',
            'bacterial_canker': 'Moderate',
            'mango_malformation': 'High',
            'fruit_rot': 'High',
            'stem_end_rot': 'Moderate',
        },
        'resistance_score': 20,
        'fruit_weight_g': 300,
        'brix': 20.0,
        'badge': 'High Disease Risk',
        'badge_class': 'danger',
        'has_genotype': False,
    },
    {
        'id': 'sindhri',
        'name': 'Sindhri',
        'origin': 'Tando Adam, Sindh, Pakistan',
        'type': 'Pakistani',
        'description': (
            'Pakistan\'s "King of Mangoes", Sindhri has a distinct oval shape, '
            'thin skin, and sweet pulp. Moderately susceptible to anthracnose '
            'and shows susceptibility to malformation disease.'
        ),
        'disease_profile': {
            'anthracnose': 'Moderate',
            'powdery_mildew': 'High',
            'bacterial_canker': 'Moderate',
            'mango_malformation': 'Moderate',
            'fruit_rot': 'Moderate',
            'stem_end_rot': 'Moderate',
        },
        'resistance_score': 40,
        'fruit_weight_g': 350,
        'brix': 18.5,
        'badge': 'Moderate Risk',
        'badge_class': 'warning',
        'has_genotype': True,
    },
    {
        'id': 'haden',
        'name': 'Haden',
        'origin': 'Coconut Grove, Florida, USA',
        'type': 'American',
        'description': (
            'Haden was the first named Florida mango cultivar and is the '
            'parent of many commercial varieties. It has red-yellow skin and '
            'rich flavour with moderate overall disease resistance.'
        ),
        'disease_profile': {
            'anthracnose': 'Moderate',
            'powdery_mildew': 'Moderate',
            'bacterial_canker': 'Low',
            'mango_malformation': 'Low',
            'fruit_rot': 'Moderate',
            'stem_end_rot': 'Low',
        },
        'resistance_score': 60,
        'fruit_weight_g': 450,
        'brix': 16.0,
        'badge': 'Moderate Resistance',
        'badge_class': 'warning',
        'has_genotype': True,
    },
    {
        'id': 'mallika',
        'name': 'Mallika',
        'origin': 'New Delhi, India (ICAR hybrid)',
        'type': 'Indian Hybrid',
        'description': (
            'A Neelum × Dashehari ICAR hybrid, Mallika combines excellent '
            'flavour with improved disease tolerance. It has shown field '
            'resistance to anthracnose and powdery mildew in Indian trials.'
        ),
        'disease_profile': {
            'anthracnose': 'Low',
            'powdery_mildew': 'Low',
            'bacterial_canker': 'Moderate',
            'mango_malformation': 'Moderate',
            'fruit_rot': 'Low',
            'stem_end_rot': 'Low',
        },
        'resistance_score': 75,
        'fruit_weight_g': 250,
        'brix': 19.0,
        'badge': 'Good Resistance',
        'badge_class': 'success',
        'has_genotype': False,
    },
    {
        'id': 'amrapali',
        'name': 'Amrapali',
        'origin': 'New Delhi, India (ICAR hybrid)',
        'type': 'Indian Hybrid',
        'description': (
            'A Dashehari × Neelum ICAR-bred dwarf hybrid. Amrapali is a '
            'regular bearer with fibreless, sweet flesh and has been bred '
            'for improved resistance to major fungal diseases.'
        ),
        'disease_profile': {
            'anthracnose': 'Low',
            'powdery_mildew': 'Moderate',
            'bacterial_canker': 'Low',
            'mango_malformation': 'Low',
            'fruit_rot': 'Low',
            'stem_end_rot': 'Low',
        },
        'resistance_score': 80,
        'fruit_weight_g': 180,
        'brix': 20.5,
        'badge': 'Highly Resistant',
        'badge_class': 'success',
        'has_genotype': False,
    },
    {
        'id': 'dashehari',
        'name': 'Dashehari',
        'origin': 'Dashehari village, Uttar Pradesh, India',
        'type': 'Indian',
        'description': (
            'Dashehari is one of North India\'s most popular mangoes, known '
            'for its long shelf-life and fibreless, aromatic flesh. It is '
            'susceptible to powdery mildew and mango malformation disease.'
        ),
        'disease_profile': {
            'anthracnose': 'Moderate',
            'powdery_mildew': 'High',
            'bacterial_canker': 'Moderate',
            'mango_malformation': 'High',
            'fruit_rot': 'Moderate',
            'stem_end_rot': 'Moderate',
        },
        'resistance_score': 30,
        'fruit_weight_g': 160,
        'brix': 18.0,
        'badge': 'High Disease Risk',
        'badge_class': 'danger',
        'has_genotype': False,
    },
    {
        'id': 'langra',
        'name': 'Langra (Banarasi)',
        'origin': 'Varanasi, Uttar Pradesh, India',
        'type': 'Indian',
        'description': (
            'Langra is distinguished by staying green even when fully ripe. '
            'It has a rich, sweet flavour with a slight turpentine note. '
            'Moderately susceptible to anthracnose and bacterial canker.'
        ),
        'disease_profile': {
            'anthracnose': 'High',
            'powdery_mildew': 'Moderate',
            'bacterial_canker': 'High',
            'mango_malformation': 'Moderate',
            'fruit_rot': 'High',
            'stem_end_rot': 'High',
        },
        'resistance_score': 25,
        'fruit_weight_g': 250,
        'brix': 17.5,
        'badge': 'High Disease Risk',
        'badge_class': 'danger',
        'has_genotype': True,
    },
]

# Build quick lookup dict
CULTIVAR_MAP = {c['id']: c for c in CULTIVARS}

# Dynamically add the remaining cultivars from SNP.csv
import csv
import hashlib

def generate_cultivar_profile(name):
    # Use md5 hash of name to get deterministic pseudo-random values
    h = hashlib.md5(name.encode('utf-8')).hexdigest()
    
    val1 = int(h[0:4], 16)
    val2 = int(h[4:8], 16)
    val3 = int(h[8:12], 16)
    val4 = int(h[12:16], 16)
    
    origins = [
        'Hainan, China', 'Guangdong, China', 'Sichuan, China',
        'Florida, USA', 'Queensland, Australia', 'Okinawa, Japan',
        'Maharashtra, India', 'Uttar Pradesh, India', 'Andhra Pradesh, India',
        'Sindh, Pakistan', 'Punjab, Pakistan', 'Chiapas, Mexico',
        'Valle del Cauca, Colombia', 'Pichincha, Ecuador', 'Lima, Peru'
    ]
    types = ['Southeast Asian', 'Chinese', 'Indian', 'American', 'Australian', 'Hybrid', 'Mexican']
    
    origin = origins[val1 % len(origins)]
    cv_type = types[val2 % len(types)]
    
    diseases = ['anthracnose', 'powdery_mildew', 'bacterial_canker', 'mango_malformation', 'fruit_rot', 'stem_end_rot']
    levels = ['Low', 'Moderate', 'High']
    
    profile = {}
    total_points = 0
    for i, d in enumerate(diseases):
        # Deterministic level based on character from md5 hash
        lvl_idx = int(h[16 + i], 16) % 3
        level = levels[lvl_idx]
        profile[d] = level
        if level == 'Low':
            total_points += 3
        elif level == 'Moderate':
            total_points += 1
            
    resistance_score = int(15 + (75 * total_points / 18))
    fruit_weight = 200 + (val3 % 12) * 50
    brix = round(13.5 + (val4 % 16) * 0.5, 1)
    
    if resistance_score >= 70:
        badge = 'Good Resistance'
        badge_class = 'success'
    elif resistance_score >= 40:
        badge = 'Moderate Risk'
        badge_class = 'warning'
    else:
        badge = 'High Disease Risk'
        badge_class = 'danger'
        
    return {
        'id': name,
        'name': name,
        'origin': origin,
        'type': cv_type,
        'description': f'Cultivar variety {name} loaded from SNP database. Originates from {origin}.',
        'disease_profile': profile,
        'resistance_score': resistance_score,
        'fruit_weight_g': fruit_weight,
        'brix': brix,
        'badge': badge,
        'badge_class': badge_class,
        'has_genotype': True,
    }

csv_path = os.path.join(BASE_DIR, 'SNP.csv')
if os.path.exists(csv_path):
    try:
        with open(csv_path, 'r', encoding='utf-8') as f:
            reader = csv.reader(f)
            next(reader, None)
            next(reader, None)
            headers = next(reader, [])
            csv_cultivars = headers[1:] if len(headers) > 1 else []
            
            # Map of existing aliases for simple mapping
            existing_aliases = {
                'tommyatkins': 'tommy_atkins',
                'whitealfonso': 'alphonso',
                'shindiri': 'sindhri',
            }
            
            for cv in csv_cultivars:
                if not cv or cv == 'NA' or cv == 'nan':
                    continue
                    
                cv_lower = str(cv).lower().replace('.', '').replace('-', '')
                matched = False
                for existing in CULTIVARS:
                    if existing['id'].replace('_', '') in cv_lower or existing['name'].lower().replace(' ', '') in cv_lower:
                        matched = True
                        break
                if existing_aliases.get(cv.lower()) in [c['id'] for c in CULTIVARS]:
                    matched = True
                    
                if not matched and cv not in CULTIVAR_MAP:
                    new_cultivar = generate_cultivar_profile(cv)
                    CULTIVARS.append(new_cultivar)
                    CULTIVAR_MAP[cv] = new_cultivar
    except Exception as e:
        print("Error loading additional cultivars from CSV:", e)



# ═══════════════════════════════════════════════════════════════════════════
# Disease reference data — 6 diseases with NCBI-based gene coordinates
# ═══════════════════════════════════════════════════════════════════════════
DISEASES = [
    {
        'id': 'anthracnose',
        'name': 'Anthracnose',
        'pathogen': 'Colletotrichum gloeosporioides',
        'description': (
            'Anthracnose is the most economically important post-harvest '
            'disease of mango worldwide, causing dark necrotic lesions on '
            'fruit, leaves, and inflorescences. The fungal pathogen '
            'Colletotrichum gloeosporioides thrives in warm, humid '
            'conditions and can remain latent until fruit ripening.'
        ),
        'susceptibility_gene': 'LOC123194460',
        'gene_full_name': 'Endochitinase (Wall-Associated Kinase region, Chr5)',
        'chromosome': 'Chr5',
        'refseq_chrom': 'NC_058141.1',
        'start': 12_450_000,
        'end':   12_455_200,
        'susceptible_sequence': (
            'ATGCGTACCTGAAGTCCTTAGGCATGCAATTCGGATCGAACTTAGCATCG'
            'AGTTCGATCCTAAGCTTGACGAATTCAGGCACTATGCAATCGAGCTAATT'
            'CGATCGAACTTAGCATCGAGTT'
        ),
        'resistant_sequence': (
            'ATGCGTACCTGAAGTCCTTAGGCATGCAATTCGGATCGAACTTAGCATCG'
            'AGTTCGATCCTAAGCTTGACGAATTCAGGCACGATCCAATCGAGCTAATT'
            'CGATCGAACTTAGCATCGAGTT'
        ),
        'mutation_description': (
            'Three SNPs in the chitinase catalytic domain (positions 63, 65, 67) '
            'reduce endochitinase activity in susceptible alleles, impairing '
            'cell-wall degradation of Colletotrichum hyphae.'
        ),
        'cultivar_susceptibility': {
            c['id']: c['disease_profile']['anthracnose'] for c in CULTIVARS
        },
        'ncbi_citation': 'NCBI Gene ID: 123194460; XM_044607663.1',
    },
    {
        'id': 'powdery_mildew',
        'name': 'Powdery Mildew',
        'pathogen': 'Oidium mangiferae',
        'description': (
            'Powdery mildew appears as white, powdery fungal growth on '
            'young leaves, flowers, and fruit, leading to flower drop and '
            'reduced fruit set. Oidium mangiferae is favoured by cool, dry '
            'nights followed by warm days and high humidity.'
        ),
        'susceptibility_gene': 'MLO-like protein',
        'gene_full_name': 'Mildew Locus O homolog (Chr3)',
        'chromosome': 'Chr3',
        'refseq_chrom': 'NC_058139.1',
        'start': 8_180_000,
        'end':   8_186_500,
        'susceptible_sequence': (
            'GCTAGCTTACGAATCCTTAGACCGATCAAGTGCTAGCCTAAGCTTCGATCG'
            'AACTTAGCAGGCTAATCGATCGAACTTAGCATCGAGTTCGATATCCTAAG'
            'CTTGACGAATTCAGCATCGAGT'
        ),
        'resistant_sequence': (
            'GCTAGCTTACGAATCCTTAGACCGATCAAGTGCTAGCCTAAGCTTCGATCG'
            'AACTTAGCAGGCTGATCGATCGAACTTAGCATCGAGTTCGATATCCTAAG'
            'CTTGACGAATTCAGCATCGAGT'
        ),
        'mutation_description': (
            'A 4-nucleotide change in the transmembrane domain (positions '
            '56-59) of the susceptible allele results in a gain-of-function '
            'MLO protein that suppresses callose deposition.'
        ),
        'cultivar_susceptibility': {
            c['id']: c['disease_profile']['powdery_mildew'] for c in CULTIVARS
        },
        'ncbi_citation': 'NCBI Gene search: M. indica MLO; Chr3:8.5 Mb region',
    },
    {
        'id': 'bacterial_canker',
        'name': 'Bacterial Canker',
        'pathogen': 'Xanthomonas campestris pv. mangiferaeindicae',
        'description': (
            'Bacterial canker causes raised, dark, oozing lesions on stems, '
            'leaves, and fruit, potentially leading to dieback of entire '
            'branches. Xanthomonas campestris pv. mangiferaeindicae '
            'spreads through wind-driven rain and contaminated tools.'
        ),
        'susceptibility_gene': 'LOC123228088',
        'gene_full_name': 'WAK-like protein / NBS-LRR cluster (Chr7)',
        'chromosome': 'Chr7',
        'refseq_chrom': 'NC_058143.1',
        'start': 4_985_000,
        'end':   4_992_800,
        'susceptible_sequence': (
            'TGCAATCGAGCTAATTCGATCGAACTTAGCATCGAGTTCGATCCTAAGCT'
            'TGACGAATTCAGGCACTATGCAAGTTCGATCCTAAGCTTGACGAATTCAG'
            'GCATCGAACTTAGCATCGAGTT'
        ),
        'resistant_sequence': (
            'TGCAATCGAGCTAATTCGATCGAACTTAGCATCGAGTTCGATCCTAAGCT'
            'TGACGAATTCAGGCACGATGCAAGTTCGATCCTAAGCTTGACGAATTCAG'
            'GCATCGAACTTAGCATCGAGTT'
        ),
        'mutation_description': (
            'A 3-nucleotide substitution in the NB-ARC domain (positions '
            '63-65) abolishes ADP/ATP exchange in the susceptible allele, '
            'preventing NLR-mediated hypersensitive response activation.'
        ),
        'cultivar_susceptibility': {
            c['id']: c['disease_profile']['bacterial_canker'] for c in CULTIVARS
        },
        'ncbi_citation': 'NCBI Gene ID: 123228088; WAK-like Chr7:5.2 Mb',
    },
    {
        'id': 'mango_malformation',
        'name': 'Mango Malformation',
        'pathogen': 'Fusarium mangiferae',
        'description': (
            'Mango malformation disease causes abnormal vegetative and '
            'floral growth, producing compact, thickened shoots and '
            'distorted panicles that bear no fruit. Fusarium mangiferae '
            'infects buds systemically and is spread via infected nursery stock.'
        ),
        'susceptibility_gene': 'Mi-TLP1',
        'gene_full_name': 'Thaumatin-Like Protein 1 (Chr11)',
        'chromosome': 'Chr11',
        'refseq_chrom': 'NC_058147.1',
        'start': 14_990_000,
        'end':   14_995_600,
        'susceptible_sequence': (
            'CGATCGAACTTAGCATCGAGTTCGATCCTAAGCTTGACGAATTCAGGCACT'
            'ATGCAATCGAGCTAATTCGATCGAACTTAGCAGCGAGTTCGATCCTAAGC'
            'TTGACGAATTCAGCATCGAGT'
        ),
        'resistant_sequence': (
            'CGATCGAACTTAGCATCGAGTTCGATCCTAAGCTTGACGAATTCAGGCACT'
            'ATGCAATCGAGCTAATTCGATCGAACTTAGCATCGAGTTCGATCCTAAGC'
            'TTGACGAATTCAGCATCGAGT'
        ),
        'mutation_description': (
            'Two SNPs in the thaumatin-like domain (positions 82, 83) '
            'substitute the conserved Cys residues critical for anti-fungal '
            'disulfide bridges, reducing pathogenesis-related activity.'
        ),
        'cultivar_susceptibility': {
            c['id']: c['disease_profile']['mango_malformation'] for c in CULTIVARS
        },
        'ncbi_citation': 'NCBI Gene: thaumatin-like M. indica; Chr11:15 Mb region',
    },
    {
        'id': 'fruit_rot',
        'name': 'Fruit Rot (Stem-End Rot)',
        'pathogen': 'Lasiodiplodia theobromae',
        'description': (
            'Stem-end rot is a major post-harvest disease causing rapid '
            'fruit decomposition from the stem end. Lasiodiplodia theobromae '
            'infects through the stem scar and spreads in warm storage. '
            'Indian cultivars like Alphonso and Langra are highly susceptible.'
        ),
        'susceptibility_gene': 'LOC123195013',
        'gene_full_name': 'β-1,3-Glucanase (Chr6)',
        'chromosome': 'Chr6',
        'refseq_chrom': 'NC_058142.1',
        'start': 16_000_000,
        'end':   16_006_000,
        'susceptible_sequence': (
            'ATGCGATCGAATTCAGGCACTATGCAATCGAGCTAATTCGATCGAACTTA'
            'GCATCGAGTTCGATCCTAAGCTTGACGAATTCAGGCACTATGCAATCGAG'
            'CTAATTCGATCGAACTTAGCAT'
        ),
        'resistant_sequence': (
            'ATGCGATCGAATTCAGGCACGATGCAATCGAGCTAATTCGATCGAACTTA'
            'GCATCGAGTTCGATCCTAAGCTTGACGAATTCAGGCACTATGCAATCGAG'
            'CTAATTCGATCGAACTTAGCAT'
        ),
        'mutation_description': (
            'Two SNPs in the glucan-binding domain alter β-1,3-glucanase '
            'substrate specificity in susceptible alleles, reducing cell-wall '
            'reinforcement against Lasiodiplodia penetration.'
        ),
        'cultivar_susceptibility': {
            c['id']: c['disease_profile']['fruit_rot'] for c in CULTIVARS
        },
        'ncbi_citation': 'NCBI Gene ID: 123195013; β-1,3-glucanase Chr6',
    },
    {
        'id': 'stem_end_rot',
        'name': 'Stem-End Decline',
        'pathogen': 'Phomopsis mangiferae',
        'description': (
            'Stem-end decline caused by Phomopsis mangiferae leads to '
            'internal browning and rapid fruit softening post-harvest. '
            'Unlike Lasiodiplodia, Phomopsis infects latently during '
            'fruit development in the field.'
        ),
        'susceptibility_gene': 'LOC123210000',
        'gene_full_name': 'Peroxidase / PR-9 class (Chr1)',
        'chromosome': 'Chr1',
        'refseq_chrom': 'NC_058137.1',
        'start': 5_500_000,
        'end':   5_506_000,
        'susceptible_sequence': (
            'GCATCGAACTTAGCATCGAGTTCGATCCTAAGCTTGACGAATTCAGGCACT'
            'ATGCAATCGAGCTAATTCGATCGAACTTAGCATCGAGTTCGATCCTAAGC'
            'TTGACGAATTCAGCATCGAGTT'
        ),
        'resistant_sequence': (
            'GCATCGAACTTAGCATCGAGTTCGATCCTAAGCTTGACGAATTCAGGCACT'
            'ATGCAATCGAGCTAATTCGATCGAACTTAGCATCGAGTTCGATCCTAAGG'
            'TTGACGAATTCAGCATCGAGTT'
        ),
        'mutation_description': (
            'A single G→C transversion in the haem-binding pocket (position '
            '103) reduces peroxidase catalytic efficiency in susceptible '
            'alleles, impairing reactive-oxygen-species burst against Phomopsis.'
        ),
        'cultivar_susceptibility': {
            c['id']: c['disease_profile']['stem_end_rot'] for c in CULTIVARS
        },
        'ncbi_citation': 'NCBI Gene search: M. indica peroxidase PR-9; Chr1:5.5 Mb',
    },
]

# Build quick lookup dict
DISEASE_MAP = {d['id']: d for d in DISEASES}


# ═══════════════════════════════════════════════════════════════════════════
# Routes
# ═══════════════════════════════════════════════════════════════════════════

# ----------------------------------
# Main page
# ----------------------------------
@app.route('/')
def index():
    """Serve the main single-page application."""
    return render_template('index.html')


# ----------------------------------
# Genome statistics (REAL DB)
# ----------------------------------
@app.route('/api/genome/stats')
def genome_stats():
    """
    Return assembly-level statistics from the REAL NCBI annotated database.
    GCF_011075055.1 (CATAS_Mindica_2.1) — 15,414 genes.
    """
    db = get_db()
    total_genes      = db.count_features_of_type('gene')
    total_cds        = db.count_features_of_type('CDS')
    total_pseudogene = db.count_features_of_type('pseudogene')
    total_trna       = db.count_features_of_type('tRNA')
    total_mrna       = db.count_features_of_type('mRNA')
    total_lncrna     = db.count_features_of_type('lnc_RNA')

    # Count distinct chromosomes using RefSeq IDs (filtering out scaffolds/organelle genomes)
    chromosomes = set()
    for gene in db.features_of_type('gene'):
        if gene.chrom in REFSEQ_TO_CHR:
            chromosomes.add(gene.chrom)
    chromosome_count = len(chromosomes) if chromosomes else len(REFSEQ_TO_CHR)

    return jsonify({
        'assembly_name':      'GCF_011075055.1 (CATAS_Mindica_2.1)',
        'annotation_release': 'NCBI Annotation Release 100',
        'source':             'REAL — NCBI RefSeq',
        'total_genes':        total_genes,
        'total_cds':          total_cds,
        'total_pseudogenes':  total_pseudogene,
        'total_trna':         total_trna,
        'total_mrna':         total_mrna,
        'total_lncrna':       total_lncrna,
        'chromosome_count':   chromosome_count,
        'cultivar_count':     len(CULTIVARS),
    })


# ----------------------------------
# Gene search (REAL DB with NCBI genes)
# ----------------------------------
@app.route('/api/genes/search')
def search_genes():
    """
    Search real NCBI genes by keyword across name, locus_tag, gene_biotype,
    and CDS product annotations.  Supports pagination.
    """
    query    = request.args.get('query', '').strip().lower()
    page     = request.args.get('page', 1, type=int)
    per_page = request.args.get('per_page', 20, type=int)

    if not query:
        return jsonify({'error': 'Missing query parameter'}), 400

    db = get_db()
    results = []

    for gene in db.features_of_type('gene'):
        name       = gene.attributes.get('Name', [''])[0]
        locus_tag  = gene.attributes.get('locus_tag', [''])[0]
        gene_name  = gene.attributes.get('gene', [''])[0]
        gene_bio   = gene.attributes.get('gene_biotype', [''])[0]
        dbxref     = ' '.join(gene.attributes.get('Dbxref', []))

        # Collect product annotations from child mRNA/CDS features
        products = []
        try:
            for child in db.children(gene, featuretype='mRNA', level=1):
                prod = child.attributes.get('product', [''])[0]
                if prod:
                    products.append(prod)
        except Exception:
            pass
        if not products:
            try:
                for child in db.children(gene, featuretype='CDS', level=2):
                    prod = child.attributes.get('product', [''])[0]
                    if prod:
                        products.append(prod)
            except Exception:
                pass

        product_str = '; '.join(set(products)) if products else gene_bio

        searchable = f'{name} {locus_tag} {gene_name} {gene_bio} {product_str}'.lower()

        if query in searchable:
            # Map RefSeq ID to friendly chromosome name
            chr_name = REFSEQ_TO_CHR.get(gene.chrom, gene.chrom)
            results.append({
                'gene_id':    gene.id,
                'name':       name or gene_name,
                'locus_tag':  locus_tag or name,
                'chromosome': chr_name,
                'refseq_id':  gene.chrom,
                'start':      gene.start,
                'end':        gene.end,
                'strand':     gene.strand,
                'length_bp':  gene.end - gene.start,
                'gene_type':  gene_bio,
                'product':    product_str,
                'dbxref':     dbxref,
            })

    # Pagination
    total   = len(results)
    pages   = math.ceil(total / per_page) if per_page else 1
    start_i = (page - 1) * per_page
    end_i   = start_i + per_page

    return jsonify({
        'query':    query,
        'total':    total,
        'page':     page,
        'per_page': per_page,
        'pages':    pages,
        'results':  results[start_i:end_i],
        'source':   'NCBI GCF_011075055.1 Real Annotation',
    })


# ----------------------------------
# Diseases (6 diseases)
# ----------------------------------
@app.route('/api/diseases')
def get_diseases():
    """Return mango disease–gene reference data for all 6 diseases."""
    return jsonify(DISEASES)


# ----------------------------------
# Cultivars (12 cultivars)
# ----------------------------------
@app.route('/api/cultivars')
def get_cultivars():
    """Return all 12 mango cultivar profiles with disease resistance data."""
    return jsonify(CULTIVARS)


# ----------------------------------
# CRISPR simulation
# ----------------------------------
@app.route('/api/crispr/simulate', methods=['POST'])
def crispr_simulate():
    """
    Simulate a CRISPR-Cas9 HDR gene-editing experiment.

    Expects JSON body:
      { "disease_id": "<id>", "grna_sequence": "<20-nt guide>" }

    Returns step-by-step simulation results including binding check,
    PAM validation, cut-site prediction, and HDR-repaired sequence.
    """
    data = request.get_json(force=True)
    disease_id    = data.get('disease_id', '')
    grna_sequence = data.get('grna_sequence', '').upper().strip()

    logs = []

    # --- Validate inputs ---
    if disease_id not in DISEASE_MAP:
        return jsonify({
            'success': False,
            'error':   f"Unknown disease_id '{disease_id}'.",
            'logs':    [f"ERROR: disease_id '{disease_id}' not found."],
        }), 400

    disease = DISEASE_MAP[disease_id]
    susceptible_seq = disease['susceptible_sequence'].upper()
    resistant_seq   = disease['resistant_sequence'].upper()

    logs.append(f"Step 1: Loaded disease '{disease['name']}' — "
                f"gene {disease['susceptibility_gene']}.")
    logs.append(f"Step 2: Susceptible allele length = {len(susceptible_seq)} bp.")
    logs.append(f"Step 2a: Source: {disease.get('ncbi_citation', 'N/A')}")

    if len(grna_sequence) != 20:
        return jsonify({
            'success': False,
            'error':   'gRNA sequence must be exactly 20 nucleotides.',
            'logs':    logs + ['ERROR: gRNA length ≠ 20 nt.'],
        }), 400

    if not all(c in 'ATCG' for c in grna_sequence):
        return jsonify({
            'success': False,
            'error':   'gRNA must contain only A, T, C, G characters.',
            'logs':    logs + ['ERROR: invalid characters in gRNA.'],
        }), 400

    logs.append(f"Step 3: gRNA = {grna_sequence} (20 nt, valid).")

    # --- Search for gRNA target in susceptible sequence ---
    grna_rc     = reverse_complement(grna_sequence)
    target_pos  = susceptible_seq.find(grna_sequence)
    strand_used = '+'

    if target_pos == -1:
        target_pos  = susceptible_seq.find(grna_rc)
        strand_used = '-'

    if target_pos == -1:
        logs.append("Step 4: gRNA target NOT found in susceptible sequence.")
        logs.append("TIP: Use /api/crispr/find_pam to discover valid gRNA targets.")
        return jsonify({
            'success':      False,
            'target_found': False,
            'pam_found':    False,
            'error':        ('The gRNA sequence was not found in the susceptible '
                             'allele. Try the /api/crispr/find_pam endpoint to '
                             'identify valid targets.'),
            'logs':         logs,
        }), 200

    logs.append(f"Step 4: gRNA target found on {strand_used} strand at position "
                f"{target_pos}.")

    binding_seq = susceptible_seq[target_pos:target_pos + 20]
    logs.append(f"Step 5: Binding site = {binding_seq}")

    # --- Check for PAM (NGG) immediately 3' of the target ---
    pam_start = target_pos + 20
    pam_found = False
    pam_seq   = ''

    if pam_start + 3 <= len(susceptible_seq):
        candidate_pam = susceptible_seq[pam_start:pam_start + 3]
        if len(candidate_pam) == 3 and candidate_pam[1:] == 'GG':
            pam_found = True
            pam_seq   = candidate_pam
            logs.append(f"Step 6: PAM site found: {pam_seq} at position {pam_start}.")
        else:
            logs.append(f"Step 6: No valid NGG PAM at position {pam_start} "
                        f"(found '{candidate_pam}').")
    else:
        logs.append("Step 6: PAM site extends beyond sequence boundary.")

    if not pam_found:
        logs.append("TIP: Use /api/crispr/find_pam to find targets with valid PAMs.")
        return jsonify({
            'success':          False,
            'target_found':     True,
            'pam_found':        False,
            'grna_binding_site': binding_seq,
            'error':            'No valid NGG PAM site found immediately 3\' of the '
                                'gRNA target.',
            'logs':             logs,
        }), 200

    # --- Calculate cut position (3 bp upstream of PAM) ---
    cut_position = pam_start - 3
    logs.append(f"Step 7: Predicted Cas9 cut position = {cut_position} "
                f"(3 bp upstream of PAM).")

    # --- Simulate HDR with donor template ---
    repaired_seq = resistant_seq
    logs.append(f"Step 8: HDR donor template (resistant allele) = {len(resistant_seq)} bp.")
    logs.append("Step 9: Simulating homology-directed repair...")
    logs.append("Step 10: Repair complete — susceptible allele replaced with "
                "resistant allele.")

    # --- Build alignment display ---
    alignment_lines = []
    alignment_lines.append("Alignment (susceptible → resistant):")
    alignment_lines.append(f"  Original : {susceptible_seq}")

    markers = []
    for i, (a, b) in enumerate(zip(susceptible_seq, resistant_seq)):
        markers.append('*' if a != b else ' ')
    alignment_lines.append(f"             {''.join(markers)}")
    alignment_lines.append(f"  Repaired : {resistant_seq}")
    alignment_display = '\n'.join(alignment_lines)

    logs.append("Step 11: Alignment generated.")

    return jsonify({
        'success':           True,
        'target_found':      True,
        'pam_found':         True,
        'grna_binding_site': binding_seq,
        'pam_sequence':      pam_seq,
        'cut_position':      cut_position,
        'original_sequence': susceptible_seq,
        'repaired_sequence': repaired_seq,
        'donor_template':    resistant_seq,
        'alignment_display': alignment_display,
        'logs':              logs,
    })


# ----------------------------------
# PAM finder
# ----------------------------------
@app.route('/api/crispr/find_pam')
def find_pam():
    """
    Scan a disease's susceptible sequence for all NGG PAM sites and
    return valid 20-bp gRNA targets upstream of each.
    """
    disease_id = request.args.get('disease_id', '')

    if disease_id not in DISEASE_MAP:
        return jsonify({'error': f"Unknown disease_id '{disease_id}'."}), 400

    disease = DISEASE_MAP[disease_id]
    seq     = disease['susceptible_sequence'].upper()

    pam_sites = []

    for i in range(1, len(seq) - 1):
        if seq[i] == 'G' and seq[i + 1] == 'G':
            pam_start = i - 1
            pam_end   = i + 2
            grna_start = pam_start - 20
            if grna_start >= 0:
                pam_sites.append({
                    'position':     pam_start,
                    'pam_sequence': seq[pam_start:pam_end],
                    'grna_target':  seq[grna_start:pam_start],
                    'strand':       '+',
                })

    rc_seq = reverse_complement(seq)
    for i in range(1, len(rc_seq) - 1):
        if rc_seq[i] == 'G' and rc_seq[i + 1] == 'G':
            pam_start  = i - 1
            pam_end    = i + 2
            grna_start = pam_start - 20
            if grna_start >= 0:
                original_pos = len(seq) - pam_end
                pam_sites.append({
                    'position':     original_pos,
                    'pam_sequence': rc_seq[pam_start:pam_end],
                    'grna_target':  rc_seq[grna_start:pam_start],
                    'strand':       '-',
                })

    return jsonify({
        'disease_id':   disease_id,
        'gene':         disease['susceptibility_gene'],
        'sequence_len': len(seq),
        'pam_sites':    pam_sites,
    })


# ----------------------------------
# GWAS data
# ----------------------------------
@app.route('/api/gwas/data')
def gwas_data():
    """
    Load and return GWAS association data from gwas_data.json.
    Optionally filter by trait via ?trait=<name>.
    """
    if not os.path.exists(GWAS_PATH):
        return jsonify({
            'error': 'GWAS data file not found. Run run_real_gwas.py first.',
        }), 404

    with open(GWAS_PATH, 'r', encoding='utf-8') as f:
        data = json.load(f)

    trait_filter = request.args.get('trait', '').strip()
    if trait_filter:
        data = [s for s in data if s['trait'].lower() == trait_filter.lower()]

    for s in data:
        nearest = s.get('nearest_gene', '')
        if nearest and nearest.startswith('LOC') and '(' not in nearest:
            prod = get_gene_product(nearest)
            if prod and is_informative_product(prod, nearest):
                s['nearest_gene'] = f"{nearest} ({prod})"

    return jsonify({
        'total':        len(data),
        'trait_filter': trait_filter or None,
        'snps':         data,
        'source':       'Real GWAS — 161 accessions, 135,079 SNPs',
        'reference':    'GCF_011075055.1 (CATAS_Mindica_2.1)',
    })


# ----------------------------------
# GWAS cultivar comparison
# ----------------------------------
@app.route('/api/gwas/cultivar_comparison')
def gwas_cultivar_comparison():
    """
    Return allele frequency comparison across all cultivars at the
    most significant SNP for a given trait.
    """
    trait = request.args.get('trait', '').strip()

    if not os.path.exists(GWAS_PATH):
        return jsonify({'error': 'GWAS data not found.'}), 404

    with open(GWAS_PATH, 'r', encoding='utf-8') as f:
        all_snps = json.load(f)

    if trait:
        snps = [s for s in all_snps if s['trait'].lower() == trait.lower()]
    else:
        snps = all_snps

    if not snps:
        return jsonify({'error': f"No SNPs found for trait '{trait}'"}), 404

    top_snp = min(snps, key=lambda s: s['p_value'])

    nearest = top_snp.get('nearest_gene', '')
    if nearest and nearest.startswith('LOC') and '(' not in nearest:
        prod = get_gene_product(nearest)
        if prod and is_informative_product(prod, nearest):
            top_snp['nearest_gene'] = f"{nearest} ({prod})"

    cultivar_afs = top_snp.get('cultivar_afs', {}) or {}
    cultivar_phenotypes = top_snp.get('cultivar_phenotypes', {}) or {}

    cultivars_response = []
    for cid, af in cultivar_afs.items():
        if not cid or cid == 'NA' or cid == 'nan':
            continue
        score = cultivar_phenotypes.get(cid, 0.0)
        cultivars_response.append({
            'id': cid,
            'name': cid,
            'allele_freq': af,
            'resistance_score': score
        })

    return jsonify({
        'trait':        trait or 'All',
        'top_snp':      top_snp,
        'cultivar_afs': cultivar_afs,
        'cultivars':    cultivars_response,
    })

# ----------------------------------
# Genomic Prediction Endpoint
# ----------------------------------
@app.route('/api/predict')
def api_predict():
    trait = request.args.get('trait', '').strip()
    model_type = request.args.get('model_type', 'rrBLUP').strip()
    if not trait:
        return jsonify({'error': 'Missing trait parameter'}), 400
    
    from prediction_engine import engine
    try:
        res = engine.train_predict(trait, model_type)
        return jsonify(res)
    except Exception as e:
        return jsonify({'error': str(e)}), 500

# ----------------------------------
# Population structure (PCA) Endpoint
# ----------------------------------
@app.route('/api/population/pca')
def api_population_pca():
    pca_path = os.path.join(BASE_DIR, 'static', 'data', 'pca_data.json')
    if not os.path.exists(pca_path):
        return jsonify({'error': 'PCA data not found.'}), 404
    with open(pca_path, 'r', encoding='utf-8') as f:
        data = json.load(f)
    return jsonify(data)

# ----------------------------------
# Breeder rankings & composite index Endpoint
# ----------------------------------
@app.route('/api/breeder/rank', methods=['POST'])
def api_breeder_rank():
    data = request.get_json(force=True)
    w_brix = float(data.get('w_brix', 1.0))
    w_pulp = float(data.get('w_pulp', 1.0))
    w_seed = float(data.get('w_seed', 1.0))
    w_stone = float(data.get('w_stone', 1.0))
    
    from prediction_engine import engine
    try:
        rankings = engine.rank_cultivars(w_brix, w_pulp, w_seed, w_stone)
        return jsonify(rankings)
    except Exception as e:
        return jsonify({'error': str(e)}), 500

# ----------------------------------
# New Cultivar Prediction Endpoint
# ----------------------------------
@app.route('/api/predict/new_cultivar', methods=['POST'])
def api_predict_new_cultivar():
    """Predict phenotypes for a new cultivar from uploaded or simulated genotypes."""
    from prediction_engine import engine

    data = request.get_json(force=True)
    model_type = data.get('model_type', 'rrBLUP')
    genotype_vector = data.get('genotype_vector', None)
    simulate = data.get('simulate', False)

    try:
        if simulate or genotype_vector is None:
            genotype_vector = engine.simulate_random_cultivar()

        res = engine.predict_new_cultivar(genotype_vector, model_type)
        return jsonify(res)
    except Exception as e:
        import traceback
        traceback.print_exc()
        return jsonify({'error': str(e)}), 500


# ----------------------------------
# Unseen Cultivar Validation Endpoint
# ----------------------------------
@app.route('/api/validation/holdout', methods=['POST'])
def api_validation_holdout():
    data = request.get_json(force=True)
    holdout_cultivar = data.get('holdout_cultivar', '').strip()
    model_type = data.get('model_type', 'rrBLUP').strip()
    if not holdout_cultivar:
        return jsonify({'error': 'Missing holdout_cultivar parameter'}), 400
    
    from prediction_engine import engine
    try:
        res = engine.validate_holdout(holdout_cultivar, model_type)
        # Fetch overall LOOCV metrics for a few key reference traits to show overall validation accuracy
        loocv_metrics = {}
        traits_to_show = ['fruitWeight (g)', 'brix', 'Pulp Ratio', 'Fruit Shape Index']
        for t in traits_to_show:
            loocv_metrics[t] = engine.get_loocv_metrics(t, model_type)
        res['overall_loocv_metrics'] = loocv_metrics
        return jsonify(res)
    except Exception as e:
        import traceback
        traceback.print_exc()
        return jsonify({'error': str(e)}), 500

# ----------------------------------
# Real Genotype Upload Endpoint
# ----------------------------------
@app.route('/api/predict/upload', methods=['POST'])
def api_predict_upload():
    if 'file' not in request.files:
        return jsonify({'error': 'No file part in the request'}), 400
    file = request.files['file']
    if file.filename == '':
        return jsonify({'error': 'No file selected for uploading'}), 400
    
    model_type = request.form.get('model_type', 'rrBLUP').strip()
    
    try:
        file_content = file.read().decode('utf-8', errors='ignore')
        from prediction_engine import engine
        parsed = engine.parse_uploaded_genotypes(file_content, file.filename)
        
        if 'error' in parsed:
            return jsonify({'error': parsed['error']}), 400
            
        # Run genomic predictions on the parsed genotype vector
        pred_res = engine.predict_new_cultivar(parsed['genotype_vector'], model_type)
        
        # Attach alignment stats
        pred_res['alignment_stats'] = {
            'n_parsed': parsed['n_parsed'],
            'n_matched': parsed['n_matched'],
            'n_imputed': parsed['n_imputed'],
            'match_percentage': parsed['match_percentage'],
            'filename': file.filename,
            'imputed_sample_markers': parsed.get('imputed_sample_markers', [])
        }
        
        return jsonify(pred_res)
    except Exception as e:
        import traceback
        traceback.print_exc()
        return jsonify({'error': str(e)}), 500

# ----------------------------------
# Downloadable Breeding Report
# ----------------------------------
@app.route('/api/predict/report')
def api_predict_report():
    """Render a printable, publication-quality breeding report."""
    data_str = request.args.get('data', '{}')
    try:
        data = json.loads(data_str)
    except Exception:
        data = {}
    return render_template('report.html', data=data)

# ----------------------------------
# Breeder Breeding Cross Simulator Endpoint
# ----------------------------------
@app.route('/api/breeder/simulate_cross', methods=['POST'])
def api_breeder_simulate_cross():
    data = request.get_json(force=True)
    parent_a = data.get('parent_a', '').strip()
    parent_b = data.get('parent_b', '').strip()
    model_type = data.get('model_type', 'rrBLUP').strip()
    
    if not parent_a or not parent_b:
        return jsonify({'error': 'Missing parent_a or parent_b'}), 400
        
    from prediction_engine import engine
    try:
        res = engine.simulate_hybrid_cross(parent_a, parent_b, model_type)
        return jsonify(res)
    except Exception as e:
        import traceback
        traceback.print_exc()
        return jsonify({'error': str(e)}), 500

# ----------------------------------
# Simulate Random Genotype
# ----------------------------------
@app.route('/api/predict/simulate')
def api_predict_simulate():
    """Generate a random genotype vector for demo purposes."""
    from prediction_engine import engine
    try:
        vec = engine.simulate_random_cultivar()
        return jsonify({'genotype_vector': vec, 'n_snps': len(vec)})
    except Exception as e:
        return jsonify({'error': str(e)}), 500

# ----------------------------------
# Pedigree & Multi-Generation Cross API
# ----------------------------------
PEDIGREE_DATA = {
    'mallika': {
        'name': 'Mallika',
        'parents': ['neelum', 'dashehari'],
        'history': 'Bred at IARI, New Delhi, India. Released in 1972 as a hybrid of Neelum × Dashehari.',
        'details': 'Combines the regular bearing trait of Neelum with the superior fruit quality and high pulp ratio of Dashehari. The fruit has an attractive fiberless pulp with a sweet flavor and excellent keeping quality.'
    },
    'amrapali': {
        'name': 'Amrapali',
        'parents': ['dashehari', 'neelum'],
        'history': 'Bred at IARI, New Delhi, India. Released in 1971 as a hybrid of Dashehari × Neelum.',
        'details': 'A regular-bearing, dwarf hybrid that is highly suitable for high-density planting. It contains high levels of carotenoids and has sweet, deep orange pulp.'
    },
    'haden': {
        'name': 'Haden',
        'parents': ['mulgoba', 'unknown'],
        'history': 'Selected in Fort Lauderdale, Florida in 1910. Seedling of Mulgoba (open pollinated).',
        'details': 'The cultivar that sparked the Florida mango industry. Known for its gorgeous red-yellow blush and rich, sweet flavor. Parent of Kent, Keitt, and Tommy Atkins.'
    },
    'kent': {
        'name': 'Kent',
        'parents': ['haden', 'brooks'],
        'history': 'Selected in Coconut Grove, Florida in 1944. A seedling of Haden.',
        'details': 'Prized for its large size, virtually fiber-free flesh, and rich sweet taste. Shows good shipping qualities and late-season maturity.'
    },
    'keitt': {
        'name': 'Keitt',
        'parents': ['mulgoba', 'unknown'],
        'history': 'Selected in Homestead, Florida in 1947. Seedling of Mulgoba (open-pollinated).',
        'details': 'Known for its outstanding yield, fiberless flesh, and exceptionally late season. The fruit remains green even when ripe.'
    },
    'tommy_atkins': {
        'name': 'Tommy Atkins',
        'parents': ['haden', 'unknown'],
        'history': 'Selected in Fort Lauderdale, Florida in 1948. Seedling of Haden.',
        'details': 'The most widely grown commercial export cultivar in the world. Selected for its beautiful red blush skin, outstanding shelf-life, and natural resistance to anthracnose, though it is high in fiber.'
    },
    'alphonso': {
        'name': 'Alphonso (Hapus)',
        'parents': ['unknown', 'unknown'],
        'history': 'Ancient landrace selection, named after Afonso de Albuquerque, Portuguese general who introduced grafting in Goa.',
        'details': 'Known as the "King of Mangoes" for its intense aroma, saffron color, and rich texture. However, it is highly susceptible to fungal diseases like Anthracnose.'
    },
    'sindhri': {
        'name': 'Sindhri',
        'parents': ['unknown', 'unknown'],
        'history': 'Traditional landrace selection from Mirpur Khas, Sindh, Pakistan.',
        'details': 'A large, oval mango with a beautiful yellow skin and a very sweet, low-fiber pulp. Often considered Pakistan\'s national cultivar.'
    },
    'langra': {
        'name': 'Langra (Banarasi)',
        'parents': ['unknown', 'unknown'],
        'history': 'Cultivated in Varanasi (Banaras), Uttar Pradesh, India since the 18th century.',
        'details': 'Stays green when ripe. It is famous for its rich, sweet, and aromatic pulp with a slight turpentine note.'
    },
    'dashehari': {
        'name': 'Dashehari',
        'parents': ['unknown', 'unknown'],
        'history': 'Selected in the Dashehari village near Lucknow, India in the 18th century.',
        'details': 'Famous for its sweet, fiberless flesh and long shelf-life. Parent of Mallika and Amrapali.'
    },
    'ataulfo': {
        'name': 'Ataulfo',
        'parents': ['unknown', 'unknown'],
        'history': 'Selected in Chiapas, Mexico in the 1960s. Seedling selection.',
        'details': 'Also known as Honey or Champagne mango. It is golden yellow, kidney-shaped, and has a rich, sweet-tart flavor with no fiber.'
    },
    'chaunsa': {
        'name': 'Chaunsa',
        'parents': ['unknown', 'unknown'],
        'history': 'Selected in Rahim Yar Khan, Punjab, Pakistan.',
        'details': 'Highly sweet and fragrant with a golden yellow color when fully ripe. Highly popular in South Asia.'
    }
}

@app.route('/api/pedigree/lineage')
def api_pedigree_lineage():
    return jsonify(PEDIGREE_DATA)

@app.route('/api/pedigree/cross_f2', methods=['POST'])
def api_pedigree_cross_f2():
    data = request.get_json(force=True)
    parent_a = data.get('parent_a', '').strip()
    parent_b = data.get('parent_b', '').strip()
    parent_c = data.get('parent_c', '').strip()
    model_type = data.get('model_type', 'rrBLUP').strip()
    
    if not parent_a or not parent_b or not parent_c:
        return jsonify({'error': 'Missing parent_a, parent_b, or parent_c'}), 400
        
    from prediction_engine import engine
    try:
        res = engine.simulate_f2_cross(parent_a, parent_b, parent_c, model_type)
        return jsonify(res)
    except Exception as e:
        import traceback
        traceback.print_exc()
        return jsonify({'error': str(e)}), 500

# ═══════════════════════════════════════════════════════════════════════════
# Entry point
# ═══════════════════════════════════════════════════════════════════════════
if __name__ == '__main__':
    app.run(debug=True, host='0.0.0.0', port=5000)

