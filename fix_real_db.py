import sqlite3
import json

# Fix the real DB with the proper gffutils meta format
db_path = r'C:\mangoproject\mango_genes_real.db'
conn = sqlite3.connect(db_path)
cur = conn.cursor()

# Clear existing bad meta
cur.execute("DELETE FROM meta")

# Insert correct dialect JSON (matching gffutils 0.14 format for GFF3)
dialect = {
    "leading semicolon": False,
    "trailing semicolon": False,
    "quoted GFF2 values": False,
    "field separator": ";",
    "semicolon in quotes": False,
    "keyval separator": "=",
    "multival separator": ",",
    "fmt": "gff3",
    "repeated keys": False,
    "order": ["ID", "Dbxref", "Name", "gbkey", "gene", "gene_biotype",
              "gene_synonym", "locus_tag", "Parent", "product",
              "protein_id", "transcript_id", "cultivar", "country"]
}

dialect_json = json.dumps(dialect)
cur.execute("INSERT INTO meta VALUES (?, ?)", (dialect_json, "0.14"))
conn.commit()

cur.execute("SELECT * FROM meta")
print("Meta set:", cur.fetchall())

conn.close()

# Now test with gffutils
import sys
sys.path.insert(0, r'C:\mangoproject\libs')
import gffutils

db = gffutils.FeatureDB(db_path)
print(f"Genes: {db.count_features_of_type('gene')}")
print(f"CDS: {db.count_features_of_type('CDS')}")
print(f"Pseudogenes: {db.count_features_of_type('pseudogene')}")
print(f"tRNA: {db.count_features_of_type('tRNA')}")

gene = next(db.features_of_type('gene'))
print(f"Sample gene: {gene.id} | chrom={gene.chrom} | start={gene.start} | end={gene.end}")
print(f"  Attributes: {dict(gene.attributes)}")
print("SUCCESS!")
