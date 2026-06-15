import sys
sys.path.insert(0, r"C:\mangoproject\libs")
import gffutils, pandas as pd, os

DB_FILE  = "mango_genes.db"
GFF_FILE = "mango.gff"

if not os.path.exists(DB_FILE):
    print("Building gene database, please wait ~1 min...")
    db = gffutils.create_db(
        GFF_FILE, dbfn=DB_FILE, force=True,
        keep_order=True, merge_strategy="merge",
        sort_attribute_values=True
    )
    print("Database ready!")
else:
    db = gffutils.FeatureDB(DB_FILE)

QUERY = input("Enter gene name or keyword to search: ")
print(f"\nSearching for '{QUERY}' in mango genome...")

records = []
for gene in db.features_of_type("gene"):
    note  = gene.attributes.get("Note",      [""])[0]
    name  = gene.attributes.get("Name",      [""])[0]
    locus = gene.attributes.get("locus_tag", [""])[0]
    records.append({
        "gene_id":    gene.id,
        "name":       name,
        "locus_tag":  locus,
        "chromosome": gene.chrom,
        "start":      gene.start,
        "end":        gene.end,
        "strand":     gene.strand,
        "length_bp":  gene.end - gene.start,
        "function":   note,
    })

df = pd.DataFrame(records)

mask = df.apply(
    lambda row: QUERY.lower() in (
        row["function"] + " " +
        row["name"]     + " " +
        row["locus_tag"]
    ).lower(),
    axis=1
)
results = df[mask].reset_index(drop=True)

if results.empty:
    print(f"\nNo genes found matching '{QUERY}'")
    print("Tips:")
    print("  Try shorter keyword  e.g. 'chitin' instead of 'chitinase precursor'")
    print("  Try gene family name e.g. 'LRR' or 'WRKY' or 'kinase'")
else:
    print(f"\n{'='*70}")
    print(f"  Found {len(results)} gene(s) matching '{QUERY}'")
    print(f"{'='*70}")
    print(f"\n{'#':<5} {'Gene ID':<20} {'Chr':<12} {'Start':>12} {'End':>12} {'Strand':<8} {'Length(bp)'}")
    print("-" * 80)
    for i, row in results.iterrows():
        print(f"{i+1:<5} {row['gene_id']:<20} {row['chromosome']:<12} {row['start']:>12} {row['end']:>12} {row['strand']:<8} {row['length_bp']}")

    print("\nFUNCTIONS:")
    for i, row in results.iterrows():
        print(f"  [{i+1}] {row['function']}")

    out_file = f"search_{QUERY.replace(' ','_')}.csv"
    results.to_csv(out_file, index=False)
    print(f"\nSaved -> {out_file}")

    if len(results) > 1:
        print("\nCHROMOSOME SUMMARY:")
        print(results.groupby("chromosome")[["gene_id"]].count().rename(
            columns={"gene_id": "count"}
        ).to_string())