import sqlite3
db = sqlite3.connect(r'C:\mangoproject\mango_genes_real.db')
cur = db.cursor()
cur.execute("SELECT name FROM sqlite_master WHERE type='table'")
tables = cur.fetchall()
print('Tables:', tables)
for t in tables:
    cur.execute(f"SELECT COUNT(*) FROM {t[0]}")
    print(f'  {t[0]}: {cur.fetchone()[0]} rows')
    cur.execute(f"PRAGMA table_info({t[0]})")
    cols = cur.fetchall()
    print(f'    Columns: {[c[1] for c in cols]}')
    cur.execute(f"SELECT * FROM {t[0]} LIMIT 2")
    rows = cur.fetchall()
    for row in rows:
        print(f'    Sample: {row}')
