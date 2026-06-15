import pypdf
import sys
import re

sys.stdout.reconfigure(encoding='utf-8')

reader = pypdf.PdfReader(r'c:\mangoproject\research_paper.pdf')
print(f"Total pages: {len(reader.pages)}")

full_text = ""
for idx, page in enumerate(reader.pages):
    full_text += page.extract_text() or ""

# Let's look for traits
print("--- Abstract & Intro Traits ---")
matches = re.findall(r'(?:trait|phenotype|property|characteristic|quality|morphological|disease|resistance|weight|length|width|thickness|brix|sugar)s?', full_text, re.IGNORECASE)
print(f"Found {len(matches)} matches of keywords.")

# Let's search for lines containing 'Table' or list of phenotypes
lines = full_text.split('\n')
print("--- Table/Trait Related Lines ---")
count = 0
for line in lines:
    if any(kw in line.lower() for kw in ['fruit weight', 'stone weight', 'seed weight', 'anthracnose', 'brix', 'powdery', 'canker', 'malformation', 'rot', 'decline', 'acidity', 'carotenoid', 'shelf life', 'flowering']):
        print(line[:120])
        count += 1
        if count > 100:
            break
