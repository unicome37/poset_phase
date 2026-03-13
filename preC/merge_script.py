import re

files = [
    'MANUSCRIPT_Section1_Abstract_Introduction.md',
    'MANUSCRIPT_Section2_Methods.md',
    'MANUSCRIPT_Section3_Results.md',
    'MANUSCRIPT_Section4_SimpsonsParadox.md',
    'MANUSCRIPT_Section5_ComponentDecomposition.md',
    'MANUSCRIPT_Section6_Discussion.md',
]
merged = []
for i, f in enumerate(files):
    with open(f, 'r', encoding='utf-8') as fh:
        content = fh.read()
    content = re.sub(r'^# Prediction C.*?Manuscript Draft.*?\n', '', content, count=1)
    if i == 0:
        content = content.replace(
            '## Hierarchy Depth Observables Predict Combinatorial Entropy',
            '# Hierarchy Depth Observables Predict Combinatorial Entropy'
        )
    merged.append(content.strip())

sep = '\n\n---\n\n'
full = sep.join(merged)

with open('MANUSCRIPT_PredictionC_Full.md', 'w', encoding='utf-8') as fh:
    fh.write(full)

lines = full.count('\n') + 1
words = len(full.split())
tables = full.count('**Table ')
refs_remaining = full.count('[ref]')
print(f'Lines: {lines}, Words: {words}, Tables: {tables}, [ref]: {refs_remaining}')
