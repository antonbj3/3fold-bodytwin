from pathlib import Path
P = Path(__file__).resolve().parent
text = (P / 'MANUSCRIPT_TEMPLATE.md').read_text()
for (token, name) in {'TABLE1': 'table1_flow', 'TABLE2': 'table2_arch_profile', 'TABLE3A': 'table3a_opposition', 'TABLE3B': 'table3b_corrected_crossbite', 'TABLE4': 'table4_report_conflicts'}.items():
    text = text.replace('[[' + token + ']]', (P / 'tables' / (name + '.md')).read_text().strip())
assert '[[' not in text
(P / 'MANUSCRIPT.md').write_text(text)
print('MANUSCRIPT.md assembled from locked evidence tables.')
