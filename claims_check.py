#!/usr/bin/env python3
"""Claims and language check for the manuscript and README.
Fails if the manuscript or README contain over-claiming phrases, if the eponym is used in
manuscript prose more than once, or if the required structural elements are missing.
Needs only the standard library. Usage: python3 claims_check.py"""
import re, sys, zipfile, os
from xml.etree import ElementTree as ET
W = '{http://schemas.openxmlformats.org/wordprocessingml/2006/main}'
HERE = os.path.dirname(os.path.abspath(__file__))
BANNED = [r'toward a vaccine', r'toward a cure', r'a cure and lower mortality', r'lower mortality from', r'reduced mortality:',
          r'saves? lives', r'\bEbola cure\b', r'will cure', r'patient benefit is', r'clinically proven']
fails = []
def check_text(label, text):
    for pat in BANNED:
        for m in re.finditer(pat, text, flags=re.I):
            fails.append(f'{label}: banned phrase "{m.group(0)}" near "{text[max(0, m.start() - 40):m.end() + 40]!r}"')
# ---- manuscript
docx = os.path.join(HERE, 'Rein_Persistence_Theory_EBOV_Nyachio_Final.docx')
root = ET.fromstring(zipfile.ZipFile(docx).read('word/document.xml'))
body = root.find(W + 'body')
prose, alltext = [], []
for el in body:
    if el.tag == W + 'p':
        t = ''.join(x.text or '' for x in el.iter(W + 't')); prose.append(t); alltext.append(t)
    elif el.tag == W + 'tbl':
        alltext.append(''.join(x.text or '' for x in el.iter(W + 't')))
full = '\n'.join(alltext)
check_text('manuscript', full)
eponym = [t for t in prose if 'Rein Persistence Theory' in t and not t.lstrip().startswith(('//', '#'))]
if len(eponym) > 1: fails.append(f'manuscript: eponym "Rein Persistence Theory" used in {len(eponym)} prose paragraphs (limit 1)')
for need in ['Table 3a', 'Table 11a', 'Table 14a', 'Table 21a', 'pre-registered', 'F6.', 'Maasai Mara University']:
    if need not in full: fails.append(f'manuscript: required element missing: {need}')
# ---- markdown notes
for path in [os.path.join(HERE, 'README.md')]:
    check_text(os.path.relpath(path, HERE), open(path, encoding='utf8').read())
if fails:
    print('\n'.join(fails)); print(f'claims check: FAIL ({len(fails)})'); sys.exit(1)
print('claims check: PASS (no banned phrases; eponym limit respected; required elements present)')
