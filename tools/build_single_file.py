#!/usr/bin/env python3
"""Inline css/ and js/ into one self-contained HTML file: dist/protein_cad.html.

    python tools/build_single_file.py

The result has no external dependencies, so it can be copied to another machine and opened by
double-clicking.
"""
import pathlib
import re

root = pathlib.Path(__file__).resolve().parent.parent
html = (root / 'index.html').read_text(encoding='utf-8')

html = re.sub(r'<link rel="stylesheet" href="([^"]+)">',
              lambda m: '<style>\n' + (root / m.group(1)).read_text(encoding='utf-8') + '\n</style>', html)
html = re.sub(r'<script src="([^"]+)"></script>',
              lambda m: '<script>\n' + (root / m.group(1)).read_text(encoding='utf-8').replace('</script', '<\\/script') + '\n</script>', html)

out = root / 'dist' / 'protein_cad.html'
out.parent.mkdir(exist_ok=True)
out.write_text(html, encoding='utf-8')
print('wrote %s (%.0f KB)' % (out, out.stat().st_size / 1024))
