#!/usr/bin/env python3
"""Serve the app folder over http://localhost (no caching), for development and the test pages.

    python tools/serve.py [--port 8770]

Opening index.html straight from disk works too; the server is only needed for tests/ (which
fetch a PDB file) and to avoid stale caches while editing.
"""
import argparse
import functools
import http.server
import pathlib

root = pathlib.Path(__file__).resolve().parent.parent


class NoCache(http.server.SimpleHTTPRequestHandler):
    def end_headers(self):
        self.send_header('Cache-Control', 'no-store')
        super().end_headers()


ap = argparse.ArgumentParser()
ap.add_argument('--port', type=int, default=8770)
args = ap.parse_args()
handler = functools.partial(NoCache, directory=str(root))
with http.server.ThreadingHTTPServer(('127.0.0.1', args.port), handler) as httpd:
    print('serving %s at http://127.0.0.1:%d/' % (root, args.port))
    httpd.serve_forever()
