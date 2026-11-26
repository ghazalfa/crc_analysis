#!/usr/bin/env python3
"""FETCH (setup): download one file listed in data/sources.csv and verify its md5.

The source is an http(s) URL or a path relative to the project directory.
Pass --md5 NA to skip the check; the computed md5 is then printed so it can be pinned.
"""
import argparse
import hashlib
import os
import shutil
import sys
import time
import urllib.request

TIMEOUT_SECONDS = 300
CHUNK_BYTES = 1 << 20


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--name", required=True, help="filename to save as")
    p.add_argument("--url", required=True, help="http(s) URL, or a path relative to --project-dir")
    p.add_argument("--md5", required=True, help="expected md5, or NA to skip the check")
    p.add_argument("--project-dir", required=True)
    return p.parse_args()


def download(url, destination, attempts=3):
    for attempt in range(1, attempts + 1):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "crc-mutation-clusters"})
            with urllib.request.urlopen(req, timeout=TIMEOUT_SECONDS) as r, open(destination, "wb") as out:
                shutil.copyfileobj(r, out, CHUNK_BYTES)
            return
        except Exception as err:  # network hiccups: retry with a short back-off
            if attempt == attempts:
                raise
            print(f"[fetch] attempt {attempt} failed ({err}); retrying", file=sys.stderr)
            time.sleep(5 * attempt)


def md5sum(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(CHUNK_BYTES), b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    a = parse_args()
    if a.url.startswith(("http://", "https://")):
        download(a.url, a.name)
    else:
        shutil.copy(os.path.join(a.project_dir, a.url), a.name)

    got = md5sum(a.name)
    if a.md5 == "NA":
        print(f"[fetch] {a.name}: md5 not pinned, computed {got}", file=sys.stderr)
    elif got != a.md5:
        sys.exit(f"[fetch] {a.name}: md5 mismatch, expected {a.md5}, got {got}")
    else:
        print(f"[fetch] {a.name}: md5 verified", file=sys.stderr)


if __name__ == "__main__":
    main()
