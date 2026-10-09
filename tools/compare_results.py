"""Compare a results folder against a reference results folder.

Used to verify that refactoring did not change any notebook output.

Usage:
    python tools/compare_results.py REFERENCE_DIR NEW_DIR [--rtol R] [--atol A]

File types:
    .csv / .tsv   loaded with pandas; same shape and columns, numeric columns
                  equal within tolerance (NaN == NaN), other columns exact.
    .png / .jpg   pixel comparison; reports max abs difference and the
                  fraction of differing pixels.
    other         byte-identical (SHA-256). PDFs and SVGs often embed
                  timestamps, so a hash mismatch there is reported as a warning.

Exit code 0 when everything matches, 1 otherwise.
"""

from __future__ import annotations

import argparse
import hashlib
import sys
from pathlib import Path

import numpy as np
import pandas as pd

TABLE_EXT = {".csv", ".tsv"}
IMAGE_EXT = {".png", ".jpg", ".jpeg"}
SOFT_HASH_EXT = {".pdf", ".svg"}


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _compare_table(ref: Path, new: Path, rtol: float, atol: float) -> str | None:
    sep = "\t" if ref.suffix == ".tsv" else ","
    a = pd.read_csv(ref, sep=sep)
    b = pd.read_csv(new, sep=sep)
    if a.shape != b.shape:
        return f"shape {a.shape} -> {b.shape}"
    if list(a.columns) != list(b.columns):
        return f"columns differ: {sorted(set(a.columns) ^ set(b.columns))}"
    problems = []
    for col in a.columns:
        x, y = a[col], b[col]
        if pd.api.types.is_numeric_dtype(x) and pd.api.types.is_numeric_dtype(y):
            xv, yv = x.to_numpy(float), y.to_numpy(float)
            close = np.isclose(xv, yv, rtol=rtol, atol=atol, equal_nan=True)
            if not close.all():
                diff = np.nanmax(np.abs(xv[~close] - yv[~close]))
                problems.append(f"{col}: {(~close).sum()} values, max |diff| {diff:.3g}")
        elif not x.astype(str).equals(y.astype(str)):
            problems.append(f"{col}: {(x.astype(str) != y.astype(str)).sum()} values")
    return "; ".join(problems) or None


def _compare_image(ref: Path, new: Path) -> str | None:
    import matplotlib.image as mpimg

    a = np.asarray(mpimg.imread(ref), dtype=float)
    b = np.asarray(mpimg.imread(new), dtype=float)
    if a.shape != b.shape:
        return f"size {a.shape} -> {b.shape}"
    diff = np.abs(a - b)
    if diff.max() == 0:
        return None
    pixel_diff = diff.reshape(diff.shape[0], diff.shape[1], -1).max(axis=2) > 0
    return f"max |diff| {diff.max():.3g}, {pixel_diff.mean():.2%} pixels differ"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("reference", type=Path)
    parser.add_argument("new", type=Path)
    parser.add_argument("--rtol", type=float, default=1e-9)
    parser.add_argument("--atol", type=float, default=1e-12)
    args = parser.parse_args()

    ref_files = {p.relative_to(args.reference) for p in args.reference.rglob("*") if p.is_file()}
    new_files = {p.relative_to(args.new) for p in args.new.rglob("*") if p.is_file()}

    n_ok = n_fail = n_warn = 0
    for rel in sorted(ref_files - new_files):
        print(f"MISSING  {rel}")
        n_fail += 1
    for rel in sorted(new_files - ref_files):
        print(f"EXTRA    {rel}")
        n_warn += 1

    for rel in sorted(ref_files & new_files):
        ref, new = args.reference / rel, args.new / rel
        ext = rel.suffix.lower()
        try:
            if ext in TABLE_EXT:
                problem = _compare_table(ref, new, args.rtol, args.atol)
            elif ext in IMAGE_EXT:
                problem = _compare_image(ref, new)
            else:
                problem = None if _sha256(ref) == _sha256(new) else "bytes differ"
        except Exception as exc:  # unreadable file: report, keep going
            problem = f"could not compare ({exc})"

        if problem is None:
            n_ok += 1
        elif ext in SOFT_HASH_EXT:
            print(f"WARN     {rel}: {problem} (may be embedded timestamp)")
            n_warn += 1
        else:
            print(f"DIFF     {rel}: {problem}")
            n_fail += 1

    print(f"\n{n_ok} identical, {n_fail} different/missing, {n_warn} warnings")
    return 0 if n_fail == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
