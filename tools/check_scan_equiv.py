"""Check that the cleaned ``data_loading.scan`` resolves the same paths as the original.

Loads ``data_loading.py`` twice: once from a git ref (default: the tag
``thesis-final-2026-10``) and once from the working tree. Both scan the
dataset and every field the new ``PatientCase`` still has is compared case by
case. Keys or fields removed by the cleanup are ignored; keys added would be
reported.

This replaces re-running every notebook after a ``data_loading`` refactor:
the notebooks only see what ``scan`` returns.

Usage:
    python tools/check_scan_equiv.py [--ref REF] [--root DATASET_ROOT]

Exit code 0 means identical, 1 means at least one difference.
"""

from __future__ import annotations

import argparse
import dataclasses
import importlib.util
import subprocess
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
MODULE = "lv-scar-segmentation/src/data_loading.py"

# PatientCase fields renamed by the cleanup: new name -> name in the ref.
RENAMED_FIELDS = {"lv_anatomy": "anatomy"}


def _load(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod  # dataclasses need the module registered
    spec.loader.exec_module(mod)
    return mod


def _load_from_ref(ref: str, tmp: Path):
    src = subprocess.run(
        ["git", "-C", str(REPO), "show", f"{ref}:{MODULE}"],
        capture_output=True, text=True, check=True, encoding="utf-8",
    ).stdout
    path = tmp / "data_loading_ref.py"
    path.write_text(src, encoding="utf-8")
    return _load(path, "data_loading_ref")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref", default="thesis-final-2026-10")
    ap.add_argument("--root", default=r"F:/RM_TEKNON_DEVELOP")
    args = ap.parse_args()

    with tempfile.TemporaryDirectory() as tmp:
        old = _load_from_ref(args.ref, Path(tmp))
        new = _load(REPO / MODULE, "data_loading_new")

        problems = []

        old_variants = getattr(old, "DE_MRI_VARIANTS", getattr(old, "_DE_MRI_VARIANTS", None))
        if tuple(new.DE_MRI_VARIANTS) != tuple(old_variants):
            problems.append(f"DE_MRI_VARIANTS: {old_variants!r} != {new.DE_MRI_VARIANTS!r}")

        print(f"Scanning with {args.ref} ...")
        old_cases = old.scan(args.root)
        print("Scanning with working tree ...")
        new_cases = new.scan(args.root)

        if len(old_cases) != len(new_cases):
            problems.append(f"case count: {len(old_cases)} != {len(new_cases)}")

        fields = [f.name for f in dataclasses.fields(new.PatientCase)]
        for o, n in zip(old_cases, new_cases):
            tag = f"{n.year}/{n.patient_id}"
            for name in fields:
                old_name = RENAMED_FIELDS.get(name, name)
                if not hasattr(o, old_name):
                    problems.append(f"{tag}: field '{name}' is new")
                    continue
                ov, nv = getattr(o, old_name), getattr(n, name)
                if isinstance(nv, dict):
                    extra = set(nv) - set(ov)
                    if extra:
                        problems.append(f"{tag}.{name}: new keys {sorted(extra)}")
                    for key in nv.keys() & ov.keys():
                        if ov[key] != nv[key]:
                            problems.append(f"{tag}.{name}[{key!r}]: {ov[key]} != {nv[key]}")
                elif ov != nv:
                    problems.append(f"{tag}.{name}: {ov} != {nv}")

    print(f"{len(new_cases)} cases, fields compared: {', '.join(fields)}")
    if problems:
        print(f"DIFFERENT ({len(problems)} problems):")
        for p in problems[:50]:
            print("  " + p)
        return 1
    print("IDENTICAL")
    return 0


if __name__ == "__main__":
    sys.exit(main())
