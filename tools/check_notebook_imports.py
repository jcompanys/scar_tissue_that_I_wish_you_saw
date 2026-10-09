"""Check that every name notebooks import from ``src`` still exists.

Parses the code cells of each notebook, collects ``from src.X import a, b``
and ``from X import a`` (for notebooks that put ``src/`` itself on
``sys.path``), imports the module from the working tree and reports names
that are missing. Catches ImportErrors caused by removing or renaming
functions without running the notebooks. Needs no dataset.

Usage:
    python tools/check_notebook_imports.py [NOTEBOOK ...]

Without arguments every ``*.ipynb`` under ``lv-scar-segmentation/notebooks``
is checked, except ``old/``. Exit code 0 means all names resolve.
"""

from __future__ import annotations

import ast
import importlib
import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
PROJECT = REPO / "lv-scar-segmentation"
SRC_MODULES = {p.stem for p in (PROJECT / "src").glob("*.py")}


def _imports(nb_path: Path):
    """Yield (module, name, cell_index) for src imports in a notebook."""
    nb = json.loads(nb_path.read_text(encoding="utf-8"))
    for i, cell in enumerate(nb["cells"]):
        if cell["cell_type"] != "code":
            continue
        # Drop IPython magics so the cell parses as Python.
        lines = [l for l in "".join(cell["source"]).splitlines()
                 if not l.lstrip().startswith(("%", "!"))]
        try:
            tree = ast.parse("\n".join(lines))
        except SyntaxError:
            continue
        for node in ast.walk(tree):
            if not isinstance(node, ast.ImportFrom) or node.module is None:
                continue
            mod = node.module
            if mod.startswith("src."):
                mod = mod[4:]
            if mod not in SRC_MODULES:
                continue
            for alias in node.names:
                yield mod, alias.name, i


def main() -> int:
    sys.path.insert(0, str(PROJECT))
    nbs = [Path(a) for a in sys.argv[1:]] or [
        p for p in sorted((PROJECT / "notebooks").rglob("*.ipynb"))
        if "old" not in p.relative_to(PROJECT / "notebooks").parts
        and ".ipynb_checkpoints" not in p.parts
    ]
    missing = 0
    for nb in nbs:
        for mod, name, cell in _imports(nb):
            module = importlib.import_module(f"src.{mod}")
            if name != "*" and not hasattr(module, name):
                print(f"MISSING  {nb.name} cell {cell}: src.{mod}.{name}")
                missing += 1
    print(f"{len(nbs)} notebooks checked, {missing} missing names")
    return 1 if missing else 0


if __name__ == "__main__":
    sys.exit(main())
