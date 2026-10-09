"""Check that edits to Python files are documentation-only.

Compares each file in the working tree against the same file at a git ref,
after stripping docstrings, comments and type annotations. If the remaining
code is identical, the edit cannot change behaviour.

Usage:
    python tools/check_ast_equiv.py [--ref REF] FILE [FILE ...]

REF defaults to the tag ``thesis-final-2026-10``. Exit code 0 means every
file is equivalent, 1 means at least one file differs.
"""

from __future__ import annotations

import argparse
import ast
import difflib
import subprocess
import sys
from pathlib import Path


class _StripDocs(ast.NodeTransformer):
    """Remove docstrings and annotations, keeping everything that runs."""

    def _strip_docstring(self, node):
        body = node.body
        if (
            body
            and isinstance(body[0], ast.Expr)
            and isinstance(body[0].value, ast.Constant)
            and isinstance(body[0].value.value, str)
        ):
            node.body = body[1:] or [ast.Pass()]
        return node

    def visit_Module(self, node):
        self.generic_visit(node)
        return self._strip_docstring(node)

    def visit_ClassDef(self, node):
        # Class-body annotations are kept: in dataclasses they define fields.
        self.generic_visit(node)
        return self._strip_docstring(node)

    def _visit_func(self, node):
        self.generic_visit(node)
        node.returns = None
        args = node.args
        for a in args.posonlyargs + args.args + args.kwonlyargs:
            a.annotation = None
        if args.vararg:
            args.vararg.annotation = None
        if args.kwarg:
            args.kwarg.annotation = None
        return self._strip_docstring(node)

    visit_FunctionDef = _visit_func
    visit_AsyncFunctionDef = _visit_func


def _normalise(source: str, filename: str) -> str:
    tree = _StripDocs().visit(ast.parse(source, filename=filename))
    ast.fix_missing_locations(tree)
    return ast.unparse(tree)


def _git_show(ref: str, path: Path) -> str:
    repo_root = subprocess.run(
        ["git", "rev-parse", "--show-toplevel"],
        capture_output=True, text=True, check=True,
    ).stdout.strip()
    rel = path.resolve().relative_to(Path(repo_root).resolve()).as_posix()
    return subprocess.run(
        ["git", "show", f"{ref}:{rel}"],
        capture_output=True, text=True, check=True, encoding="utf-8",
    ).stdout


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("files", nargs="+", type=Path)
    parser.add_argument("--ref", default="thesis-final-2026-10")
    args = parser.parse_args()

    all_ok = True
    for path in args.files:
        old = _normalise(_git_show(args.ref, path), f"{args.ref}:{path}")
        new = _normalise(path.read_text(encoding="utf-8"), str(path))
        if old == new:
            print(f"OK    {path}  (documentation-only vs {args.ref})")
            continue
        all_ok = False
        print(f"DIFF  {path}  (code differs from {args.ref})")
        sys.stdout.writelines(
            difflib.unified_diff(
                old.splitlines(keepends=True), new.splitlines(keepends=True),
                fromfile=f"{args.ref}:{path}", tofile=str(path), n=2,
            )
        )
    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main())
