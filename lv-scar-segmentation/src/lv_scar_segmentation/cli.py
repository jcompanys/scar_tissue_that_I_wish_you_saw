"""Small CLI for dataset structure scanning."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from .data_loading import scan_dataset


def main() -> int:
    parser = argparse.ArgumentParser(description="Scan LV scar dataset folder structure")
    parser.add_argument(
        "--dataset-root",
        type=Path,
        required=True,
        help="Path like F:/RM_TEKNON_DEVELOP",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=None,
        help="Optional output JSON report path",
    )
    args = parser.parse_args()

    report = scan_dataset(args.dataset_root)
    payload = report.to_dict()

    if args.out is not None:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    else:
        print(json.dumps(payload, indent=2))

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
