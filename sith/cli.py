from __future__ import annotations

import argparse
import sys

from . import __version__
from .pipeline import run


def main() -> None:
    parser = argparse.ArgumentParser(
        prog="sith",
        description="SITH 2.0 (S-PLUS Integrated Tool for Homogeneous Photometry)",
    )
    parser.add_argument("field", help="S-PLUS field, e.g. SGR107")
    parser.add_argument("config", help="JEDIS.conf")
    parser.add_argument("--check", action="store_true", help="validate inputs/dependencies and stop")
    parser.add_argument("--force", action="store_true", help="recreate existing products")
    parser.add_argument("--version", action="version", version=f"SITH {__version__}")
    args = parser.parse_args()
    try:
        run(args.field, args.config, force=args.force, check_only=args.check)
    except KeyboardInterrupt:
        raise SystemExit(130)
    except Exception as exc:
        print(f"\nSITH 2.0 ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
