"""Command line entry point: ``python -m pepdesign.cli`` or ``pepdesign``."""

from __future__ import annotations

import argparse
from pathlib import Path

from pepdesign import __version__


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="pepdesign",
        description="De novo peptide binder design with non-circular evaluation",
    )
    parser.add_argument("--version", action="version", version=f"pepdesign {__version__}")
    parser.add_argument("--data-dir", type=Path, default=Path("data"))
    parser.add_argument("--results-dir", type=Path, default=Path("results"))
    sub = parser.add_subparsers(dest="command", required=True)

    targets = sub.add_parser("targets", help="prepare targets, hotspots and held-out binders")
    targets.add_argument("--target", required=True)

    generate = sub.add_parser("generate", help="backbone generation and sequence design")
    generate.add_argument("--target", required=True)
    generate.add_argument("--generator", default="rfdiffusion")
    generate.add_argument(
        "--scorer",
        default="boltz2",
        help="must differ from the generator; enforced, not merely advised",
    )
    generate.add_argument("--n-designs", type=int, default=1000)

    sub.add_parser("controls", help="scrambled and length-matched controls")
    sub.add_parser("evaluate", help="separation of designs, controls and known binders")

    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    raise SystemExit(f"'{args.command}' is not implemented yet; see README milestones")


if __name__ == "__main__":
    raise SystemExit(main())
