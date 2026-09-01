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

    analysis = sub.add_parser(
        "analysis", help="peptides, controls, sequence-level scoring and separation"
    )
    analysis.add_argument("--search-limit", type=int, default=500)
    analysis.add_argument("--max-peptides", type=int, default=120)
    analysis.add_argument("--seed", type=int, default=0)

    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)

    if args.command == "analysis":
        from pepdesign.controls import CONTROL_KINDS, composition_distance, make_controls
        from pepdesign.evaluate import build_findings, write
        from pepdesign.score import EsmScorer
        from pepdesign.targets import build

        peptides, dropped = build(args.data_dir / "peptides.json", limit=args.search_limit)
        peptides = peptides[: args.max_peptides]
        pairs = [(p.entity_id, p.sequence) for p in peptides]
        print(f"{len(pairs)} peptides (dropped {dropped})")

        controls = make_controls(pairs, kinds=CONTROL_KINDS, seed=args.seed)
        scrambled = [c.sequence for c in controls if c.kind == "scrambled"]
        distance = composition_distance([s for _, s in pairs], scrambled)
        print(
            f"{len(controls)} controls; composition distance real vs scrambled {distance:.6f}"
        )

        items = [(pid, seq, "real") for pid, seq in pairs]
        items += [(c.control_id, c.sequence, c.kind) for c in controls]
        scored = EsmScorer().score_many(items)

        findings = build_findings(
            scored,
            notes={
                "peptide_source": "RCSB short chains in multi-protein structures",
                "dropped": dropped,
                "composition_distance_real_vs_scrambled": round(distance, 8),
            },
        )
        write(findings, args.results_dir / "findings.json")
        for separation in findings["separations"]:
            print(
                f"  vs {separation['control_kind']:<20} AUC {separation['auc']:.3f} "
                f"[{separation['auc_ci_low']:.3f}, {separation['auc_ci_high']:.3f}]"
            )
        return 0

    raise SystemExit(
        f"'{args.command}' is not implemented: generation and structure scoring need a GPU. "
        "Run 'analysis' for the sequence-level experiment that does not."
    )


if __name__ == "__main__":
    raise SystemExit(main())
