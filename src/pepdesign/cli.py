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

    targets = sub.add_parser("targets", help="retrieve and cache the peptide set from RCSB")
    targets.add_argument("--search-limit", type=int, default=500)
    targets.add_argument(
        "--hotspots",
        action="store_true",
        help="also define target hotspots from a predicted complex; needs a GPU",
    )

    generate = sub.add_parser("generate", help="backbone generation and sequence design")
    generate.add_argument("--target", required=True)
    generate.add_argument("--generator", default="rfdiffusion")
    generate.add_argument(
        "--scorer",
        default="boltz2",
        help="must differ from the generator; enforced, not merely advised",
    )
    generate.add_argument("--n-designs", type=int, default=1000)

    controls = sub.add_parser(
        "controls", help="write the scrambled and length-matched control sets"
    )
    controls.add_argument("--search-limit", type=int, default=500)
    controls.add_argument("--max-peptides", type=int, default=None, help=MAX_PEPTIDES_HELP)
    controls.add_argument("--seed", type=int, default=0)

    sub.add_parser("evaluate", help="print the separation table from an existing run")

    analysis = sub.add_parser(
        "analysis", help="peptides, controls, sequence-level scoring and separation"
    )
    analysis.add_argument("--search-limit", type=int, default=500)
    analysis.add_argument("--max-peptides", type=int, default=None, help=MAX_PEPTIDES_HELP)
    analysis.add_argument("--seed", type=int, default=0)
    analysis.add_argument(
        "--force",
        action="store_true",
        help="overwrite findings even if they came from a larger peptide set",
    )

    return parser


#: A default of "all of them" rather than a number, because a default that silently scored
#: part of the cached set is how the published run came to cover 120 of 190 peptides.
MAX_PEPTIDES_HELP = "how many cached peptides to score; the default is all of them"

#: Subcommands that genuinely need hardware this project has never had.
GPU_GATED = {
    "generate": "RFdiffusion / ProteinMPNN generation",
}

#: The one step inside an otherwise CPU-only subcommand that needs the same hardware.
HOTSPOTS_NEED_GPU = (
    "hotspot definition needs a predicted complex, which needs a GPU this project has not "
    "had access to. Drop --hotspots to retrieve and cache the peptide set, which is pure "
    "network and CPU and is the step 'pepdesign analysis' depends on."
)


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        return _dispatch(args)
    except RuntimeError as exc:
        raise SystemExit(str(exc)) from exc
    except OSError as exc:
        raise SystemExit(f"could not reach RCSB or the model hub: {exc}") from exc


def _load_peptides(args: argparse.Namespace):
    from pepdesign.targets import build

    peptides, dropped = build(args.data_dir / "peptides.json", limit=args.search_limit)
    if not peptides:
        raise SystemExit(
            "no usable peptides were retrieved from RCSB. "
            "Check network access, or raise --search-limit."
        )
    if args.max_peptides is None:
        return peptides, dropped
    return peptides[: args.max_peptides], dropped


def _dispatch(args: argparse.Namespace) -> int:
    if args.command in GPU_GATED:
        raise SystemExit(
            f"'{args.command}' is not implemented: {GPU_GATED[args.command]} needs a GPU "
            "this project has not had access to.\n"
            "Run 'pepdesign analysis' for the sequence-level experiment, which does not."
        )

    if args.command == "targets":
        # Retrieval is network plus CPU and is the code that built data/peptides.json.
        # Only the hotspot step needs a GPU, so only the hotspot step is refused.
        if args.hotspots:
            raise SystemExit(HOTSPOTS_NEED_GPU)
        from pepdesign.targets import build

        path = args.data_dir / "peptides.json"
        cached = path.exists()
        peptides, dropped = build(path, limit=args.search_limit)
        source = "cache" if cached else "RCSB"
        print(f"{len(peptides)} peptides from {source} -> {path} (dropped {dropped})")
        if cached:
            print(
                "This is the committed set the published AUCs came from. Delete the file to "
                "re-retrieve, which returns a different set because the PDB grows."
            )
        return 0

    if args.command == "controls":
        # Implemented and CPU-only. This used to be refused with a GPU message.
        from pepdesign.controls import CONTROL_KINDS, composition_distance, make_controls

        peptides, _ = _load_peptides(args)
        pairs = [(p.entity_id, p.sequence) for p in peptides]
        controls = make_controls(pairs, kinds=CONTROL_KINDS, seed=args.seed)
        scrambled = [c.sequence for c in controls if c.kind == "scrambled"]
        distance = composition_distance([s for _, s in pairs], scrambled)

        args.results_dir.mkdir(parents=True, exist_ok=True)
        out = args.results_dir / "controls.json"
        import json as _json

        out.write_text(
            _json.dumps(
                {
                    "n_peptides": len(pairs),
                    "composition_distance_real_vs_scrambled": round(distance, 8),
                    "controls": [
                        {
                            "id": c.control_id,
                            "kind": c.kind,
                            "sequence": c.sequence,
                            "derived_from": c.derived_from,
                        }
                        for c in controls
                    ],
                },
                indent=1,
            )
        )
        print(f"{len(controls)} controls for {len(pairs)} peptides -> {out}")
        print(f"composition distance real vs scrambled: {distance:.6f} (must be 0 by design)")
        return 0

    if args.command == "evaluate":
        # Implemented and CPU-only: reads what analysis already produced.
        import json as _json

        path = args.results_dir / "findings.json"
        if not path.exists():
            raise SystemExit(f"no findings at {path}. Run 'pepdesign analysis' first.")
        findings = _json.loads(path.read_text())
        print(f"filter: {findings['filter']}")
        # The threshold aims at 5% of controls and a finite sample lands beside it, so the
        # admitted fraction is printed next to the recall it buys rather than assumed.
        print(
            f"{'control':<22} {'AUC':>7} {'95% CI':>18} {'d':>7} {'ctrl in':>8} {'recall':>8}"
        )
        for sep in findings["separations"]:
            ci = f"[{sep['auc_ci_low']:.3f}, {sep['auc_ci_high']:.3f}]"
            admitted = sep.get("control_admitted_at_threshold")
            admitted_text = "?" if admitted is None else f"{admitted:.3f}"
            print(
                f"{sep['control_kind']:<22} {sep['auc']:>7.3f} {ci:>18} "
                f"{sep['effect_size']:>7.2f} {admitted_text:>8} "
                f"{sep['real_recall_at_threshold']:>8.3f}"
            )
        print(f"\nnot run: {findings['not_run']}")
        return 0

    if args.command == "analysis":
        # A 4-peptide smoke run silently replacing a 120-peptide result is how a published
        # number quietly becomes wrong. Refuse unless asked twice.
        import json as _json

        existing = args.results_dir / "findings.json"
        if existing.exists() and not args.force:
            previous = _json.loads(existing.read_text())
            previous_n = previous.get("populations", {}).get("real", {}).get("n", 0)
            limit = args.max_peptides
            if limit is not None and previous_n > limit:
                raise SystemExit(
                    f"{existing} holds a run over {previous_n} peptides and this run would "
                    f"use {limit}. Refusing to overwrite a larger result.\n"
                    "Pass --force, or raise --max-peptides."
                )

        from pepdesign.controls import CONTROL_KINDS, composition_distance, make_controls
        from pepdesign.evaluate import build_findings, write, write_scores
        from pepdesign.score import EsmScorer
        from pepdesign.targets import SNAPSHOT_RELEASE_DATE

        peptides, dropped = _load_peptides(args)
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
                "peptide_set": str(args.data_dir / "peptides.json"),
                "peptide_release_cutoff": SNAPSHOT_RELEASE_DATE,
                "n_peptides_scored": len(pairs),
                "dropped": dropped,
                "composition_distance_real_vs_scrambled": round(distance, 8),
            },
        )
        write(findings, args.results_dir / "findings.json")
        scores_path = write_scores(scored, args.results_dir / "scores.csv")
        print(f"per-sequence scores -> {scores_path}")
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
