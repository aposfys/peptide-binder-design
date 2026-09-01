"""The CLI surface: what is reachable, what is refused, and what cannot be clobbered."""

from __future__ import annotations

import json

import pytest

from pepdesign.cli import GPU_GATED, build_parser, main


@pytest.mark.parametrize("command", sorted(GPU_GATED))
def test_gpu_gated_commands_say_why(command):
    with pytest.raises(SystemExit) as excinfo:
        main([command, "--target", "X"] if command != "controls" else [command])
    message = str(excinfo.value)
    assert "GPU" in message
    assert "analysis" in message


def test_controls_is_reachable_not_gpu_gated():
    """controls is implemented and CPU-only; it used to be refused with a GPU message."""
    assert "controls" not in GPU_GATED


def test_evaluate_without_findings_says_what_to_run(tmp_path):
    with pytest.raises(SystemExit) as excinfo:
        main(["--results-dir", str(tmp_path), "evaluate"])
    assert "pepdesign analysis" in str(excinfo.value)


def test_a_small_run_cannot_overwrite_a_larger_one(tmp_path):
    """A 4-peptide smoke run silently replacing a 120-peptide result is how a published
    number quietly becomes wrong. This happened during development."""
    (tmp_path / "findings.json").write_text(json.dumps({"populations": {"real": {"n": 120}}}))
    with pytest.raises(SystemExit) as excinfo:
        main(["--results-dir", str(tmp_path), "analysis", "--max-peptides", "4"])
    message = str(excinfo.value)
    assert "Refusing to overwrite" in message
    assert "--force" in message


def test_subcommands_parse():
    parser = build_parser()
    assert parser.parse_args(["analysis"]).command == "analysis"
    assert parser.parse_args(["controls"]).command == "controls"
