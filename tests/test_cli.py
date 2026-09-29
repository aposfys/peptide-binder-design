"""The CLI surface: what is reachable, what is refused, and what cannot be clobbered."""

from __future__ import annotations

import json

import pytest

from pepdesign.cli import GPU_GATED, build_parser, main


@pytest.mark.parametrize("command", sorted(GPU_GATED))
def test_gpu_gated_commands_say_why(command):
    with pytest.raises(SystemExit) as excinfo:
        main([command, "--target", "X"])
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


def test_targets_is_reachable_because_retrieval_needs_no_gpu():
    """Retrieval is network plus CPU and is the code that built data/peptides.json."""
    assert "targets" not in GPU_GATED


def test_targets_reads_the_committed_set_and_says_so(tmp_path, capsys):
    data = tmp_path / "peptides.json"
    data.write_text(
        json.dumps(
            {
                "dropped": {},
                "peptides": [
                    {"entity_id": "1ABC_1", "sequence": "ACDEFGHIK", "description": "x"}
                ],
            }
        )
    )
    assert main(["--data-dir", str(tmp_path), "targets"]) == 0
    out = capsys.readouterr().out
    assert "1 peptides from cache" in out
    assert "published AUCs came from" in out


def test_only_the_hotspot_step_of_targets_is_gpu_gated(tmp_path):
    with pytest.raises(SystemExit) as excinfo:
        main(["--data-dir", str(tmp_path), "targets", "--hotspots"])
    message = str(excinfo.value)
    assert "GPU" in message
    assert "Drop --hotspots" in message


def test_max_peptides_defaults_to_the_whole_cached_set():
    """A numeric default is how the published run came to cover 120 of 190 peptides."""
    parser = build_parser()
    assert parser.parse_args(["analysis"]).max_peptides is None
    assert parser.parse_args(["controls"]).max_peptides is None
