"""The peptide query, which has to be pinned or the numbers computed over it cannot be checked."""

from __future__ import annotations

import json

from pepdesign import targets


def test_query_is_bounded_by_a_release_date():
    payload = targets.search_payload(limit=10, released_on_or_before="2026-09-01")
    text = json.dumps(payload)
    assert "rcsb_accession_info.initial_release_date" in text
    assert "2026-09-01" in text
    assert payload["request_options"]["sort"], "an unsorted query returns a moving set"


def test_identifiers_come_back_sorted(monkeypatch):
    monkeypatch.setattr(
        targets,
        "_post",
        lambda url, payload: {
            "result_set": [{"identifier": i} for i in ["9ZZZ_1", "1AAA_2", "5MMM_1"]]
        },
    )
    assert targets.search_entity_ids(limit=3) == ["1AAA_2", "5MMM_1", "9ZZZ_1"]


def test_the_committed_set_is_preferred_over_the_network(tmp_path, monkeypatch):
    """The cache is the published set, so it must win over a fresh retrieval."""

    def explode(*args, **kwargs):
        raise AssertionError("build must not hit the network when the cache exists")

    monkeypatch.setattr(targets, "search_entity_ids", explode)
    path = tmp_path / "peptides.json"
    path.write_text(
        json.dumps(
            {
                "dropped": {"too_short": 1},
                "peptides": [
                    {"entity_id": "1ABC_1", "sequence": "ACDEFGHIK", "description": "x"}
                ],
            }
        )
    )
    peptides, dropped = targets.build(path)
    assert [p.entity_id for p in peptides] == ["1ABC_1"]
    assert dropped == {"too_short": 1}
