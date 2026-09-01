"""Known experimental peptide binders, from the PDB.

**What this set is, precisely:** short protein chains (8-30 residues, standard amino acids
only) that appear in a solved structure alongside at least one other protein entity. Those
are peptides observed bound to a protein partner, which is the closest thing to a
"validated binder" that can be assembled without a wet lab.

**What it is not:** a curated binder set. Some short chains in multi-protein structures are
subunits of a complex rather than ligands, and a few are crystallisation tags. The
population is a proxy, it is noisy in a known direction, and every result here should be
read with that in mind rather than after it has been forgotten.
"""

from __future__ import annotations

import json
import time
import urllib.parse
import urllib.request
from dataclasses import dataclass
from pathlib import Path

RCSB_SEARCH = "https://search.rcsb.org/rcsbsearch/v2/query"
RCSB_GRAPHQL = "https://data.rcsb.org/graphql"

STANDARD_AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")

MIN_LENGTH = 8
MAX_LENGTH = 30


@dataclass(frozen=True)
class Peptide:
    """One peptide chain observed bound in a complex."""

    entity_id: str
    sequence: str
    description: str

    @property
    def length(self) -> int:
        return len(self.sequence)


def _post(url: str, payload: dict, attempts: int = 4) -> dict:
    body = json.dumps(payload).encode("utf-8")
    last: Exception | None = None
    for attempt in range(attempts):
        try:
            request = urllib.request.Request(
                url, data=body, headers={"Content-Type": "application/json"}
            )
            with urllib.request.urlopen(request, timeout=120) as response:
                return json.loads(response.read().decode("utf-8"))
        except Exception as exc:
            last = exc
            time.sleep(2**attempt)
    raise RuntimeError(f"RCSB search failed after {attempts} attempts") from last


def _get(url: str, attempts: int = 4) -> dict:
    last: Exception | None = None
    for attempt in range(attempts):
        try:
            with urllib.request.urlopen(url, timeout=120) as response:
                return json.loads(response.read().decode("utf-8"))
        except Exception as exc:
            last = exc
            time.sleep(2**attempt)
    raise RuntimeError(f"RCSB GraphQL failed after {attempts} attempts") from last


def search_entity_ids(*, limit: int = 400) -> list[str]:
    """Short protein entities that share their structure with another protein entity."""
    payload = {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "entity_poly.rcsb_entity_polymer_type",
                        "operator": "exact_match",
                        "value": "Protein",
                    },
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_polymer_entity.formula_weight",
                        "operator": "less",
                        "value": 3.5,
                    },
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_entry_info.polymer_entity_count_protein",
                        "operator": "greater",
                        "value": 1,
                    },
                },
            ],
        },
        "return_type": "polymer_entity",
        "request_options": {
            "paginate": {"start": 0, "rows": limit},
            "results_content_type": ["experimental"],
        },
    }
    response = _post(RCSB_SEARCH, payload)
    return [row["identifier"] for row in response.get("result_set", [])]


def fetch_sequences(entity_ids: list[str], *, batch: int = 50) -> list[Peptide]:
    """Sequences and descriptions for a list of polymer entities."""
    peptides: list[Peptide] = []
    for start in range(0, len(entity_ids), batch):
        chunk = entity_ids[start : start + batch]
        ids = ",".join(f'"{value}"' for value in chunk)
        query = f"""{{
          polymer_entities(entity_ids:[{ids}]) {{
            rcsb_id
            entity_poly {{ pdbx_seq_one_letter_code_can }}
            rcsb_polymer_entity {{ pdbx_description }}
          }}
        }}"""
        response = _get(f"{RCSB_GRAPHQL}?{urllib.parse.urlencode({'query': query})}")
        for entity in response.get("data", {}).get("polymer_entities") or []:
            if entity is None:
                continue
            sequence = (entity["entity_poly"] or {}).get("pdbx_seq_one_letter_code_can") or ""
            sequence = sequence.strip().upper()
            description = (entity["rcsb_polymer_entity"] or {}).get("pdbx_description") or ""
            peptides.append(
                Peptide(
                    entity_id=entity["rcsb_id"], sequence=sequence, description=description
                )
            )
    return peptides


def usable(peptides: list[Peptide]) -> tuple[list[Peptide], dict[str, int]]:
    """Keep peptides in range and free of non-standard residues, counting what was dropped.

    Non-standard residues arrive as ``X``. They are dropped rather than substituted,
    because substituting them invents sequence that was never observed -- and a language
    model scoring an invented residue would report a confidence about nothing.
    """
    kept: list[Peptide] = []
    dropped = {"too_short": 0, "too_long": 0, "non_standard": 0, "duplicate": 0}
    seen: set[str] = set()

    for peptide in peptides:
        if len(peptide.sequence) < MIN_LENGTH:
            dropped["too_short"] += 1
        elif len(peptide.sequence) > MAX_LENGTH:
            dropped["too_long"] += 1
        elif set(peptide.sequence) - STANDARD_AMINO_ACIDS:
            dropped["non_standard"] += 1
        elif peptide.sequence in seen:
            # The PDB is full of the same peptide solved many times. Keeping duplicates
            # would weight those sequences and shrink the effective sample.
            dropped["duplicate"] += 1
        else:
            seen.add(peptide.sequence)
            kept.append(peptide)
    return kept, dropped


def build(out_path: Path, *, limit: int = 400) -> tuple[list[Peptide], dict[str, int]]:
    """Fetch and filter the peptide set, cached."""
    if out_path.exists():
        stored = json.loads(out_path.read_text())
        return [Peptide(**row) for row in stored["peptides"]], stored["dropped"]

    entity_ids = search_entity_ids(limit=limit)
    peptides, dropped = usable(fetch_sequences(entity_ids))

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(
        json.dumps(
            {
                "source": "RCSB PDB: protein entities < 3.5 kDa in structures with >1 protein entity",
                "requested": limit,
                "kept": len(peptides),
                "dropped": dropped,
                "peptides": [
                    {
                        "entity_id": p.entity_id,
                        "sequence": p.sequence,
                        "description": p.description,
                    }
                    for p in peptides
                ],
            },
            indent=1,
        )
    )
    return peptides, dropped
