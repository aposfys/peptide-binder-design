"""What happens when the scorer's backends are not installed."""

from __future__ import annotations

import sys

import pytest

from pepdesign.score import import_backends


def test_missing_backends_name_the_extra_to_install(monkeypatch):
    """The default install pulls in nothing, so this is the common first failure."""
    monkeypatch.setitem(sys.modules, "torch", None)
    monkeypatch.setitem(sys.modules, "transformers", None)
    with pytest.raises(SystemExit) as excinfo:
        import_backends()
    message = str(excinfo.value)
    assert "[score]" in message
    assert "controls" in message and "evaluate" in message
