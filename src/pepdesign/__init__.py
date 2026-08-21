"""De novo peptide binder design, with the evaluation taken more seriously than the design.

The package separates *generation* from *scoring* at the module level and refuses to let
them share a model, because a pipeline whose generator and scorer are the same network
measures agreement rather than binding.
"""

__version__ = "0.1.0"

#: Models used to propose backbones and sequences.
GENERATOR_MODELS = frozenset({"rfdiffusion", "rfpeptides", "proteinmpnn"})

#: Models used to score proposals. Deliberately disjoint from GENERATOR_MODELS; the test
#: suite asserts the disjointness so the constraint cannot erode over time.
SCORER_MODELS = frozenset({"boltz2", "af2-multimer", "chai1"})

__all__ = ["GENERATOR_MODELS", "SCORER_MODELS", "__version__"]
