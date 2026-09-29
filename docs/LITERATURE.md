# Where this sits in the literature

Numbers in brackets point at [REFERENCES.md](REFERENCES.md).

The general principle, that how you construct negatives drives apparent performance, is not
new, and this repository is not the first to say it. It is best established in TCR-pMHC
specificity prediction. Building negatives by mispairing peptides with other TCRs from the
same dataset lets a model score a pair by recognising a TCR it has seen before rather than by
learning recognition, and every model tested drops substantially once the negatives come from
outside the dataset, which is the argument STAPLER makes [9]. That the construction technique
is itself the bias, and that apparent performance can collapse to chance on unseen epitopes,
has been argued directly [10], and work since then builds harder negatives rather than
shuffled ones [11]. Decoy selection has its own literature in virtual screening.

What I could not find published is this test applied to protein-language-model scoring of
peptide binders. The peptide-binder design literature, PepMLM and target-conditioned masked
language modelling [6], DiffPepBuilder [7], and contrastive target-conditioned design, reports
performance against controls that differ from the positives in composition and length.
Composition-preserving scrambles are not standard practice there, and the result in
[results/RESULTS.md](../results/RESULTS.md) is what happens when you use them. A filter that
reads as working at AUC 0.684 against a composition-matched null drops to 0.569 once the
control keeps each peptide's residue census fixed, and its recall at a matched control
admission rate drops by a factor of 4.3.

There is also a known confound in the same direction worth naming. Protein language models
transfer unevenly to peptide-length sequences, so a score calibrated on protein-length input
is already on uncertain ground before the control question arises. The scorer here is ESM-2
[8], and the generation stack this repository was built around is RFdiffusion [1] with
ProteinMPNN [4].
