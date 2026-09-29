# References

Every work named in the README, `ANALYSIS.md`, `docs/DESIGN.md` and `docs/LITERATURE.md`.
Method names in brackets are the names in common use, which in several cases do not appear
in the published title.

## Design and generation

1. Watson, J. L., Juergens, D., Bennett, N. R., et al. De novo design of protein structure
   and function with RFdiffusion. *Nature* 620, 1089-1100 (2023).
   doi:10.1038/s41586-023-06415-8 [RFdiffusion]
2. Bennett, N. R., Watson, J. L., Ragotte, R. J., et al. Atomically accurate de novo design
   of antibodies with RFdiffusion. *Nature* 649, 183-193 (2025).
   doi:10.1038/s41586-025-09721-5
3. Rettie, S. A., Campbell, K. V., Bera, A. K., et al. Accurate de novo design of
   high-affinity protein-binding macrocycles using deep learning. *Nature Chemical Biology*
   21, 1948-1956 (2025). doi:10.1038/s41589-025-01929-w [RFpeptides]
4. Dauparas, J., Anishchenko, I., Bennett, N., et al. Robust deep learning-based protein
   sequence design using ProteinMPNN. *Science* 378, 49-56 (2022).
   doi:10.1126/science.add2187
5. Yang, W., Wang, S., Lee, G. R., et al. The past, present and future of de novo protein
   design. *Nature* 652, 1139-1152 (2026). doi:10.1038/s41586-026-10328-7

## Peptide binder design

6. Chen, L. T., Dumas, M., Watson, R., et al. Target sequence-conditioned design of peptide
   binders using masked language modeling. *Nature Biotechnology* 44, 1002-1010 (2026).
   doi:10.1038/s41587-025-02761-2 [PepMLM, named in the earlier preprint
   arXiv:2310.03842]
7. Wang, F., Wang, Y., Feng, L., Zhang, C. & Lai, L. Target-specific de novo peptide binder
   design with DiffPepBuilder. *Journal of Chemical Information and Modeling* 64, 9135-9149
   (2024). doi:10.1021/acs.jcim.4c00975

## The scorer used here

8. Lin, Z., Akin, H., Rao, R., et al. Evolutionary-scale prediction of atomic-level protein
   structure with a language model. *Science* 379, 1123-1130 (2023).
   doi:10.1126/science.ade2574 [ESM-2, the checkpoint family this repository scores with]

## Negative-control construction, which is the prior art for the experiment here

9. Kwee, B. P. Y., Messemaker, M., Marcus, E., et al. STAPLER: efficient learning of
   TCR-peptide specificity prediction from full-length TCR-peptide data. *bioRxiv*
   2023.04.25.538237 (2023). doi:10.1101/2023.04.25.538237. Preprint, no peer-reviewed
   version located.
10. Dens, C., Laukens, K., Bittremieux, W. & Meysman, P. The pitfalls of negative data bias
    for the T-cell epitope specificity challenge. *Nature Machine Intelligence* 5, 1060-1062
    (2023). doi:10.1038/s42256-023-00727-0. Answered by Gao et al., *Nature Machine
    Intelligence* 5, 1063-1065 (2023), doi:10.1038/s42256-023-00725-2.
11. Mi, X., Zhu, J., Dai, Z., et al. Mitigating negative data bias to enhance TCR-epitope
    binding and residue interaction prediction. *Briefings in Bioinformatics* 27, bbag418
    (2026). doi:10.1093/bib/bbag418
12. Moris, P., De Pauw, J., Postovskaya, A., et al. Current challenges for unseen-epitope
    TCR interaction prediction and a new perspective derived from image classification.
    *Briefings in Bioinformatics* 22, bbaa318 (2021). doi:10.1093/bib/bbaa318
