# Wound and PDAC fibroblast analysis

Lightweight reproduction of selected results from *A shared matrix program and distinct contractile states in human wound and pancreatic ductal adenocarcinoma fibroblasts*.

## Run

From the repository directory, use base R. The release was checked with R 4.5.1.

```sh
Rscript reproduce.R
```

The script writes nine CSV tables to `results/`. Pass an output directory as the first argument to save them elsewhere.

## Results reproduced

| Manuscript panel | Steps and input level | Expected result |
|---|---|---|
| 1D–E | Intersect the saved W2 and F1 top-200 gene lists and count the NABA annotations by division. | 74 shared genes, including 15 core and 12 associated matrisome genes. |
| 3A–D | Summarize saved estimates from 12 primary PDAC sections and the fitted PDAC-minus-D30 curves over shared ECM score ranges. | Median F1 Moran's I of 0.808. Positive within-section changes in all 12 sections and positive fitted D30 contrasts in both cohorts. |
| 3E | Recompute Spearman correlations within each of the 16 external patients from observed Contractile10 scores and saved predictions. | Median 0.764 for the elastic net and 0.529 for the baseline. Performance improved in all 16 patients. |
| 4D–E | Rank measurable genes using discovery regulatory potential and calculate external patient-level means across the fixed target sets. Input values are high-minus-low log₂ CPM. | Positive mean responses in 16/16 patients for TGFB1 and 15/16 for TGFB2. |

## Input provenance

- `Wound_W2_top200.csv`, `PDAC_F1_top200.csv` and `matrisome_annotations.csv` contain the retained cNMF rankings and NABA gene annotation classes. The single-cell discovery inputs are [GSE241132](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE241132) and [GSE212966](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE212966).
- `pdac_section_*.csv` contains section-level results from the 12 primary CytAssist sections of [GSE278687](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE278687). `*_D30_contrast_curve.csv` contains saved fitted comparisons with D30 skin-wound fibroblasts.
- `CRA001160_fixed_predictions.csv` contains one anonymized score-and-prediction row per eligible fibroblast from [CRA001160](https://ngdc.cncb.ac.cn/gsa/browse/CRA001160). The patient codes are retained to calculate patient-wise correlations.
- `CRA001160_fixed_TGFB_targets.csv` contains the discovery-fixed targets, their discovery regulatory potentials and the external patient-level expression differences. Blank or `NA` values identify targets unavailable in the external evaluation.

The reproduction starts from these versioned intermediate outputs. The article's Supplementary Methods document the upstream preprocessing, cNMF, spatial modeling, elastic-net fitting and NicheNet analyses.
