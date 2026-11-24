# MetaDrugSig: leave-one-study-out (LOSO) validation of τ²

Code and outputs for the LOSO analysis in the MetaDrugSig manuscript (*Database*, Oxford), Figure 6 and Supplementary Tables.

**Question:** is τ² estimated from the *other* studies informative about how a gene behaves in a study that was held out (interval coverage, direction and magnitude)?

## Design
- **Held-out unit:** one GEO series (GSE), with all its comparisons for that drug pooled by fixed-effect inverse-variance weighting, because comparisons within a series share controls.
- **Eligible drugs:** Tier 1, not stratum C, at least 3 GSEs, so at least 2 remain for training (351 drugs; 326 retained after quality control).
- **Folds:** every GSE, capped at 10 randomly chosen GSEs per drug (`--max_folds`, seed 42).
- **Training:** REML consensus effect, SE and τ² from the remaining series, computed per gene (genes need at least 2 training series). The vectorised REML (Fisher scoring) was checked against `metafor::rma` (maximum |Δτ²| < 5e-5 on simulated data).
- **SE** recovered from the edgeR F statistic as |logFC|/√F, as in the portal.
- **Outcomes in the held-out study:** (1) coverage of the 95% prediction interval θ̂ ± 1.96·√(SE² + τ² + SE²ₕ) versus the same interval without τ²; (2) same-direction rate and fraction of the consensus effect retained, stratified by training τ² and by |consensus log₂FC|; (3) MDS-Broad and MDS-Robust gene sets; (4) top-N gene ranking by single-study |z|, fixed-effect |z|, REML |z|, REML |log₂FC| and the continuous MDS-Robust score |θ̂|/√(τ² + SE²).
- **Uncertainty:** percentile bootstrap over drugs (500 replicates).
- Main figure and tables use genes with at least 5 training comparisons (panel C uses all genes).

## Inputs
Both are in the MetaDrugSig Zenodo record (https://doi.org/10.5281/zenodo.19432765):
`DEG_results_v3.tar.gz` (per-study edgeR results, read in one pass without unpacking) and `drug_metadata_v3.csv`.

## Run
```bash
# 1. per-drug LOSO tables (Python 3, numpy, pandas; scipy optional). Restartable.
python3 loso_stream.py --tar DEG_results_v3.tar.gz --meta drug_metadata_v3.csv --out results
#    --selftest runs a simulated check; --max_drugs 15 runs a pilot; --n_sub 0 skips the optional gene subsample table

# 2. summary tables and figure (R: data.table, ggplot2, patchwork)
Rscript loso_summarize.R --res results --out loso_report
```
To reproduce the figure from the saved per-drug outputs without rerunning step 1:
```bash
mkdir -p results/per_drug && tar xzf loso_per_drug_outputs.tgz -C results/per_drug
mv results/per_drug/drugs_skipped.csv results/
Rscript loso_summarize.R --res results --out loso_report
```

## Files
| File | Contents |
|---|---|
| `loso_stream.py` | Streaming LOSO analysis (produces the manuscript results) |
| `loso_summarize.R` | Bootstrap summaries, tables and Figure 6 |
| `loso_functions.R` | R reference implementation of the same model (used to cross-check the Python REML) |
| `loso_per_drug_outputs.tgz` | Per-drug output tables (326 drugs) |
| `results_tables/` | Supplementary Tables S1–S4 as CSV, and drugs skipped |
| `figures/` | Figure 6 |
