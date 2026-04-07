# MetaDrugSig

**A heterogeneity-aware transcriptomic drug signature database**

REML random-effects meta-analysis of 1,092 drugs across 6,663 independent RNA-seq comparisons from 3,683 GEO studies, with gene-level uncertainty quantification and the MDS signature framework for heterogeneity-aware gene selection.

---

## Database

Full database available on Zenodo:  
**https://doi.org/10.5281/zenodo.19154955**

### Gene-by-drug matrices (63,474 genes × 1,092 drugs)

| File | Description |
|------|-------------|
| `drug_signatures_logFC.csv` | Consensus log₂ fold changes |
| `drug_signatures_FDR.csv` | Benjamini-Hochberg adjusted p-values |
| `drug_signatures_SE.csv` | Standard errors |
| `drug_signatures_tau2.csv` | Between-study variance τ² (Tier 1 only) |

### Metadata files

| File | Description |
|------|-------------|
| `drug_metadata_v3.csv` | Drug annotations, ATC codes, tier designation |
| `drug_comparisons_v3.csv` | Per-comparison study metadata |
| `comparison_metadata_v3.csv` | Tissue, dose, duration annotations per comparison |
| `spearman_similarity_v3.rds` | Drug-drug Spearman similarity matrix |
| `DEG_results_v3.tar.gz` | Per-study edgeR differential expression results |

---

## Portal

Interactive portal: **[URL]**

### Run locally

```r
# Install dependencies
source("requirements.R")

# Place data files in MetaDrugSig_app/data/ (download from Zenodo)
# Then run:
setwd("MetaDrugSig_app/")
shiny::runApp("app.R")
```

---

## MDS Signature Framework

Three gene selection filters of increasing stringency:

| Filter | Criterion | Use case |
|--------|-----------|----------|
| **MDS-Broad** | FDR < 0.05, \|logFC\| > 0.5 | Exploratory analyses |
| **MDS-Robust** | 95% prediction interval (logFC ± 1.96√(τ²+SE²)) excludes zero | Directionally reliable genes across any new study |
| **MDS-Context** | Filter studies by tissue/dose/duration before recomputing REML | Context-specific drug effects |

---

## Tiers

| Tier | Criterion | N drugs | Method |
|------|-----------|---------|--------|
| Tier 1 | ≥2 independent comparisons | 556 | REML random-effects meta-analysis |
| Tier 2 | Single comparison | 536 | edgeR differential expression |

---

## Citation

Joseph DF, Zhang X, Grover S, Talukdar N, Palmer Z, Koirala S, Leschiner D, Shankar R, Chen B.  
**MetaDrugSig: A heterogeneity-aware transcriptomic drug signature database.**  
*Database (Oxford)*, 2026.  
Zenodo: https://doi.org/10.5281/zenodo.19154955

---

## Contact

Dimitri F. Joseph — josephd9@msu.edu  
Bin Chen Lab, Michigan State University  
Department of Pharmacology and Toxicology
