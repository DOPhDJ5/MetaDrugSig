# ============================================================
# generate_fig5_rges_pathway.R
# Figure 5: RGES pathway validation
# Known drug classes rank significantly toward expected pathway extremes
# Panels A-C: TNFa/NF-kB, Cell Cycle, Cholesterol Homeostasis
# ============================================================

library(data.table)
library(ggplot2)
library(patchwork)
library(msigdbr)

BASE    <- "~/Desktop/MetaDrugSig_local"
OUT_DIR <- "~/Desktop/MetaDrugSig_local/figures"
dir.create(OUT_DIR, showWarnings=FALSE)

# ============================================================
# LOAD DATA
# ============================================================
cat("Loading data...\n")
logfc_dt <- fread(file.path(BASE, "drug_signatures_logFC.csv"))
genes    <- logfc_dt$gene
logfc_dt[, gene := NULL]
logfc_mat <- as.matrix(logfc_dt)
rownames(logfc_mat) <- genes

meta <- fread(file.path(BASE, "drug_metadata_v3.csv"))
meta <- meta[has_signature == TRUE & tier == 1]
tier1_ids <- meta$DrugBank_ID

cat("Tier 1 drugs:", length(tier1_ids), "\n")

# ============================================================
# RGES FUNCTION — canonical Chen Lab / OCTAD implementation
# rank(-1 * logFC) so rank 1 = most upregulated
# Negative RGES = drug REVERSES query = candidate therapeutic
# ============================================================
cmap_score_new <- function(sig_up, sig_down, drug_signature, gamma=1) {
  num_genes <- nrow(drug_signature)
  drug_signature[,"rank"] <- rank(-1 * drug_signature[,"rank"], ties.method="random")
  up_tags_rank   <- merge(drug_signature, sig_up,   by.x="ids", by.y=1)
  down_tags_rank <- merge(drug_signature, sig_down, by.x="ids", by.y=1)
  up_pos   <- sort(up_tags_rank$rank)
  down_pos <- sort(down_tags_rank$rank)
  ks <- function(pos, n) {
    a <- max(sapply(seq_len(n), function(j) j/n - pos[j]/num_genes))
    b <- max(sapply(seq_len(n), function(j) pos[j]/num_genes - (j-1)/n))
    if (a > b) a else -b
  }
  ks_up   <- if (length(up_pos)   > 1) ks(up_pos,   length(up_pos))   else 0
  ks_down <- if (length(down_pos) > 1) ks(down_pos, length(down_pos)) else 0
  if (ks_up == 0 & ks_down != 0) return(-ks_down)
  if (ks_up != 0 & ks_down == 0) return(gamma * ks_up)
  return(gamma * ks_up - ks_down)
}

compute_rges_all <- function(sig_up_genes, sig_down_genes=character(0)) {
  sig_up   <- data.frame(gene=sig_up_genes)
  sig_down <- data.frame(gene=sig_down_genes)
  ids <- intersect(tier1_ids, colnames(logfc_mat))
  scores <- vapply(ids, function(id) {
    v <- logfc_mat[, id]; v <- v[!is.na(v)]
    tryCatch(
      cmap_score_new(sig_up, sig_down, data.frame(ids=names(v), rank=v)),
      error=function(e) NA_real_)
  }, numeric(1))
  data.frame(DrugBank_ID=ids, RGES=scores) |>
    merge(as.data.frame(meta[, .(DrugBank_ID, DrugBank_Name,
                                  ATC_L2_Name, ATC_L4_Name)]),
          by="DrugBank_ID") |>
    subset(!is.na(RGES))
}

# ============================================================
# LOAD MSIGDB PATHWAYS
# ============================================================
cat("Loading MSigDB...\n")
msig_h <- msigdbr(species="Homo sapiens", collection="H")
get_pathway <- function(name) {
  msig_h$gene_symbol[msig_h$gs_name == name]
}

tnfa_genes   <- get_pathway("HALLMARK_TNFA_SIGNALING_VIA_NFKB")
g2m_genes    <- get_pathway("HALLMARK_G2M_CHECKPOINT")
chol_genes   <- get_pathway("HALLMARK_CHOLESTEROL_HOMEOSTASIS")

cat("TNFa pathway genes in DB:",  sum(tnfa_genes %in% rownames(logfc_mat)), "\n")
cat("G2M pathway genes in DB:",   sum(g2m_genes  %in% rownames(logfc_mat)), "\n")
cat("Chol pathway genes in DB:",  sum(chol_genes %in% rownames(logfc_mat)), "\n")

# ============================================================
# COMPUTE RGES FOR EACH PATHWAY
# ============================================================
cat("Computing RGES — TNFa/NF-kB...\n")
rges_tnfa <- compute_rges_all(tnfa_genes)

cat("Computing RGES — Cell Cycle (G2M)...\n")
rges_g2m  <- compute_rges_all(g2m_genes)

cat("Computing RGES — Cholesterol Homeostasis...\n")
rges_chol <- compute_rges_all(chol_genes)

# ============================================================
# DEFINE DRUG CLASSES
# ============================================================
# Corticosteroids: ATC H02AB
corticosteroid_ids <- meta$DrugBank_ID[
  !is.na(meta$ATC_L4_Name) &
  grepl("corticosteroid|glucocorticoid", meta$ATC_L4_Name, ignore.case=TRUE)]

# CDK4/6 inhibitors: palbociclib, ribociclib, abemaciclib
cdk_names <- c("Palbociclib","Ribociclib","Abemaciclib","Trilaciclib")
cdk_ids   <- meta$DrugBank_ID[meta$DrugBank_Name %in% cdk_names]

# Statins: ATC C10AA
statin_ids <- meta$DrugBank_ID[
  !is.na(meta$ATC_L4_Name) &
  grepl("statin|HMG CoA reductase", meta$ATC_L4_Name, ignore.case=TRUE)]

cat("Corticosteroids:", length(corticosteroid_ids), "\n")
cat("CDK4/6 inhibitors:", length(cdk_ids), "\n")
cat("Statins:", length(statin_ids), "\n")

# ============================================================
# PREPARE PLOT DATA
# ============================================================
prep_plot_data <- function(rges_df, class_ids, class_label) {
  rges_df$group <- ifelse(rges_df$DrugBank_ID %in% class_ids,
                           class_label, "All other drugs")
  rges_df$group <- factor(rges_df$group,
                           levels=c("All other drugs", class_label))
  rges_df
}

pd_tnfa <- prep_plot_data(rges_tnfa, corticosteroid_ids, "Corticosteroids")
pd_g2m  <- prep_plot_data(rges_g2m,  cdk_ids,            "CDK4/6 inhibitors")
pd_chol <- prep_plot_data(rges_chol, statin_ids,         "Statins")

# ============================================================
# WILCOXON TESTS
# ============================================================
wt_tnfa <- wilcox.test(
  RGES ~ group, data=pd_tnfa,
  alternative="greater")  # corticosteroids more negative
wt_g2m  <- wilcox.test(
  RGES ~ group, data=pd_g2m,
  alternative="greater")  # CDK inhibitors more negative
wt_chol <- wilcox.test(
  RGES ~ group, data=pd_chol,
  alternative="greater")  # statins more negative

fmt_p <- function(p) {
  if (p < 0.001) sprintf("p = %.5f", p)
  else sprintf("p = %.4f", p)
}

cat("\nWilcoxon results:\n")
cat("TNFa/NF-kB (corticosteroids):", fmt_p(wt_tnfa$p.value), "\n")
cat("Cell Cycle (CDK4/6):",          fmt_p(wt_g2m$p.value),  "\n")
cat("Cholesterol (statins):",         fmt_p(wt_chol$p.value), "\n")

# ============================================================
# PLOT THEME
# ============================================================
theme_fig5 <- function() {
  theme_classic(base_size=10) +
  theme(
    axis.title.x    = element_blank(),
    axis.text.x     = element_text(size=9, color="black"),
    axis.text.y     = element_text(size=8),
    axis.title.y    = element_text(size=9),
    plot.title      = element_text(size=10, face="bold", hjust=0.5),
    plot.subtitle   = element_text(size=9,  hjust=0.5, color="#333333"),
    legend.position = "none",
    panel.grid.major.y = element_line(color="grey93", linewidth=0.3)
  )
}

GREY  <- "#555555"
RED   <- "#C0392B"
JITTER_W <- 0.15

make_panel <- function(pd, class_label, pathway_title, pval) {
  n_class <- sum(pd$group == class_label)
  n_other <- sum(pd$group == "All other drugs")

  ggplot(pd, aes(x=group, y=RGES)) +
    geom_jitter(aes(color=group), width=JITTER_W, size=0.7, alpha=0.4,
                show.legend=FALSE) +
    geom_boxplot(aes(color=group), width=0.45, outlier.shape=NA,
                 fill=NA, linewidth=0.7, show.legend=FALSE) +
    scale_color_manual(values=c("All other drugs"=GREY, setNames(RED, class_label))) +
    scale_x_discrete(labels=function(x) {
      ifelse(x=="All other drugs",
             paste0("All other drugs\n(n=",n_other,")"),
             paste0(x,"\n(n=",n_class,")"))
    }) +
    labs(y="RGES Score",
         title=pathway_title,
         subtitle=paste0("Wilcoxon ", fmt_p(pval))) +
    theme_fig5()
}

pA <- make_panel(pd_tnfa, "Corticosteroids",   "TNF\u03b1 / NF-\u03baB Signaling", wt_tnfa$p.value)
pB <- make_panel(pd_g2m,  "CDK4/6 inhibitors", "Cell Cycle",                         wt_g2m$p.value)
pC <- make_panel(pd_chol, "Statins",            "Cholesterol Homeostasis",            wt_chol$p.value)

# ============================================================
# COMBINE AND SAVE
# ============================================================
fig5 <- pA + pB + pC +
  plot_layout(ncol=3) +
  plot_annotation(
    title   = "Known drug classes rank significantly toward expected pathway extremes",
    theme   = theme(
      plot.title = element_text(size=11, face="bold", hjust=0.5,
                                margin=margin(b=8))
    )
  )

out_pdf <- file.path(OUT_DIR, "Fig5_rges_pathway_validation.pdf")
out_png <- file.path(OUT_DIR, "Fig5_rges_pathway_validation.png")

ggsave(out_pdf, fig5, width=10, height=4, useDingbats=FALSE)
ggsave(out_png, fig5, width=10, height=4, dpi=300)

cat("\nSaved:\n", out_pdf, "\n", out_png, "\n")

# ============================================================
# SAVE RGES RESULTS FOR SUPPLEMENTARY
# ============================================================
rges_tnfa$pathway <- "HALLMARK_TNFA_SIGNALING_VIA_NFKB"
rges_g2m$pathway  <- "HALLMARK_G2M_CHECKPOINT"
rges_chol$pathway <- "HALLMARK_CHOLESTEROL_HOMEOSTASIS"

rges_all <- rbind(rges_tnfa, rges_g2m, rges_chol)
fwrite(rges_all, file.path(OUT_DIR, "Fig5_rges_scores_all_pathways.csv"))
cat("RGES scores saved.\n")
