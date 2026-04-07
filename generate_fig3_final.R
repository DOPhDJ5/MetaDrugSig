# ============================================================
# MetaDrugSig Figure 3 — Three-Way AUROC Benchmarking
# Panel A: ROC curve (Pyrimidine analogues)
# Panel B: Grouped bar chart (7 classes, 3 databases)
# Panel C: Scatter (mean studies vs REML AUROC, all Tier 1 L4)
# ============================================================

library(data.table)
library(ggplot2)
library(patchwork)
library(pROC)

BASE   <- "~/Desktop/MetaDrugSig_local"
OUTDIR <- file.path(BASE, "figures")

# ============================================================
# COLORS — consistent across all panels
# ============================================================
db_colors <- c(
  "MetaDrugSig REML" = "#B2182B",
  "Single-study"     = "#F4A582",
  "CREEDS"           = "#2166AC"
)

# ============================================================
# LOAD DATA
# ============================================================
cat("Loading data...\n")
results <- fread(file.path(OUTDIR, "threeway_auroc_hpcc_results.csv"))
scatter <- fread(file.path(OUTDIR, "class_auroc_vs_studies.csv"))
sim     <- readRDS(file.path(BASE, "spearman_similarity_v3.rds"))
meta    <- fread(file.path(BASE, "drug_metadata_v3.csv"))
meta_t1 <- meta[tier == 1]

# ============================================================
# CLASSES TO INCLUDE IN PANEL B
# 7 mechanistically defined classes
# ============================================================
keep_classes <- c(
  "Pyrimidine analogues",
  "BCR-ABL tyrosine kinase inhibitors",
  "Retinoids for cancer treatment",
  "ANTINEOPLASTIC AGENTS",
  "ENDOCRINE THERAPY",
  "IMMUNOSUPPRESSANTS",
  "SEX HORMONES AND MODULATORS OF THE GENITAL SYSTEM"
)

# Short display labels
class_labels <- c(
  "Pyrimidine analogues"                              = "Pyrimidine analogues",
  "BCR-ABL tyrosine kinase inhibitors"                = "BCR-ABL inhibitors",
  "Retinoids for cancer treatment"                    = "Retinoids",
  "ANTINEOPLASTIC AGENTS"                             = "Antineoplastic agents",
  "ENDOCRINE THERAPY"                                 = "Endocrine therapy",
  "IMMUNOSUPPRESSANTS"                                = "Immunosuppressants",
  "SEX HORMONES AND MODULATORS OF THE GENITAL SYSTEM" = "Sex hormones"
)

# ============================================================
# PANEL A — ROC curve for Pyrimidine analogues
# ============================================================
cat("Building Panel A — ROC curve...\n")

# Get pyrimidine analogue drug IDs
pyr_ids <- meta_t1$DrugBank_ID[
  !is.na(meta_t1$ATC_L4_Name) &
  meta_t1$ATC_L4_Name == "Pyrimidine analogues" &
  meta_t1$DrugBank_ID %in% rownames(sim)
]
bg_ids <- rownames(sim)

cat("Pyrimidine analogue drugs:", length(pyr_ids), "\n")
cat("Drug names:", paste(meta_t1$DrugBank_Name[meta_t1$DrugBank_ID %in% pyr_ids], collapse=", "), "\n")

# MetaDrugSig ROC — pool all query drugs
all_scores_mds <- c()
all_labels_mds <- c()
for (q in pyr_ids) {
  scores <- sim[q, bg_ids[bg_ids != q]]
  labels <- as.integer(names(scores) %in% pyr_ids)
  all_scores_mds <- c(all_scores_mds, scores)
  all_labels_mds <- c(all_labels_mds, labels)
}
roc_mds  <- roc(all_labels_mds, all_scores_mds, quiet=TRUE)
auc_mds  <- round(as.numeric(auc(roc_mds)), 3)
cat("MetaDrugSig AUC:", auc_mds, "\n")

# Build ROC data frame for plotting
roc_df <- data.frame(
  FPR      = 1 - roc_mds$specificities,
  TPR      = roc_mds$sensitivities,
  Database = "MetaDrugSig REML"
)

# Get mean AUROC from results for annotation
pyr_row <- results[class == "Pyrimidine analogues"]
ss_auc  <- pyr_row$ss_auroc_mean
ss_sd   <- pyr_row$ss_auroc_sd
cr_auc  <- pyr_row$creeds_auroc_mean
cr_sd   <- pyr_row$creeds_auroc_sd

panel_A <- ggplot(roc_df, aes(x=FPR, y=TPR)) +
  geom_abline(intercept=0, slope=1, linetype="dotted",
              color="gray60", linewidth=0.6) +
  geom_line(color="#B2182B", linewidth=1.5) +
  annotate("segment", x=0.35, xend=0.35, y=0.05, yend=0.38,
           color="gray85", linewidth=0.3) +
  annotate("text", x=0.36, y=0.36,
           label=sprintf("MetaDrugSig REML  AUC = %.3f", auc_mds),
           color="#B2182B", size=3.4, hjust=0, fontface="bold") +
  annotate("text", x=0.36, y=0.27,
           label=sprintf("Single-study           AUC = %.3f ± %.3f", ss_auc, ss_sd),
           color="#D6604D", size=3.1, hjust=0) +
  annotate("text", x=0.36, y=0.18,
           label=sprintf("CREEDS                   AUC = %.3f ± %.3f", cr_auc, cr_sd),
           color="#2166AC", size=3.1, hjust=0) +
  annotate("text", x=0.36, y=0.09,
           label="(single-study and CREEDS: mean ± SD, 100 permutations)",
           color="gray50", size=2.6, hjust=0, fontface="italic") +
  scale_x_continuous(limits=c(0,1), expand=c(0.01,0),
                     breaks=seq(0,1,0.25)) +
  scale_y_continuous(limits=c(0,1), expand=c(0.01,0),
                     breaks=seq(0,1,0.25)) +
  coord_equal() +
  labs(title="Pyrimidine analogues classification",
       x="False positive rate",
       y="True positive rate") +
  theme_classic(base_size=12) +
  theme(plot.title   = element_text(face="bold", hjust=0.5, size=11),
        panel.grid   = element_blank(),
        axis.line    = element_line(linewidth=0.5))

# ============================================================
# PANEL B — Grouped bar chart
# ============================================================
cat("Building Panel B — grouped bar chart...\n")

plot_data <- results[class %in% keep_classes]
plot_data[, class_label := class_labels[class]]

# Order by MDS AUROC ascending for horizontal bar chart
class_order_b <- plot_data$class_label[order(plot_data$mds_auroc)]
plot_data[, class_label := factor(class_label, levels=class_order_b)]

# Reshape to long format
long_data <- rbind(
  data.frame(
    class_label = plot_data$class_label,
    database    = "MetaDrugSig REML",
    auroc       = plot_data$mds_auroc,
    ymin        = plot_data$mds_auroc,
    ymax        = plot_data$mds_auroc
  ),
  data.frame(
    class_label = plot_data$class_label,
    database    = "Single-study",
    auroc       = plot_data$ss_auroc_mean,
    ymin        = plot_data$ss_auroc_mean - plot_data$ss_auroc_sd,
    ymax        = plot_data$ss_auroc_mean + plot_data$ss_auroc_sd
  ),
  data.frame(
    class_label = plot_data$class_label,
    database    = "CREEDS",
    auroc       = plot_data$creeds_auroc_mean,
    ymin        = plot_data$creeds_auroc_mean - plot_data$creeds_auroc_sd,
    ymax        = plot_data$creeds_auroc_mean + plot_data$creeds_auroc_sd
  )
)
long_data$class_label <- factor(long_data$class_label, levels=class_order_b)
long_data$database    <- factor(long_data$database,
                                 levels=c("MetaDrugSig REML","Single-study","CREEDS"))

panel_B <- ggplot(long_data,
                  aes(x=class_label, y=auroc, fill=database)) +
  geom_col(position=position_dodge(0.75), width=0.70) +
  geom_errorbar(aes(ymin=pmax(ymin,0), ymax=pmin(ymax,1)),
                position=position_dodge(0.75),
                width=0.22, linewidth=0.45, color="gray30") +
  geom_hline(yintercept=0.5, linetype="dashed",
             color="gray40", linewidth=0.5) +
  annotate("text", x=0.6, y=0.513,
           label="random", color="gray40", size=2.8, hjust=0) +
  scale_fill_manual(values=db_colors, name=NULL) +
  scale_y_continuous(limits=c(0, 1.05),
                     breaks=seq(0,1,0.2),
                     expand=c(0,0)) +
  coord_flip() +
  labs(title="Drug class prediction by transcriptomic\nsignature similarity",
       x=NULL, y="AUROC") +
  theme_classic(base_size=11) +
  theme(legend.position    = "bottom",
        legend.text        = element_text(size=9),
        legend.key.size    = unit(0.4,"cm"),
        plot.title         = element_text(face="bold", size=10, hjust=0.5),
        axis.text.y        = element_text(size=9.5),
        panel.grid.major.x = element_line(color="gray93", linewidth=0.3),
        panel.grid.minor   = element_blank())

# ============================================================
# PANEL C — Scatter: mean studies per drug vs REML AUROC
# All Tier 1 ATC L4 classes
# ============================================================
cat("Building Panel C — scatter plot...\n")

# Exclude catch-all categories
exclude_cls <- c(
  "Other antineoplastic agents",
  "Other protein kinase inhibitors",
  "Other immunosuppressants",
  "Other cytotoxic antibiotics",
  "Other antiepileptics",
  "Other nervous system drugs"
)
scatter_clean <- scatter[!class %in% exclude_cls & !is.na(auroc)]

# Label key classes
label_classes <- c(
  "Anthracyclines and related substances",
  "Histone deacetylase (HDAC) inhibitors",
  "Cyclin-dependent kinase (CDK) inhibitors",
  "Proteasome inhibitors",
  "Epidermal growth factor receptor (EGFR) tyrosine kinase inhibitors",
  "Mitogen-activated protein kinase (MEK) inhibitors",
  "Pyrimidine analogues"
)

# Short labels for scatter
scatter_clean[, label := ifelse(class %in% label_classes,
  fcase(
    class == "Anthracyclines and related substances",                              "Anthracyclines",
    class == "Histone deacetylase (HDAC) inhibitors",                              "HDAC inhibitors",
    class == "Cyclin-dependent kinase (CDK) inhibitors",                           "CDK inhibitors",
    class == "Proteasome inhibitors",                                               "Proteasome\ninhibitors",
    class == "Epidermal growth factor receptor (EGFR) tyrosine kinase inhibitors", "EGFR inhibitors",
    class == "Mitogen-activated protein kinase (MEK) inhibitors",                  "MEK inhibitors",
    class == "Pyrimidine analogues",                                                "Pyrimidine\nanalogues",
    default = ""
  ), "")]

scatter_clean[, highlighted := class %in% label_classes]

# Spearman correlation
cor_val <- cor(scatter_clean$mean_studies, scatter_clean$auroc,
               method="spearman", use="complete.obs")
cat("Spearman rho (mean studies vs AUROC):", round(cor_val, 3), "\n")

panel_C <- ggplot(scatter_clean,
                  aes(x=mean_studies, y=auroc,
                      size=n_drugs, color=highlighted)) +
  geom_point(alpha=0.75) +
  geom_smooth(aes(group=1), method="lm", se=TRUE,
              color="gray40", linewidth=0.8,
              linetype="dashed", fill="gray88",
              inherit.aes=FALSE,
              data=scatter_clean,
              mapping=aes(x=mean_studies, y=auroc)) +
  ggrepel::geom_text_repel(
    data=scatter_clean[highlighted==TRUE],
    aes(label=label),
    size=2.8, color="gray20",
    segment.color="gray60",
    segment.size=0.3,
    box.padding=0.4,
    max.overlaps=20,
    inherit.aes=TRUE
  ) +
  scale_color_manual(values=c("FALSE"="gray75","TRUE"="#B2182B"),
                     guide="none") +
  scale_size_continuous(name="Drugs\nin class", range=c(2,7),
                        breaks=c(3,5,8,12)) +
  annotate("text",
           x=max(scatter_clean$mean_studies, na.rm=TRUE)*0.62,
           y=min(scatter_clean$auroc, na.rm=TRUE)+0.01,
           label=paste0("Spearman ρ = ", round(cor_val, 2)),
           size=3.4, color="gray30", hjust=0) +
  scale_x_continuous(expand=c(0.05,0)) +
  scale_y_continuous(limits=c(0.5,1.02), breaks=seq(0.5,1.0,0.1)) +
  labs(title="Meta-analytic performance\nby drug study depth",
       x="Mean GEO studies per drug (within class)",
       y="REML AUROC") +
  theme_classic(base_size=11) +
  theme(plot.title    = element_text(face="bold", size=10, hjust=0.5),
        legend.position = "right",
        legend.text   = element_text(size=8.5),
        panel.grid.major = element_line(color="gray93", linewidth=0.3))

# ============================================================
# SAVE INDIVIDUAL PANELS
# ============================================================
cat("Saving individual panels...\n")

ggsave(file.path(OUTDIR, "Fig3A_ROC_pyrimidine.pdf"),
       panel_A, width=6, height=5.5, useDingbats=FALSE)
ggsave(file.path(OUTDIR, "Fig3A_ROC_pyrimidine.png"),
       panel_A, width=6, height=5.5, dpi=300)

ggsave(file.path(OUTDIR, "Fig3B_bar_chart.pdf"),
       panel_B, width=7, height=5, useDingbats=FALSE)
ggsave(file.path(OUTDIR, "Fig3B_bar_chart.png"),
       panel_B, width=7, height=5, dpi=300)

ggsave(file.path(OUTDIR, "Fig3C_scatter_studies_auroc.pdf"),
       panel_C, width=7, height=5.5, useDingbats=FALSE)
ggsave(file.path(OUTDIR, "Fig3C_scatter_studies_auroc.png"),
       panel_C, width=7, height=5.5, dpi=300)

# ============================================================
# COMBINED FIGURE 3
# Layout: A (ROC) | B (bar chart) on top row
#         C (scatter) on bottom spanning full width
# ============================================================
cat("Combining Figure 3...\n")

fig3 <- (panel_A | panel_B) / panel_C +
  plot_layout(heights=c(1.1, 1)) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(plot.tag = element_text(face="bold", size=14))
  )

ggsave(file.path(OUTDIR, "Fig3_threeway_AUROC.pdf"),
       fig3, width=13, height=11, useDingbats=FALSE)
ggsave(file.path(OUTDIR, "Fig3_threeway_AUROC.png"),
       fig3, width=13, height=11, dpi=300)

cat("\n=== Figure 3 saved ===\n")
cat("Individual panels: Fig3A, Fig3B, Fig3C\n")
cat("Combined: Fig3_threeway_AUROC\n")
cat("Output directory:", OUTDIR, "\n")

# ============================================================
# PRINT MANUSCRIPT NUMBERS
# ============================================================
cat("\n=== MANUSCRIPT NUMBERS FOR RESULTS ===\n")
cat("Background: 90 matched drugs\n")
cat("Classes shown in Figure 3B:\n")
for (cls in keep_classes) {
  row <- results[class == cls]
  lbl <- class_labels[cls]
  cat(sprintf("  %-25s MDS=%.3f | SS=%.3f±%.3f | CREEDS=%.3f±%.3f\n",
              lbl,
              row$mds_auroc, row$ss_auroc_mean, row$ss_auroc_sd,
              row$creeds_auroc_mean, row$creeds_auroc_sd))
}
cat(sprintf("\nSpearman rho (studies vs AUROC): %.3f\n", cor_val))
cat(sprintf("AUROC range (Panel C): %.3f - %.3f\n",
            min(scatter_clean$auroc, na.rm=TRUE),
            max(scatter_clean$auroc, na.rm=TRUE)))










panel_C <- ggplot(scatter_clean,
                  aes(x=mean_studies, y=auroc,
                      size=n_drugs, color=highlighted)) +
  geom_point(alpha=0.75) +
  geom_smooth(method="lm", se=TRUE,
              color="gray40", linewidth=0.8,
              linetype="dashed", fill="gray88") +
  ggrepel::geom_text_repel(
    data=scatter_clean[highlighted==TRUE],
    aes(label=label),
    size=2.8, color="gray20",
    segment.color="gray60",
    segment.size=0.3,
    box.padding=0.4,
    max.overlaps=20
  ) +
  scale_color_manual(values=c("FALSE"="gray75","TRUE"="#B2182B"),
                     guide="none") +
  scale_size_continuous(name="Drugs\nin class", range=c(2,7),
                        breaks=c(3,5,8,12)) +
  annotate("text",
           x=max(scatter_clean$mean_studies, na.rm=TRUE)*0.62,
           y=min(scatter_clean$auroc, na.rm=TRUE)+0.01,
           label=paste0("Spearman \u03c1 = ", round(cor_val, 2)),
           size=3.4, color="gray30", hjust=0) +
  scale_x_continuous(expand=c(0.05,0)) +
  scale_y_continuous(limits=c(0.5,1.02), breaks=seq(0.5,1.0,0.1)) +
  labs(title="Meta-analytic performance\nby drug study depth",
       x="Mean GEO studies per drug (within class)",
       y="REML AUROC") +
  theme_classic(base_size=11) +
  theme(plot.title       = element_text(face="bold", size=10, hjust=0.5),
        legend.position  = "right",
        legend.text      = element_text(size=8.5),
        panel.grid.major = element_line(color="gray93", linewidth=0.3))

# Save Panel C
ggsave(file.path(OUTDIR, "Fig3C_scatter_studies_auroc.pdf"),
       panel_C, width=7, height=5.5, useDingbats=FALSE)
ggsave(file.path(OUTDIR, "Fig3C_scatter_studies_auroc.png"),
       panel_C, width=7, height=5.5, dpi=300)

# Combined figure
fig3 <- (panel_A | panel_B) / panel_C +
  plot_layout(heights=c(1.1, 1)) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(plot.tag = element_text(face="bold", size=14))
  )

ggsave(file.path(OUTDIR, "Fig3_threeway_AUROC.pdf"),
       fig3, width=13, height=11, useDingbats=FALSE)
ggsave(file.path(OUTDIR, "Fig3_threeway_AUROC.png"),
       fig3, width=13, height=11, dpi=300)
cat("Figure 3 complete.\n")


label=paste0("Spearman rho = ", round(cor_val, 2)),

ggplot(scatter_clean,
       aes(x=mean_studies, y=auroc, color=highlighted)) +
  geom_point(aes(size=n_drugs), alpha=0.75) +

  panel_C <- ggplot(scatter_clean,
                    aes(x=mean_studies, y=auroc, color=highlighted)) +
  geom_point(aes(size=n_drugs), alpha=0.75) +
  geom_smooth(method="lm", se=TRUE,
              color="gray40", linewidth=0.8,
              linetype="dashed", fill="gray88") +
  ggrepel::geom_text_repel(
    data=scatter_clean[highlighted==TRUE],
    aes(label=label),
    size=2.8, color="gray20",
    segment.color="gray60",
    segment.size=0.3,
    box.padding=0.4,
    max.overlaps=20
  ) +
  scale_color_manual(values=c("FALSE"="gray75","TRUE"="#B2182B"),
                     guide="none") +
  scale_size_continuous(name="Drugs\nin class", range=c(2,7),
                        breaks=c(3,5,8,12)) +
  annotate("text",
           x=max(scatter_clean$mean_studies, na.rm=TRUE)*0.62,
           y=min(scatter_clean$auroc, na.rm=TRUE)+0.01,
           label=paste0("Spearman rho = ", round(cor_val, 2)),
           size=3.4, color="gray30", hjust=0) +
  scale_x_continuous(expand=c(0.05,0)) +
  scale_y_continuous(limits=c(0.5,1.02), breaks=seq(0.5,1.0,0.1)) +
  labs(title="Meta-analytic performance\nby drug study depth",
       x="Mean GEO studies per drug (within class)",
       y="REML AUROC") +
  theme_classic(base_size=11) +
  theme(plot.title       = element_text(face="bold", size=10, hjust=0.5),
        legend.position  = "right",
        legend.text      = element_text(size=8.5),
        panel.grid.major = element_line(color="gray93", linewidth=0.3))

ggsave(file.path(OUTDIR, "Fig3C_scatter_studies_auroc.png"),
       panel_C, width=7, height=5.5, dpi=300)

fig3 <- (panel_A | panel_B) / panel_C +
  plot_layout(heights=c(1.1, 1)) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(plot.tag = element_text(face="bold", size=14))
  )

ggsave(file.path(OUTDIR, "Fig3_threeway_AUROC.png"),
       fig3, width=13, height=11, dpi=300)
cat("Done.\n")


# Use mean of per-drug AUROCs instead of pooled
pyr_aurocs <- sapply(pyr_ids, function(q) {
  scores <- sim[q, bg_ids[bg_ids != q]]
  labels <- as.integer(names(scores) %in% pyr_ids)
  as.numeric(pROC::roc(labels, scores, quiet=TRUE)$auc)
})
cat("Per-drug AUROCs:", round(pyr_aurocs, 3), "\n")
cat("Mean AUROC:", round(mean(pyr_aurocs), 3), "\n")

# Use pooled ROC for curve shape but annotate with mean AUROC
auc_mds <- round(mean(pyr_aurocs), 3)
