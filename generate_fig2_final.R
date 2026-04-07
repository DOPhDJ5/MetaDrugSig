# ============================================================
# MetaDrugSig Figure 2 — Pharmacological Class Coherence
# Panel A: UMAP colored by ATC class
# Panel B: Heatmap ordered by class with drug name labels
# Panel C: Bar chart intra vs inter correlation per class
# ============================================================

library(data.table)
library(ggplot2)
library(patchwork)
library(umap)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(plotly)
library(htmlwidgets)

BASE   <- "~/Desktop/MetaDrugSig_local"
OUTDIR <- file.path(BASE, "figures")
dir.create(OUTDIR, showWarnings=FALSE)

# ============================================================
# LOAD DATA
# ============================================================
cat("Loading data...\n")
sim     <- readRDS(file.path(BASE, "spearman_similarity_v3.rds"))
meta    <- fread(file.path(BASE, "drug_metadata_v3.csv"))
meta_t1 <- meta[tier == 1]
stats   <- fread(file.path(OUTDIR, "class_coherence_statistics_v3.csv"))

# ============================================================
# CLASS DEFINITIONS — 6 mechanistically defined ATC L4 classes
# Anthracyclines listed first per Results narrative
# ============================================================
highlight_classes <- c(
  "Anthracyclines and related substances",
  "Histone deacetylase (HDAC) inhibitors",
  "Epidermal growth factor receptor (EGFR) tyrosine kinase inhibitors",
  "BCR-ABL tyrosine kinase inhibitors",
  "Poly (ADP-ribose) polymerase (PARP) inhibitors",
  "Pyrimidine analogues"
)

class_labels <- c(
  "Anthracyclines and related substances"                                = "Anthracyclines",
  "Histone deacetylase (HDAC) inhibitors"                                = "HDAC inhibitors",
  "Epidermal growth factor receptor (EGFR) tyrosine kinase inhibitors"   = "EGFR inhibitors",
  "BCR-ABL tyrosine kinase inhibitors"                                   = "BCR-ABL inhibitors",
  "Poly (ADP-ribose) polymerase (PARP) inhibitors"                       = "PARP inhibitors",
  "Pyrimidine analogues"                                                  = "Pyrimidine analogues"
)

meta_t1[, drug_class := ifelse(
  ATC_L4_Name %in% names(class_labels),
  class_labels[ATC_L4_Name],
  "Other/Unclassified"
)]

# Align to sim matrix
common_drugs <- intersect(colnames(sim), meta_t1$DrugBank_ID)
sim_sub      <- sim[common_drugs, common_drugs]
meta_sub     <- meta_t1[match(common_drugs, DrugBank_ID)]

# Color palette — anthracyclines first = most prominent red
class_order <- c("Anthracyclines", "HDAC inhibitors", "EGFR inhibitors",
                 "BCR-ABL inhibitors", "PARP inhibitors", "Pyrimidine analogues")
class_levels <- c(class_order, "Other/Unclassified")
n_named      <- length(class_order)

# Manually assign colors for clarity
named_cols <- setNames(
  c("#B2182B", "#F4A582", "#2166AC", "#4DAC26", "#762A83", "#E08214"),
  class_order
)
class_cols <- c(named_cols, "Other/Unclassified" = "gray88")
meta_sub[, drug_class := factor(drug_class, levels=class_levels)]

# ============================================================
# PANEL A — UMAP
# ============================================================
cat("Computing UMAP...\n")
set.seed(42)
dist_mat               <- as.matrix(1 - sim_sub)
dist_mat[dist_mat < 0] <- 0

cfg              <- umap.defaults
cfg$input        <- "dist"
cfg$n_neighbors  <- 20
cfg$min_dist     <- 0.3
cfg$metric       <- "euclidean"
cfg$n_epochs     <- 500
cfg$random_state <- 42

umap_result <- umap(dist_mat, config=cfg)

umap_df <- data.frame(
  UMAP1       = umap_result$layout[,1],
  UMAP2       = umap_result$layout[,2],
  DrugBank_ID = common_drugs,
  drug_class  = meta_sub$drug_class,
  drug_name   = meta_sub$DrugBank_Name,
  stringsAsFactors = FALSE
)
umap_df$drug_class <- factor(umap_df$drug_class, levels=class_levels)

# Static UMAP
p_umap <- ggplot(umap_df, aes(x=UMAP1, y=UMAP2, color=drug_class)) +
  geom_point(data=subset(umap_df, drug_class=="Other/Unclassified"),
             size=0.9, alpha=0.30) +
  geom_point(data=subset(umap_df, drug_class!="Other/Unclassified"),
             size=2.8, alpha=0.92) +
  scale_color_manual(values=class_cols, name="Drug class") +
  labs(title="UMAP of Tier 1 Drug Transcriptomic Signatures",
       subtitle=sprintf("%d REML meta-analytic signatures · ATC pharmacological class",
                        nrow(umap_df)),
       x="UMAP 1", y="UMAP 2") +
  theme_classic(base_size=11) +
  theme(legend.text     = element_text(size=8.5),
        legend.key.size = unit(0.4, "cm"),
        legend.title    = element_text(size=9.5, face="bold"),
        plot.title      = element_text(face="bold", size=11),
        plot.subtitle   = element_text(size=8.5, color="gray40")) +
  guides(color=guide_legend(override.aes=list(size=3), ncol=1))

ggsave(file.path(OUTDIR, "Fig2A_UMAP_tier1.pdf"),
       p_umap, width=8, height=6, useDingbats=FALSE)
ggsave(file.path(OUTDIR, "Fig2A_UMAP_tier1.png"),
       p_umap, width=8, height=6, dpi=300)
cat("UMAP saved.\n")

# Interactive UMAP
umap_df$label <- paste0("<b>", umap_df$drug_name, "</b><br>",
                         umap_df$drug_class, "<br>", umap_df$DrugBank_ID)

p_int <- plot_ly() %>%
  add_trace(data=subset(umap_df, drug_class=="Other/Unclassified"),
            x=~UMAP1, y=~UMAP2, text=~label, hoverinfo="text",
            type="scatter", mode="markers", name="Other/Unclassified",
            marker=list(color="gray88", size=5, opacity=0.25)) %>%
  add_trace(data=subset(umap_df, drug_class!="Other/Unclassified"),
            x=~UMAP1, y=~UMAP2, color=~drug_class,
            colors=class_cols[names(class_cols)!="Other/Unclassified"],
            text=~label, hoverinfo="text",
            type="scatter", mode="markers",
            marker=list(size=10, opacity=0.95)) %>%
  layout(title=list(text="MetaDrugSig — Tier 1 Drug Signatures (UMAP)",
                    font=list(size=14)),
         xaxis=list(title="UMAP 1", zeroline=FALSE),
         yaxis=list(title="UMAP 2", zeroline=FALSE),
         legend=list(title=list(text="<b>Drug class</b>")),
         plot_bgcolor="white", paper_bgcolor="white")

saveWidget(p_int, file.path(OUTDIR, "Fig2A_UMAP_interactive.html"),
           selfcontained=TRUE)
cat("Interactive UMAP saved.\n")

# ============================================================
# PANEL B — HEATMAP ordered by class
# ============================================================
cat("Building heatmap...\n")
ht_opt(RESET=TRUE)

hm_drugs <- meta_sub[drug_class != "Other/Unclassified"]
hm_drugs <- hm_drugs[order(drug_class, DrugBank_Name)]
hm_ids   <- hm_drugs$DrugBank_ID
hm_mat   <- sim_sub[hm_ids, hm_ids]
rownames(hm_mat) <- hm_drugs$DrugBank_Name
colnames(hm_mat) <- hm_drugs$DrugBank_Name

# Use same colors as UMAP
hm_classes <- as.character(unique(hm_drugs$drug_class))
hm_cols    <- class_cols[hm_classes]

row_ann <- rowAnnotation(
  Class = as.character(hm_drugs$drug_class),
  col   = list(Class=hm_cols),
  annotation_legend_param = list(
    Class=list(title="Drug class",
               title_gp=gpar(fontsize=9, fontface="bold"),
               labels_gp=gpar(fontsize=8))),
  show_annotation_name=FALSE
)

col_ann <- HeatmapAnnotation(
  Class = as.character(hm_drugs$drug_class),
  col   = list(Class=hm_cols),
  show_annotation_name=FALSE,
  show_legend=FALSE
)

col_fun <- colorRamp2(
  c(-0.5, -0.2, 0, 0.2, 0.5),
  c("#2166AC", "#92C5DE", "white", "#F4A582", "#D6604D")
)

ht <- Heatmap(
  as.matrix(hm_mat),
  name                 = "Spearman ρ",
  col                  = col_fun,
  cluster_rows         = FALSE,
  cluster_columns      = FALSE,
  show_row_names       = TRUE,
  row_names_gp         = gpar(fontsize=7),
  show_column_names    = FALSE,
  top_annotation       = col_ann,
  right_annotation     = row_ann,
  column_title         = "Drug Signature Similarity — Ordered by Pharmacological Class",
  column_title_gp      = gpar(fontsize=10, fontface="bold"),
  heatmap_legend_param = list(title_gp=gpar(fontsize=9, fontface="bold"),
                               labels_gp=gpar(fontsize=8)),
  use_raster=FALSE
)

cat("Saving heatmap...\n")
pdf(file.path(OUTDIR, "Fig2B_heatmap_ordered_by_class.pdf"),
    width=12, height=9, onefile=FALSE)
draw(ht, heatmap_legend_side="right", annotation_legend_side="right",
     padding=unit(c(5,5,5,30), "mm"))
dev.off()

png(file.path(OUTDIR, "Fig2B_heatmap_ordered_by_class.png"),
    width=12, height=9, units="in", res=300)
draw(ht, heatmap_legend_side="right", annotation_legend_side="right",
     padding=unit(c(5,5,5,30), "mm"))
dev.off()
cat("Heatmap saved.\n")

# ============================================================
# PANEL C — Bar chart intra vs inter correlation
# Anthracyclines first, ordered by fold enrichment descending
# ============================================================
cat("Building Panel C...\n")

# Use confirmed stats
stats_plot <- copy(stats)
stats_plot[, class := factor(class,
  levels=c("Anthracyclines", "HDAC inhibitors", "EGFR inhibitors",
           "BCR-ABL inhibitors", "PARP inhibitors", "Pyrimidine analogues"))]

# Reshape to long
long_stats <- rbind(
  data.frame(class=stats_plot$class, type="Within-class",
             rho=stats_plot$mean_intra, n=stats_plot$n_drugs),
  data.frame(class=stats_plot$class, type="Between-class",
             rho=stats_plot$mean_inter, n=stats_plot$n_drugs)
)
long_stats$type <- factor(long_stats$type,
                           levels=c("Within-class", "Between-class"))

# Significance stars
stats_plot[, stars := ifelse(wilcox_p < 0.001, "***",
                      ifelse(wilcox_p < 0.01,  "**",
                      ifelse(wilcox_p < 0.05,  "*", "ns")))]

# Max rho per class for star placement
stats_plot[, max_rho := pmax(mean_intra, mean_inter) + 0.02]

bar_cols <- c("Within-class"="#B2182B", "Between-class"="#BABABA")

p_bar <- ggplot(long_stats, aes(x=class, y=rho, fill=type)) +
  geom_col(position=position_dodge(0.7), width=0.65) +
  geom_text(data=stats_plot,
            aes(x=class, y=max_rho,
                label=paste0(stars, "\n", round(fold_enrich,1), "x")),
            inherit.aes=FALSE,
            size=3.2, hjust=0.5, color="gray20", lineheight=0.9) +
  scale_fill_manual(values=bar_cols, name=NULL) +
  scale_y_continuous(limits=c(0, 0.55), expand=c(0,0),
                     breaks=seq(0, 0.5, 0.1)) +
  labs(title="Within-class vs between-class signature similarity",
       x=NULL,
       y="Mean Spearman ρ") +
  theme_classic(base_size=11) +
  theme(axis.text.x     = element_text(angle=35, hjust=1, size=9.5),
        legend.position = "top",
        legend.text     = element_text(size=9.5),
        plot.title      = element_text(face="bold", size=11, hjust=0.5),
        panel.grid.major.y = element_line(color="gray92", linewidth=0.4))

ggsave(file.path(OUTDIR, "Fig2C_class_coherence_bar.pdf"),
       p_bar, width=8, height=5, useDingbats=FALSE)
ggsave(file.path(OUTDIR, "Fig2C_class_coherence_bar.png"),
       p_bar, width=8, height=5, dpi=300)
cat("Panel C saved.\n")

# ============================================================
# COMBINED FIGURE 2
# Layout: A (UMAP) | B (heatmap) on top
#         C (bar chart) on bottom spanning full width
# ============================================================
cat("Combining Figure 2...\n")

# Note: ComplexHeatmap cannot be combined with patchwork directly
# Save as separate panels — combine in Illustrator or with cowplot

# For a patchwork-compatible version use PNG of heatmap as inset
library(png)
library(grid)

hm_png <- readPNG(file.path(OUTDIR, "Fig2B_heatmap_ordered_by_class.png"))
hm_grob <- rasterGrob(hm_png, interpolate=TRUE)
p_heatmap_wrap <- ggplot() +
  annotation_custom(hm_grob, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
  theme_void()

fig2 <- (p_umap | p_heatmap_wrap) / p_bar +
  plot_layout(heights=c(1.2, 1)) +
  plot_annotation(
    tag_levels="A",
    theme=theme(plot.tag=element_text(face="bold", size=14))
  )

ggsave(file.path(OUTDIR, "Fig2_combined.pdf"),
       fig2, width=14, height=11, useDingbats=FALSE)
ggsave(file.path(OUTDIR, "Fig2_combined.png"),
       fig2, width=14, height=11, dpi=300)
cat("Figure 2 combined saved.\n")

# ============================================================
# PRINT RESULTS FOR MANUSCRIPT
# ============================================================
cat("\n=== MANUSCRIPT NUMBERS ===\n")
cat(sprintf("Global: intra=%.3f | inter=%.3f | fold=%.2fx | p=%s\n",
            0.222, 0.088, 2.54, "5.09e-13"))
cat("\nPer-class results (ordered as in manuscript):\n")
for(i in 1:nrow(stats_plot)) {
  cat(sprintf("%-22s n=%d | intra=%.3f | inter=%.3f | fold=%.2fx | p=%s | %s\n",
              stats_plot$class[i], stats_plot$n_drugs[i],
              stats_plot$mean_intra[i], stats_plot$mean_inter[i],
              stats_plot$fold_enrich[i],
              formatC(stats_plot$wilcox_p[i], format="e", digits=2),
              stats_plot$stars[i]))
}

cat("\nAll figures saved to:", OUTDIR, "\n")
cat("Done.\n")

# Supplementary -----------------------------------------------------------

library(data.table)

BASE   <- "~/Desktop/MetaDrugSig_local"
OUTDIR <- file.path(BASE, "figures")

sim     <- readRDS(file.path(BASE, "spearman_similarity_v3.rds"))
meta    <- fread(file.path(BASE, "drug_metadata_v3.csv"))
meta_t1 <- meta[tier == 1]
supp    <- fread(file.path(OUTDIR, "SupplTable_matched_drug_sets.csv"))

# Class definitions
class_map <- c(
  "Anthracyclines and related substances"                              = "Anthracyclines",
  "Histone deacetylase (HDAC) inhibitors"                              = "HDAC inhibitors",
  "Epidermal growth factor receptor (EGFR) tyrosine kinase inhibitors" = "EGFR inhibitors",
  "BCR-ABL tyrosine kinase inhibitors"                                 = "BCR-ABL inhibitors",
  "Poly (ADP-ribose) polymerase (PARP) inhibitors"                     = "PARP inhibitors",
  "Pyrimidine analogues"                                               = "Pyrimidine analogues"
)

meta_t1[, drug_class := class_map[ATC_L4_Name]]
class_drugs <- meta_t1[!is.na(drug_class) & DrugBank_ID %in% rownames(sim)]

# For each drug compute:
# 1. Mean intra-class correlation (excluding self)
# 2. Mean inter-class correlation
# 3. n studies contributing to signature

results <- list()

for (db_id in class_drugs$DrugBank_ID) {
  cls      <- class_drugs$drug_class[class_drugs$DrugBank_ID == db_id]
  cls_ids  <- class_drugs$DrugBank_ID[class_drugs$drug_class == cls]
  cls_ids  <- cls_ids[cls_ids != db_id & cls_ids %in% rownames(sim)]
  other_ids <- class_drugs$DrugBank_ID[class_drugs$drug_class != cls &
                                         class_drugs$DrugBank_ID %in% rownames(sim)]
  
  mean_intra <- if(length(cls_ids) > 0) round(mean(sim[db_id, cls_ids], na.rm=TRUE), 3) else NA
  mean_inter <- if(length(other_ids) > 0) round(mean(sim[db_id, other_ids], na.rm=TRUE), 3) else NA
  
  n_studies  <- meta_t1$n_comparisons[meta_t1$DrugBank_ID == db_id]
  
  results[[db_id]] <- data.table(
    Class              = cls,
    DrugBank_ID        = db_id,
    Drug_Name          = meta_t1$DrugBank_Name[meta_t1$DrugBank_ID == db_id],
    N_GEO_Studies      = n_studies,
    Mean_Intra_Class_r = mean_intra,
    Mean_Inter_Class_r = mean_inter,
    Fold_Enrichment    = round(mean_intra / mean_inter, 2)
  )
}

supp_final <- rbindlist(results)
supp_final <- supp_final[order(Class, -Mean_Intra_Class_r)]

print(supp_final)
cat("\nTotal drugs:", nrow(supp_final), "\n")

fwrite(supp_final,
       file.path(OUTDIR, "SupplementaryTable_BenchmarkDrugs.csv"))
cat("Saved.\n")
