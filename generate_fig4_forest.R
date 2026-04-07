# ============================================================
# MetaDrugSig Figure 4 — Doxorubicin Forest Plots
# Panel A: SNRPEP4 — stable gene (I² = 0%, n=27 studies)
# Panel B: HAND1   — tissue-specific gene (I² = 99.6%, n=36 studies)
# Panel C: HAND1 meta-regression — tissue explains heterogeneity
# ============================================================

library(data.table)
library(ggplot2)
library(metafor)
library(patchwork)

BASE     <- "~/Desktop/Michigan State University/Chen Lab/Project_OCTAD2_drug_signature/Manuscript3/MetaDrugSig_app/data"
OUTDIR   <- "~/Desktop/MetaDrugSig_local/figures"
DRUG_ID  <- "DB00997"
DRUG_NM  <- "Doxorubicin"

comp_meta <- fread(file.path(BASE, "comparison_metadata_v3.csv"))

# ============================================================
# LOAD PER-GSE DATA FOR A GENE
# ============================================================
load_per_gse <- function(drug_name, gene_name) {
  drug_dir <- file.path(BASE, "per_GSE", tolower(drug_name))
  if (!dir.exists(drug_dir)) {
    dirs  <- list.dirs(file.path(BASE,"per_GSE"), recursive=FALSE, full.names=FALSE)
    match <- dirs[grepl(tolower(drug_name), tolower(dirs))]
    if (length(match)==0) { cat("Directory not found for", drug_name, "\n"); return(NULL) }
    drug_dir <- file.path(BASE, "per_GSE", match[1])
  }
  files <- list.files(drug_dir, pattern="\\.csv$", full.names=TRUE)
  cat("Files found:", length(files), "\n")
  if (length(files)==0) return(NULL)

  results <- list()
  for (f in files) {
    fname     <- basename(f)
    gse_id    <- sub("DEGs_","",sub("_Cell_Line.*\\.csv$","",fname))
    cell_line <- gsub("_"," ",sub(".*Cell_Line__","",sub("\\.csv$","",fname)))
    d <- tryCatch(fread(f), error=function(e) NULL)
    if (is.null(d)||nrow(d)==0) next
    gr <- d[!is.na(d$gene) & d$gene==gene_name,]
    if (nrow(gr)==0) next
    gr  <- gr[1,]
    lfc <- as.numeric(gr$logFC)
    fst <- as.numeric(gr$F)
    fdr <- as.numeric(gr$FDR)
    if (is.na(lfc)||is.na(fst)||fst<=0) next
    se <- abs(lfc)/sqrt(fst)
    if (is.na(se)||se<=0||se>10) next
    cm      <- comp_meta[GSE_ID==gse_id & DrugBank_ID==DRUG_ID,]
    tissue  <- if (nrow(cm)>0&&!is.na(cm$tissue_category[1])) cm$tissue_category[1] else "Unknown"
    dose    <- if (nrow(cm)>0&&!is.na(cm$dose[1]))     cm$dose[1]     else NA_character_
    dur     <- if (nrow(cm)>0&&!is.na(cm$duration[1])) cm$duration[1] else NA_character_
    results[[fname]] <- data.frame(
      gse_id=gse_id, cell_line=cell_line, tissue=tissue,
      dose_raw=ifelse(is.na(dose),"NA",dose),
      duration_raw=ifelse(is.na(dur),"NA",dur),
      logFC=lfc, SE=se, FDR=fdr,
      label=paste0(gse_id," | ",cell_line),
      stringsAsFactors=FALSE)
  }
  if (length(results)==0) return(NULL)
  do.call(rbind, results)
}

# ============================================================
# BUILD FOREST PLOT — color by tissue
# ============================================================
build_forest <- function(gene_data, gene_name, drug_name,
                         color_by_tissue=FALSE, title_suffix="") {
  if (is.null(gene_data)||nrow(gene_data)<3) {
    cat("Insufficient data for", gene_name, "\n"); return(NULL)
  }
  meta_fit <- tryCatch(
    rma(yi=logFC, sei=SE, data=gene_data, method="REML"),
    error=function(e) { cat("rma failed:", e$message, "\n"); NULL })
  if (is.null(meta_fit)) return(NULL)

  pooled <- list(
    est   = as.numeric(meta_fit$beta),
    ci_lb = meta_fit$ci.lb,
    ci_ub = meta_fit$ci.ub,
    pval  = meta_fit$pval,
    I2    = meta_fit$I2,
    tau2  = meta_fit$tau2
  )

  gene_data$weight     <- 1/(gene_data$SE^2 + meta_fit$tau2)
  gene_data$weight_pct <- 100*gene_data$weight/sum(gene_data$weight)
  gene_data$ci_lower   <- gene_data$logFC - 1.96*gene_data$SE
  gene_data$ci_upper   <- gene_data$logFC + 1.96*gene_data$SE
  gene_data$significant <- !is.na(gene_data$FDR) & gene_data$FDR < 0.05
  gene_data <- gene_data[order(gene_data$logFC),]
  gene_data$row_id     <- seq_len(nrow(gene_data))
  gene_data$effect_str <- sprintf("%.2f [%.2f, %.2f]",
                                   gene_data$logFC,
                                   gene_data$ci_lower,
                                   gene_data$ci_upper)
  gene_data$weight_str <- sprintf("%.1f%%", gene_data$weight_pct)
  gene_data$is_summary <- FALSE

  n_studies <- nrow(gene_data); y_max <- n_studies

  summary_row <- data.frame(
    gse_id="", cell_line="", tissue="Pooled",
    dose_raw="", duration_raw="",
    logFC=pooled$est, SE=as.numeric(meta_fit$se), FDR=pooled$pval,
    label="Overall (RE Model)",
    weight=NA_real_, weight_pct=100,
    ci_lower=pooled$ci_lb, ci_upper=pooled$ci_ub,
    significant=pooled$pval<0.05, row_id=0L,
    effect_str=sprintf("%.2f [%.2f, %.2f]",
                        pooled$est,pooled$ci_lb,pooled$ci_ub),
    weight_str="", is_summary=TRUE,
    stringsAsFactors=FALSE)

  keep    <- names(summary_row)
  plot_df <- rbind(gene_data[,intersect(keep,names(gene_data))], summary_row)
  plot_df <- plot_df[order(-plot_df$row_id),]

  x_min <- min(c(plot_df$ci_lower,-0.3), na.rm=TRUE) - 0.2
  x_max <- max(c(plot_df$ci_upper,  0.3), na.rm=TRUE) + 0.4

  # Tissue color palette
  tissue_cols <- c(
    "Heart"         = "#B2182B",
    "Breast"        = "#E08214",
    "Lung"          = "#FDB863",
    "Skin"          = "#B2ABD2",
    "Intestine"     = "#5AAE61",
    "Kidney"        = "#1B7837",
    "Brain/CNS"     = "#4393C3",
    "Blood/Immune"  = "#2166AC",
    "Liver"         = "#8C510A",
    "Adrenal gland" = "#DFC27D",
    "Cervix"        = "#DE77AE",
    "Eye/Retina"    = "#7FBC41",
    "Pooled"        = "black",
    "Unknown"       = "gray60"
  )

  point_color <- if (color_by_tissue) {
    plot_df$tissue_col <- tissue_cols[plot_df$tissue]
    plot_df$tissue_col[is.na(plot_df$tissue_col)] <- "gray60"
    "tissue"
  } else {
    "steelblue4"
  }

  # Left panel — study labels
  p_left <- ggplot(plot_df, aes(y=row_id)) +
    geom_text(aes(x=0, label=label), hjust=0, size=2.2,
              fontface=ifelse(plot_df$is_summary,"bold","plain")) +
    annotate("text", x=0, y=y_max+1.5,
             label="Study | Cell line", hjust=0, fontface="bold", size=2.6) +
    scale_x_continuous(limits=c(-0.5,10), expand=c(0,0)) +
    scale_y_continuous(limits=c(-2, y_max+2), expand=c(0,0)) +
    theme_void() + theme(plot.margin=margin(5,0,5,5))

  # Forest panel
  if (color_by_tissue) {
    study_df <- plot_df[!plot_df$is_summary,]
    study_df$tissue_col <- tissue_cols[study_df$tissue]
    study_df$tissue_col[is.na(study_df$tissue_col)] <- "gray60"

    p_forest <- ggplot(plot_df, aes(y=row_id)) +
      geom_vline(xintercept=0, linetype="dashed", color="gray50", linewidth=0.5) +
      geom_segment(data=plot_df[!plot_df$is_summary,],
                   aes(x=ci_lower,xend=ci_upper,y=row_id,yend=row_id),
                   color="gray50", linewidth=0.5) +
      geom_point(data=study_df,
                 aes(x=logFC, size=weight_pct),
                 color=study_df$tissue_col, shape=15) +
      geom_polygon(
        data=data.frame(
          x=c(pooled$ci_lb,pooled$est,pooled$ci_ub,pooled$est),
          y=c(0,0.38,0,-0.38)),
        aes(x=x,y=y), fill="black") +
      scale_size_continuous(range=c(1.5,5), guide="none") +
      scale_y_continuous(breaks=NULL,
                         expand=expansion(add=c(0.5,0.5)),
                         limits=c(-2,y_max+2)) +
      scale_x_continuous(limits=c(x_min,x_max)) +
      labs(x="log\u2082 fold change", y=NULL) +
      theme_classic(base_size=10) +
      theme(axis.line.y=element_blank(),
            axis.ticks.y=element_blank(),
            plot.margin=margin(5,5,5,5))

    # Tissue legend
    present_tissues <- unique(study_df$tissue[study_df$tissue!="Unknown"])
    legend_df <- data.frame(
      tissue=present_tissues,
      color=tissue_cols[present_tissues],
      x=1, y=seq_along(present_tissues))
    legend_df <- legend_df[!is.na(legend_df$color),]

    p_legend <- ggplot(legend_df, aes(y=y)) +
      geom_point(aes(x=0), color=legend_df$color, shape=15, size=3) +
      geom_text(aes(x=0.15, label=tissue), hjust=0, size=2.8) +
      scale_x_continuous(limits=c(-0.1,2), expand=c(0,0)) +
      scale_y_continuous(limits=c(0, max(legend_df$y)+1), expand=c(0,0)) +
      theme_void() +
      labs(title="Tissue") +
      theme(plot.title=element_text(size=8,face="bold"))

  } else {
    p_forest <- ggplot(plot_df, aes(y=row_id)) +
      geom_vline(xintercept=0, linetype="dashed", color="gray50", linewidth=0.5) +
      geom_segment(data=plot_df[!plot_df$is_summary,],
                   aes(x=ci_lower,xend=ci_upper,y=row_id,yend=row_id),
                   color="gray40", linewidth=0.5) +
      geom_point(data=plot_df[!plot_df$is_summary,],
                 aes(x=logFC, size=weight_pct, shape=significant),
                 color="steelblue4") +
      geom_polygon(
        data=data.frame(
          x=c(pooled$ci_lb,pooled$est,pooled$ci_ub,pooled$est),
          y=c(0,0.38,0,-0.38)),
        aes(x=x,y=y), fill="black") +
      scale_size_continuous(range=c(1.5,5), guide="none") +
      scale_shape_manual(values=c("TRUE"=15,"FALSE"=0),
                         labels=c("TRUE"="FDR<0.05","FALSE"="ns"),
                         name=NULL) +
      scale_y_continuous(breaks=NULL,
                         expand=expansion(add=c(0.5,0.5)),
                         limits=c(-2,y_max+2)) +
      scale_x_continuous(limits=c(x_min,x_max)) +
      labs(x="log\u2082 fold change", y=NULL) +
      theme_classic(base_size=10) +
      theme(axis.line.y=element_blank(),
            axis.ticks.y=element_blank(),
            legend.position="bottom",
            plot.margin=margin(5,5,5,5))
  }

  # Right panel — effect sizes
  p_right <- ggplot(plot_df, aes(y=row_id)) +
    geom_text(aes(x=0,   label=effect_str), hjust=0, size=2.2,
              fontface=ifelse(plot_df$is_summary,"bold","plain")) +
    geom_text(aes(x=2.7, label=weight_str), hjust=1, size=2.2) +
    annotate("text", x=0,   y=y_max+1.5,
             label="logFC [95% CI]", hjust=0, fontface="bold", size=2.6) +
    annotate("text", x=2.7, y=y_max+1.5,
             label="Weight", hjust=1, fontface="bold", size=2.6) +
    scale_x_continuous(limits=c(-0.5,3.2), expand=c(0,0)) +
    scale_y_continuous(limits=c(-2,y_max+2), expand=c(0,0)) +
    theme_void() + theme(plot.margin=margin(5,5,5,0))

  subtitle <- sprintf(
    "n=%d studies | I\u00b2=%.1f%% | \u03c4\u00b2=%.3f | Pooled logFC=%.3f [%.3f, %.3f] p=%.2e",
    n_studies, pooled$I2, pooled$tau2,
    pooled$est, pooled$ci_lb, pooled$ci_ub, pooled$pval)

  if (color_by_tissue) {
    combined <- (p_left + p_forest + p_right + p_legend) +
      plot_layout(widths=c(3,1.5,1.1,0.8))
  } else {
    combined <- (p_left + p_forest + p_right) +
      plot_layout(widths=c(3,1.5,1.1))
  }

  combined + plot_annotation(
    title    = sprintf("%s effect on %s%s", drug_name, gene_name, title_suffix),
    subtitle = subtitle,
    caption  = "SE approximated from edgeR F-statistic. Square size proportional to study weight.",
    theme=theme(
      plot.title    = element_text(face="bold", size=12),
      plot.subtitle = element_text(size=8, color="gray40"),
      plot.caption  = element_text(size=7, hjust=0, color="gray40")))
}

# ============================================================
# LOAD DATA
# ============================================================
cat("Loading TOP2B data...\n")
gd_top2b <- load_per_gse(DRUG_NM, "TOP2B")
cat("TOP2B studies found:", if(!is.null(gd_top2b)) nrow(gd_top2b) else 0, "\n")

cat("\nLoading HAND1 data...\n")
gd_hand1 <- load_per_gse(DRUG_NM, "HAND1")
cat("HAND1 studies found:", if(!is.null(gd_hand1)) nrow(gd_hand1) else 0, "\n")

# ============================================================
# GENERATE FIGURES
# ============================================================

# Panel A — TOP2B (known dox target, color by tissue to show cardiac enrichment)
cat("\nGenerating Panel A — TOP2B...\n")
p_top2b <- build_forest(gd_top2b, "TOP2B", DRUG_NM,
                         color_by_tissue=TRUE,
                         title_suffix=" — canonical cardiotoxicity target (I\u00b2 = 97.3%)")

# Panel B — HAND1 (cardiac-specific, color by tissue)
cat("Generating Panel B — HAND1...\n")
p_hand1 <- build_forest(gd_hand1, "HAND1", DRUG_NM,
                         color_by_tissue=TRUE,
                         title_suffix=" — cardiac-specific transcription factor (I\u00b2 = 99.6%)")

# ============================================================
# META-REGRESSION — HAND1 tissue moderator
# ============================================================
cat("\nRunning HAND1 meta-regression...\n")
if (!is.null(gd_hand1) && nrow(gd_hand1) >= 10) {
  null_fit <- rma(yi=logFC, sei=SE, data=gd_hand1, method="REML")
  mod_fit  <- rma(yi=logFC, sei=SE, data=gd_hand1,
                  mods=~tissue, method="REML")

  cat("Null model: I2=", round(null_fit$I2,1),
      "tau2=", round(null_fit$tau2,3), "\n")
  cat("Tissue model: I2=", round(mod_fit$I2,1),
      "tau2=", round(mod_fit$tau2,3), "\n")
  r2 <- max(0, (null_fit$tau2 - mod_fit$tau2)/null_fit$tau2)*100
  delta_aic <- null_fit$fit.stats["AIC","REML"] - mod_fit$fit.stats["AIC","REML"]
  cat("R2:", round(r2,1), "%\n")
  cat("ΔAIC:", round(delta_aic,2), "\n")
  cat("QM p:", mod_fit$QMp, "\n")

  # Meta-regression coefficient plot (Panel C)
  coef_df <- as.data.frame(coef(summary(mod_fit)))
  coef_df$tissue <- gsub("^tissue","",rownames(coef_df))
  coef_df$tissue[1] <- paste0(coef_df$tissue[1]," (ref)")
  coef_df <- coef_df[order(coef_df$estimate),]
  coef_df$tissue <- factor(coef_df$tissue, levels=coef_df$tissue)

  p_metareg <- ggplot(coef_df, aes(x=estimate, y=tissue)) +
    geom_vline(xintercept=0, linetype="dashed", color="gray50", linewidth=0.5) +
    geom_errorbarh(aes(xmin=ci.lb, xmax=ci.ub),
                   height=0.3, color="gray40", linewidth=0.6) +
    geom_point(aes(color=estimate > 0), size=3) +
    scale_color_manual(values=c("TRUE"="#B2182B","FALSE"="#2166AC"), guide="none") +
    labs(
      title = "HAND1 meta-regression — tissue moderator",
      subtitle = sprintf("R\u00b2=%.1f%% | \u0394AIC=%.1f | QM p=%.2e | I\u00b2 null=%.1f%% \u2192 residual=%.1f%%",
                         r2, delta_aic, mod_fit$QMp,
                         null_fit$I2, mod_fit$I2),
      x = "Tissue effect estimate (log\u2082 FC relative to reference)",
      y = NULL) +
    theme_classic(base_size=10) +
    theme(plot.title=element_text(face="bold",size=11),
          plot.subtitle=element_text(size=8,color="gray40"),
          axis.text.y=element_text(size=9))
}

# ============================================================
# SAVE ALL PANELS
# ============================================================
cat("\nSaving figures...\n")

if (!is.null(p_snrpep4)) {
  n_a <- if(!is.null(gd_snrpep4)) nrow(gd_snrpep4) else 10
  ggsave(file.path(OUTDIR, "Fig4A_SNRPEP4_forest.png"),
         p_snrpep4,
         width=13, height=max(5, n_a*0.28+3), dpi=300)
  ggsave(file.path(OUTDIR, "Fig4A_SNRPEP4_forest.pdf"),
         p_snrpep4,
         width=13, height=max(5, n_a*0.28+3), useDingbats=FALSE)
  cat("Panel A saved.\n")
}

if (!is.null(p_hand1)) {
  n_b <- if(!is.null(gd_hand1)) nrow(gd_hand1) else 10
  ggsave(file.path(OUTDIR, "Fig4B_HAND1_forest.png"),
         p_hand1,
         width=14, height=max(5, n_b*0.28+3), dpi=300)
  ggsave(file.path(OUTDIR, "Fig4B_HAND1_forest.pdf"),
         p_hand1,
         width=14, height=max(5, n_b*0.28+3), useDingbats=FALSE)
  cat("Panel B saved.\n")
}

if (exists("p_metareg")) {
  ggsave(file.path(OUTDIR, "Fig4C_HAND1_metareg.png"),
         p_metareg, width=8, height=6, dpi=300)
  ggsave(file.path(OUTDIR, "Fig4C_HAND1_metareg.pdf"),
         p_metareg, width=8, height=6, useDingbats=FALSE)
  cat("Panel C saved.\n")
}

cat("\nAll figures saved to:", OUTDIR, "\n")
