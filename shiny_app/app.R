library(shiny)
library(data.table)
library(dplyr)
library(ggplot2)
library(plotly)
library(metafor)
library(patchwork)
library(msigdbr)

# ============================================================
# PATHS
# ============================================================
LOGFC_PATH  <- "../data/drug_signatures_logFC.csv"
FDR_PATH    <- "../data/drug_signatures_FDR.csv"
SE_PATH     <- "../data/drug_signatures_SE.csv"
META_PATH   <- "../data/drug_metadata_v3.csv"
COMP_META   <- "../data/comparison_metadata_v3.csv"
PERGSE_ROOT <- "../data/per_GSE"
ZENODO_URL  <- "https://doi.org/10.5281/zenodo.19154955"

MIN_STUDIES_TOTAL    <- 3
MIN_STUDIES_METAREG  <- 10

# ============================================================
# LOAD DATA
# ============================================================
cat("Loading metadata...\n")
meta <- fread(META_PATH)
meta <- meta[has_signature == TRUE]

drug_choices <- setNames(meta$DrugBank_ID,
                         paste0(meta$DrugBank_Name, " (", meta$DrugBank_ID, ")"))
drug_choices <- drug_choices[order(names(drug_choices))]

cat("Loading comparison metadata...\n")
comp_meta <- fread(COMP_META)

cat("Loading logFC matrix...\n")
logfc_dt  <- fread(LOGFC_PATH)
genes     <- logfc_dt$gene
logfc_dt[, gene := NULL]
logfc_mat <- as.matrix(logfc_dt)
rownames(logfc_mat) <- genes

cat("Loading FDR matrix...\n")
fdr_dt  <- fread(FDR_PATH)
fdr_dt[, gene := NULL]
fdr_mat <- as.matrix(fdr_dt)
rownames(fdr_mat) <- genes

cat("Loading SE matrix...\n")
se_dt  <- fread(SE_PATH)
se_dt[, gene := NULL]
se_mat <- as.matrix(se_dt)
rownames(se_mat) <- genes

tier1_ids <- meta[tier == 1, DrugBank_ID]
cat(sprintf("Total drugs: %d | Tier 1: %d | Tier 2: %d\n",
            nrow(meta), sum(meta$tier==1), sum(meta$tier==2)))

# ============================================================
# LOAD MSIGDB
# ============================================================
cat("Loading MSigDB gene sets...\n")
msig_h     <- msigdbr(species="Homo sapiens", collection="H")
msig_kegg  <- msigdbr(species="Homo sapiens", collection="C2", subcollection="CP:KEGG_MEDICUS")
msig_react <- msigdbr(species="Homo sapiens", collection="C2", subcollection="CP:REACTOME")
msig_gobp  <- msigdbr(species="Homo sapiens", collection="C5", subcollection="GO:BP")

build_pathway_list <- function(df) split(df$gene_symbol, df$gs_name)

pathway_db <- list(
  "Hallmark (H)"          = build_pathway_list(msig_h),
  "KEGG"                  = build_pathway_list(msig_kegg),
  "Reactome"              = build_pathway_list(msig_react),
  "GO Biological Process" = build_pathway_list(msig_gobp)
)
collection_choices <- names(pathway_db)
cat("Ready.\n")

# ============================================================
# HELPERS
# ============================================================
get_drug_sig <- function(drug_id) {
  if (!drug_id %in% colnames(logfc_mat)) return(NULL)
  lfc <- logfc_mat[, drug_id]
  fdr <- if (drug_id %in% colnames(fdr_mat)) fdr_mat[, drug_id] else rep(NA_real_, length(lfc))
  se  <- if (drug_id %in% colnames(se_mat))  se_mat[, drug_id]  else rep(NA_real_, length(lfc))
  df  <- data.frame(gene=names(lfc), logFC=lfc, FDR=fdr, SE=se)
  df[!is.na(df$logFC), ]
}

get_drug_study_summary <- function(drug_id) {
  # Returns tissue/dose/duration breakdown for a drug from comparison metadata
  drug_name <- meta$DrugBank_Name[meta$DrugBank_ID == drug_id][1]
  cm <- comp_meta[DrugBank_ID == drug_id]
  if (nrow(cm) == 0) return(NULL)
  list(
    n_total   = nrow(cm),
    n_gses    = length(unique(cm$GSE_ID)),
    tier      = meta$tier[meta$DrugBank_ID == drug_id][1],
    tissue    = cm[!is.na(tissue_category) & tissue_category != "",
                   .N, by=tissue_category][order(-N)],
    dose      = cm[!is.na(dose_category) & dose_category != "",
                   .N, by=dose_category][order(-N)],
    duration  = cm[!is.na(duration_category) & duration_category != "",
                   .N, by=duration_category][order(-N)],
    drug_name = drug_name
  )
}

load_per_gse <- function(drug_name, gene_name) {
  drug_dir <- file.path(PERGSE_ROOT, tolower(drug_name))
  if (!dir.exists(drug_dir)) {
    dirs  <- list.dirs(PERGSE_ROOT, recursive=FALSE, full.names=FALSE)
    match <- dirs[grepl(tolower(drug_name), tolower(dirs))]
    if (length(match)==0) return(NULL)
    drug_dir <- file.path(PERGSE_ROOT, match[1])
  }
  files <- list.files(drug_dir, pattern="\\.csv$", full.names=TRUE)
  if (length(files)==0) return(NULL)
  drug_id <- meta$DrugBank_ID[tolower(meta$DrugBank_Name)==tolower(drug_name)][1]
  results <- list()
  for (f in files) {
    fname     <- basename(f)
    gse_id    <- sub("DEGs_","",sub("_Cell_Line.*\\.csv$","",fname))
    cell_line <- gsub("_"," ",sub(".*Cell_Line__","",sub("\\.csv$","",fname)))
    d <- tryCatch(fread(f), error=function(e) NULL)
    if (is.null(d)||nrow(d)==0) next
    gr <- d[!is.na(d$gene) & d$gene==gene_name,]
    if (nrow(gr)==0) next
    gr <- gr[1,]
    lfc <- as.numeric(gr$logFC); fst <- as.numeric(gr$F); fdr <- as.numeric(gr$FDR)
    if (is.na(lfc)||is.na(fst)||fst<=0) next
    se <- abs(lfc)/sqrt(fst)
    if (is.na(se)||se<=0||se>10) next
    cm       <- comp_meta[GSE_ID==gse_id & DrugBank_ID==drug_id,]
    tissue   <- if (nrow(cm)>0&&!is.na(cm$tissue[1]))   cm$tissue[1]   else "Unknown"
    dose     <- if (nrow(cm)>0&&!is.na(cm$dose[1]))     cm$dose[1]     else NA_character_
    duration <- if (nrow(cm)>0&&!is.na(cm$duration[1])) cm$duration[1] else NA_character_
    lineage  <- if (nrow(cm)>0&&!is.na(cm$lineage[1]))  cm$lineage[1]  else "Unknown"
    dose_cat <- tryCatch({
      d_num <- as.numeric(gsub("[^0-9.]","",dose))
      if (is.na(d_num)) "Unknown" else if (d_num<0.1) "Low" else if (d_num<=1) "Medium" else "High"
    }, error=function(e) "Unknown")
    dur_cat <- tryCatch({
      d_h <- as.numeric(gsub("[^0-9.]","",duration))
      if (!is.na(d_h)&&grepl("day|Day",duration)) d_h <- d_h*24
      if (is.na(d_h)) "Unknown" else if (d_h<=24) "Acute (<=24h)" else "Chronic (>24h)"
    }, error=function(e) "Unknown")
    results[[fname]] <- data.frame(
      gse_id=gse_id, cell_line=cell_line, tissue=tissue, lineage=lineage,
      dose_raw=ifelse(is.na(dose),"NA",dose), dose_cat=dose_cat,
      duration_raw=ifelse(is.na(duration),"NA",duration), duration_cat=dur_cat,
      logFC=lfc, SE=se, FDR=fdr,
      label=paste0(gse_id," | ",cell_line," | ",
                   ifelse(is.na(dose),"NA",dose)," | ",
                   ifelse(is.na(duration),"NA",duration)),
      stringsAsFactors=FALSE)
  }
  if (length(results)==0) return(NULL)
  do.call(rbind, results)
}

run_meta_regression <- function(gene_data, moderators) {
  n <- nrow(gene_data); results <- list()
  null_fit <- tryCatch(rma(yi=logFC,sei=SE,data=gene_data,method="REML"),error=function(e) NULL)
  if (is.null(null_fit)) return(list(error="Meta-analysis failed on null model."))
  results$null <- list(I2=null_fit$I2,tau2=null_fit$tau2,AIC=null_fit$fit.stats["AIC","REML"])
  results$moderators <- list()
  for (mod in moderators) {
    col <- switch(mod,"Tissue"="tissue","Dose"="dose_cat","Duration"="duration_cat","Lineage"="lineage",NULL)
    if (is.null(col)||!col %in% names(gene_data)) next
    vals <- gene_data[[col]]; vals[is.na(vals)|vals=="NA"|vals=="Unknown"] <- "Unknown"
    gene_data[[col]] <- vals
    n_levels <- length(unique(vals[vals!="Unknown"]))
    n_per_level <- table(vals); warnings <- character(0)
    if (n<MIN_STUDIES_METAREG)
      warnings <- c(warnings,sprintf("Only %d studies — meta-regression requires >= %d.",n,MIN_STUDIES_METAREG))
    if (n_levels<2) {
      warnings <- c(warnings,sprintf("%s has only %d informative level(s).",mod,n_levels))
      results$moderators[[mod]] <- list(warnings=warnings,skipped=TRUE); next
    }
    sparse_levels <- names(n_per_level)[n_per_level<3]
    if (length(sparse_levels)>0)
      warnings <- c(warnings,sprintf("Sparse levels (n<3): %s.",paste(sparse_levels,collapse=", ")))
    mod_fit <- tryCatch(
      rma(yi=logFC,sei=SE,data=gene_data,mods=as.formula(paste0("~ ",col)),method="REML"),
      error=function(e) NULL)
    if (is.null(mod_fit)) {
      results$moderators[[mod]] <- list(warnings=c(warnings,"Model fitting failed."),skipped=TRUE); next
    }
    r2 <- max(0,(null_fit$tau2-mod_fit$tau2)/null_fit$tau2)*100
    delta_aic <- null_fit$fit.stats["AIC","REML"]-mod_fit$fit.stats["AIC","REML"]
    coef_df <- data.frame(
      Level=rownames(coef(summary(mod_fit))),
      Estimate=round(coef(summary(mod_fit))[,"estimate"],3),
      SE=round(coef(summary(mod_fit))[,"se"],3),
      CI_lower=round(coef(summary(mod_fit))[,"ci.lb"],3),
      CI_upper=round(coef(summary(mod_fit))[,"ci.ub"],3),
      p_value=formatC(coef(summary(mod_fit))[,"pval"],format="e",digits=2),
      stringsAsFactors=FALSE)
    coef_df$Level <- gsub(paste0("^",col),"",coef_df$Level)
    coef_df$Level[1] <- paste0(coef_df$Level[1]," (reference)")
    interpretation <- if (delta_aic>=2&&mod_fit$QMp<0.05)
      paste0(mod," explains heterogeneity (ΔAIC=",round(delta_aic,1),", R²=",round(r2,1),"%, p=",formatC(mod_fit$QMp,format="e",digits=2),").")
    else if (delta_aic>=2) paste0(mod," improves fit (ΔAIC=",round(delta_aic,1),") but QM test not significant.")
    else paste0(mod," does not explain heterogeneity (ΔAIC=",round(delta_aic,1),", R²=",round(r2,1),"%).")
    results$moderators[[mod]] <- list(
      warnings=warnings,skipped=FALSE,R2=round(r2,1),delta_AIC=round(delta_aic,2),
      I2_null=round(null_fit$I2,1),I2_residual=round(mod_fit$I2,1),
      QM_stat=round(mod_fit$QM,2),QM_df=mod_fit$m,
      QM_pval=formatC(mod_fit$QMp,format="e",digits=2),
      coef_df=coef_df,interpretation=interpretation,
      n_levels=n_levels,level_counts=as.data.frame(n_per_level))
  }
  results
}

build_forest_plot <- function(gene_data, gene_name, drug_name) {
  if (is.null(gene_data)||nrow(gene_data)<MIN_STUDIES_TOTAL) return(NULL)
  meta_fit <- tryCatch(rma(yi=logFC,sei=SE,data=gene_data,method="REML"),error=function(e) NULL)
  if (is.null(meta_fit)) return(NULL)
  pooled <- list(est=as.numeric(meta_fit$beta),ci_lb=meta_fit$ci.lb,ci_ub=meta_fit$ci.ub,
                 pval=meta_fit$pval,I2=meta_fit$I2,tau2=meta_fit$tau2)
  gene_data$weight     <- 1/(gene_data$SE^2+meta_fit$tau2)
  gene_data$weight_pct <- 100*gene_data$weight/sum(gene_data$weight)
  gene_data$ci_lower   <- gene_data$logFC-1.96*gene_data$SE
  gene_data$ci_upper   <- gene_data$logFC+1.96*gene_data$SE
  gene_data$significant <- !is.na(gene_data$FDR)&gene_data$FDR<0.05
  gene_data <- gene_data[order(gene_data$logFC),]
  gene_data$row_id     <- seq_len(nrow(gene_data))
  gene_data$effect_str <- sprintf("%.2f [%.2f, %.2f]",gene_data$logFC,gene_data$ci_lower,gene_data$ci_upper)
  gene_data$weight_str <- sprintf("%.1f%%",gene_data$weight_pct)
  gene_data$is_summary <- FALSE
  n_studies <- nrow(gene_data); y_max <- n_studies
  summary_row <- data.frame(
    gse_id="",cell_line="",tissue="Pooled",lineage="",
    dose_raw="",dose_cat="",duration_raw="",duration_cat="",
    logFC=pooled$est,SE=as.numeric(meta_fit$se),FDR=pooled$pval,
    label="Overall (RE Model)",weight=NA_real_,weight_pct=100,
    ci_lower=pooled$ci_lb,ci_upper=pooled$ci_ub,
    significant=pooled$pval<0.05,row_id=0L,
    effect_str=sprintf("%.2f [%.2f, %.2f]",pooled$est,pooled$ci_lb,pooled$ci_ub),
    weight_str="",is_summary=TRUE,stringsAsFactors=FALSE)
  keep <- names(summary_row)
  plot_df <- rbind(gene_data[,intersect(keep,names(gene_data))],summary_row)
  plot_df  <- plot_df[order(-plot_df$row_id),]
  x_min <- min(c(plot_df$ci_lower,-0.5),na.rm=TRUE)-0.3
  x_max <- max(c(plot_df$ci_upper, 0.5),na.rm=TRUE)+0.5
  p_left <- ggplot(plot_df,aes(y=row_id))+
    geom_text(aes(x=0,label=label),hjust=0,size=2.3,fontface=ifelse(plot_df$is_summary,"bold","plain"))+
    annotate("text",x=0,y=y_max+1.3,label="Study | Cell line | Dose | Duration",hjust=0,fontface="bold",size=2.7)+
    scale_x_continuous(limits=c(-0.5,10),expand=c(0,0))+
    scale_y_continuous(limits=c(-2,y_max+2),expand=c(0,0))+
    theme_void()+theme(plot.margin=margin(5,0,5,5))
  p_forest <- ggplot(plot_df,aes(y=row_id))+
    geom_vline(xintercept=0,linetype="dashed",color="gray50",linewidth=0.5)+
    geom_segment(data=plot_df[!plot_df$is_summary,],
                 aes(x=ci_lower,xend=ci_upper,y=row_id,yend=row_id),color="gray40",linewidth=0.6)+
    geom_point(data=plot_df[!plot_df$is_summary,],
               aes(x=logFC,size=weight_pct,shape=significant),color="steelblue4")+
    geom_polygon(data=data.frame(x=c(pooled$ci_lb,pooled$est,pooled$ci_ub,pooled$est),y=c(0,0.38,0,-0.38)),
                 aes(x=x,y=y),fill="black")+
    scale_size_continuous(range=c(1.5,5),guide="none")+
    scale_shape_manual(values=c("TRUE"=15,"FALSE"=0),labels=c("TRUE"="FDR<0.05","FALSE"="ns"),name=NULL)+
    scale_y_continuous(breaks=NULL,expand=expansion(add=c(0.5,0.5)),limits=c(-2,y_max+2))+
    scale_x_continuous(limits=c(x_min,x_max))+
    labs(x="log2 fold change",y=NULL)+
    theme_classic(base_size=10)+
    theme(axis.line.y=element_blank(),axis.ticks.y=element_blank(),
          legend.position="bottom",plot.margin=margin(5,5,5,5))
  p_right <- ggplot(plot_df,aes(y=row_id))+
    geom_text(aes(x=0,label=effect_str),hjust=0,size=2.3,fontface=ifelse(plot_df$is_summary,"bold","plain"))+
    geom_text(aes(x=2.7,label=weight_str),hjust=1,size=2.3)+
    annotate("text",x=0,  y=y_max+1.3,label="logFC [95% CI]",hjust=0,fontface="bold",size=2.7)+
    annotate("text",x=2.7,y=y_max+1.3,label="Weight",       hjust=1,fontface="bold",size=2.7)+
    scale_x_continuous(limits=c(-0.5,3.2),expand=c(0,0))+
    scale_y_continuous(limits=c(-2,y_max+2),expand=c(0,0))+
    theme_void()+theme(plot.margin=margin(5,5,5,0))
  p_left+p_forest+p_right+
    plot_layout(widths=c(3,1.5,1.1))+
    plot_annotation(
      title=sprintf("%s effect on %s",drug_name,gene_name),
      subtitle=sprintf("n=%d | I2=%.1f%% | tau2=%.3f | Pooled: %.3f [%.3f,%.3f] p=%.2e",
                       n_studies,pooled$I2,pooled$tau2,pooled$est,pooled$ci_lb,pooled$ci_ub,pooled$pval),
      caption="SE approximated from edgeR F-statistic. Square size proportional to study weight.",
      theme=theme(plot.title=element_text(face="bold",size=12),
                  plot.subtitle=element_text(size=8,color="gray40"),
                  plot.caption=element_text(size=7,hjust=0,color="gray40")))
}

# ============================================================
# RGES — canonical cmap_score_new (Chen Lab / OCTAD)
# Source: github.com/Bin-Chen-Lab/RGES/runRGESExample.R
# rank(-1 * logFC) so rank 1 = most upregulated
# Negative RGES = drug REVERSES query = candidate therapeutic
# Positive  RGES = drug MIMICS query
# ============================================================
cmap_score_new <- function(sig_up, sig_down, drug_signature, gamma=1) {
  num_genes <- nrow(drug_signature)
  drug_signature[,"rank"] <- rank(-1 * drug_signature[,"rank"], ties.method="random")
  up_tags_rank   <- merge(drug_signature, sig_up,   by.x="ids", by.y=1)
  down_tags_rank <- merge(drug_signature, sig_down, by.x="ids", by.y=1)
  up_tags_position   <- sort(up_tags_rank$rank)
  down_tags_position <- sort(down_tags_rank$rank)
  ks <- function(tags_position, num_tags) {
    a <- max(sapply(seq_len(num_tags), function(j) j/num_tags - tags_position[j]/num_genes))
    b <- max(sapply(seq_len(num_tags), function(j) tags_position[j]/num_genes - (j-1)/num_tags))
    if (a > b) a else -b
  }
  ks_up   <- if (length(up_tags_position)   > 1) ks(up_tags_position,   length(up_tags_position))   else 0
  ks_down <- if (length(down_tags_position) > 1) ks(down_tags_position, length(down_tags_position)) else 0
  if (ks_up == 0 & ks_down != 0) return(-ks_down)
  if (ks_up != 0 & ks_down == 0) return(gamma * ks_up)
  return(gamma * ks_up - ks_down)
}

# ============================================================
# UI
# ============================================================
ui <- navbarPage(
  title="MetaDrugSig",

  # ---- DRUG EXPLORER ----
  tabPanel("Drug Explorer",
    sidebarLayout(
      sidebarPanel(width=3,
        h4("Select drug"),
        selectizeInput("drug_select","Drug:",choices=NULL,
                       options=list(placeholder="Type drug name...")),
        hr(),
        h4("Forest plot"),
        p(style="font-size:11px;color:gray;",
          "Click a gene on the volcano, or type directly:"),
        selectizeInput("gene_direct","Gene:",choices=NULL,
                       options=list(placeholder="e.g. TOP2B, CDKN1A...")),
        actionButton("run_forest","Generate forest plot",
                     class="btn-primary btn-sm",width="100%"),
        hr(),
        tags$details(
          tags$summary(style="cursor:pointer;font-weight:bold;color:#2c3e50;",
                       "\u25b6 Test heterogeneity moderators"),
          br(),
          p(style="font-size:11px;color:gray;",
            "REML meta-regression testing whether tissue, dose, or duration explain I\u00b2. Requires \u2265 10 studies."),
          checkboxGroupInput("metareg_mods","Moderators:",
                             choices=c("Tissue","Dose","Duration","Lineage"),
                             selected=c("Tissue","Dose","Duration")),
          actionButton("run_metareg","Run meta-regression",
                       class="btn-warning btn-sm",width="100%")
        ),
        hr(),
        numericInput("fdr_thresh","FDR threshold (volcano):",
                     value=0.05,min=0,max=1,step=0.01),
        hr(),
        uiOutput("drug_info_box"),
        hr(),
        downloadButton("dl_sig","Download signature"),
        br(),br(),
        downloadButton("dl_forest","Download forest PDF")
      ),
      mainPanel(width=9,
        tabsetPanel(id="main_tabs",

          # ---- TAB 1: STUDY OVERVIEW (NEW — shown first) ----
          tabPanel("Study overview",
            br(),
            uiOutput("study_overview_header"),
            br(),
            fluidRow(
              # Tissue distribution
              column(6,
                div(style="background:#f8f9fa;padding:12px;border-radius:6px;",
                  h5(style="margin:0 0 8px 0;color:#2c3e50;",
                     icon("flask"), " Tissue distribution"),
                  plotOutput("tissue_bar", height="260px")
                )
              ),
              # Dose + Duration donuts
              column(6,
                div(style="background:#f8f9fa;padding:12px;border-radius:6px;margin-bottom:10px;",
                  h5(style="margin:0 0 8px 0;color:#2c3e50;",
                     icon("pills"), " Treatment dose"),
                  plotOutput("dose_donut", height="120px")
                ),
                div(style="background:#f8f9fa;padding:12px;border-radius:6px;",
                  h5(style="margin:0 0 8px 0;color:#2c3e50;",
                     icon("clock"), " Treatment duration"),
                  plotOutput("duration_donut", height="120px")
                )
              )
            ),
            br(),
            # Gene filter options
            div(style="background:#eef4ff;padding:12px;border-radius:6px;",
              h5(style="margin:0 0 10px 0;color:#2c3e50;",
                 icon("filter"), " Signature gene filters"),
              fluidRow(
                column(4,
                  selectInput("gene_filter_type","Gene selection:",
                    choices=c(
                      "All significant genes (FDR < threshold)" = "all",
                      "Consensus genes (low I\u00b2 < 30%)"     = "consensus",
                      "Context-dependent genes (I\u00b2 \u2265 75%)" = "context"
                    ), selected="all")
                ),
                column(4,
                  selectInput("context_filter","Context filter:",
                    choices=c("None — use full consensus" = "none",
                              "Tissue type" = "tissue",
                              "Dose category" = "dose",
                              "Duration category" = "duration"),
                    selected="none")
                ),
                column(4,
                  uiOutput("context_value_ui")
                )
              ),
              p(style="font-size:11px;color:gray;margin-top:6px;",
                "Gene filter applies to volcano plot and top genes tables. Context filter selects studies with matching metadata before recomputing the consensus.")
            )
          ),

          # ---- TAB 2: VOLCANO ----
          tabPanel("Volcano plot",
            br(),
            p(style="font-size:11px;color:gray;",
              "Hover for gene details. Click to generate forest plot."),
            plotlyOutput("volcano_plot", height="500px")
          ),

          # ---- TAB 3: FOREST PLOT ----
          tabPanel("Forest plot",
            br(),
            uiOutput("forest_status"),
            plotOutput("forest_plot", height="auto", width="100%")
          ),

          # ---- TAB 4: HETEROGENEITY MODERATORS ----
          tabPanel("Heterogeneity moderators",
            br(),
            uiOutput("metareg_output")
          ),

          # ---- TAB 5: TOP GENES ----
          tabPanel("Top genes",
            br(),
            fluidRow(
              column(6, h4("Top upregulated"),  DT::DTOutput("tbl_up")),
              column(6, h4("Top downregulated"), DT::DTOutput("tbl_down"))
            )
          )
        )
      )
    )
  ),

  # ---- RGES REPURPOSING ----
  tabPanel("RGES Repurposing",
    sidebarLayout(
      sidebarPanel(width=3,
        h4("Disease gene signature"),
        div(style="background:#fff8e1;padding:8px;border-radius:4px;margin-bottom:10px;",
          p(style="margin:0;font-size:11px;color:#b8860b;font-weight:bold;",
            "\u26a0 Convention: Negative RGES = reversal (candidate therapeutic). Positive RGES = mimicry.")),
        hr(),
        div(style="background:#eef4ff;padding:10px;border-radius:6px;margin-bottom:12px;",
          h5(style="margin:0 0 8px 0;color:#2c3e50;","MSigDB pathway"),
          selectInput("msig_collection","Collection:",
                      choices=c("\u2014 Select \u2014",collection_choices),
                      selected="\u2014 Select \u2014"),
          selectizeInput("msig_pathway","Pathway:",choices=NULL,
                         options=list(placeholder="Type to search...")),
          uiOutput("msig_preview"),
          selectInput("msig_direction","Load genes as:",
                      choices=c("Disease upregulated (UP in disease)"="up",
                                "Disease downregulated (DOWN in disease)"="down"),
                      selected="up"),
          actionButton("load_msig","Load pathway genes",
                       class="btn-info btn-sm",width="100%")
        ),
        p(style="text-align:center;color:gray;font-size:11px;",
          "\u2014 or enter genes manually \u2014"),
        br(),
        textAreaInput("sig_up",
          HTML("<b>Disease upregulated</b> (genes UP in disease):"),
          rows=4, placeholder="TP53\nMYC\n..."),
        textAreaInput("sig_down",
          HTML("<b>Disease downregulated</b> (genes DOWN in disease):"),
          rows=4, placeholder="KLHL24\nPRKN\n..."),
        hr(),
        div(style="background:#f0f0f0;padding:10px;border-radius:6px;margin-bottom:10px;",
          radioButtons("rges_tier","Drug signature tier:",
                       choices=c("Tier 1 only \u2014 REML meta-analytic (recommended)"="tier1",
                                 "All drugs \u2014 Tier 1 + Tier 2"="all"),
                       selected="tier1"),
          uiOutput("tier_count_label")
        ),
        hr(),
        actionButton("run_rges","Run RGES",class="btn-primary btn-block",width="100%"),
        hr(),
        downloadButton("dl_rges","Download full results CSV")
      ),
      mainPanel(width=9,
        conditionalPanel(condition="input.run_rges==0",
          div(style="padding:60px;text-align:center;color:gray;",
            h4("Enter a disease gene signature and click Run RGES"),
            div(style="background:#f8f8f8;padding:20px;border-radius:8px;display:inline-block;text-align:left;margin-top:16px;max-width:520px;",
              p(style="margin:6px 0;font-size:13px;",
                strong("\u2212 Negative RGES"),
                " = drug reverses the disease signature = candidate therapeutic"),
              p(style="margin:6px 0;font-size:13px;",
                strong("+ Positive RGES"),
                " = drug mimics the disease signature"),
              hr(),
              p(style="margin:4px 0;font-size:11px;color:#555;",
                strong("Tier 1 (recommended):"),
                " REML meta-analytic consensus signatures from \u22652 independent studies."),
              p(style="margin:4px 0;font-size:11px;color:#555;",
                strong("Tier 2:"),
                " Single-study edgeR signatures (noisier, use with caution)."),
              hr(),
              p(style="margin:4px 0;font-size:12px;color:#888;",
                strong("Example 1:"),
                " Load HALLMARK_G2M_CHECKPOINT as \u2018Disease upregulated\u2019",
                " \u2192 negative RGES = anti-proliferative drugs",
                " (CDK inhibitors, MEK inhibitors, microtubule agents)."),
              p(style="margin:4px 0;font-size:12px;color:#888;",
                strong("Example 2:"),
                " Enter proteostasis genes suppressed in doxorubicin cardiotoxicity as \u2018Disease downregulated\u2019",
                " \u2192 negative RGES = candidate cardioprotective drugs.")
            )
          )
        ),
        conditionalPanel(condition="input.run_rges>0",
          tabsetPanel(
            tabPanel("Top & Bottom 5",
              br(),
              uiOutput("rges_query_summary"),
              br(),
              fluidRow(
                column(6,
                  div(style="background:#eaf4ea;padding:12px;border-radius:6px;",
                    h4(style="color:#1a7a1a;margin-top:0;",
                       "\u2212 Top 5 Reversers (most negative RGES)"),
                    p(style="font-size:11px;color:gray;",
                      "Strongest transcriptional reversal \u2014 candidate therapeutics"),
                    DT::DTOutput("rges_top5"))),
                column(6,
                  div(style="background:#fdf0f0;padding:12px;border-radius:6px;",
                    h4(style="color:#a01010;margin-top:0;",
                       "+ Top 5 Mimics (most positive RGES)"),
                    p(style="font-size:11px;color:gray;",
                      "Strongest transcriptional concordance \u2014 potential mimics"),
                    DT::DTOutput("rges_bot5")))
              )
            ),
            tabPanel("Full results",
              br(),
              p(style="font-size:11px;color:gray;",
                "All drugs shown. Click RGES column header to sort."),
              DT::DTOutput("rges_table")),
            tabPanel("Distribution",
              br(),
              plotOutput("rges_plot",height="320px"),
              br(),
              p(style="font-size:11px;color:gray;",
                "Negative RGES = reversal. Positive RGES = mimicry. Red dashed line at 0. Blue dotted at \u22120.3 (strong reversal threshold)."))
          )
        )
      )
    )
  ),

  # ---- DOWNLOAD ----
  tabPanel("Download",
    fluidPage(
      h3("Download MetaDrugSig"), br(),
      fluidRow(
        column(4, wellPanel(
          h4("Full database"),
          tags$ul(
            tags$li("drug_signatures_logFC.csv"),
            tags$li("drug_signatures_FDR.csv"),
            tags$li("drug_signatures_SE.csv"),
            tags$li("drug_metadata_v3.csv"),
            tags$li("comparison_metadata_v3.csv"),
            tags$li("spearman_similarity_v3.rds"),
            tags$li("DEG_results_v3.tar.gz")),
          a(href=ZENODO_URL, target="_blank",
            class="btn btn-primary", "Zenodo"))),
        column(4, wellPanel(
          h4("Single drug"),
          selectizeInput("dl_drug","Drug:",choices=NULL,
                         options=list(placeholder="Type drug name...")),
          downloadButton("dl_single","Download CSV"))),
        column(4, wellPanel(
          h4("Metadata"),
          p("1,092 drugs | Tier 1 & 2 | ATC codes"),
          downloadButton("dl_meta","Download metadata")))
      ), br(),
      h4("Citation"),
      wellPanel(
        p("Joseph DF et al. MetaDrugSig: A meta-analytic transcriptomic drug signature database.",
          em("Database (Oxford)"), ", 2026."),
        p(a(href=ZENODO_URL, target="_blank", ZENODO_URL)))
    )
  )
)

# ============================================================
# SERVER
# ============================================================
server <- function(input, output, session) {

  updateSelectizeInput(session,"drug_select",choices=drug_choices,server=TRUE)
  updateSelectizeInput(session,"dl_drug",choices=drug_choices,server=TRUE)

  observeEvent(input$drug_select, {
    req(input$drug_select)
    sig <- get_drug_sig(input$drug_select)
    if (!is.null(sig))
      updateSelectizeInput(session,"gene_direct",
                           choices=sig$gene[order(-abs(sig$logFC))],server=TRUE)
  })

  output$drug_info_box <- renderUI({
    req(input$drug_select)
    m <- meta[DrugBank_ID==input$drug_select]
    if (nrow(m)==0) return(NULL)
    tier_badge_col <- if (m$tier==1) "#28a745" else "#fd7e14"
    div(style="background:#f0f4ff;padding:10px;border-radius:6px;",
      h5(m$DrugBank_Name, style="margin:0 0 6px 0;"),
      p(style="margin:0;font-size:11px;color:#555;",
        span(style=paste0("background:",tier_badge_col,
                          ";color:white;padding:2px 6px;border-radius:8px;",
                          "font-size:10px;font-weight:bold;margin-right:6px;"),
             paste0("Tier ",m$tier)),
        strong("Studies:"), m$n_comparisons, " | ",
        strong("ATC:"),
        ifelse(is.na(m$ATC_L2_Name)||m$ATC_L2_Name=="","N/A",m$ATC_L2_Name)))
  })

  drug_sig_rv <- reactive({
    req(input$drug_select)
    get_drug_sig(input$drug_select)
  })

  study_summary_rv <- reactive({
    req(input$drug_select)
    get_drug_study_summary(input$drug_select)
  })

  # ---- STUDY OVERVIEW ----
  output$study_overview_header <- renderUI({
    s <- study_summary_rv(); req(!is.null(s))
    tier_col <- if (s$tier==1) "#28a745" else "#fd7e14"
    tier_lbl <- if (s$tier==1) "Tier 1 — REML meta-analytic" else "Tier 2 — single-study edgeR"
    div(style="background:#f0f4ff;padding:14px;border-radius:8px;",
      fluidRow(
        column(3,
          div(style="text-align:center;",
            div(style=paste0("font-size:36px;font-weight:bold;color:",tier_col,";"),
                s$n_gses),
            div(style="font-size:12px;color:gray;","GEO studies"))),
        column(3,
          div(style="text-align:center;",
            div(style="font-size:36px;font-weight:bold;color:#2c3e50;",
                s$n_total),
            div(style="font-size:12px;color:gray;","comparisons"))),
        column(3,
          div(style="text-align:center;",
            div(style=paste0("font-size:14px;font-weight:bold;color:white;",
                             "background:",tier_col,";padding:6px 10px;",
                             "border-radius:6px;display:inline-block;"),
                tier_lbl))),
        column(3,
          div(style="text-align:center;",
            div(style="font-size:14px;color:#555;font-weight:bold;",
                s$drug_name),
            div(style="font-size:11px;color:gray;",input$drug_select)))
      )
    )
  })

  output$tissue_bar <- renderPlot({
    s <- study_summary_rv(); req(!is.null(s) && nrow(s$tissue) > 0)
    df <- s$tissue
    # Top 12 + Other
    if (nrow(df) > 12) {
      top <- df[1:12]
      other_n <- sum(df[13:nrow(df)]$N)
      df <- rbind(top, data.table(tissue_category="Other", N=other_n))
    }
    df[, tissue_category := factor(tissue_category, levels=rev(tissue_category))]
    ggplot(df, aes(x=tissue_category, y=N)) +
      geom_col(fill="#4393C3", width=0.7) +
      geom_text(aes(label=N), hjust=-0.15, size=3, color="gray30") +
      scale_y_continuous(expand=c(0,0),
                         limits=c(0, max(df$N)*1.15)) +
      coord_flip() +
      labs(x=NULL, y="Studies") +
      theme_classic(base_size=10) +
      theme(axis.text.y=element_text(size=9),
            panel.grid.major.x=element_line(color="gray93",linewidth=0.3))
  })

  output$dose_donut <- renderPlot({
    s <- study_summary_rv(); req(!is.null(s) && nrow(s$dose) > 0)
    df <- s$dose
    total <- sum(df$N)
    df[, pct := round(N/total*100)]
    df[, label := paste0(dose_category, " (", pct, "%)")]
    df[, ymax := cumsum(N)]
    df[, ymin := c(0, head(cumsum(N),-1))]
    cols <- c("low"="#4DAC26","medium"="#F1B300","high"="#B2182B","Unknown"="#CCCCCC")
    df[, color := cols[dose_category]]
    ggplot(df, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=2.8, fill=color)) +
      geom_rect() +
      geom_text(aes(x=4.6, y=(ymax+ymin)/2, label=label),
                size=2.8, hjust=0) +
      scale_fill_identity() +
      coord_polar(theta="y") +
      xlim(c(0,5)) +
      theme_void(base_size=9)
  })

  output$duration_donut <- renderPlot({
    s <- study_summary_rv(); req(!is.null(s) && nrow(s$duration) > 0)
    df <- s$duration
    total <- sum(df$N)
    df[, pct := round(N/total*100)]
    df[, label := paste0(duration_category, " (", pct, "%)")]
    df[, ymax := cumsum(N)]
    df[, ymin := c(0, head(cumsum(N),-1))]
    cols <- c("acute"="#B2182B","chronic"="#2166AC","Unknown"="#CCCCCC")
    df[, color := cols[duration_category]]
    ggplot(df, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=2.8, fill=color)) +
      geom_rect() +
      geom_text(aes(x=4.6, y=(ymax+ymin)/2, label=label),
                size=2.8, hjust=0) +
      scale_fill_identity() +
      coord_polar(theta="y") +
      xlim(c(0,5)) +
      theme_void(base_size=9)
  })

  output$context_value_ui <- renderUI({
    req(input$context_filter != "none")
    drug_id <- input$drug_select; req(!is.null(drug_id))
    cm <- comp_meta[DrugBank_ID == drug_id]
    choices <- switch(input$context_filter,
      "tissue"   = sort(unique(cm$tissue_category[!is.na(cm$tissue_category)])),
      "dose"     = sort(unique(cm$dose_category[!is.na(cm$dose_category)])),
      "duration" = sort(unique(cm$duration_category[!is.na(cm$duration_category)])),
      character(0))
    if (length(choices)==0) return(NULL)
    selectInput("context_value", "Context value:", choices=choices)
  })

  # ---- Tier count label ----
  output$tier_count_label <- renderUI({
    n <- if (input$rges_tier=="tier1") length(tier1_ids) else nrow(meta)
    p(style="font-size:10px;color:gray;margin:4px 0 0 0;",
      sprintf("Scoring against %d drugs", n))
  })

  # ---- VOLCANO ----
  output$volcano_plot <- renderPlotly({
    sig <- drug_sig_rv(); req(!is.null(sig)&&nrow(sig)>0)
    thresh <- input$fdr_thresh
    sig$sig_cat <- "NS"
    sig$sig_cat[!is.na(sig$FDR)&sig$FDR<thresh&sig$logFC>0] <- "Up"
    sig$sig_cat[!is.na(sig$FDR)&sig$FDR<thresh&sig$logFC<0] <- "Down"
    sig$neglog10fdr <- -log10(sig$FDR+1e-300)
    sig$hover_text  <- paste0(
      "<b>",sig$gene,"</b><br>logFC: ",round(sig$logFC,3),
      "<br>FDR: ",formatC(sig$FDR,format="e",digits=2),
      "<br>SE: ",round(sig$SE,3))
    plot_ly(sig,x=~logFC,y=~neglog10fdr,type="scatter",mode="markers",
            color=~sig_cat,
            colors=c("Up"="#B2182B","Down"="#2166AC","NS"="#B4B4B4"),
            text=~hover_text,hoverinfo="text",
            marker=list(size=5,opacity=0.7),source="volcano") %>%
      add_segments(x=0,xend=0,y=0,yend=max(sig$neglog10fdr,na.rm=TRUE)+1,
                   line=list(color="gray",dash="dash",width=1),showlegend=FALSE) %>%
      add_segments(x=min(sig$logFC,na.rm=TRUE),xend=max(sig$logFC,na.rm=TRUE),
                   y=-log10(thresh),yend=-log10(thresh),
                   line=list(color="gray",dash="dash",width=1),showlegend=FALSE) %>%
      layout(title=list(text=paste("Volcano:",input$drug_select),font=list(size=13)),
             xaxis=list(title="log2 fold change",zeroline=FALSE),
             yaxis=list(title="-log10(FDR)"),
             legend=list(orientation="h",y=-0.15),
             clickmode="event+select") %>%
      event_register("plotly_click")
  })

  clicked_gene  <- reactiveVal(NULL)
  forest_data_rv <- reactiveVal(NULL)
  forest_gene_rv <- reactiveVal(NULL)

  observeEvent(event_data("plotly_click",source="volcano"), {
    click <- event_data("plotly_click",source="volcano"); req(!is.null(click))
    sig   <- drug_sig_rv(); req(!is.null(sig))
    nearest <- sig$gene[which.min(
      abs(sig$logFC-click$x)+abs((-log10(sig$FDR+1e-300))-click$y))]
    clicked_gene(nearest)
    updateSelectizeInput(session,"gene_direct",selected=nearest)
    updateTabsetPanel(session,"main_tabs",selected="Forest plot")
    showNotification(paste("Generating forest plot for",nearest),type="message")
  })

  observeEvent(list(clicked_gene(),input$run_forest), {
    gene <- if (!is.null(clicked_gene())&&input$run_forest==0)
      clicked_gene() else input$gene_direct
    req(!is.null(gene)&&gene!="")
    drug_name <- meta$DrugBank_Name[meta$DrugBank_ID==input$drug_select][1]
    req(!is.na(drug_name))
    withProgress(message=paste("Loading studies for",gene,"..."), {
      gd <- load_per_gse(drug_name,gene)
      forest_data_rv(gd); forest_gene_rv(gene)
    })
    metareg_rv(NULL)
  }, ignoreInit=TRUE)

  output$forest_status <- renderUI({
    gd <- forest_data_rv(); gene <- forest_gene_rv()
    if (is.null(gene))
      return(div(style="padding:40px;text-align:center;color:gray;",
        h4("Click a gene on the volcano or use direct input")))
    if (is.null(gd)||nrow(gd)<MIN_STUDIES_TOTAL)
      return(div(style="padding:20px;color:darkorange;",
        icon("exclamation-triangle"),
        p(paste("Insufficient per-study data for",gene,
                "— need at least",MIN_STUDIES_TOTAL,"studies."))))
    NULL
  })

  output$forest_plot <- renderPlot({
    gd <- forest_data_rv(); gene <- forest_gene_rv()
    req(!is.null(gd)&&nrow(gd)>=MIN_STUDIES_TOTAL&&!is.null(gene))
    drug_name <- meta$DrugBank_Name[meta$DrugBank_ID==input$drug_select][1]
    build_forest_plot(gd,gene,drug_name)
  }, height=function() {
    gd <- forest_data_rv()
    if (is.null(gd)) 400 else max(400,as.integer(nrow(gd)*22+120))
  })

  metareg_rv <- reactiveVal(NULL)
  observeEvent(input$run_metareg, {
    gd <- forest_data_rv(); gene <- forest_gene_rv()
    req(!is.null(gd)&&!is.null(gene)&&length(input$metareg_mods)>0)
    updateTabsetPanel(session,"main_tabs",selected="Heterogeneity moderators")
    withProgress(message="Running meta-regression...",{
      res <- run_meta_regression(gd,input$metareg_mods); metareg_rv(res)
    })
  })

  output$metareg_output <- renderUI({
    res <- metareg_rv(); gene <- forest_gene_rv(); gd <- forest_data_rv()
    if (is.null(res))
      return(div(style="padding:40px;text-align:center;color:gray;",
        h4("Generate a forest plot first, then click 'Run meta-regression'")))
    if (!is.null(res$error))
      return(div(style="padding:20px;color:red;",p(res$error)))
    n_studies <- if (!is.null(gd)) nrow(gd) else "?"
    tagList(
      div(style="background:#f8f9fa;padding:12px;border-radius:6px;margin-bottom:16px;",
        h4(style="margin:0 0 6px 0;",
           paste("Meta-regression:",gene,"\u2014",n_studies,"studies")),
        p(style="margin:0;font-size:12px;color:#555;",strong("Null model: "),
          sprintf("I\u00b2 = %.1f%% | tau\u00b2 = %.3f | AIC = %.1f",
                  res$null$I2,res$null$tau2,res$null$AIC))),
      lapply(names(res$moderators), function(mod) {
        r <- res$moderators[[mod]]
        warn_ui <- if (length(r$warnings)>0)
          div(style="background:#fff3cd;padding:8px;border-radius:4px;margin-bottom:8px;",
              lapply(r$warnings,function(w)
                p(style="margin:2px 0;font-size:11px;",icon("exclamation-triangle")," ",w)))
        else NULL
        if (isTRUE(r$skipped))
          return(div(style="border:1px solid #dee2e6;border-radius:6px;padding:12px;margin-bottom:12px;",
            h5(mod,style="margin:0 0 8px 0;"),warn_ui,
            p(style="color:gray;font-size:12px;","Moderator skipped \u2014 see warnings.")))
        badge_color <- if (r$delta_AIC>=2&&as.numeric(r$QM_pval)<0.05) "#28a745"
                       else if (r$delta_AIC>=2) "#fd7e14" else "#6c757d"
        div(style="border:1px solid #dee2e6;border-radius:6px;padding:14px;margin-bottom:14px;",
          div(style="display:flex;justify-content:space-between;align-items:center;margin-bottom:10px;",
            h5(mod,style="margin:0;"),
            span(style=paste0("background:",badge_color,
                              ";color:white;padding:3px 8px;border-radius:12px;",
                              "font-size:11px;font-weight:bold;"),
              if (r$delta_AIC>=2&&as.numeric(r$QM_pval)<0.05) "SIGNIFICANT"
              else if (r$delta_AIC>=2) "MARGINAL" else "NOT SIGNIFICANT")),
          warn_ui,
          div(style="display:grid;grid-template-columns:repeat(4,1fr);gap:8px;margin-bottom:10px;",
            div(style="background:#f8f9fa;padding:8px;border-radius:4px;text-align:center;",
              div(style="font-size:20px;font-weight:bold;color:#2c3e50;",paste0(r$R2,"%")),
              div(style="font-size:10px;color:gray;","R\u00b2")),
            div(style="background:#f8f9fa;padding:8px;border-radius:4px;text-align:center;",
              div(style="font-size:20px;font-weight:bold;color:#2c3e50;",round(r$delta_AIC,1)),
              div(style="font-size:10px;color:gray;","\u0394AIC")),
            div(style="background:#f8f9fa;padding:8px;border-radius:4px;text-align:center;",
              div(style="font-size:20px;font-weight:bold;color:#2c3e50;",paste0(r$I2_residual,"%")),
              div(style="font-size:10px;color:gray;",
                  paste0("Residual I\u00b2 (was ",r$I2_null,"%)"))),
            div(style="background:#f8f9fa;padding:8px;border-radius:4px;text-align:center;",
              div(style="font-size:20px;font-weight:bold;color:#2c3e50;",r$QM_pval),
              div(style="font-size:10px;color:gray;","QM p-value"))),
          p(style="font-size:12px;font-style:italic;color:#495057;margin-bottom:10px;",
            r$interpretation),
          h6("Per-level coefficients:"),
          div(style="font-size:11px;",
              renderTable(r$coef_df,striped=TRUE,bordered=TRUE,spacing="xs",digits=3)),
          h6(style="margin-top:8px;","Study counts per level:"),
          div(style="font-size:11px;color:gray;",
              paste(paste0(r$level_counts$vals," (n=",r$level_counts$Freq,")"),
                    collapse=" | ")))
      }))
  })

  output$tbl_up <- DT::renderDT({
    sig <- drug_sig_rv(); req(!is.null(sig))
    top <- head(sig[sig$logFC>0,][order(-sig[sig$logFC>0,]$logFC),],50)
    top$logFC <- round(top$logFC,3); top$SE <- round(top$SE,3)
    top$FDR <- formatC(top$FDR,format="e",digits=2)
    DT::datatable(top[,c("gene","logFC","SE","FDR")],
                  options=list(pageLength=10,dom="tp"),rownames=FALSE)
  })
  output$tbl_down <- DT::renderDT({
    sig <- drug_sig_rv(); req(!is.null(sig))
    top <- head(sig[sig$logFC<0,][order(sig[sig$logFC<0,]$logFC),],50)
    top$logFC <- round(top$logFC,3); top$SE <- round(top$SE,3)
    top$FDR <- formatC(top$FDR,format="e",digits=2)
    DT::datatable(top[,c("gene","logFC","SE","FDR")],
                  options=list(pageLength=10,dom="tp"),rownames=FALSE)
  })

  output$dl_sig <- downloadHandler(
    filename=function() paste0(input$drug_select,"_signature.csv"),
    content=function(file)
      write.csv(drug_sig_rv()[order(-drug_sig_rv()$logFC),],file,row.names=FALSE))
  output$dl_forest <- downloadHandler(
    filename=function() paste0(input$drug_select,"_",forest_gene_rv(),"_forest.pdf"),
    content=function(file) {
      gd <- forest_data_rv(); gene <- forest_gene_rv()
      req(!is.null(gd)&&!is.null(gene))
      drug_name <- meta$DrugBank_Name[meta$DrugBank_ID==input$drug_select][1]
      p <- build_forest_plot(gd,gene,drug_name)
      ggplot2::ggsave(file,p,width=13,
                      height=max(6,0.30*nrow(gd)+3),useDingbats=FALSE)
    })
  output$dl_single <- downloadHandler(
    filename=function() paste0(input$dl_drug,"_signature.csv"),
    content=function(file) {
      sig <- get_drug_sig(input$dl_drug)
      write.csv(sig[order(-sig$logFC),],file,row.names=FALSE)
    })
  output$dl_meta <- downloadHandler(
    filename=function() "drug_metadata_v3.csv",
    content=function(file) fwrite(meta,file))

  # ---- MSIGDB ----
  observeEvent(input$msig_collection, {
    req(input$msig_collection != "\u2014 Select \u2014")
    updateSelectizeInput(session,"msig_pathway",
                         choices=names(pathway_db[[input$msig_collection]]),
                         server=TRUE,
                         options=list(placeholder="Type to search pathways..."))
  })
  output$msig_preview <- renderUI({
    req(input$msig_collection != "\u2014 Select \u2014", input$msig_pathway != "")
    gp <- pathway_db[[input$msig_collection]][[input$msig_pathway]]
    if (is.null(gp)) return(NULL)
    n_total <- length(gp); n_in_db <- sum(gp %in% rownames(logfc_mat))
    div(style="font-size:11px;color:#555;margin:6px 0;",
      p(style="margin:2px 0;",
        sprintf("Pathway: %d total genes, %d in MetaDrugSig",n_total,n_in_db)),
      p(style="margin:2px 0;color:gray;",
        paste("Preview:",paste(head(gp,5),collapse=", "),
              if(n_total>5)"..."else"")))
  })
  observeEvent(input$sig_up, {
    updateSelectizeInput(session,"msig_pathway",selected="")
  }, ignoreInit=TRUE)
  observeEvent(input$sig_down, {
    updateSelectizeInput(session,"msig_pathway",selected="")
  }, ignoreInit=TRUE)
  observeEvent(input$load_msig, {
    req(input$msig_collection != "\u2014 Select \u2014", input$msig_pathway != "")
    gp <- pathway_db[[input$msig_collection]][[input$msig_pathway]]
    req(!is.null(gp))
    gf <- gp[gp %in% rownames(logfc_mat)]
    if (input$msig_direction=="up") {
      updateTextAreaInput(session,"sig_up",value=paste(gf,collapse="\n"))
      updateTextAreaInput(session,"sig_down",value="")
    } else {
      updateTextAreaInput(session,"sig_up",value="")
      updateTextAreaInput(session,"sig_down",value=paste(gf,collapse="\n"))
    }
    showNotification(
      sprintf("Loaded %d/%d genes from %s as '%s'.",
              length(gf),length(gp),input$msig_pathway,
              if(input$msig_direction=="up")
                "disease upregulated" else "disease downregulated"),
      type="message",duration=5)
  })

  # ============================================================
  # RGES COMPUTATION
  # ============================================================
  rges_rv <- eventReactive(input$run_rges, {
    up   <- trimws(unlist(strsplit(input$sig_up,   "\n"))); up   <- up[up   != ""]
    down <- trimws(unlist(strsplit(input$sig_down, "\n"))); down <- down[down != ""]
    validate(need(length(up)+length(down)>0,"Enter at least one gene."))

    sig_up   <- data.frame(gene=up)
    sig_down <- data.frame(gene=down)

    ids <- if (input$rges_tier=="tier1")
      intersect(colnames(logfc_mat), tier1_ids) else colnames(logfc_mat)

    withProgress(
      message=sprintf("Computing RGES across %d drugs...", length(ids)),
      value=0, {
        scores <- vapply(seq_along(ids), function(i) {
          incProgress(1/length(ids))
          v <- logfc_mat[, ids[i]]; v <- v[!is.na(v)]
          tryCatch(
            cmap_score_new(sig_up, sig_down,
                           data.frame(ids=names(v), rank=v)),
            error=function(e) NA_real_)
        }, numeric(1))
      })

    data.frame(DrugBank_ID=ids, RGES=round(scores,4)) %>%
      left_join(as.data.frame(
        meta[, .(DrugBank_ID, DrugBank_Name, tier,
                 n_comparisons, ATC_L2_Name, ATC_L4_Name)]),
        by="DrugBank_ID") %>%
      arrange(RGES) %>%
      filter(!is.na(RGES))
  })

  output$rges_top5 <- DT::renderDT({
    top5 <- head(rges_rv(),5); top5$RGES <- round(top5$RGES,4)
    DT::datatable(
      top5[,c("DrugBank_Name","RGES","tier","n_comparisons","ATC_L4_Name")],
      options=list(dom="t",pageLength=5),rownames=FALSE,
      colnames=c("Drug","RGES","Tier","Studies","ATC Class")) %>%
      DT::formatStyle("RGES",color="#1a7a1a",fontWeight="bold")
  })
  output$rges_bot5 <- DT::renderDT({
    bot5 <- tail(rges_rv(),5); bot5$RGES <- round(bot5$RGES,4)
    DT::datatable(
      bot5[,c("DrugBank_Name","RGES","tier","n_comparisons","ATC_L4_Name")],
      options=list(dom="t",pageLength=5),rownames=FALSE,
      colnames=c("Drug","RGES","Tier","Studies","ATC Class")) %>%
      DT::formatStyle("RGES",color="#a01010",fontWeight="bold")
  })

  output$rges_query_summary <- renderUI({
    df   <- rges_rv()
    up   <- trimws(unlist(strsplit(input$sig_up,  "\n"))); up   <- up[up  !=""]
    down <- trimws(unlist(strsplit(input$sig_down,"\n"))); down <- down[down!=""]
    n_neg <- sum(df$RGES < 0,na.rm=TRUE)
    n_str <- sum(df$RGES < -0.3,na.rm=TRUE)
    tier_label <- if (input$rges_tier=="tier1")
      "Tier 1 only (REML)" else "All drugs (Tier 1 + 2)"
    pathway_genes <- if (!is.null(input$msig_pathway)&&input$msig_pathway!=""&&
                         !is.null(input$msig_collection)&&
                         input$msig_collection!="\u2014 Select \u2014")
      pathway_db[[input$msig_collection]][[input$msig_pathway]] else character(0)
    genes_from_pathway <- length(pathway_genes)>0 &&
      all(c(up,down) %in% pathway_genes)
    pathway_label <- if (genes_from_pathway) input$msig_pathway
      else if (length(up)+length(down)>0)
        sprintf("Custom gene list (%d genes total)",length(up)+length(down))
      else "No genes entered"
    div(style="background:#f8f9fa;padding:12px;border-radius:6px;",
      h5(style="margin:0 0 8px 0;","Query summary"),
      p(style="margin:2px 0;font-size:12px;",strong("Query: "),pathway_label),
      p(style="margin:2px 0;font-size:12px;",
        strong("Genes: "),length(up)," disease upregulated | ",
        length(down)," disease downregulated"),
      p(style="margin:2px 0;font-size:12px;",strong("Tier: "),tier_label),
      p(style="margin:2px 0;font-size:12px;",strong("Drugs scored: "),nrow(df)),
      p(style="margin:2px 0;font-size:12px;",
        strong("Reversers (RGES < 0): "),
        sprintf("%d (%.1f%%)",n_neg,100*n_neg/nrow(df))),
      p(style="margin:2px 0;font-size:12px;",
        strong("Strong reversers (RGES < \u22120.3): "),
        sprintf("%d (%.1f%%)",n_str,100*n_str/nrow(df))))
  })

  output$rges_table <- DT::renderDT({
    DT::datatable(rges_rv(),
      options=list(pageLength=25,dom="lfrtip",
                   order=list(list(1,"asc"))),
      rownames=FALSE) %>%
      DT::formatStyle("RGES",
                      color=DT::styleInterval(0,c("#2166AC","#B2182B")))
  })

  output$rges_plot <- renderPlot({
    df <- rges_rv()
    ggplot(df,aes(x=RGES))+
      geom_histogram(bins=60,fill="gray70",color="white")+
      geom_vline(xintercept=0,linetype="dashed",color="red",linewidth=0.8)+
      geom_vline(xintercept=-0.3,linetype="dotted",color="steelblue",linewidth=0.8)+
      annotate("text",x=-0.32,y=Inf,
               label="Strong reversal\nthreshold (\u22120.3)",
               hjust=1,vjust=1.5,size=3,color="steelblue")+
      labs(x="RGES score",y="Number of drugs",
           title=sprintf("RGES distribution \u2014 %d drugs",nrow(df)))+
      theme_minimal(base_size=12)
  })

  output$dl_rges <- downloadHandler(
    filename=function() {
      query <- if (!is.null(input$msig_pathway)&&input$msig_pathway!="")
        gsub("[^A-Za-z0-9_]","_",input$msig_pathway) else "custom"
      tier  <- if (input$rges_tier=="tier1") "_tier1" else "_all"
      paste0("MetaDrugSig_RGES_",query,tier,".csv")
    },
    content=function(file) write.csv(rges_rv(),file,row.names=FALSE))
}

shinyApp(ui=ui, server=server)
