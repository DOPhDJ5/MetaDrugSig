#!/usr/bin/env Rscript
# ==============================================================================
# MetaDrugSig LOSO validation - combine per-drug results, tables, figure
#   Rscript loso_summarize.R --res loso_out --out loso_report [--B 500] [--N 500] [--min_k_bin 5-9]
# All uncertainty = percentile bootstrap resampling DRUGS (drugs, not genes, are the independent units).
# ==============================================================================
script_dir <- dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(script_dir, "loso_functions.R"))
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })

args <- commandArgs(TRUE)
getopt <- function(key, default) { i <- match(paste0("--", key), args); if (is.na(i)) default else args[i + 1] }
res_dir <- getopt("res", "loso_out"); out_dir <- getopt("out", "loso_report")
B <- as.integer(getopt("B", 500)); HEAD_N <- as.integer(getopt("N", 500))
min_k_bin <- getopt("min_k_bin", "5-9")     # training comparisons per gene used for the main figure
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (dir.exists(file.path(res_dir, "per_drug"))) {
  # output of loso_stream.py (per-drug CSVs)
  pd_dir <- file.path(res_dir, "per_drug")
  stack <- function(nm) rbindlist(lapply(list.files(pd_dir, paste0("__", nm, "\\.csv$"), full.names = TRUE), fread))
  bin <- stack("bin"); sets <- stack("sets"); topn <- stack("topn"); folds <- stack("fold")
  if (!nrow(folds)) stop("no per-drug results in ", pd_dir)
  sk <- file.path(res_dir, "drugs_skipped.csv")
  skipped <- if (file.exists(sk)) fread(sk) else data.table()
  L <- unique(folds$drug)
} else {
  files <- list.files(res_dir, "^drug_.*\\.rds$", full.names = TRUE)
  if (!length(files)) stop("no drug_*.rds files in ", res_dir)
  L <- lapply(files, readRDS)
  is_skip <- vapply(L, function(x) !is.null(x$skipped), TRUE)
  skipped <- rbindlist(lapply(L[is_skip], function(x) data.table(drug_id = x$drug_id, drug_name = x$drug_name, why = x$skipped)))
  L <- L[!is_skip]
  stack <- function(nm) rbindlist(lapply(L, function(x) { t <- copy(x[[nm]]); t[, drug := x$drug_id]; t }))
  bin <- stack("bin"); sets <- stack("sets"); topn <- stack("topn"); folds <- stack("fold")
}
bin[, `:=`(tau_bin = factor(tau_bin, TAU_LEVELS), th_bin = factor(th_bin, TH_LEVELS), k_bin = factor(k_bin, K_LEVELS))]
if (nrow(skipped)) fwrite(skipped, file.path(out_dir, "drugs_skipped.csv"))
cat(sprintf("%d drugs analysed, %d skipped, %d folds\n", length(L), nrow(skipped), nrow(folds)))

# ---- bootstrap over drugs for ratio statistics ---------------------------------
boot_ratio <- function(dt, by, num, den, B = 500, seed = 1) {
  d <- copy(dt)
  d[, grp := do.call(paste, c(.SD, sep = "\r")), .SDcols = by]
  g <- d[, .(num = sum(get(num)), den = sum(get(den))), by = .(drug, grp)]
  ud <- unique(g$drug); ug <- unique(g$grp)
  N <- matrix(0, length(ud), length(ug)); D <- N
  N[cbind(match(g$drug, ud), match(g$grp, ug))] <- g$num
  D[cbind(match(g$drug, ud), match(g$grp, ug))] <- g$den
  est <- colSums(N) / colSums(D)
  set.seed(seed)
  bs <- matrix(replicate(B, { i <- sample.int(length(ud), replace = TRUE)
                              colSums(N[i, , drop = FALSE]) / colSums(D[i, , drop = FALSE]) }), nrow = length(ug))
  key <- unique(d[, c("grp", by), with = FALSE])
  out <- data.table(grp = ug, est = est, lo = apply(bs, 1, quantile, .025, na.rm = TRUE),
                    hi = apply(bs, 1, quantile, .975, na.rm = TRUE), denom = colSums(D), n_drugs = length(ud))
  merge(key, out, by = "grp")[, grp := NULL][]
}

# ---- T1: predictive-interval calibration by tau^2 --------------------------------
cal <- function(b) rbind(boot_ratio(b, "tau_bin", "cover_pi", "n", B)[, interval := "95% prediction interval (uses tau^2)"],
                         boot_ratio(b, "tau_bin", "cover_ci", "n", B)[, interval := "95% confidence interval (ignores tau^2)"])
bin_main <- bin[as.character(k_bin) %in% K_LEVELS[match(min_k_bin, K_LEVELS):length(K_LEVELS)]]
if (!nrow(bin_main)) { warning("no genes in k_bin >= ", min_k_bin, "; using all"); bin_main <- bin }
T1 <- cal(bin_main); T1_all <- cal(bin)
fwrite(T1, file.path(out_dir, "T1_calibration_by_tau2_mainK.csv")); fwrite(T1_all, file.path(out_dir, "T1b_calibration_by_tau2_allK.csv"))
T1k <- rbind(boot_ratio(bin, "k_bin", "cover_pi", "n", B)[, interval := "PI"], boot_ratio(bin, "k_bin", "cover_ci", "n", B)[, interval := "CI"])
fwrite(T1k, file.path(out_dir, "T1c_calibration_by_n_training_comparisons.csv"))

# ---- T2: replication by tau^2 at FIXED |consensus effect| ------------------------
T2s <- boot_ratio(bin_main, c("tau_bin", "th_bin"), "same", "n", B)
T2m <- boot_ratio(bin_main, c("tau_bin", "th_bin"), "aligned", "abs_theta", B)
fwrite(T2s, file.path(out_dir, "T2_same_direction_by_tau2_and_effect.csv"))
fwrite(T2m, file.path(out_dir, "T2b_magnitude_retained_by_tau2_and_effect.csv"))

# ---- T3: MDS gene sets ---------------------------------------------------------------
sets[, one := n]
T3 <- rbindlist(lapply(c("same", "sig_same"), function(m) boot_ratio(sets, "set", m, "n", B)[, metric := m]))
pooled <- sets[, .(n = sum(n), sx = sum(sx), sy = sum(sy), sxx = sum(sxx), syy = sum(syy), sxy = sum(sxy), aligned = sum(aligned)), by = set]
pooled[, `:=`(slope_heldout_on_consensus = (n * sxy - sx * sy) / (n * sxx - sx^2),
              pearson_r = (n * sxy - sx * sy) / sqrt((n * sxx - sx^2) * (n * syy - sy^2)),
              mean_aligned_heldout_effect = aligned / n)]
T3w <- dcast(T3, set ~ metric, value.var = c("est", "lo", "hi"))
T3w <- merge(T3w, pooled[, .(set, n_gene_folds = n, slope_heldout_on_consensus, pearson_r, mean_aligned_heldout_effect)], by = "set")
fwrite(T3w, file.path(out_dir, "T3_MDS_gene_sets.csv"))

# ---- T4: head-to-head ranking at matched N --------------------------------------------
topn[, key := paste(drug, fold)]
broad_m <- c("broad_low_tau2", "broad_low_I2"); allm <- unique(topn$method)
T4 <- list(); T4p <- list()
for (Nn in sort(unique(topn$N))) for (grp in c("core", "with_broad")) {
  meths <- if (grp == "core") setdiff(allm, broad_m) else allm
  t <- topn[N == Nn & method %in% meths]
  ok <- t[, .(ok = .N == length(meths) && all(n_sel == Nn)), by = key][ok == TRUE, key]
  t <- t[key %in% ok]; if (!nrow(t)) next
  t[, one := 1]
  for (mt in c("same", "sig_same", "aligned")) {
    t[, m := get(mt)]
    T4[[length(T4) + 1]] <- boot_ratio(t, "method", "m", "one", B)[, `:=`(N = Nn, comparison_set = grp, metric = mt, n_matched_folds = length(ok))]
    for (ref in c("fixed_effect_z", "single_study_z", "random_effects_z")) {
      w <- dcast(t, drug + fold ~ method, value.var = "m")
      for (mm in setdiff(meths, ref)) {
        dd <- w[, .(drug, diff = get(mm) - get(ref), one = 1)]
        T4p[[length(T4p) + 1]] <- boot_ratio(dd, rep("x", 1)[0] , "diff", "one", B)[, `:=`(method = mm, minus = ref, N = Nn, comparison_set = grp, metric = mt)]
      }
    }
  }
}
T4 <- rbindlist(T4); T4p <- rbindlist(T4p, fill = TRUE)
fwrite(T4, file.path(out_dir, "T4_ranking_head_to_head.csv")); fwrite(T4p, file.path(out_dir, "T4b_paired_differences.csv"))

# ---- Figure ------------------------------------------------------------------------------
pal <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9")
th <- theme_classic(base_size = 11) + theme(legend.position = "bottom", legend.title = element_blank(),
                                              axis.text.x = element_text(angle = 35, hjust = 1))
pA <- ggplot(T1, aes(tau_bin, est, colour = interval, group = interval)) +
  geom_hline(yintercept = 0.95, linetype = 2, colour = "grey50") +
  geom_pointrange(aes(ymin = lo, ymax = hi), position = position_dodge(.3), size = .3) + geom_line(position = position_dodge(.3)) +
  scale_colour_manual(values = pal[1:2]) + coord_cartesian(ylim = c(0, 1)) + guides(colour = guide_legend(nrow = 2)) +
  labs(x = expression("Training "*hat(tau)^2*" (log"[2]*"FC"^2*")"), y = "Held-out estimates inside interval", tag = "A") + th
pB <- ggplot(T2s, aes(tau_bin, est, colour = th_bin, group = th_bin)) +
  geom_pointrange(aes(ymin = lo, ymax = hi), size = .25, position = position_dodge(.25)) + geom_line(position = position_dodge(.25)) +
  geom_hline(yintercept = .5, linetype = 3, colour = "grey60") + scale_colour_manual(values = pal) +
  labs(x = expression("Training "*hat(tau)^2), y = "Same direction in held-out study", tag = "B", colour = "|Consensus log2FC|") +
  th + theme(legend.title = element_text()) +
  guides(colour = guide_legend(nrow = 1, title.position = "top"))
pD <- ggplot(T2m, aes(tau_bin, est, colour = th_bin, group = th_bin)) +
  geom_pointrange(aes(ymin = lo, ymax = hi), size = .25, position = position_dodge(.25)) + geom_line(position = position_dodge(.25)) +
  geom_hline(yintercept = 1, linetype = 2, colour = "grey50") + scale_colour_manual(values = pal) +
  labs(x = expression("Training "*hat(tau)^2), y = "Held-out effect / consensus effect", tag = "D") + th + theme(legend.position = "none")
lab <- c(single_study_z = "Single study\n(|z|)", fixed_effect_z = "Fixed-effect\nmeta (|z|)", random_effects_z = "REML meta\n(|z|)",
         effect_size = "REML meta\n(|log2FC|)", robust_score = "REML meta\n(PI score)")
t4 <- T4[comparison_set == "core" & metric == "same" & N == HEAD_N]
if (!nrow(t4)) t4 <- T4[comparison_set == "core" & metric == "same" & N == max(N)]
t4[, method := factor(method, names(lab), lab)]
pC <- ggplot(t4, aes(method, est)) + geom_pointrange(aes(ymin = lo, ymax = hi), colour = pal[1]) +
  geom_hline(yintercept = .5, linetype = 3, colour = "grey60") +
  labs(x = NULL, y = sprintf("Same direction in held-out study\n(top %d genes per drug)", t4$N[1]), tag = "C") + th
fig <- (pA | pB) / (pC | pD)
ggsave(file.path(out_dir, "Figure_LOSO_tau2_validation.pdf"), fig, width = 10, height = 8.5)
ggsave(file.path(out_dir, "Figure_LOSO_tau2_validation.png"), fig, width = 10, height = 8.5, dpi = 300)

# ---- plain-language headline numbers --------------------------------------------------------
f3 <- function(x) sprintf("%.3f", x)
lines <- c(
  sprintf("Drugs analysed: %d (skipped %d); LOSO folds: %d; training-comparison filter for main figure: k >= %s", length(L), nrow(skipped), nrow(folds), min_k_bin),
  sprintf("Median held-out genes evaluated per fold: %d", as.integer(median(folds$n_eval))),
  "", "Calibration (main-k subset): share of held-out estimates inside the interval",
  capture.output(print(dcast(T1, tau_bin ~ interval, value.var = "est")[, lapply(.SD, function(v) if (is.numeric(v)) round(v, 3) else v)])),
  "", "Same-direction rate by tau^2 at fixed |effect| (see T2 csv)",
  capture.output(print(dcast(T2s, tau_bin ~ th_bin, value.var = "est")[, lapply(.SD, function(v) if (is.numeric(v)) round(v, 3) else v)])),
  "", sprintf("Head-to-head same-direction rate, top %d genes/drug", t4$N[1]),
  capture.output(print(t4[, .(method = gsub("\n", " ", method), est = round(est, 3), lo = round(lo, 3), hi = round(hi, 3))])))
writeLines(lines, file.path(out_dir, "RESULTS_SUMMARY.txt")); cat(paste(lines, collapse = "\n"), "\n")
cat("\nWrote tables + figure to", out_dir, "\n")
