# ==============================================================================
# MetaDrugSig  -  leave-one-study-out (LOSO) validation of tau^2
# Shared functions.  Sourced by loso_run.R, loso_summarize.R, loso_selftest.R
#
# Question (Reviewer 1): does tau^2, estimated from the *other* studies, tell us
# which genes will reproduce (direction + magnitude) in a *held-out* study?
# ==============================================================================
suppressPackageStartupMessages(library(data.table))

TAU_LEVELS <- c("0", "(0,.01]", "(.01,.03]", "(.03,.1]", "(.1,.3]", "(.3,1]", ">1")
TH_LEVELS  <- c("[0,.25)", "[.25,.5)", "[.5,1)", "[1,Inf)")   # |consensus log2FC| bins
K_LEVELS   <- c("2", "3-4", "5-9", "10+")

# ---- 1. Vectorised REML across genes (Fisher scoring) -------------------------
# Y, V : genes x comparisons matrices (effect, sampling variance); NA = not measured
# Same model as metafor::rma(method = "REML"), no Knapp-Hartung adjustment.
.reml_core <- function(Y, V, maxit = 200L, tol = 1e-9) {
  M  <- is.finite(Y) & is.finite(V) & V > 0
  Y[!M] <- 0; V[!M] <- 1
  k  <- rowSums(M)
  ok <- k >= 2L
  # DerSimonian-Laird start
  wf  <- M / V
  S1f <- rowSums(wf); S2f <- rowSums(wf^2)
  thf <- rowSums(wf * Y) / S1f
  Q   <- rowSums(wf * (Y - thf)^2)
  tau2 <- pmax(0, (Q - (k - 1)) / (S1f - S2f / S1f))
  tau2[!is.finite(tau2)] <- 0
  active <- ok
  for (it in seq_len(maxit)) {
    W  <- M / (V + tau2)
    S1 <- rowSums(W); S2 <- rowSums(W^2); S3 <- rowSums(W^3)
    th <- rowSums(W * Y) / S1
    R2 <- rowSums(W^2 * (Y - th)^2)
    score <- 0.5 * (R2 - (S1 - S2 / S1))
    info  <- 0.5 * (S2 - 2 * S3 / S1 + S2^2 / S1^2)
    step  <- score / info
    step[!is.finite(step)] <- 0
    new   <- pmin(pmax(0, tau2 + step), 1e3)
    delta <- abs(new - tau2)
    tau2  <- ifelse(active, new, tau2)
    active <- active & delta > tol
    if (!any(active)) break
  }
  W  <- M / (V + tau2)
  S1 <- rowSums(W)
  th <- rowSums(W * Y) / S1
  s2 <- (k - 1) * S1f / (S1f^2 - S2f)                 # Higgins & Thompson typical variance
  I2 <- 100 * tau2 / (tau2 + s2)
  thf[!ok] <- NA
  th[!ok] <- NA; tau2[!ok] <- NA; I2[!ok] <- NA
  list(theta = th, se = ifelse(ok, 1 / sqrt(S1), NA_real_), tau2 = tau2, I2 = I2, k = k,
       theta_fe = thf, se_fe = ifelse(ok, 1 / sqrt(S1f), NA_real_),   # fixed-effect (tau^2 ignored)
       converged = !active)
}

fast_reml <- function(Y, V, chunk = 8000L, ...) {
  G <- nrow(Y)
  starts <- seq(1L, G, by = chunk)
  res <- lapply(starts, function(s) {
    ix <- s:min(G, s + chunk - 1L)
    .reml_core(Y[ix, , drop = FALSE], V[ix, , drop = FALSE], ...)
  })
  nms <- names(res[[1]])
  setNames(lapply(nms, function(nm) unlist(lapply(res, `[[`, nm), use.names = FALSE)), nms)
}

# ---- 2. Read all per-comparison DEG files for one drug ------------------------
# Layout (same as the portal): <pergse_dir>/<drug folder>/DEGs_<GSE>_Cell_Line__<x>.csv
# columns used: gene, logFC, F.   SE is recovered as |logFC| / sqrt(F)  (QL F = t^2),
# exactly as in the portal (app_test.R).
read_drug_comparisons <- function(drug_dir, max_se = 10) {
  files <- list.files(drug_dir, pattern = "\\.csv$", full.names = TRUE)
  if (!length(files)) return(NULL)
  one <- lapply(files, function(f) {
    d <- tryCatch(fread(f, select = c("gene", "logFC", "F")), error = function(e) NULL)
    if (is.null(d) || !nrow(d)) return(NULL)
    d <- d[!is.na(gene) & !duplicated(gene)]
    d[, se := abs(logFC) / sqrt(F)]
    d <- d[is.finite(logFC) & is.finite(F) & F > 0 & is.finite(se) & se > 0 & se <= max_se]
    if (!nrow(d)) return(NULL)
    d[, .(gene, y = logFC, v = se^2)]
  })
  gse <- sub("^DEGs_", "", sub("_Cell_Line.*\\.csv$", "", basename(files)))
  keep <- !vapply(one, is.null, TRUE)
  one <- one[keep]; gse <- gse[keep]; files <- files[keep]
  if (!length(one)) return(NULL)
  genes <- sort(unique(unlist(lapply(one, `[[`, "gene"), use.names = FALSE)))
  Y <- matrix(NA_real_, length(genes), length(one)); V <- Y
  for (j in seq_along(one)) {
    i <- match(one[[j]]$gene, genes)
    Y[i, j] <- one[[j]]$y; V[i, j] <- one[[j]]$v
  }
  list(Y = Y, V = V, gse = gse, genes = genes, files = basename(files))
}

# ---- 3. One LOSO fold ----------------------------------------------------------
# Train: every comparison NOT in the held-out GSE  -> REML consensus, tau^2, I^2
# Test : held-out GSE (comparisons within that GSE pooled by fixed-effect IVW,
#        because they share controls/lab and are not independent evidence)
.one_fold <- function(Y, V, gse, genes, g, n_sub, topN) {
  held <- gse == g
  Yt <- Y[, !held, drop = FALSE]; Vt <- V[, !held, drop = FALSE]; gt <- gse[!held]
  Yh <- Y[,  held, drop = FALSE]; Vh <- V[,  held, drop = FALSE]
  okh <- is.finite(Yh) & is.finite(Vh) & Vh > 0
  Wh <- ifelse(okh, 1 / Vh, 0); Yh[!okh] <- 0
  sw <- rowSums(Wh); y_h <- rowSums(Wh * Yh) / sw; se_h <- 1 / sqrt(sw)

  fit <- fast_reml(Yt, Vt)
  # number of distinct *GSEs* (not comparisons) that measured each gene in training
  ng <- Reduce(`+`, lapply(unique(gt), function(u)
    as.integer(rowSums(is.finite(Yt[, gt == u, drop = FALSE])) > 0)))
  ok <- is.finite(fit$tau2) & ng >= 2L & sw > 0 & is.finite(y_h)
  if (sum(ok) < 100L) return(NULL)

  # baseline mimicking single-study resources: ONE random training comparison
  j  <- sample.int(ncol(Yt), 1L)
  z1 <- abs(Yt[, j]) / sqrt(Vt[, j]); z1[!is.finite(z1)] <- NA_real_
  d <- data.table(gene = genes[ok], theta = fit$theta[ok], se = fit$se[ok],
                  tau2 = fit$tau2[ok], I2 = fit$I2[ok], k = fit$k[ok],
                  theta_fe = fit$theta_fe[ok], se_fe = fit$se_fe[ok], z1 = z1[ok],
                  y = y_h[ok], se_h = se_h[ok])
  d[, z := theta / se]
  d[, z_fe := theta_fe / se_fe]
  d[, fdr := p.adjust(2 * pnorm(-abs(z)), "BH")]
  d[, broad := fdr < 0.05 & abs(theta) > 0.5]                       # MDS-Broad
  d[, robust := abs(theta) > 1.96 * sqrt(tau2 + se^2)]              # MDS-Robust (PI excludes 0)
  # outcomes in the held-out study
  d[, cover_pi := abs(y - theta) <= 1.96 * sqrt(se_h^2 + se^2 + tau2)]  # predictive interval for a NEW study's estimate
  d[, cover_ci := abs(y - theta) <= 1.96 * sqrt(se_h^2 + se^2)]         # same, ignoring tau^2
  d[, same := sign(y) == sign(theta)]
  d[, sig_h := abs(y / se_h) > 1.96]
  d[, sig_same := sig_h & same]
  d[, aligned := y * sign(theta)]
  d[, tau_bin := factor(as.character(cut(tau2, c(-Inf, 0, .01, .03, .1, .3, 1, Inf), labels = TAU_LEVELS)), TAU_LEVELS)]
  d[, th_bin  := factor(as.character(cut(abs(theta), c(-Inf, .25, .5, 1, Inf), labels = TH_LEVELS, right = FALSE)), TH_LEVELS)]
  d[, k_bin   := factor(as.character(cut(k, c(1, 2, 4, 9, Inf), labels = K_LEVELS)), K_LEVELS)]

  bin <- d[, .(n = .N, cover_pi = sum(cover_pi), cover_ci = sum(cover_ci), same = sum(same),
               sig_h = sum(sig_h), sig_same = sum(sig_same),
               aligned = sum(aligned), abs_theta = sum(abs(theta))),
           by = .(tau_bin, th_bin, k_bin)]

  sumset <- function(m, nm) {
    x <- d[m]; if (!nrow(x)) return(NULL)
    x[, .(set = nm, n = .N, same = sum(same), sig_h = sum(sig_h), sig_same = sum(sig_same),
          aligned = sum(aligned), sx = sum(theta), sy = sum(y), sxx = sum(theta^2),
          syy = sum(y^2), sxy = sum(theta * y))]
  }
  sets <- rbindlist(list(
    sumset(rep(TRUE, nrow(d)), "all_evaluable"),
    sumset(d$broad, "Broad"), sumset(d$robust, "Robust"),
    sumset(d$broad & !d$robust, "Broad_not_Robust"),
    sumset(d$robust & !d$broad, "Robust_not_Broad"),
    sumset(d$broad & d$robust, "Broad_and_Robust"),
    sumset(!d$broad, "not_Broad")))

  # Head-to-head at matched N: which ranking picks genes that replicate best?
  sc <- list(
    single_study_z = list(pool = is.finite(d$z1), s = d$z1),                       # what single-study resources offer
    fixed_effect_z = list(pool = rep(TRUE, nrow(d)), s = abs(d$z_fe)),              # naive meta-analysis, tau^2 ignored
    random_effects_z = list(pool = rep(TRUE, nrow(d)), s = abs(d$z)),               # REML z (tau^2 enters via SE)
    effect_size    = list(pool = rep(TRUE, nrow(d)), s = abs(d$theta)),
    robust_score   = list(pool = rep(TRUE, nrow(d)), s = abs(d$theta) / sqrt(d$tau2 + d$se^2)),  # continuous MDS-Robust
    broad_low_tau2 = list(pool = d$broad, s = -d$tau2),
    broad_low_I2   = list(pool = d$broad, s = -d$I2))
  topn <- rbindlist(lapply(names(sc), function(nm) {
    idx <- which(sc[[nm]]$pool & is.finite(sc[[nm]]$s))
    ord <- idx[order(-sc[[nm]]$s[idx], -abs(d$z[idx]))]
    rbindlist(lapply(topN, function(N) {
      sel <- d[head(ord, N)]
      if (!nrow(sel)) return(NULL)
      data.table(method = nm, N = N, n_sel = nrow(sel), same = mean(sel$same),
                 sig_same = mean(sel$sig_same), aligned = mean(sel$aligned))
    }))
  }))

  sub <- d[sample.int(.N, min(.N, n_sub)),
           .(gene, theta, se, tau2, I2, k, theta_fe, se_fe, y, se_h, z, broad, robust)]
  for (x in list(bin, sets, topn, sub)) x[, fold := g]
  list(bin = bin, sets = sets, topn = topn, sub = sub,
       fold = data.table(fold = g, n_train_comp = sum(!held), n_heldout_comp = sum(held),
                         n_eval = nrow(d), n_broad = sum(d$broad), n_robust = sum(d$robust),
                         median_tau2 = median(d$tau2)))
}

loso_one_drug <- function(Y, V, gse, genes, max_folds = 10L, n_sub = 2000L,
                          topN = c(100L, 500L, 2000L)) {
  ug <- unique(gse)
  if (length(ug) < 3L) return(NULL)            # need >=2 training GSEs after holding one out
  folds <- if (length(ug) > max_folds) sample(ug, max_folds) else ug
  out <- lapply(folds, function(g) .one_fold(Y, V, gse, genes, g, n_sub, topN))
  out <- out[!vapply(out, is.null, TRUE)]
  if (!length(out)) return(NULL)
  setNames(lapply(c("bin", "sets", "topn", "sub", "fold"),
                  function(nm) rbindlist(lapply(out, `[[`, nm))),
           c("bin", "sets", "topn", "sub", "fold"))
}

# ---- 4. Guard: make sure the fast REML == metafor ------------------------------
check_reml_against_metafor <- function(Y, V, n = 300L, tol = 5e-4) {
  stopifnot(requireNamespace("metafor", quietly = TRUE))
  fit <- fast_reml(Y, V)
  idx <- which(rowSums(is.finite(Y) & is.finite(V) & V > 0) >= 2L)
  idx <- head(sample(idx), n)
  mf <- t(vapply(idx, function(i) {
    m <- is.finite(Y[i, ]) & is.finite(V[i, ]) & V[i, ] > 0
    r <- tryCatch(metafor::rma(yi = Y[i, m], vi = V[i, m], method = "REML"), error = function(e) NULL)
    if (is.null(r)) c(NA_real_, NA_real_) else c(r$tau2, as.numeric(r$beta))
  }, numeric(2)))
  d_tau <- abs(fit$tau2[idx] - mf[, 1]); d_th <- abs(fit$theta[idx] - mf[, 2])
  list(max_abs_diff_tau2 = max(d_tau, na.rm = TRUE), max_abs_diff_theta = max(d_th, na.rm = TRUE),
       n_checked = sum(is.finite(d_tau)), pass = all(d_tau < tol, na.rm = TRUE) && all(d_th < tol, na.rm = TRUE))
}
