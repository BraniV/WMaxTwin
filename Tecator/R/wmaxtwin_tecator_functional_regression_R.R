
################################################################################
# WMaxTwin for Tecator Functional Regression, R version
#
# This script reproduces the Tecator functional-regression example used in the
# WMaxTwin paper.  It is intentionally self-contained and uses only base R plus
# the recommended stats and utils packages that ship with R.
#
# The workflow is:
#   1. Load the StatLib Tecator data.
#   2. Interpolate each 100-channel spectrum to 128 grid points.
#   3. Standardize spectra using calibration-sample statistics only.
#   4. Compute an orthonormal Haar representation of each spectrum.
#   5. Construct raw MaxTwin and nested WMaxTwin pairwise distance matrices.
#   6. Compare Random, SPlit, Twinning, DUPLEX, MaxTwin, and WMaxTwin splits.
#   7. Select jmax and ridge penalty by internal validation inside calibration.
#   8. Refit on all 172 calibration samples and evaluate on the 43-sample
#      external test set.
#
# Default run: R_REPS = 8.  Increase to 50 or more for a more stable Monte
# Carlo comparison.
################################################################################

DATA_URL <- "https://lib.stat.cmu.edu/datasets/tecator"
OUTDIR <- file.path(getwd(), "tecator_wmaxtwin_outputs_R")
DATA_PATH <- file.path(getwd(), "tecator_statlib.txt")

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------
# Data handling
# -----------------------------------------------------------------------------

download_tecator <- function(path = DATA_PATH) {
  if (file.exists(path) && file.info(path)$size > 100000) {
    return(path)
  }
  message("Downloading Tecator data from ", DATA_URL)
  utils::download.file(DATA_URL, path, quiet = FALSE, mode = "wb")
  path
}

load_tecator <- function(path = DATA_PATH) {
  text <- paste(readLines(path, warn = FALSE), collapse = "\n")
  hits <- gregexpr("extrapolation_examples=25", text, fixed = TRUE)[[1]]
  if (length(hits) < 2 || hits[1] < 0) {
    stop("Could not locate Tecator data block.")
  }
  sub <- substring(text, hits[2])
  nl <- regexpr("\n", sub, fixed = TRUE)[1]
  start <- hits[2] + nl - 1
  nums <- scan(text = substring(text, start + 1), quiet = TRUE)
  needed <- 240 * 125
  if (length(nums) < needed) {
    stop(sprintf("Parsed only %d numeric values; expected at least %d.", length(nums), needed))
  }
  arr <- matrix(nums[seq_len(needed)], nrow = 240, ncol = 125, byrow = TRUE)
  X <- arr[, 1:100, drop = FALSE]
  endpoints <- arr[, 123:125, drop = FALSE]
  colnames(endpoints) <- c("moisture", "fat", "protein")
  wavelengths <- seq(850, 1050, length.out = 100)
  list(X = X, endpoints = endpoints, wavelengths = wavelengths)
}

interp_to_power2 <- function(X, wavelengths, m = 128) {
  grid <- seq(min(wavelengths), max(wavelengths), length.out = m)
  Xnew <- t(apply(X, 1, function(row) stats::approx(wavelengths, row, xout = grid, rule = 2)$y))
  list(X = Xnew, grid = grid)
}

standardize_fit <- function(X) {
  mu <- colMeans(X)
  sig <- apply(X, 2, stats::sd)
  sig[sig < 1e-12] <- 1
  list(center = mu, scale = sig)
}

standardize_apply <- function(X, sc) {
  sweep(sweep(X, 2, sc$center, "-"), 2, sc$scale, "/")
}

# -----------------------------------------------------------------------------
# Haar transform
# -----------------------------------------------------------------------------

haar_transform_matrix <- function(n) {
  J <- log2(n)
  if (abs(J - round(J)) > 1e-12) stop("n must be a power of two")
  J <- as.integer(round(J))
  rows <- list(rep(1 / sqrt(n), n))
  for (lev in seq_len(J)) {
    block_len <- n / (2^(lev - 1))
    half <- block_len / 2
    nblocks <- 2^(lev - 1)
    for (b in seq_len(nblocks)) {
      v <- rep(0, n)
      start <- (b - 1) * block_len + 1
      v[start:(start + half - 1)] <- 1 / sqrt(block_len)
      v[(start + half):(start + block_len - 1)] <- -1 / sqrt(block_len)
      rows[[length(rows) + 1]] <- v
    }
  }
  do.call(rbind, rows)
}

haar_coefficients <- function(X) {
  n_grid <- ncol(X)
  H <- haar_transform_matrix(n_grid)
  C <- X %*% t(H)
  J <- as.integer(round(log2(n_grid)))
  blocks <- list(scaling = 1L)
  idx <- 2L
  for (j in seq_len(J)) {
    L <- 2^(j - 1)
    blocks[[as.character(j)]] <- idx:(idx + L - 1)
    idx <- idx + L
  }
  list(C = C, blocks = blocks, H = H)
}

# -----------------------------------------------------------------------------
# Geometry and scale weights
# -----------------------------------------------------------------------------

pairwise_sq_dists <- function(Z) {
  G <- tcrossprod(Z)
  s <- rowSums(Z * Z)
  D2 <- outer(s, s, "+") - 2 * G
  D2[D2 < 0] <- 0
  D2
}

normalize_distance <- function(D2) {
  vals <- D2[upper.tri(D2)]
  pos <- vals[vals > 0]
  med <- if (length(pos) > 0) stats::median(pos) else 1
  D2 / med
}

scale_relevance_weights <- function(C_cal, y_cal, blocks) {
  y <- as.numeric(scale(y_cal, center = TRUE, scale = TRUE))
  w <- numeric(0)
  for (nm in names(blocks)) {
    if (nm == "scaling") next
    idx <- blocks[[nm]]
    Z <- C_cal[, idx, drop = FALSE]
    Zc <- sweep(Z, 2, colMeans(Z), "-")
    zs <- apply(Zc, 2, stats::sd)
    ok <- zs > 1e-12
    if (!any(ok)) {
      w[nm] <- 0
    } else {
      corr <- as.vector(crossprod(Zc[, ok, drop = FALSE], y)) / ((length(y) - 1) * zs[ok])
      w[nm] <- mean(corr^2)
    }
  }
  total <- sum(w)
  if (total <= 0) {
    w[] <- 1 / length(w)
  } else {
    w <- w / total
  }
  w
}

wavelet_weighted_distance <- function(C, blocks, weights) {
  n <- nrow(C)
  D <- matrix(0, n, n)
  for (nm in names(weights)) {
    idx <- blocks[[nm]]
    Z <- C[, idx, drop = FALSE]
    Zc <- sweep(Z, 2, colMeans(Z), "-")
    sdv <- apply(Zc, 2, stats::sd)
    sdv[sdv < 1e-8] <- 1
    Zs <- sweep(Zc, 2, sdv, "/")
    D <- D + as.numeric(weights[nm]) * normalize_distance(pairwise_sq_dists(Zs))
  }
  D
}

# -----------------------------------------------------------------------------
# Split constructors
# -----------------------------------------------------------------------------

random_split <- function(n, m_val) {
  sort(sample.int(n, size = m_val, replace = FALSE))
}

support_split <- function(D, m_val) {
  n <- nrow(D)
  Dj <- D + matrix(rnorm(n * n, sd = 1e-10), n, n)
  meanD <- rowMeans(Dj)
  selected <- integer(0)
  remaining <- seq_len(n)
  sum_mean <- 0
  sum_pair <- 0
  for (step in 0:(m_val - 1)) {
    best_i <- NA_integer_
    best_obj <- Inf
    m_new <- step + 1
    for (i in remaining) {
      add_pair <- if (step == 0) 0 else 2 * sum(Dj[i, selected])
      obj <- 2 * (sum_mean + meanD[i]) / m_new - (sum_pair + add_pair) / (m_new^2)
      if (obj < best_obj) {
        best_obj <- obj
        best_i <- i
      }
    }
    if (step > 0) sum_pair <- sum_pair + 2 * sum(Dj[best_i, selected])
    selected <- c(selected, best_i)
    remaining <- setdiff(remaining, best_i)
    sum_mean <- sum_mean + meanD[best_i]
  }
  selected
}

duplex_split <- function(D, m_val) {
  n <- nrow(D)
  Dj <- D + matrix(rnorm(n * n, sd = 1e-10), n, n)
  U <- Dj
  U[lower.tri(U, diag = TRUE)] <- -Inf
  pos <- which(U == max(U), arr.ind = TRUE)[1, ]
  selected <- as.integer(c(pos[1], pos[2]))
  remaining <- setdiff(seq_len(n), selected)
  while (length(selected) < m_val) {
    cand <- remaining
    mind <- apply(Dj[cand, selected, drop = FALSE], 1, min)
    chosen <- cand[which.max(mind)]
    selected <- c(selected, chosen)
    remaining <- setdiff(remaining, chosen)
  }
  selected
}

twinning_split <- function(D, m_val) {
  n <- nrow(D)
  anchors <- support_split(D, m_val)
  Dj <- D + matrix(rnorm(n * n, sd = 1e-10), n, n)
  unused <- seq_len(n)
  val <- integer(0)
  for (a0 in anchors) {
    if (length(val) >= m_val) break
    a <- as.integer(a0)
    if (!(a %in% unused) || length(unused) < 2) next
    candidates <- setdiff(unused, a)
    b <- candidates[which.min(Dj[a, candidates])]
    val <- c(val, if (runif(1) < 0.5) a else b)
    unused <- setdiff(unused, c(a, b))
  }
  while (length(val) < m_val && length(unused) >= 2) {
    u <- unused
    subD <- Dj[u, u, drop = FALSE]
    diag(subD) <- Inf
    pos <- which(subD == min(subD), arr.ind = TRUE)[1, ]
    a <- u[pos[1]]; b <- u[pos[2]]
    val <- c(val, if (runif(1) < 0.5) a else b)
    unused <- setdiff(unused, c(a, b))
  }
  if (length(val) < m_val) {
    extra <- sample(unused, size = m_val - length(val), replace = FALSE)
    val <- c(val, extra)
  }
  as.integer(val)
}

maxtwin_split <- function(D, m_val) {
  n <- nrow(D)
  Dj <- D + matrix(rnorm(n * n, sd = 1e-10), n, n)
  unused <- seq_len(n)
  val <- integer(0)
  while (length(val) < m_val && length(unused) >= 2) {
    u <- unused
    subD <- Dj[u, u, drop = FALSE]
    diag(subD) <- Inf
    pos <- which(subD == min(subD), arr.ind = TRUE)[1, ]
    a <- u[pos[1]]; b <- u[pos[2]]
    val <- c(val, if (runif(1) < 0.5) a else b)
    unused <- setdiff(unused, c(a, b))
  }
  if (length(val) < m_val) {
    extra <- sample(unused, size = m_val - length(val), replace = FALSE)
    val <- c(val, extra)
  }
  as.integer(val)
}

# -----------------------------------------------------------------------------
# Wavelet-ridge model selection and final test evaluation
# -----------------------------------------------------------------------------

features_by_jmax <- function(C, blocks, jmax) {
  idx <- 1L
  for (j in seq_len(jmax)) idx <- c(idx, blocks[[as.character(j)]])
  C[, idx, drop = FALSE]
}

ridge_fit <- function(X, y, alpha) {
  X_aug <- cbind(Intercept = 1, X)
  p <- ncol(X_aug)
  P <- diag(c(0, rep(alpha, p - 1)), p, p)
  A <- crossprod(X_aug) + P
  b <- crossprod(X_aug, y)
  beta <- tryCatch(as.numeric(solve(A, b)), error = function(e) as.numeric(qr.solve(A, b)))
  list(beta = beta)
}

ridge_predict <- function(fit, X) {
  X_aug <- cbind(Intercept = 1, X)
  as.numeric(X_aug %*% fit$beta)
}

fit_select_refit <- function(C_cal, y_cal, C_test, y_test, blocks, train_idx, val_idx,
                             alphas, jmax_grid) {
  best_vmse <- Inf
  best_j <- NA_integer_
  best_alpha <- NA_real_
  records <- data.frame(jmax = integer(), alpha = numeric(), val_mse = numeric())
  for (jmax in jmax_grid) {
    F <- features_by_jmax(C_cal, blocks, jmax)
    Xtr <- F[train_idx, , drop = FALSE]
    Xva <- F[val_idx, , drop = FALSE]
    ytr <- y_cal[train_idx]
    yva <- y_cal[val_idx]
    sc <- standardize_fit(Xtr)
    Xtr_s <- standardize_apply(Xtr, sc)
    Xva_s <- standardize_apply(Xva, sc)
    for (alpha in alphas) {
      fit <- ridge_fit(Xtr_s, ytr, alpha)
      pred <- ridge_predict(fit, Xva_s)
      vmse <- mean((yva - pred)^2)
      records <- rbind(records, data.frame(jmax = jmax, alpha = alpha, val_mse = vmse))
      if (vmse < best_vmse) {
        best_vmse <- vmse
        best_j <- jmax
        best_alpha <- alpha
      }
    }
  }
  Fcal <- features_by_jmax(C_cal, blocks, best_j)
  Ftest <- features_by_jmax(C_test, blocks, best_j)
  sc <- standardize_fit(Fcal)
  fit <- ridge_fit(standardize_apply(Fcal, sc), y_cal, best_alpha)
  pred_test <- ridge_predict(fit, standardize_apply(Ftest, sc))
  test_mse <- mean((y_test - pred_test)^2)
  list(jmax = best_j,
       alpha = best_alpha,
       val_mse = best_vmse,
       test_mse = test_mse,
       test_rmse = sqrt(test_mse),
       pred_test = pred_test,
       val_curve = records)
}

method_order <- function(results = NULL) {
  c("Random", "SPlit", "Twinning", "DUPLEX", "MaxTwin",
    "WMaxTwin gamma=0.10", "WMaxTwin gamma=0.25",
    "WMaxTwin gamma=0.50", "WMaxTwin gamma=0.75",
    "WMaxTwin equal gamma=0.25")
}

# -----------------------------------------------------------------------------
# Figures and summary tables
# -----------------------------------------------------------------------------

summarize_results <- function(results) {
  ord <- method_order()
  rows <- list()
  for (m in ord[ord %in% unique(results$method)]) {
    sub <- results[results$method == m, ]
    rows[[length(rows) + 1]] <- data.frame(
      method = m,
      mean_jmax = mean(sub$jmax),
      sd_jmax = stats::sd(sub$jmax),
      median_alpha = stats::median(sub$alpha),
      mean_val_mse = mean(sub$val_mse),
      sd_val_mse = stats::sd(sub$val_mse),
      mean_test_rmse = mean(sub$test_rmse),
      sd_test_rmse = stats::sd(sub$test_rmse),
      mean_test_mse = mean(sub$test_mse),
      sd_test_mse = stats::sd(sub$test_mse)
    )
  }
  ans <- do.call(rbind, rows)
  ans[is.na(ans)] <- 0
  ans
}

make_figures <- function(X_cal_raw, wl, y_cal, weights, results, val_curves_to_plot) {
  ord <- method_order()
  ord <- ord[ord %in% unique(results$method)]

  png(file.path(OUTDIR, "fig1_tecator_spectra.png"), width = 1400, height = 850, res = 180)
  matplot(wl, t(X_cal_raw[seq(1, nrow(X_cal_raw), by = 4), ]), type = "l", lty = 1,
          col = rgb(0, 0, 0, 0.25), xlab = "Wavelength (nm)", ylab = "Absorbance",
          main = "Tecator absorbance spectra, calibration samples")
  dev.off()

  qs <- quantile(y_cal, probs = seq(0, 1, 0.25), type = 7)
  qs <- unique(qs)
  qg <- cut(y_cal, breaks = qs, include.lowest = TRUE, labels = FALSE)
  png(file.path(OUTDIR, "fig2_quartile_mean_spectra.png"), width = 1400, height = 850, res = 180)
  plot(wl, colMeans(X_cal_raw[qg == 1, , drop = FALSE]), type = "l", lwd = 2,
       xlab = "Wavelength (nm)", ylab = "Mean absorbance",
       ylim = range(X_cal_raw), main = "Mean spectra by fat-content quartile")
  for (g in sort(unique(qg))[-1]) lines(wl, colMeans(X_cal_raw[qg == g, , drop = FALSE]), lwd = 2, lty = g)
  legend("topright", legend = paste("fat quartile", sort(unique(qg))), lty = sort(unique(qg)), lwd = 2, bty = "n", cex = 0.8)
  dev.off()

  png(file.path(OUTDIR, "fig3_response_scale_weights.png"), width = 1100, height = 750, res = 180)
  barplot(as.numeric(weights), names.arg = names(weights), xlab = "Wavelet scale j (coarse to fine)",
          ylab = "Normalized response-relevance weight", main = "Calibration-only scale weights for WMaxTwin")
  dev.off()

  if (length(val_curves_to_plot) > 0) {
    vc <- do.call(rbind, val_curves_to_plot)
    png(file.path(OUTDIR, "fig4_validation_curves_rep0.png"), width = 1450, height = 850, res = 180)
    plot(NA, xlim = range(vc$jmax), ylim = range(vc$val_mse), xlab = "Maximum included wavelet scale jmax",
         ylab = "Best validation MSE over ridge penalties", main = "Example validation curves, first repetition")
    labs <- unique(vc$method)
    for (i in seq_along(labs)) {
      sub <- vc[vc$method == labs[i], ]
      agg <- aggregate(val_mse ~ jmax, data = sub, FUN = min)
      lines(agg$jmax, agg$val_mse, type = "b", pch = i, lty = i, lwd = 1.8)
    }
    legend("topright", legend = labs, pch = seq_along(labs), lty = seq_along(labs), lwd = 1.8, bty = "n", cex = 0.8)
    dev.off()
  }

  data <- lapply(ord, function(m) results$test_rmse[results$method == m])
  png(file.path(OUTDIR, "fig5_test_rmse_boxplot.png"), width = 1600, height = 900, res = 180)
  boxplot(data, names = ord, las = 2, ylab = "External test RMSE: fat percent", main = "Tecator external test error across repeated internal splits")
  points(seq_along(data), sapply(data, mean), pch = 19)
  dev.off()

  data <- lapply(ord, function(m) results$jmax[results$method == m])
  png(file.path(OUTDIR, "fig6_selected_jmax_boxplot.png"), width = 1600, height = 900, res = 180)
  boxplot(data, names = ord, las = 2, ylab = "Selected jmax", main = "Selected wavelet-resolution cutoff")
  points(seq_along(data), sapply(data, mean), pch = 19)
  dev.off()

  data <- lapply(ord, function(m) log10(results$alpha[results$method == m]))
  png(file.path(OUTDIR, "fig7_selected_alpha_boxplot.png"), width = 1600, height = 900, res = 180)
  boxplot(data, names = ord, las = 2, ylab = "log10 selected ridge penalty", main = "Selected shrinkage penalty across repeated internal splits")
  points(seq_along(data), sapply(data, mean), pch = 19)
  dev.off()

  labs <- c("MaxTwin", sprintf("WMaxTwin gamma=%.2f", c(0.10, 0.25, 0.50, 0.75)))
  labs <- labs[labs %in% unique(results$method)]
  xvals <- c(0.0, 0.10, 0.25, 0.50, 0.75)[seq_along(labs)]
  y <- sapply(labs, function(m) mean(results$test_rmse[results$method == m]))
  e <- sapply(labs, function(m) stats::sd(results$test_rmse[results$method == m]) / sqrt(sum(results$method == m)))
  png(file.path(OUTDIR, "fig8_gamma_path_rmse.png"), width = 1150, height = 800, res = 180)
  plot(xvals, y, type = "b", pch = 19, xlab = "Nested WMaxTwin mixing parameter gamma",
       ylab = "Mean external test RMSE", main = "Nested path: MaxTwin is gamma=0")
  arrows(xvals, y - e, xvals, y + e, angle = 90, code = 3, length = 0.05)
  dev.off()
}

write_latex_table <- function(summary, weights) {
  bs <- "\\\\"
  lines <- c("% Generated by wmaxtwin_tecator_functional_regression_R.R", "",
             "\\begin{tabular}{lrrrrr}",
             "\\hline",
             paste0("split & mean $\\widehat{j}_{\\max}$ & median $\\widehat\\lambda$ & mean val. MSE & mean test RMSE & sd test RMSE ", bs),
             "\\hline")
  for (i in seq_len(nrow(summary))) {
    r <- summary[i, ]
    lines <- c(lines, paste0(sprintf("%s & %.2f & %.3g & %.4f & %.4f & %.4f",
                                     r$method, r$mean_jmax, r$median_alpha,
                                     r$mean_val_mse, r$mean_test_rmse, r$sd_test_rmse),
                              " ", bs))
  }
  lines <- c(lines, "\\hline", "\\end{tabular}")
  writeLines(lines, file.path(OUTDIR, "tecator_wmaxtwin_table_R.tex"))

  wlines <- c("scale_j,response_weight", sprintf("%s,%.12g", names(weights), as.numeric(weights)))
  writeLines(wlines, file.path(OUTDIR, "tecator_scale_weights_R.csv"))
}

# -----------------------------------------------------------------------------
# Main study
# -----------------------------------------------------------------------------

run_study <- function(R = 8, seed = 20260619) {
  download_tecator(DATA_PATH)
  dat <- load_tecator(DATA_PATH)
  X <- dat$X
  endpoints <- dat$endpoints
  wavelengths <- dat$wavelengths
  y_fat <- endpoints[, "fat"]

  ip <- interp_to_power2(X, wavelengths, m = 128)
  X128 <- ip$X
  wl128 <- ip$grid

  X_cal_raw <- X128[1:172, , drop = FALSE]
  y_cal <- y_fat[1:172]
  X_test_raw <- X128[173:215, , drop = FALSE]
  y_test <- y_fat[173:215]

  sc_curve <- standardize_fit(X_cal_raw)
  X_cal <- standardize_apply(X_cal_raw, sc_curve)
  X_test <- standardize_apply(X_test_raw, sc_curve)

  hc <- haar_coefficients(X_cal)
  C_cal <- hc$C
  blocks <- hc$blocks
  C_test <- X_test %*% t(hc$H)

  D0 <- normalize_distance(pairwise_sq_dists(X_cal))
  weights <- scale_relevance_weights(C_cal, y_cal, blocks)
  Dw <- wavelet_weighted_distance(C_cal, blocks, weights)
  equal_weights <- weights
  equal_weights[] <- 1 / length(equal_weights)
  Dw_equal <- wavelet_weighted_distance(C_cal, blocks, equal_weights)

  gamma_path <- c(0, 0.10, 0.25, 0.50, 0.75)
  D_gamma <- list()
  for (g in gamma_path) {
    key <- sprintf("%.2f", g)
    D_gamma[[key]] <- normalize_distance((1 - g) * D0 + g * Dw)
  }
  D_gamma_equal <- normalize_distance(0.75 * D0 + 0.25 * Dw_equal)

  n <- nrow(X_cal)
  m_val <- 43
  all_idx <- seq_len(n)
  alphas <- 10^seq(-4, 4, length.out = 13)
  jmax_grid <- seq_len(round(log2(ncol(X_cal))))

  result_rows <- list()
  val_curves_to_plot <- list()

  for (r in 0:(R - 1)) {
    set.seed(seed + r)
    methods <- list()
    methods[["Random"]] <- random_split(n, m_val)
    methods[["SPlit"]] <- support_split(D0, m_val)
    methods[["Twinning"]] <- twinning_split(D0, m_val)
    methods[["DUPLEX"]] <- duplex_split(D0, m_val)
    methods[["MaxTwin"]] <- maxtwin_split(D0, m_val)
    for (g in gamma_path[-1]) {
      nm <- sprintf("WMaxTwin gamma=%.2f", g)
      methods[[nm]] <- maxtwin_split(D_gamma[[sprintf("%.2f", g)]], m_val)
    }
    methods[["WMaxTwin equal gamma=0.25"]] <- maxtwin_split(D_gamma_equal, m_val)

    for (nm in names(methods)) {
      val_idx <- sort(unique(as.integer(methods[[nm]])))
      if (length(val_idx) < m_val) {
        missing <- setdiff(all_idx, val_idx)
        extra <- sample(missing, size = m_val - length(val_idx), replace = FALSE)
        val_idx <- sort(c(val_idx, extra))
      }
      train_idx <- setdiff(all_idx, val_idx)
      res <- fit_select_refit(C_cal, y_cal, C_test, y_test, blocks,
                              train_idx, val_idx, alphas, jmax_grid)
      result_rows[[length(result_rows) + 1]] <- data.frame(
        rep = r,
        method = nm,
        jmax = res$jmax,
        alpha = res$alpha,
        val_mse = res$val_mse,
        test_mse = res$test_mse,
        test_rmse = res$test_rmse
      )
      if (r == 0 && nm %in% c("Random", "MaxTwin", "WMaxTwin gamma=0.25", "WMaxTwin gamma=0.50")) {
        vc <- res$val_curve
        vc$method <- nm
        val_curves_to_plot[[length(val_curves_to_plot) + 1]] <- vc
      }
    }
  }

  results <- do.call(rbind, result_rows)
  utils::write.csv(results, file.path(OUTDIR, "tecator_split_results_R.csv"), row.names = FALSE)
  scale_table <- data.frame(scale_j = as.integer(names(weights)),
                            response_weight = as.numeric(weights),
                            equal_weight = as.numeric(equal_weights))
  utils::write.csv(scale_table, file.path(OUTDIR, "tecator_scale_weights_R.csv"), row.names = FALSE)

  summary <- summarize_results(results)
  utils::write.csv(summary, file.path(OUTDIR, "tecator_split_summary_R.csv"), row.names = FALSE)
  make_figures(X_cal_raw, wl128, y_cal, weights, results, val_curves_to_plot)
  write_latex_table(summary, weights)

  list(results = results, summary = summary, weights = weights, outdir = OUTDIR)
}

if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  R_REPS <- 8
  if (length(args) >= 1) R_REPS <- as.integer(args[1])
  ans <- run_study(R = R_REPS, seed = 20260619)
  cat("\nCalibration-only scale weights:\n")
  print(data.frame(scale_j = names(ans$weights), weight = as.numeric(ans$weights)), row.names = FALSE)
  cat("\nSummary:\n")
  print(ans$summary, row.names = FALSE)
  cat("\nOutputs written to:\n", ans$outdir, "\n")
}
