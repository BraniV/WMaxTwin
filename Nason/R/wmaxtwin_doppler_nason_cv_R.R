
############################################################
# WMaxTwin Doppler/Nason-type wavelet cross-validation in R
#
# This is a self-contained R translation of the Python notebook
# for Example II in the WMaxTwin paper.
#
# Main design point:
#   Nason's even/odd split is grid-compatible.
#   MaxTwin/WMaxTwin splits are generally not dyadic subgrids.
#   Therefore the fit is a masked wavelet penalized least-squares
#   problem on the original grid, not an ordinary DWT of retained
#   observations.
#
# Dependencies: base R only.  No Wavelet Toolbox/package is needed.
# The periodic DWT/IDWT for sym4 and db2 is implemented directly.
############################################################

options(stringsAsFactors = FALSE)

############################################################
# User settings
############################################################

# If this file is sourced from the R Markdown/Jupyter notebook, set
# RUN_MAIN <- FALSE before source().  If run directly by Rscript,
# the script executes the default simulation.
if (!exists("RUN_MAIN", inherits = FALSE)) RUN_MAIN <- TRUE

N      <- 1024
SNR    <- 7
SIGMA  <- 1
WAVELET <- "sym4"
J0     <- 3
NREP   <- 18
SEED   <- 2026
N_ITER <- 30
WEIGHT_MODE <- "paper"  # "paper" uses explicit omega_j; "script" reproduces the original Python split weights
OUTPUT_DIR <- "outputs_R"

############################################################
# Orthogonal wavelet filters, PyWavelets periodic convention
############################################################

wavelet_filters <- function(wavelet = "sym4") {
  if (wavelet == "sym4") {
    h <- c(
      -0.07576571478927333,
      -0.02963552764599851,
       0.49761866763201545,
       0.8037387518059161,
       0.29785779560527736,
      -0.09921954357684722,
      -0.012603967262037833,
       0.0322231006040427
    )
    g <- c(
      -0.0322231006040427,
      -0.012603967262037833,
       0.09921954357684722,
       0.29785779560527736,
      -0.8037387518059161,
       0.49761866763201545,
       0.02963552764599851,
      -0.07576571478927333
    )
  } else if (wavelet == "db2") {
    h <- c(
      -0.12940952255126034,
       0.2241438680420134,
       0.8365163037378079,
       0.48296291314453416
    )
    g <- c(
      -0.48296291314453416,
       0.8365163037378079,
      -0.2241438680420134,
      -0.12940952255126034
    )
  } else {
    stop("Unknown wavelet. Use 'sym4' or 'db2'.")
  }
  list(h = h, g = g)
}

periodic_idx <- function(n, L) {
  m <- n / 2
  outer(2 * (0:(m - 1)), 0:(L - 1), function(a, b) ((a + b) %% n) + 1)
}

dwt_step <- function(x, h, g) {
  n <- length(x)
  idx <- periodic_idx(n, length(h))
  X <- matrix(x[as.vector(idx)], nrow = nrow(idx), ncol = ncol(idx))
  list(a = as.vector(X %*% h), d = as.vector(X %*% g))
}

idwt_step <- function(a, d, h, g) {
  m <- length(a)
  n <- 2 * m
  idx <- periodic_idx(n, length(h))
  vals <- outer(a, h, "*") + outer(d, g, "*")
  as.numeric(tabulate(as.vector(idx), weights = as.vector(vals), nbins = n))
}

wavedec_periodic <- function(x, j0 = 3, wavelet = "sym4") {
  filt <- wavelet_filters(wavelet)
  h <- filt$h
  g <- filt$g
  n <- length(x)
  J <- log2(n)
  if (abs(J - round(J)) > 1e-12) stop("n must be a power of two")
  J <- as.integer(round(J))
  if (!(1 <= j0 && j0 < J)) stop("j0 must satisfy 1 <= j0 < J")

  a <- as.numeric(x)
  details_by_level <- list()
  for (lev in seq(J - 1, j0, by = -1)) {
    st <- dwt_step(a, h, g)
    a <- st$a
    details_by_level[[as.character(lev)]] <- st$d
  }
  details <- lapply(j0:(J - 1), function(j) details_by_level[[as.character(j)]])
  names(details) <- as.character(j0:(J - 1))
  list(a = a, details = details)
}

waverec_periodic <- function(a, details, j0 = 3, wavelet = "sym4") {
  filt <- wavelet_filters(wavelet)
  h <- filt$h
  g <- filt$g
  cur <- as.numeric(a)
  for (r in seq_along(details)) {
    cur <- idwt_step(cur, details[[r]], h, g)
  }
  cur
}

soft <- function(x, lambda) {
  sign(x) * pmax(abs(x) - lambda, 0)
}

shrink_full <- function(y, lams, j0 = 3, wavelet = "sym4") {
  dec <- wavedec_periodic(y, j0, wavelet)
  details2 <- vector("list", length(dec$details))
  names(details2) <- names(dec$details)
  for (r in seq_along(dec$details)) {
    j <- as.integer(names(dec$details)[r])
    details2[[r]] <- soft(dec$details[[r]], lams[as.character(j)])
  }
  waverec_periodic(dec$a, details2, j0, wavelet)
}

############################################################
# Doppler signal and shrinkage utilities
############################################################

doppler_signal <- function(n, snr = 7, sigma = 1) {
  t <- ((1:n) - 0.5) / n
  f <- sqrt(t * (1 - t)) * sin(2 * pi * 1.05 / (t + 0.05))
  f <- f - mean(f)
  f <- f * sqrt(snr * sigma^2 / mean(f^2))
  list(t = t, f = f)
}

mse <- function(a, b) {
  mean((a - b)^2)
}

sure_threshold_detail <- function(d, sigma = 1) {
  x <- sort(abs(d))
  n <- length(x)
  s2 <- sigma^2
  csum <- cumsum(x^2)
  k <- 1:n
  risks <- n * s2 + csum + (n - k) * x^2 - 2 * s2 * k
  x[which.min(risks)]
}

level_sure_lams <- function(y, j0 = 3, wavelet = "sym4", sigma = 1) {
  J <- as.integer(log2(length(y)))
  dec <- wavedec_periodic(y, j0, wavelet)
  vals <- sapply(seq_along(dec$details), function(r) {
    sure_threshold_detail(dec$details[[r]], sigma)
  })
  names(vals) <- as.character(j0:(J - 1))
  vals
}

fixed_lams <- function(lambda, j0, J) {
  vals <- rep(lambda, J - j0)
  names(vals) <- as.character(j0:(J - 1))
  vals
}

ramp_lams <- function(slope, j0, J, power = 1.5) {
  js <- j0:(J - 1)
  vals <- slope * ((js - j0) / (J - 1 - j0))^power
  names(vals) <- as.character(js)
  vals
}

############################################################
# Masked wavelet lasso by FISTA
############################################################

fit_masked_wavelet <- function(y, mask, lams, j0 = 3, wavelet = "sym4", n_iter = 35) {
  dec <- wavedec_periodic(y * mask, j0, wavelet)
  a <- dec$a
  details <- dec$details
  za <- a
  zd <- details
  tk <- 1

  for (iter in 1:n_iter) {
    fz <- waverec_periodic(za, zd, j0, wavelet)
    resid <- mask * (fz - y)
    grad <- wavedec_periodic(resid, j0, wavelet)

    a_new <- za - grad$a
    d_new <- vector("list", length(zd))
    names(d_new) <- names(zd)
    for (r in seq_along(zd)) {
      j <- as.integer(names(zd)[r])
      d_new[[r]] <- soft(zd[[r]] - grad$details[[r]], lams[as.character(j)])
    }

    t_new <- (1 + sqrt(1 + 4 * tk^2)) / 2
    momentum <- (tk - 1) / t_new
    za <- a_new + momentum * (a_new - a)
    zd <- lapply(seq_along(d_new), function(r) d_new[[r]] + momentum * (d_new[[r]] - details[[r]]))
    names(zd) <- names(d_new)
    a <- a_new
    details <- d_new
    tk <- t_new
  }

  waverec_periodic(a, details, j0, wavelet)
}

############################################################
# MaxTwin and WMaxTwin split construction
############################################################

detail_contributions <- function(y, j0 = 3, wavelet = "sym4") {
  dec <- wavedec_periodic(y, j0, wavelet)
  zero_a <- rep(0, length(dec$a))
  out <- list()
  for (r in seq_along(dec$details)) {
    tmp <- lapply(dec$details, function(d) rep(0, length(d)))
    names(tmp) <- names(dec$details)
    tmp[[r]] <- dec$details[[r]]
    j <- names(dec$details)[r]
    out[[j]] <- waverec_periodic(zero_a, tmp, j0, wavelet)
  }
  out
}

wmaxtwin_scale_weights <- function(j0, J, mode = "paper") {
  js <- j0:(J - 1)
  if (mode == "paper" && j0 == 3 && J == 10) {
    vals <- c(0.015, 0.036, 0.074, 0.124, 0.182, 0.248, 0.321)
    names(vals) <- as.character(3:9)
    return(vals)
  }
  if (mode == "script") {
    vals <- numeric(length(js))
    for (ii in seq_along(js)) {
      j <- js[ii]
      if (j >= J - 3) {
        rel <- j - (J - 3) + 1
        vals[ii] <- c(1.0, 1.8, 2.6)[rel]
      } else {
        vals[ii] <- 0.10 * ((j - j0 + 1) / (J - j0))^2
      }
    }
    names(vals) <- as.character(js)
    return(vals)
  }
  # Generic fallback: normalized high-frequency ramp with a small floor.
  eta <- 0.05
  vals <- eta + ((js - j0) / (J - 1 - j0))^1.5
  vals <- vals / sum(vals)
  names(vals) <- as.character(js)
  vals
}

max_or_wmax_split <- function(y, j0 = 3, wavelet = "sym4", kind = "wmax", weight_mode = "paper") {
  n <- length(y)
  t <- ((1:n) - 0.5) / n

  if (kind == "nason") {
    train <- seq(1, n, by = 2)
    val <- seq(2, n, by = 2)
    return(list(train = train, val = val))
  }

  if (kind == "max") {
    train <- integer(0)
    val <- integer(0)
    pair_starts <- seq(1, n, by = 2)
    for (q in seq_along(pair_starts)) {
      i <- pair_starts[q]
      j <- i + 1
      if (q %% 2 == 1) {
        train <- c(train, i)
        val <- c(val, j)
      } else {
        train <- c(train, j)
        val <- c(val, i)
      }
    }
    return(list(train = train, val = val))
  }

  cont <- detail_contributions(y, j0, wavelet)
  J <- as.integer(log2(n))
  weights <- wmaxtwin_scale_weights(j0, J, mode = weight_mode)

  feats <- list(scale(as.numeric(t))[, 1])
  for (j in j0:(J - 1)) {
    q <- abs(cont[[as.character(j)]])
    q <- as.numeric(scale(q))
    q[is.na(q)] <- 0
    wt <- weights[as.character(j)]
    feats[[length(feats) + 1]] <- sqrt(wt) * q
  }
  X <- do.call(cbind, feats)

  activity <- rep(0, n)
  for (j in max(j0, J - 3):(J - 1)) {
    activity <- activity + abs(cont[[as.character(j)]])
  }
  ord <- order(activity, decreasing = TRUE)
  remaining <- rep(TRUE, n)
  pairs <- matrix(integer(0), ncol = 2)

  for (i in ord) {
    if (!remaining[i]) next
    remaining[i] <- FALSE
    rem <- which(remaining)
    if (length(rem) == 0) break
    Xi <- matrix(X[i, ], nrow = length(rem), ncol = ncol(X), byrow = TRUE)
    dist <- rowSums((X[rem, , drop = FALSE] - Xi)^2)
    jj <- rem[which.min(dist)]
    remaining[jj] <- FALSE
    pairs <- rbind(pairs, c(i, jj))
  }

  bal <- rep(0, ncol(X))
  train <- integer(0)
  val <- integer(0)
  for (q in 1:nrow(pairs)) {
    i <- pairs[q, 1]
    j <- pairs[q, 2]
    diff <- X[i, ] - X[j, ]
    if (sqrt(sum((bal + diff)^2)) <= sqrt(sum((bal - diff)^2))) {
      train <- c(train, i)
      val <- c(val, j)
      bal <- bal + diff
    } else {
      train <- c(train, j)
      val <- c(val, i)
      bal <- bal - diff
    }
  }
  list(train = train, val = val)
}

############################################################
# Cross-validation selectors
############################################################

cv_select_fixed <- function(y, train, val, j0 = 3, wavelet = "sym4", n_iter = 30) {
  n <- length(y)
  J <- as.integer(log2(n))
  mask_t <- rep(0, n); mask_t[train] <- 1
  mask_v <- rep(0, n); mask_v[val] <- 1

  rows <- data.frame(score = numeric(0), lambda = numeric(0))
  for (lambda in seq(0.2, 2.4, length.out = 12)) {
    lams <- fixed_lams(lambda, j0, J)
    pred_t <- fit_masked_wavelet(y, mask_t, lams, j0, wavelet, n_iter)
    pred_v <- fit_masked_wavelet(y, mask_v, lams, j0, wavelet, n_iter)
    score <- mean((pred_t[val] - y[val])^2) + mean((pred_v[train] - y[train])^2)
    rows <- rbind(rows, data.frame(score = score, lambda = lambda))
  }
  rows <- rows[order(rows$score), ]
  list(lams = fixed_lams(rows$lambda[1], j0, J), chosen = rows[1, ], rows = rows)
}

cv_select_ramp <- function(y, train, val, j0 = 3, wavelet = "sym4", n_iter = 30, flat_tol = 1.015) {
  n <- length(y)
  J <- as.integer(log2(n))
  mask_t <- rep(0, n); mask_t[train] <- 1
  mask_v <- rep(0, n); mask_v[val] <- 1

  rows <- data.frame(score = numeric(0), slope = numeric(0))
  for (slope in seq(0.6, 3.4, length.out = 15)) {
    lams <- ramp_lams(slope, j0, J, power = 1.5)
    pred_t <- fit_masked_wavelet(y, mask_t, lams, j0, wavelet, n_iter)
    pred_v <- fit_masked_wavelet(y, mask_v, lams, j0, wavelet, n_iter)
    score <- mean((pred_t[val] - y[val])^2) + mean((pred_v[train] - y[train])^2)
    rows <- rbind(rows, data.frame(score = score, slope = slope))
  }
  rows <- rows[order(rows$score), ]
  min_score <- rows$score[1]
  cand <- rows[rows$score <= flat_tol * min_score, ]
  chosen <- cand[which.max(cand$slope), ]
  list(lams = ramp_lams(chosen$slope, j0, J, power = 1.5), chosen = chosen, rows = rows)
}

oracle_ramp <- function(y, f, j0 = 3, wavelet = "sym4") {
  J <- as.integer(log2(length(y)))
  best_err <- Inf
  best_slope <- NA_real_
  for (slope in seq(0.4, 4.0, length.out = 37)) {
    lams <- ramp_lams(slope, j0, J, power = 1.5)
    err <- mse(shrink_full(y, lams, j0, wavelet), f)
    if (err < best_err) {
      best_err <- err
      best_slope <- slope
    }
  }
  c(error = best_err, slope = best_slope)
}

############################################################
# Plotting helpers
############################################################

save_plots <- function(saved, summary_df, output_dir) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  png(file.path(output_dir, "fig01_one_replication.png"), width = 1100, height = 650, res = 130)
  par(mar = c(4, 4, 2, 1))
  plot(saved$t, saved$y, type = "l", col = "grey70", lwd = 1, xlab = "t", ylab = "signal", main = "One Doppler replication")
  lines(saved$t, saved$f, lwd = 2)
  lines(saved$t, saved$f_sure, lwd = 1.5, lty = 2)
  lines(saved$t, saved$f_nason, lwd = 1.5, lty = 3)
  lines(saved$t, saved$f_maxtwin, lwd = 1.5, lty = 4)
  lines(saved$t, saved$f_wmax, lwd = 2, lty = 1)
  legend("topright", legend = c("noisy", "truth", "LevelSure", "Nason CV", "MaxTwin fixed", "WMaxTwin ramp"),
         lty = c(1, 1, 2, 3, 4, 1), lwd = c(1, 2, 1.5, 1.5, 1.5, 2), col = c("grey70", rep("black", 5)), bty = "n")
  dev.off()

  png(file.path(output_dir, "fig02_feature_weights.png"), width = 850, height = 550, res = 130)
  par(mar = c(4, 4, 2, 1))
  barplot(saved$weights, xlab = "level j", ylab = expression(omega[j]), main = "WMaxTwin feature weights")
  dev.off()

  png(file.path(output_dir, "fig03_split_mask.png"), width = 1100, height = 500, res = 130)
  par(mar = c(4, 4, 2, 1))
  train_flag <- rep(0, length(saved$t)); train_flag[saved$train_w] <- 1
  plot(saved$t, train_flag, pch = 19, cex = 0.55, xlab = "t", ylab = "training indicator", main = "WMaxTwin training-validation mask")
  dev.off()

  png(file.path(output_dir, "fig04_thresholds.png"), width = 850, height = 550, res = 130)
  par(mar = c(4, 4, 2, 1))
  js <- as.integer(names(saved$l_wmax))
  yr <- range(c(saved$l_wmax, saved$l_sure))
  plot(js, saved$l_wmax, type = "b", pch = 19, ylim = yr, xlab = "level j", ylab = expression(lambda[j]), main = "Selected thresholds")
  lines(js, saved$l_sure, type = "b", pch = 17, lty = 2)
  legend("topleft", legend = c("WMaxTwin ramp", "Levelwise SURE"), pch = c(19, 17), lty = c(1, 2), bty = "n")
  dev.off()

  png(file.path(output_dir, "fig05_cv_curve.png"), width = 850, height = 550, res = 130)
  par(mar = c(4, 4, 2, 1))
  rows <- saved$rows_w[order(saved$rows_w$slope), ]
  plot(rows$slope, rows$score, type = "b", pch = 19, xlab = "ramp slope s", ylab = "symmetric CV score", main = "Ramp-slope cross-validation")
  abline(v = saved$chosen_w$slope, lty = 2)
  dev.off()

  png(file.path(output_dir, "fig06_amse_bar.png"), width = 950, height = 600, res = 130)
  par(mar = c(8, 4, 2, 1))
  vals <- summary_df$AMSE
  names(vals) <- summary_df$Method
  barplot(vals, las = 2, ylab = "AMSE", main = "Monte Carlo AMSE comparison")
  dev.off()
}

############################################################
# Main simulation
############################################################

run_simulation <- function(n = 1024, snr = 7, sigma = 1, wavelet = "sym4",
                           j0 = 3, nrep = 18, seed = 2026, n_iter = 30,
                           weight_mode = "paper", output_dir = "outputs_R") {
  J <- as.integer(log2(n))
  if (2^J != n) stop("n must be a power of two")
  if (j0 != 3) warning("This example was tuned for j0 = 3.")
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  sig <- doppler_signal(n, snr, sigma)
  t <- sig$t
  f <- sig$f
  set.seed(seed)

  records <- data.frame()
  saved <- NULL

  for (r in 1:nrep) {
    cat("replication", r, "of", nrep, "\n")
    y <- f + rnorm(n, mean = 0, sd = sigma)

    l_sure <- level_sure_lams(y, j0, wavelet, sigma)
    f_sure <- shrink_full(y, l_sure, j0, wavelet)

    sp_n <- max_or_wmax_split(y, j0, wavelet, kind = "nason", weight_mode = weight_mode)
    cv_n <- cv_select_fixed(y, sp_n$train, sp_n$val, j0, wavelet, n_iter)
    f_n <- shrink_full(y, cv_n$lams, j0, wavelet)

    sp_m <- max_or_wmax_split(y, j0, wavelet, kind = "max", weight_mode = weight_mode)
    cv_m <- cv_select_fixed(y, sp_m$train, sp_m$val, j0, wavelet, n_iter)
    f_m <- shrink_full(y, cv_m$lams, j0, wavelet)

    sp_w <- max_or_wmax_split(y, j0, wavelet, kind = "wmax", weight_mode = weight_mode)
    cv_w <- cv_select_ramp(y, sp_w$train, sp_w$val, j0, wavelet, n_iter)
    f_w <- shrink_full(y, cv_w$lams, j0, wavelet)

    or <- oracle_ramp(y, f, j0, wavelet)

    records <- rbind(records, data.frame(
      rep = r,
      Noisy = mse(y, f),
      LevelSure = mse(f_sure, f),
      NasonEvenOddCV = mse(f_n, f),
      MaxTwinFixedCV = mse(f_m, f),
      WMaxTwinRampCV = mse(f_w, f),
      OracleRamp = as.numeric(or["error"]),
      WMaxTwinSlope = cv_w$chosen$slope,
      OracleSlope = as.numeric(or["slope"])
    ))

    if (is.null(saved)) {
      saved <- list(
        t = t, f = f, y = y,
        f_sure = f_sure, f_nason = f_n, f_maxtwin = f_m, f_wmax = f_w,
        l_sure = l_sure, l_wmax = cv_w$lams,
        rows_w = cv_w$rows, chosen_w = cv_w$chosen,
        train_w = sp_w$train, val_w = sp_w$val,
        weights = wmaxtwin_scale_weights(j0, J, weight_mode)
      )
    }
  }

  method_names <- c("Noisy", "LevelSure", "NasonEvenOddCV", "MaxTwinFixedCV", "WMaxTwinRampCV", "OracleRamp")
  summary_df <- data.frame(
    Method = method_names,
    AMSE = sapply(method_names, function(m) mean(records[[m]])),
    SD = sapply(method_names, function(m) sd(records[[m]])),
    SE = sapply(method_names, function(m) sd(records[[m]]) / sqrt(nrep))
  )
  rownames(summary_df) <- NULL

  write.csv(records, file.path(output_dir, "doppler_replicates_R.csv"), row.names = FALSE)
  write.csv(summary_df, file.path(output_dir, "doppler_amse_summary_R.csv"), row.names = FALSE)
  save_plots(saved, summary_df, output_dir)

  list(summary = summary_df, records = records, saved = saved)
}

############################################################
# Execute when run as a script, but not when sourced from R Markdown
############################################################

if (isTRUE(RUN_MAIN)) {
  result <- run_simulation(
    n = N, snr = SNR, sigma = SIGMA, wavelet = WAVELET,
    j0 = J0, nrep = NREP, seed = SEED, n_iter = N_ITER,
    weight_mode = WEIGHT_MODE, output_dir = OUTPUT_DIR
  )
  print(result$summary)
}
