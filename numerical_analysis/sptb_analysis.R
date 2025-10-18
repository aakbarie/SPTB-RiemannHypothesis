suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(purrr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(readr)
  library(stringr)
  library(scales)
  library(MASS)   # for robust rlm; set cfg$use_huber = FALSE to avoid
})

# ============================================================
# CONFIG
# ============================================================
cfg <- list(
  # Odlyzko zeros (one positive γ per line)
  zeros_path       = "data/zeros_1.txt",
  zeros_n_max      = 300000L,
  
  # SPTB core parameters
  sigma            = 0.60,
  alpha            = 1.0,
  kappa            = 2*pi,            # variance runs only
  c_lambda         = 1.0,             # λ = c_lambda / log(T)^2
  
  # Precision mode (use "double" for speed; switch to "mpfr" for extreme T or η)
  precision        = "double",        # "double" | "mpfr"
  mpfr_bits        = 192L,            # if precision="mpfr"
  
  # Sampling/grid
  m_min            = 12L,             # target samples per block (variance)
  min_nt           = 4000L,
  max_nt           = 30000L,          # variance runs
  
  # Zero thinning
  zeros_cap        = 6000L,
  zeros_pick       = "by_gamma_threshold",    # "uniform" | "head" | "by_gamma_threshold"
  gamma_max        = 200,
  
  # Chunking for Hσ accumulation
  chunk            = 300L,
  
  # VARIANCE sweep
  variance_T_grid  = c(2e2, 5e2, 1e3, 2e3, 5e3, 1e4),
  
  # BIAS sweep (grid of {η, T})
  bias_T_list      = c(1e4, 2e4, 5e4),
  bias_eta_list    = c(1e-4, 5e-4, 1e-3),
  
  # Synthetic zero frequency
  bias_gamma_mode  = "auto",          # "auto" or "fixed"
  gamma0_factor    = 1.5,             # γ0 ≈ gamma0_factor / Δ if auto
  bias_gamma0      = 100,             # used if mode="fixed"
  
  # Standardization of Hσ (based on baseline only)
  standardize      = "center",        # "center" | "zscore" | "none"
  
  # --- Bias-run stability knobs ---
  target_samples_per_block = 8,       # aim ≥ this many samples per block
  kappa_bias_min   = 2*pi,            # clamp κ for bias
  kappa_bias_max   = 60*pi,
  max_nt_bias      = 120000L,         # allow denser grid for bias runs
  use_huber        = TRUE,            # robust tail slope via MASS::rlm
  
  # Output
  out_dir          = "numerical_analysis/sptb_out",
  save_plots       = TRUE,
  save_tables      = TRUE,
  seed             = 7
)

dir.create(cfg$out_dir, showWarnings = FALSE, recursive = TRUE)
set.seed(cfg$seed)

# ============================================================
# 0) Odlyzko zeros loader (γ>0) with ± symmetry via weight=2
# ============================================================
load_zeros_odlyzko <- function(path, n_max = 1e6) {
  g <- scan(path, what = numeric(), quiet = TRUE)
  g <- g[g != 0]
  g <- g[seq_len(min(length(g), n_max))]
  data.table(beta = 0.5, gamma = g, idx = seq_along(g), weight = 2.0)
}
zeros <- load_zeros_odlyzko(cfg$zeros_path, cfg$zeros_n_max)
setDT(zeros)

# ============================================================
# 1) Utilities
# ============================================================
choose_nt <- function(Tmax, Delta, m_min = 12L, min_nt = 4000L, max_nt = 20000L) {
  n_blocks <- ceiling(Tmax / Delta)
  nt <- as.integer(n_blocks * m_min)
  nt <- max(nt, min_nt)
  nt <- min(nt, max_nt)
  nt
}

thin_zeros <- function(zeros, cap = 5000L,
                       method = c("head","uniform","by_gamma_threshold"),
                       gamma_max = NULL) {
  method <- match.arg(method)
  if (nrow(zeros) <= cap || method == "head") {
    return(zeros[1:min(nrow(zeros), cap)])
  }
  if (method == "uniform") {
    idx <- sort(sample.int(nrow(zeros), cap))
    return(zeros[idx])
  }
  if (method == "by_gamma_threshold") {
    stopifnot(!is.null(gamma_max))
    z <- zeros[abs(gamma) <= gamma_max]
    if (nrow(z) > cap) z <- z[1:cap]
    return(z)
  }
}

# Optional: contribution-aware thinning (can replace thin_zeros in bias runs)
thin_by_weight <- function(z, cap = 6000L, sigma = 0.6, alpha = 1, Tref = 2e4) {
  if (nrow(z) <= cap) return(z)
  w <- (sqrt(z$beta^2 + z$gamma^2))^(-alpha) * exp((z$beta - sigma) * Tref)
  idx <- order(w, decreasing = TRUE)[1:cap]
  z[idx]
}

trapz <- function(x, y) {
  if (length(x) < 2L) return(0)
  sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
}

# Compensated cumulative sum (Kahan) to reduce drift
kahan_cumsum <- function(x) {
  s <- 0; c <- 0; out <- numeric(length(x))
  for (i in seq_along(x)) {
    y <- x[i] - c
    t <- s + y
    c <- (t - s) - y
    s <- t
    out[i] <- s
  }
  out
}

# ============================================================
# 2) Optional MPFR helpers
# ============================================================
if (cfg$precision == "mpfr") {
  suppressPackageStartupMessages(library(Rmpfr))
}
to_mpfr <- function(x) if (cfg$precision == "mpfr") Rmpfr::mpfr(x, cfg$mpfr_bits) else x
mp_cos  <- function(x) if (cfg$precision == "mpfr") asNumeric(Rmpfr::cos(x)) else cos(x)
mp_exp  <- function(x) if (cfg$precision == "mpfr") asNumeric(Rmpfr::exp(x)) else exp(x)

# ============================================================
# 3) H_sigma builder (chunked, with weights, mpfr-capable)
#    H(t) = Σ |ρ|^{-α} e^{(β−σ)t} cos(γt) · weight
# ============================================================
build_H_sigma <- function(zeros, sigma = 0.6, alpha = 1,
                          Tmax = 1e3, n_t = 4000, chunk = 250L) {
  z <- copy(zeros)
  if (!"weight" %in% names(z)) z[, weight := 1.0]
  
  t <- seq(0, Tmax, length.out = n_t)
  H <- numeric(length(t))
  if (nrow(z) == 0L) return(tibble(t = t, H = H))
  
  eta <- (z$beta - sigma)
  amp <- (sqrt(z$beta^2 + z$gamma^2))^(-alpha)
  
  # mpfr conversions only if requested (keeps doubles fast)
  t_m   <- to_mpfr(t)
  eta_m <- to_mpfr(eta)
  gam_m <- to_mpfr(z$gamma)
  w_m   <- to_mpfr(amp * z$weight)
  
  i <- 1L
  while (i <= nrow(z)) {
    j <- min(i + chunk - 1L, nrow(z))
    Ej <- mp_exp(outer(t_m, eta_m[i:j], `*`))
    Cj <- mp_cos(outer(t_m, gam_m[i:j], `*`))
    Hj <- (Ej * Cj) %*% as.matrix(w_m[i:j], ncol = 1)
    H  <- H + as.numeric(Hj)
    i <- j + 1L
  }
  tibble(t = t, H = H)
}

# ============================================================
# 4) Blockwise affine fit + action (with CENTRAL differences)
# ============================================================
fit_affine_block <- function(t_block, f_block, ridge = 1e-12) {
  n <- length(t_block)
  if (n < 3L || var(t_block) == 0) return(list(a = mean(f_block), b = 0))
  tc <- t_block - mean(t_block)
  X  <- cbind(1, tc)
  XtX <- crossprod(X)
  beta <- tryCatch(
    solve(XtX + diag(ridge, 2L), crossprod(X, f_block)),
    error = function(e) {
      coef <- qr.coef(qr(X), f_block)
      if (anyNA(coef)) c(mean(f_block), 0) else coef
    }
  )
  a <- beta[1] - beta[2] * mean(t_block)
  b <- beta[2]
  list(a = a, b = b)
}

central_diff <- function(x, t) {
  n <- length(x); out <- numeric(n)
  if (n == 1) return(0*x)
  if (n == 2) return(c( (x[2]-x[1])/(t[2]-t[1]), (x[2]-x[1])/(t[2]-t[1]) ))
  out[2:(n-1)] <- (x[3:n] - x[1:(n-2)]) / (t[3:n] - t[1:(n-2)])
  out[1] <- (x[2] - x[1]) / (t[2] - t[1])
  out[n] <- (x[n] - x[n-1]) / (t[n] - t[n-1])
  out
}

blockwise_functionals <- function(series, Delta, lambda) {
  t <- series$t; f <- series$H
  blk <- floor(t / Delta)
  df  <- tibble(t = t, f = f, blk = blk)
  
  res <- df %>%
    group_by(blk) %>%
    group_modify(~{
      tb <- .x$t; fb <- .x$f
      n  <- length(tb)
      if (n < 3L || var(tb) == 0) return(tibble(int_eps2 = 0, int_deps2 = 0))
      
      ab  <- fit_affine_block(tb, fb)
      S   <- ab$a + ab$b * tb
      eps <- fb - S
      
      # CENTRAL differences (stable), then trapezoid norms
      dS   <- central_diff(S,  tb)
      df_  <- central_diff(fb, tb)
      deps <- df_ - dS
      
      tibble(
        int_eps2  = trapz(tb, eps^2),
        int_deps2 = trapz(tb, deps^2)
      )
    }) %>% ungroup() %>%
    mutate(
      block_start = (blk) * Delta,
      block_end   = (blk + 1) * Delta,
      total_block = int_eps2 + lambda * int_deps2
    ) %>%
    mutate(
      cum_lambda_d = kahan_cumsum(lambda * int_deps2),
      cum_total    = kahan_cumsum(total_block)
    )
  
  list(
    per_block = res,
    totals = tibble(
      blocks   = n_distinct(blk),
      sum_eps2 = sum(res$int_eps2),
      sum_d2   = sum(res$int_deps2),
      F_lambda = sum(res$total_block)
    )
  )
}

# ============================================================
# 5) Off-line zero injector (weight=1 for the synthetic zero)
# ============================================================
inject_offline_zero <- function(zeros, sigma, eta = 5e-5, gamma0 = 100) {
  stopifnot(all(c("beta","gamma") %in% names(zeros)))
  fake <- data.table(beta = sigma + eta, gamma = gamma0, weight = 1.0)
  for (m in setdiff(names(zeros), names(fake))) fake[, (m) := NA]
  setcolorder(fake, names(zeros))
  rbindlist(list(zeros, fake), use.names = TRUE, fill = TRUE)
}

# ============================================================
# 6) Variance sweep (polynomial regime sanity)
# ============================================================
run_variance_regime <- function(zeros, sigma, alpha, T_grid, kappa, c_lambda,
                                m_min, min_nt, max_nt,
                                zeros_cap, zeros_pick, gamma_max, chunk) {
  z_use <- thin_zeros(zeros, cap = zeros_cap, method = zeros_pick, gamma_max = gamma_max)
  map_dfr(T_grid, function(Tmax) {
    Delta  <- kappa / log(Tmax)
    lambda <- c_lambda / (log(Tmax)^2)
    n_t    <- choose_nt(Tmax, Delta, m_min = m_min, min_nt = min_nt, max_nt = max_nt)
    
    H  <- build_H_sigma(z_use, sigma, alpha, Tmax, n_t, chunk = chunk)
    H$H <- H$H - mean(H$H)  # center for cross-T comparability
    
    F  <- blockwise_functionals(H, Delta, lambda)
    tibble(
      Tmax   = Tmax,
      Delta  = Delta,
      lambda = lambda,
      blocks = F$totals$blocks,
      sum_eps2 = F$totals$sum_eps2,
      sum_d2   = F$totals$sum_d2,
      F_lambda = F$totals$F_lambda,
      T_logT_loglogT = Tmax * log(Tmax) * log(log(Tmax)),
      T_log2T        = Tmax * (log(Tmax)^2)
    )
  })
}

# ============================================================
# 7) Adaptive κ for bias runs (ensure ≥ target samples/block)
#     avg samples/block ≈ n_t * κ / (T log T)
#     choose κ = clamp(target * T log T / n_t, [κ_min, κ_max])
# ============================================================
choose_kappa_bias <- function(Tmax, n_t, target, kmin, kmax) {
  k <- target * Tmax * log(Tmax) / n_t
  min(kmax, max(kmin, k))
}

# ============================================================
# 8) Single bias run {η, T} with stable TOTAL-Δ differencing
# ============================================================
run_bias_once <- function(zeros, sigma, alpha, Tmax, eta, kappa, c_lambda,
                          m_min, min_nt, max_nt_bias,
                          zeros_cap, zeros_pick, gamma_max,
                          gamma_mode, gamma0, gamma0_factor,
                          standardize, chunk, target_spb, kmin, kmax) {
  
  # provisional Δ, pick n_t, then adapt κ to hit target samples/block
  Delta0  <- kappa / log(Tmax)
  n_t0    <- choose_nt(Tmax, Delta0, m_min = m_min, min_nt = min_nt, max_nt = max_nt_bias)
  kappa_b <- choose_kappa_bias(Tmax, n_t0, target_spb, kmin, kmax)
  Delta   <- kappa_b / log(Tmax)
  n_t     <- choose_nt(Tmax, Delta, m_min = m_min, min_nt = min_nt, max_nt = max_nt_bias)
  
  lambda <- c_lambda / (log(Tmax)^2)
  
  # Synthetic zero frequency: auto + Nyquist guard (γ0 < 0.9 * π / Δt)
  if (gamma_mode == "auto") gamma0 <- gamma0_factor / Delta
  dt  <- Tmax / max(1, n_t - 1)
  nyq <- pi / dt
  gamma0 <- min(gamma0, 0.9 * nyq)
  
  # Thin zeros (keep your method; swap to thin_by_weight() if desired)
  z_use <- thin_zeros(zeros, cap = zeros_cap, method = zeros_pick, gamma_max = gamma_max)
  
  # Baseline series
  H0 <- build_H_sigma(z_use, sigma, alpha, Tmax, n_t, chunk = chunk)
  mu0 <- mean(H0$H); sd0 <- sd(H0$H)
  trf <- switch(standardize,
                center = function(x) x - mu0,
                none   = function(x) x,
                zscore = function(x) if (sd0 > 0) (x - mu0)/sd0 else (x - mu0))
  H0$H <- trf(H0$H)
  
  # Biased series (inject off-line zero)
  z1 <- inject_offline_zero(z_use, sigma, eta = eta, gamma0 = gamma0)
  H1 <- build_H_sigma(z1, sigma, alpha, Tmax, n_t, chunk = chunk)
  H1$H <- trf(H1$H)
  
  P0 <- blockwise_functionals(H0, Delta, lambda)$per_block %>% arrange(blk)
  P1 <- blockwise_functionals(H1, Delta, lambda)$per_block %>% arrange(blk)
  
  # Per-block TOTAL difference (NO CLIPPING)
  inc_diff_total <- P1$total_block - P0$total_block
  cum_diff_total <- kahan_cumsum(inc_diff_total)
  
  out_blocks <- tibble(
    blk = P0$blk, block_start = P0$block_start,
    inc_total_diff = inc_diff_total,
    cum_total_diff = cum_diff_total
  )
  
  list(
    meta = tibble(Tmax=Tmax, eta=eta, Delta=Delta, kappa_bias=kappa_b,
                  n_t=n_t, blocks = max(P0$blk)+1, gamma0=gamma0, lambda=lambda),
    blocks = out_blocks
  )
}

# ============================================================
# 9) Bias sweep across (η, T) with tail-slope on last third
# ============================================================
run_bias_sweep <- function(zeros, cfg) {
  grid <- expand.grid(eta = cfg$bias_eta_list, Tmax = cfg$bias_T_list)
  rows <- list(); all_blocks <- list()
  for (i in seq_len(nrow(grid))) {
    eta  <- grid$eta[i]; Tmax <- grid$Tmax[i]
    res <- run_bias_once(
      zeros, cfg$sigma, cfg$alpha, Tmax, eta, cfg$kappa, cfg$c_lambda,
      cfg$m_min, cfg$min_nt, cfg$max_nt_bias,
      cfg$zeros_cap, cfg$zeros_pick, cfg$gamma_max,
      cfg$bias_gamma_mode, cfg$bias_gamma0, cfg$gamma0_factor,
      cfg$standardize, cfg$chunk,
      cfg$target_samples_per_block, cfg$kappa_bias_min, cfg$kappa_bias_max
    )
    blocks <- res$blocks %>% mutate(eta=eta, Tmax=Tmax)
    meta   <- res$meta
    
    # Tail slope on TOTAL-Δ cumulative (last third of blocks)
    tail_mask <- blocks$block_start > (max(blocks$block_start) * 2/3)
    eps <- 1e-18
    
    emp_rate <- tryCatch({
      if (!any(tail_mask)) NA_real_ else {
        if (cfg$use_huber) {
          fit <- MASS::rlm(log(pmax(cum_total_diff, eps)) ~ block_start,
                           data = blocks[tail_mask, ], psi = MASS::psi.huber)
          coef(fit)[2]
        } else {
          coef(lm(log(pmax(cum_total_diff, eps)) ~ block_start,
                  data = blocks[tail_mask, ]))[2]
        }
      }
    }, error = function(e) NA_real_)
    
    rows[[i]] <- cbind(meta, tibble(emp_tail_rate = as.numeric(emp_rate), theory_2eta = 2*eta))
    all_blocks[[i]] <- blocks
  }
  list(summary = bind_rows(rows), blocks = bind_rows(all_blocks))
}

# ============================================================
# 10) RUN (variance + bias)
# ============================================================
if (cfg$precision == "mpfr") message(sprintf("Running with MPFR precision %d bits", cfg$mpfr_bits))

var_tbl <- run_variance_regime(
  zeros,
  sigma      = cfg$sigma,
  alpha      = cfg$alpha,
  T_grid     = cfg$variance_T_grid,
  kappa      = cfg$kappa,
  c_lambda   = cfg$c_lambda,
  m_min      = cfg$m_min,
  min_nt     = cfg$min_nt,
  max_nt     = cfg$max_nt,
  zeros_cap  = cfg$zeros_cap,
  zeros_pick = cfg$zeros_pick,
  gamma_max  = cfg$gamma_max,
  chunk      = cfg$chunk
)

fit_var1 <- lm(F_lambda ~ 0 + T_logT_loglogT, data = var_tbl)
fit_var2 <- lm(F_lambda ~ 0 + T_log2T,        data = var_tbl)

bias_res <- run_bias_sweep(zeros, cfg)

# ============================================================
# 11) PLOTS
# ============================================================
p_var1 <- ggplot(var_tbl, aes(T_logT_loglogT, F_lambda)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = expression(F[lambda] ~ " vs " ~ T%*%log(T)%*%log*log(T)),
       x = expression(T%.%log(T)%.%log*log(T)), y = expression(F[lambda])) +
  theme_minimal(base_size = 13)

p_var2 <- ggplot(var_tbl, aes(T_log2T, F_lambda)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = expression(F[lambda] ~ " vs " ~ T%.%log(T)^2),
       x = expression(T%.%log(T)^2), y = expression(F[lambda])) +
  theme_minimal(base_size = 13)

p_bias_grid <- ggplot(bias_res$blocks,
                      aes(block_start, pmax(cum_total_diff, 1e-18))) +
  geom_line(linewidth = 0.8, color = "steelblue") +
  scale_y_log10() +
  facet_grid(rows = vars(paste0("eta=", format(eta, digits = 3))),
             cols = vars(paste0("T=", label_number(accuracy=1)(Tmax)))) +
  labs(title = "TOTAL action Δ (β>σ − baseline), semilog",
       x = "t (block start)",
       y = expression(paste("Δ cum ( ‖ε‖^2 + ", lambda, "‖∂", t, "ε‖^2 )"))) +
  theme_minimal(base_size = 12)

p_bias_summary <- ggplot(bias_res$summary,
                         aes(theory_2eta, emp_tail_rate,
                             label = paste0("T=", label_number(accuracy=1)(Tmax)))) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_point(size = 2) +
  geom_text_repel(size = 3, min.segment.length = 0) +
  labs(title = "Empirical tail slope (TOTAL Δ) vs theory (2η)",
       x = "Theory 2η", y = "Empirical tail slope") +
  theme_minimal(base_size = 13)

print(p_var1); print(p_var2); print(p_bias_grid); print(p_bias_summary)

if (isTRUE(cfg$save_plots)) {
  ggsave(file.path(cfg$out_dir, "variance_T_logT_loglogT.png"), p_var1, width = 7.5, height = 5.2, dpi = 140)
  ggsave(file.path(cfg$out_dir, "variance_T_log2T.png"),        p_var2, width = 7.5, height = 5.2, dpi = 140)
  ggsave(file.path(cfg$out_dir, "bias_total_diff_grid.png"),    p_bias_grid, width = 10, height = 7, dpi = 140)
  ggsave(file.path(cfg$out_dir, "bias_tail_slope_vs_theory.png"), p_bias_summary, width = 7.5, height = 5.2, dpi = 140)
}

if (isTRUE(cfg$save_tables)) {
  write_csv(var_tbl,                 file.path(cfg$out_dir, "variance_table.csv"))
  write_csv(bias_res$summary,        file.path(cfg$out_dir, "bias_summary.csv"))
  write_csv(bias_res$blocks,         file.path(cfg$out_dir, "bias_blocks.csv"))
}

# ============================================================
# 12) CONSOLE SUMMARY
# ============================================================
cat("\n=== VARIANCE SUMMARY ===\n")
cat("Slope vs T log T log log T: ", round(coef(fit_var1), 6), "\n")
cat("Slope vs T log^2 T:         ", round(coef(fit_var2), 6), "\n")

cat("\n=== BIAS SWEEP (tail slopes on TOTAL Δ) ===\n")
print(
  bias_res$summary %>%
    mutate(across(c(emp_tail_rate, theory_2eta, Delta, lambda, gamma0, kappa_bias, n_t, blocks),
                  ~format(., digits = 6)))
)
