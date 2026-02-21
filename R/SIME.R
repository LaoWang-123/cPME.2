# SIME code
# Minimal SIME core loop:
# - pme1: fixed reference pme, provides f1_fun = pme1$embedding_map
# - init2: initialization list for pme2 side: centers, theta_hat, parameterization
# - lambda, eta: given scalars
# Returns: a pme-like object for updated f2 (plus history)

sime_fit_minimal <- function(
    pme1,
    init2,
    eta,
    lambda,
    # loop control
    epsilon = 0.05,
    max_iter = 100,
    SSD_ratio_threshold = 5,
    verbose = TRUE
) {
  # ----------------------------
  # 0) Extract
  # ----------------------------
  f1_fun <- pme1$embedding_map

  X2 <- as.matrix(init2$centers)                # I x 3
  U2 <- as.matrix(init2$parameterization)       # I x 2
  I  <- nrow(X2)

  theta <- init2$theta_hat
  if (is.null(theta)) theta <- rep(1, I)
  theta <- as.numeric(theta)

  # ----------------------------
  # 1) Build anchors (fixed)
  #    Ua fixed = init2 parameters
  # ----------------------------
  Ua <- U2
  Xa <- t(apply(Ua, 1, function(u) as.numeric(f1_fun(as.numeric(u)))))  # I x 3

  # combined dataset (data + anchor)
  X_all <- rbind(X2, Xa)    # 2I x 3

  # weights:
  # - data term uses theta_hat
  # - anchor penalty uses constant eta per anchor point
  w_all <- c(theta, rep(eta, I))
  W_all <- diag(w_all)

  # ----------------------------
  # 2) Embedding builder (same form as PME)
  # ----------------------------
  build_embedding <- function(spline_coefs, params_all) {
    I_all <- nrow(params_all)  # should be 2I
    d <- 2
    function(u) {
      u <- as.numeric(u)
      as.vector(
        (t(spline_coefs[1:I_all, , drop = FALSE]) %*% etaFunc(u, params_all, 4 - d)) +
          (t(spline_coefs[(I_all + 1):(I_all + d + 1), , drop = FALSE]) %*%
             matrix(c(1, u), ncol = 1))
      )
    }
  }

  # ----------------------------
  # 3) Initialize f2 with current U2 (anchors fixed)
  # ----------------------------
  params_all <- rbind(U2, Ua)  # 2I x 2
  I_all <- nrow(params_all)

  spline_coefs <- calc_coefficients(X_all, params_all, W_all, lambda)
  f2_embedding <- build_embedding(spline_coefs, params_all)
  U2 <- calc_params(f2_embedding, X2, U2, f_input = "vector")

  # Need to consider whether to use those to calculate SSD
  SSD <- calc_SSD(f2_fun, X2, U2)
  hist <- data.frame(iter = 0, SSD = SSD, SSD_ratio = NA_real_)

  # ----------------------------
  # 4) Main loop: update f2 coefs, then update ONLY U2 (project X2)
  # ----------------------------
  count <- 1
  SSD_ratio <- 10 * epsilon

  while ((SSD_ratio > epsilon) &&
         (SSD_ratio <= SSD_ratio_threshold) &&
         (count <= max_iter)) {

    SSD_prev <- SSD
    f_prev <- f2_embedding
    U2_prev <- U2
    coefs_prev <- spline_coefs

    # (A) update embedding (Ua fixed)
    params_all <- rbind(U2, Ua)
    spline_coefs <- calc_coefficients(X_all, params_all, W_all, lambda)
    f2_embedding <- build_embedding(spline_coefs, params_all)

    # (B) update ONLY pme2 params (anchors fixed)
    U2 <- calc_params(f2_embedding, X2, U2, f_input = "vector")

    # recompute SSD (on combined set)
    SSD <- calc_SSD(f2_embedding, X2, U2)

    SSD_ratio <- abs(SSD - SSD_prev) / SSD_prev
    hist <- rbind(hist, data.frame(iter = count, SSD = SSD, SSD_ratio = SSD_ratio))

    # same “blow-up rollback” pattern as PME
    if (SSD_ratio > SSD_ratio_threshold) {
      f2_embedding <- f_prev
      U2 <- U2_prev
      spline_coefs <- coefs_prev
      SSD <- SSD_prev
    }

    if (verbose) {
      if (exists("print_SSD", mode = "function")) {
        print_SSD(lambda, SSD, SSD_ratio, count)
      } else {
        message(sprintf("[SIME] iter=%d SSD=%.6g ratio=%.3g", count, SSD, SSD_ratio))
      }
    }

    count <- count + 1
  }

  # ----------------------------
  # 5) Return minimal “pme-like” result
  # ----------------------------
  list(
    embedding_map = f2_embedding,
    params_opt = U2,
    centers = X2,
    theta_hat = theta,
    anchors = list(Ua = Ua, Xa = Xa, eta = eta),
    tuning = lambda,
    kernel_coefs = spline_coefs[1:I_all, , drop = FALSE],
    polynomial_coefs = spline_coefs[(I_all + 1):(I_all + 3), , drop = FALSE], # d=2 => +3 rows
    history = hist
  )
}
