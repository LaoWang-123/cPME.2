# SIME code

#' SIME with anchor penalty + lambda screening
#'
#' @param f1_fun f1 function with vector input (u -> R^3)
#' @param init2 list: $centers (I x 3), $parameterization (I x 2),
#'   $theta_hat (optional), $km (required for calc_msd)
#' @param data2 raw data matrix (n x D) used in calc_msd (PME criterion)
#' @param eta anchor weight (penalty strength)
#' @param lambda tuning vector (can be length 1 or longer, like PME exp(-15:5))
#' @param d intrinsic dimension (default 2 for surface)
#' @param epsilon,max_iter,SSD_ratio_threshold same as PME loop controls
#' @param verbose print_SSD like PME
#'
#' @returns list: optimal fit + full tuning path (MSD, coefs, parameterization, embeddings)
#' @export
SIME <- function(
    f1_fun,
    init2,
    data2,
    eta = 0.05,
    lambda = exp(-15:5),
    d = 2,
    epsilon = 1,
    max_iter = 100,
    SSD_ratio_threshold = 5,
    verbose = FALSE
) {

  if (!exists("calc_msd", mode = "function"))
    stop("[SIMEpme] calc_msd() not found. Load PME utilities first.")
  if (is.null(init2$km))
    stop("[SIMEpme] init2$km is required for calc_msd().")

  data2 <- as.matrix(data2)
  D <- ncol(data2)

  # ----------------------------
  # Extract init2
  # ----------------------------
  X2 <- as.matrix(init2$centers)                # I x 3
  U2_init <- as.matrix(init2$parameterization)  # I x 2
  I <- nrow(X2)

  theta <- init2$theta_hat
  if (is.null(theta)) theta <- rep(1, I)
  theta <- as.numeric(theta)

  # ----------------------------
  # Build anchors (fixed)
  # ----------------------------
  Ua <- U2_init
  Xa <- t(apply(Ua, 1, function(u) as.numeric(f1_fun(as.numeric(u)))))  # I x 3

  X_all <- rbind(X2, Xa)                 # 2I x 3
  weights_all <- diag(c(theta, rep(eta/I, I)))
  I_all <- nrow(X_all)                   # 2I

  # storage like PME
  mse <- vector()
  coefs <- list()
  parameterization <- list()
  embeddings <- list()

  # ----------------------------
  # PME-style tuning loop over lambda
  # ----------------------------
  for (tuning_idx in 1:length(lambda)) {

    # params refers to *data params* (I x 2), like PME
    params <- U2_init

    # full params includes anchors (fixed)
    params_all <- rbind(params, Ua)

    spline_coefs <- calc_coefficients(
      X_all,
      params_all,
      weights_all,
      lambda[tuning_idx]
    )

    # embedding uses current spline_coefs + current params_all (like PME uses params)
    f_embedding <- function(parameters) {
      as.vector(
        (t(spline_coefs[1:I_all, , drop = FALSE]) %*% etaFunc(parameters, params_all, 4 - d)) +
          (t(spline_coefs[(I_all + 1):(I_all + d + 1), , drop = FALSE]) %*%
             matrix(c(1, parameters), ncol = 1))
      )
    }

    f0 <- f_embedding

    # project ONLY centers X2 -> update params (anchors fixed)
    params <- calc_params(f_embedding, X2, params, "vector")
    params_all <- rbind(params, Ua)

    SSD <- calc_SSD(f_embedding, X2, params)

    count <- 1
    SSD_ratio <- 10 * epsilon

    # PME-style inner loop
    while ((SSD_ratio > epsilon) &
           (SSD_ratio <= SSD_ratio_threshold) &
           (count <= (max_iter - 1))) {

      SSD_prev <- SSD
      f0 <- f_embedding
      params_prev <- params
      coefs_prev <- spline_coefs
      params_all_prev <- params_all

      # update coefs with UPDATED params_all (crucial PME pattern)
      spline_coefs <- calc_coefficients(
        X_all,
        params_all,
        weights_all,
        lambda[tuning_idx]
      )

      f_embedding <- function(parameters) {
        as.vector(
          (t(spline_coefs[1:I_all, , drop = FALSE]) %*% etaFunc(parameters, params_all, 4 - d)) +
            (t(spline_coefs[(I_all + 1):(I_all + d + 1), , drop = FALSE]) %*%
               matrix(c(1, parameters), ncol = 1))
        )
      }

      # update ONLY params (I x 2), anchors fixed
      params <- calc_params(f_embedding, X2, params, "vector")
      params_all <- rbind(params, Ua)

      SSD <- calc_SSD(f_embedding, X2, params)
      SSD_ratio <- abs(SSD - SSD_prev) / SSD_prev
      count <- count + 1

      # rollback like PME
      if (SSD_ratio > SSD_ratio_threshold) {
        f_embedding <- f0
        params <- params_prev
        params_all <- params_all_prev
        spline_coefs <- coefs_prev
        SSD <- SSD_prev
      }

      if (verbose == TRUE) {
        if (exists("print_SSD", mode = "function")) {
          print_SSD(lambda[tuning_idx], SSD, SSD_ratio, count)
        } else {
          message(sprintf("[SIME] lambda=%g SSD=%g ratio=%g iter=%d",
                          lambda[tuning_idx], SSD, SSD_ratio, count))
        }
      }
    }

    # PME criterion: MSD on raw data
    mse[tuning_idx] <- calc_msd(data2, init2$km, f_embedding, params, D, d) # Maybe need to check later

    if (verbose == TRUE) {
      message(sprintf("When lambda = %s, MSD = %s.",
                      as.character(lambda[tuning_idx]),
                      as.character(mse[tuning_idx])))
    }

    # store path objects like PME
    coefs[[tuning_idx]] <- spline_coefs
    parameterization[[tuning_idx]] <- params

    embeddings[[tuning_idx]] <- function(parameters) {
        as.vector(
          (t(coefs[[tuning_idx]][1:I_all, , drop = FALSE]) %*% etaFunc(parameters, rbind(parameterization[[tuning_idx]], Ua), 4 - d)) +
            (t(coefs[[tuning_idx]][(I_all + 1):(I_all + d + 1), , drop = FALSE]) %*%
               matrix(c(1, parameters), ncol = 1))
        )
      }


    # PME early stop: last 4 MSD nondecreasing => break
    if (tuning_idx >= 4) {
      if (!is.unsorted(mse[(tuning_idx - 3):tuning_idx])) {
        break
      }
    }
  }

  optimal_idx <- min(which(mse == min(mse)))

  coefs_opt <- coefs[[optimal_idx]]
  params_opt <- parameterization[[optimal_idx]]

  embedding_opt <- function(parameters) {
    as.vector(
      (t(coefs_opt[1:I_all, , drop = FALSE]) %*% etaFunc(parameters, rbind(params_opt, Ua), 4 - d)) +
        (t(coefs_opt[(I_all + 1):(I_all + d + 1), , drop = FALSE]) %*%
           matrix(c(1, parameters), ncol = 1))
    )
  }

  # return pme-like + sime-specific extras
  list(
    embedding_map = embedding_opt,
    params_opt = params_opt,
    centers = X2,
    theta_hat = theta,
    anchors = list(Ua = Ua, Xa = Xa, eta = eta),
    knots = init2$km,
    tuning = lambda[optimal_idx],
    MSD = mse,
    coefs = coefs,
    parameterization = parameterization,
    tuning_vec = lambda,
    embeddings = embeddings
  )
}
