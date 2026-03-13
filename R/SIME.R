# SIME code

#' SIME with anchor penalty + lambda screening
#' This code is run with given eta and chose lambda automatically
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



############################################################
## Validation projection MSD
############################################################
calc_validation_proj_msd <- function(val_data, fit_obj) {
  val_data <- as.matrix(val_data)
  if (nrow(val_data) == 0) return(NA_real_)
  embedding_map <- fit_obj$embedding_map

  # calc_params with pca initialization
  init_params_val <- pme_initial_guess(X = val_data, d = 2, method = "pca")

  val_params_opt <- calc_params(
    f = embedding_map,
    X = val_data,
    init_params = init_params_val,
    f_input = "vector"
  )

  val_proj <- t(apply(val_params_opt, 1, function(u) {
    as.numeric(embedding_map(as.numeric(u)))
  }))

  mean(rowSums((val_data - val_proj)^2))
}


#' Cross-validation selection of the structural weight \eqn{\eta} for SIME
#'
#' Performs K-fold cross-validation to select the optimal structural weight
#' \eqn{\eta} for the SIME (Structure-Informed Manifold Estimation) model.
#' For each candidate \eqn{\eta}, the procedure:
#'
#' \enumerate{
#' \item Splits the data into K folds.
#' \item For each fold, initializes the PME manifold using \code{initialize_pme()}.
#' \item On the first fold, runs \code{SIME()} with the full \code{lambda} grid
#'       to determine the optimal smoothing parameter.
#' \item Uses this selected \code{lambda} for the remaining folds to reduce
#'       computation time.
#' \item Computes validation projection MSD (mean squared distance) between
#'       validation data and the fitted manifold.
#' }
#'
#' The \eqn{\eta} with the smallest average validation MSD is selected, and the
#' final SIME model is refit on the full dataset using the selected
#' \eqn{\eta} and corresponding \code{lambda}.
#'
#' @param f1_fun Function representing the reference manifold (typically
#'   obtained from MRI PME). Takes a parameter vector \eqn{u} and returns a
#'   point in the ambient space.
#'
#' @param f2_fun Function representing the current PET manifold estimate (eta=0) used
#'   to update the parameterization of initialization centers.
#'
#' @param data2 A numeric matrix of observations used to fit the PET manifold.
#'   Rows correspond to data points and columns correspond to ambient
#'   coordinates.
#'
#' @param eta_vec Numeric vector of candidate \eqn{\eta} values controlling the
#'   strength of structural anchoring toward \code{f1_fun}.
#'
#' @param lambda Numeric vector of smoothing parameters used by \code{SIME()}.
#'   The first fold selects the optimal value.
#'
#' @param K Integer. Number of cross-validation folds.
#'
#' @param d Intrinsic dimension of the manifold (default = 2).
#'
#' @param epsilon Convergence threshold for the inner SIME optimization loop.
#'
#' @param max_iter Maximum number of iterations allowed in the SIME fitting
#'   procedure.
#'
#' @param SSD_ratio_threshold Threshold used for stability control during
#'   the SIME optimization iterations.
#'
#' @param init_args_f2 A list of arguments passed to \code{initialize_pme()}
#'   when generating fold-specific manifold initializations.
#'
#' @param seed Optional random seed used for reproducible fold assignment.
#'
#' @param verbose Logical. If TRUE, prints progress messages during
#'   cross-validation.
#'
#' @return A list containing:
#' \describe{
#'   \item{best_eta}{Selected value of \eqn{\eta}.}
#'   \item{best_lambda}{Selected smoothing parameter.}
#'   \item{eta_mean_msd}{Average validation MSD for each candidate \eqn{\eta}.}
#'   \item{eta_fold_msd}{Fold-wise validation MSD values.}
#'   \item{eta_lambda_star}{Selected \code{lambda} for each \eqn{\eta}.}
#'   \item{final_fit}{Final SIME model fitted on the full dataset.}
#' }
#'
#' @details
#' Initialization is performed once per fold using \code{initialize_pme()} and
#' reused across candidate \eqn{\eta} values to avoid redundant computation.
#'
#' Validation error is measured using projection mean squared distance (MSD),
#' computed by projecting validation points onto the fitted manifold.
#'
#' @examples
#' \dontrun{
#' cv_fit <- SIME_cv(
#'   f1_fun = f1_fun,
#'   f2_fun = f2_fun,
#'   data2 = data2,
#'   eta_vec = c(0.01, 0.05, 0.1, 0.5),
#'   lambda = exp(-15:5),
#'   K = 5
#' )
#'
#' cv_fit$best_eta
#' cv_fit$best_lambda
#' }
#'
#' @export
SIME_cv <- function(
    f1_fun,
    f2_fun,
    data2,
    eta_vec = c(exp(-15:-2),0.1,exp(-1:5)),
    lambda = exp(-15:5),
    K = 5,
    d = 2,
    epsilon = 1,
    max_iter = 100,
    SSD_ratio_threshold = 5,
    init_args_f2 = list(
      min_clusters = 10,
      alpha = 0.01,
      max_clusters = 100,
      algorithm = "isomap",
      rescale = FALSE,
      component_type = "centers",
      subsample_size = 5,
      d=2
    ),
    seed = NULL,
    verbose = TRUE
) {
  data2 <- as.matrix(data2)
  n <- nrow(data2)

  if (length(eta_vec) < 1) stop("eta_vec must contain at least one value.")
  if (K < 2) stop("K must be at least 2.")
  if (n < K) stop("n must be >= K.")

  if (!is.null(seed)) set.seed(seed)

  ##########################################################
  ## Make folds once
  ##########################################################
  perm <- sample(seq_len(n))
  fold_id <- rep(seq_len(K), length.out = n)
  folds <- split(perm, fold_id)

  ##########################################################
  ## Precompute train/val split and init2 for each fold once
  ##########################################################
  train_data_list <- vector("list", K)
  val_data_list   <- vector("list", K)
  init2_list      <- vector("list", K)

  for (k in seq_len(K)) {
    val_idx <- folds[[k]]
    train_idx <- setdiff(seq_len(n), val_idx)

    train_data_list[[k]] <- data2[train_idx, , drop = FALSE]
    val_data_list[[k]]   <- data2[val_idx, , drop = FALSE]

    if (verbose) {
      message("==================================================")
      message(sprintf("Precomputing init2 for fold %d", k))
      message("==================================================")
    }

    init2_k <- do.call(
      initialize_pme,
      c(list(x = train_data_list[[k]]),
        init_args_f2)
    )

    init2_k$parameterization <- calc_params(
      f = f2_fun,
      X = init2_k$centers,
      init_params = init2_k$parameterization,
      f_input = "uv"
    )

    init2_list[[k]] <- init2_k
  }

  ##########################################################
  ## Storage over eta
  ##########################################################
  eta_mean_msd <- rep(NA_real_, length(eta_vec))
  eta_fold_msd <- vector("list", length(eta_vec))
  eta_lambda_star <- rep(NA_real_, length(eta_vec))

  ##########################################################
  ## Loop over eta
  ##########################################################
  for (e in seq_along(eta_vec)) {
    eta_now <- eta_vec[e]

    if (verbose) {
      message("==================================================")
      message(sprintf("Evaluating eta = %g", eta_now))
      message("==================================================")
    }

    fold_msd <- rep(NA_real_, K)

    #######################################################
    ## Fold 1: full lambda vector, let SIME choose lambda
    #######################################################
    if (verbose) {
      message(sprintf("eta = %g, fold 1: running full lambda grid", eta_now))
    }

    fit1 <- SIME(
      f1_fun = f1_fun,
      init2 = init2_list[[1]],
      data2 = train_data_list[[1]],
      eta = eta_now,
      lambda = lambda,
      d = d,
      epsilon = epsilon,
      max_iter = max_iter,
      SSD_ratio_threshold = SSD_ratio_threshold,
      verbose = verbose
    )

    lambda_star <- fit1$tuning
    fold_msd[1] <- calc_validation_proj_msd(val_data_list[[1]], fit1)

    eta_lambda_star[e] <- lambda_star

    if (verbose) {
      message(sprintf(
        "eta = %g, fold 1 selected lambda = %g, val MSD = %g",
        eta_now, lambda_star, fold_msd[1]
      ))
    }

    #######################################################
    ## Folds 2..K: fixed lambda_star
    #######################################################
    if (K >= 2) {
      for (k in 2:K) {
        if (verbose) {
          message(sprintf(
            "eta = %g, fold %d: running fixed lambda = %g",
            eta_now, k, lambda_star
          ))
        }

        fit_k <- SIME(
          f1_fun = f1_fun,
          init2 = init2_list[[k]],
          data2 = train_data_list[[k]],
          eta = eta_now,
          lambda = lambda_star,
          d = d,
          epsilon = epsilon,
          max_iter = max_iter,
          SSD_ratio_threshold = SSD_ratio_threshold,
          verbose = verbose
        )

        fold_msd[k] <- calc_validation_proj_msd(val_data_list[[k]], fit_k)

        if (verbose) {
          message(sprintf(
            "eta = %g, fold %d: val MSD = %g",
            eta_now, k, fold_msd[k]
          ))
        }
      }
    }

    eta_fold_msd[[e]] <- fold_msd
    eta_mean_msd[e] <- mean(fold_msd, na.rm = TRUE)

    if (verbose) {
      message(sprintf(
        "eta = %g, mean CV MSD = %g",
        eta_now, eta_mean_msd[e]
      ))
    }

    #######################################################
    ## Early stop over eta: last 4 mean CV MSD nondecreasing
    #######################################################
    if (e >= 4) {
      if (!is.unsorted(eta_mean_msd[(e - 3):e])) {
        if (verbose) {
          message("Early stopping on eta: last 4 mean CV MSD values are nondecreasing.")
        }
        break
      }
    }
  }

  ##########################################################
  ## Select best eta
  ##########################################################
  best_eta_idx <- min(which(eta_mean_msd == min(eta_mean_msd, na.rm = TRUE)))
  best_eta <- eta_vec[best_eta_idx]
  best_lambda <- eta_lambda_star[best_eta_idx]

  if (verbose) {
    message("==================================================")
    message(sprintf("Best eta = %g", best_eta))
    message(sprintf("Corresponding lambda = %g", best_lambda))
    message("==================================================")
  }

  ##########################################################
  ## Final refit on full data
  ##########################################################
  init2_full <- do.call(
    initialize_pme,
    c(list(x = data2),
      init_args_f2)
  )

  init2_full$parameterization <- calc_params(
    f = f2_fun,
    X = init2_full$centers,
    init_params = init2_full$parameterization,
    f_input = "uv"
  )

  final_fit <- SIME(
    f1_fun = f1_fun,
    init2 = init2_full,
    data2 = data2,
    eta = best_eta,
    lambda = best_lambda,
    d = d,
    epsilon = epsilon,
    max_iter = max_iter,
    SSD_ratio_threshold = SSD_ratio_threshold,
    verbose = verbose
  )

  list(
    folds = folds,
    train_data_list = train_data_list,
    val_data_list = val_data_list,
    init2_list = init2_list,
    eta_vec = eta_vec,
    lambda_grid = lambda,
    eta_mean_msd = eta_mean_msd,
    eta_fold_msd = eta_fold_msd,
    eta_lambda_star = eta_lambda_star,
    best_eta = best_eta,
    SIME_final = final_fit
  )
}
