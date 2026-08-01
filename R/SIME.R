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

#' Select eta by threshold rule with early stopping
#'
#' Fit SIME along an increasing eta grid. Use the baseline MSD from f2_fun
#' as a reference scale, define threshold = c * baseline_msd, and select
#' the largest eta whose fitted SIME MSD does not exceed the threshold.
#'
#' Early stopping rule:
#' stop as soon as the fitted MSD exceeds the threshold.
#'
#' @param f1_fun Reference manifold function (e.g. MRI PME embedding map).
#' @param f2_fun Baseline target manifold function (e.g. PET PME embedding map).
#' @param data2 Raw target data matrix (n x D).
#' @param c Positive multiplier for threshold. threshold = c * baseline_msd.
#' @param eta_vec Candidate eta values. Will be sorted increasingly.
#' @param lambda Candidate lambda values for SIME().
#' @param d Intrinsic dimension.
#' @param epsilon,max_iter,SSD_ratio_threshold Same as in SIME().
#' @param init_args_f2 Arguments passed to initialize_pme().
#' @param seed Optional random seed.
#' @param verbose Logical.
#'
#' @return A list containing the baseline MSD, threshold, eta path,
#'   selected eta, and selected SIME fit.
#' @export
SIME_select <- function(
    f1_fun,
    f2_fun,
    data2,
    c = 1.3,
    eta_vec = c(exp(-15:-2), 0.1, exp(-1:5)),
    lambda = exp(-15:5),
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
      d = 2
    ),
    seed = NULL,
    verbose = TRUE
) {


  data2 <- as.matrix(data2)
  n <- nrow(data2)
  D <- ncol(data2)


  eta_vec <- sort(unique(as.numeric(eta_vec)))   # increasing eta
  lambda  <- as.numeric(lambda)

  if (!is.null(seed)) set.seed(seed)

  # ----------------------------
  # initialize once on full data
  # ----------------------------
  if (verbose) {
    message("==================================================")
    message("Initializing PME on full data")
    message("==================================================")
  }

  init2 <- do.call(
    initialize_pme,
    c(list(x = data2), init_args_f2)
  )

  # align center parameterization to f2_fun
  init2$parameterization <- calc_params(
    f = f2_fun,
    X = init2$centers,
    init_params = init2$parameterization,
    f_input = "vector"
  )

  # ----------------------------
  # baseline MSD from f2_fun
  # ----------------------------
  baseline_msd <- calc_msd(x = data2, km = init2$km,f = f2_fun,t = init2$parameterization,D = D,d = d)

  threshold <- c * baseline_msd

  if (verbose) {
    message(sprintf("Baseline MSD = %g", baseline_msd))
    message(sprintf("Threshold = %g * %g = %g", c, baseline_msd, threshold))
  }

  # ----------------------------
  # storage
  # ----------------------------
  eta_fit_msd <- rep(NA_real_, length(eta_vec))
  eta_lambda_star <- rep(NA_real_, length(eta_vec))
  eta_fit_list <- vector("list", length(eta_vec))
  is_admissible <- rep(FALSE, length(eta_vec))

  selected_idx <- NA_integer_

  # ----------------------------
  # loop over eta, early stop when MSD > threshold
  # ----------------------------
  for (e in seq_along(eta_vec)) {
    eta_now <- eta_vec[e]

    fit_e <- SIME(
      f1_fun = f1_fun,
      init2 = init2,
      data2 = data2,
      eta = eta_now,
      lambda = lambda,
      d = d,
      epsilon = epsilon,
      max_iter = max_iter,
      SSD_ratio_threshold = SSD_ratio_threshold,
      verbose = FALSE
    )

    idx_star <- match(fit_e$tuning, fit_e$tuning_vec)
    if (is.na(idx_star)) idx_star <- which.min(fit_e$MSD)

    fit_msd_now <- fit_e$MSD[idx_star]

    eta_fit_msd[e] <- fit_msd_now
    eta_lambda_star[e] <- fit_e$tuning
    eta_fit_list[[e]] <- fit_e

    if (fit_msd_now <= threshold) {
      is_admissible[e] <- TRUE
      selected_idx <- e
    } else {
      is_admissible[e] <- FALSE

      if (verbose) {
        message(sprintf(
          "eta = %g, lambda* = %g, MSD = %g > threshold = %g. Early stopping.",
          eta_now, fit_e$tuning, fit_msd_now, threshold
        ))
      }
      break
    }

    if (verbose) {
      message(sprintf(
        "eta = %g, lambda* = %g, MSD = %g <= threshold = %g",
        eta_now, fit_e$tuning, fit_msd_now, threshold
      ))
    }
  }

  # ----------------------------
  # prepare output
  # ----------------------------

  selected_eta <- eta_vec[selected_idx]
  selected_fit <- eta_fit_list[[selected_idx]]

  if (verbose) {
    message("==================================================")
    message(sprintf("Selected eta = %g", selected_eta))
    message(sprintf("Selected lambda = %g", eta_lambda_star[selected_idx]))
    message(sprintf("Selected MSD = %g", eta_fit_msd[selected_idx]))
    message("==================================================")
  }

  list(
    rule = "largest eta with SIME MSD <= c * baseline MSD; early stop at first exceedance",
    baseline_msd = baseline_msd,
    threshold = threshold,
    c = c,
    eta_fit_msd = eta_fit_msd,
    eta_lambda_star = eta_lambda_star,
    admissible_eta = eta_vec[which(is_admissible)],
    selected_eta = selected_eta,
    selected_fit = selected_fit
  )
}


#' Select eta for one or more threshold multipliers along one SIME path
#'
#' This is an experimental multi-\code{c} wrapper for \code{SIME_select()}.
#' For a single \code{c}, it returns the same core fields as \code{SIME_select()}
#' and adds optional shrinkage metrics. For multiple \code{c} values, it fits
#' each eta at most once and reuses the fitted path when moving from smaller to
#' larger thresholds.
#'
#' @param f1_fun Reference manifold function.
#' @param f2_fun Baseline target manifold function.
#' @param data2 Raw target data matrix (n x D).
#' @param c Positive threshold multiplier(s). Values are evaluated increasingly.
#' @param eta_vec Candidate eta values. Will be sorted increasingly.
#' @param lambda Candidate lambda values for \code{SIME()}.
#' @param d Intrinsic dimension.
#' @param epsilon,max_iter,SSD_ratio_threshold Same as in \code{SIME()}.
#' @param init_args_f2 Arguments passed to \code{initialize_pme()}.
#' @param seed Optional random seed.
#' @param output_mode Output size mode. \code{"compact"} keeps only selected
#'   unique eta fits; \code{"full"} keeps every evaluated eta fit.
#' @param shrinkage_grid Parameter grid used for shrinkage. If \code{NULL},
#'   a 40 by 40 disk grid from \code{make_uv_grid()} is used.
#' @param shrinkage_seed Seed passed to \code{calc_correspondence_msd()}.
#' @param verbose Logical.
#'
#' @return A list with \code{summary}, \code{eta_path}, and \code{fits}. In
#'   compact mode, \code{fits} stores only selected unique eta fits. In full
#'   mode, \code{fits} stores every evaluated eta fit.
#' @export
SIME_select_path <- function(
    f1_fun,
    f2_fun,
    data2,
    c = 1.3,
    eta_vec = c(exp(-15:-2), 0.1, exp(-1:5)),
    lambda = exp(-15:5),
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
      d = 2
    ),
    seed = NULL,
    output_mode = c("compact", "full"),
    shrinkage_grid = NULL,
    shrinkage_seed = 123,
    verbose = TRUE
) {

  output_mode <- match.arg(output_mode)

  data2 <- as.matrix(data2)
  D <- ncol(data2)

  c_input <- as.numeric(c)
  if (any(!is.finite(c_input)) || any(c_input <= 0)) {
    stop("c must contain positive finite values.")
  }
  c_values <- sort(unique(c_input))

  eta_vec <- sort(unique(as.numeric(eta_vec)))
  lambda <- as.numeric(lambda)

  if (!is.null(seed)) set.seed(seed)

  if (is.null(shrinkage_grid)) {
    shrinkage_grid <- make_uv_grid(n_u = 40, n_v = 40, grid_type = "disk")
  }

  compute_saved_shrinkage <- function(fits) {
    empty <- data.frame(
      fit_ref = character(0),
      sime_to_f1_msd = numeric(0),
      sime_to_f2_msd = numeric(0),
      f2_to_f1_msd = numeric(0),
      shrinkage = numeric(0),
      stringsAsFactors = FALSE
    )
    if (length(fits) == 0L) {
      return(empty)
    }

    f1_data <- generate_surface_data(
      f = f1_fun,
      noise_sd = 0,
      seed = shrinkage_seed,
      uv_grid = shrinkage_grid
    )
    f2_data <- generate_surface_data(
      f = f2_fun,
      noise_sd = 0,
      seed = shrinkage_seed,
      uv_grid = shrinkage_grid
    )

    if (!all(dim(f1_data$XYZ) == dim(f2_data$XYZ))) {
      stop("Reference and baseline manifolds do not have matching output dimensions.")
    }

    f2_to_f1_msd <- mean(rowSums((f2_data$XYZ - f1_data$XYZ)^2))

    rows <- lapply(names(fits), function(fit_ref) {
      fit <- fits[[fit_ref]]
      sime_data <- generate_surface_data(
        f = fit$embedding_map,
        noise_sd = 0,
        seed = shrinkage_seed,
        uv_grid = shrinkage_grid
      )
      if (!all(dim(sime_data$XYZ) == dim(f1_data$XYZ))) {
        stop("SIME and reference manifolds do not have matching output dimensions.")
      }

      sime_to_f1_msd <- mean(rowSums((sime_data$XYZ - f1_data$XYZ)^2))
      sime_to_f2_msd <- mean(rowSums((sime_data$XYZ - f2_data$XYZ)^2))
      shrinkage_ratio <- if (isTRUE(all.equal(f2_to_f1_msd, 0))) {
        NA_real_
      } else {
        1 - sime_to_f1_msd / f2_to_f1_msd
      }

      data.frame(
        fit_ref = fit_ref,
        sime_to_f1_msd = sime_to_f1_msd,
        sime_to_f2_msd = sime_to_f2_msd,
        f2_to_f1_msd = f2_to_f1_msd,
        shrinkage = shrinkage_ratio,
        stringsAsFactors = FALSE
      )
    })

    do.call(rbind, rows)
  }

  if (verbose) {
    message("==================================================")
    message("Initializing PME on full data")
    message("==================================================")
  }

  init2 <- do.call(
    initialize_pme,
    c(list(x = data2), init_args_f2)
  )

  init2$parameterization <- calc_params(
    f = f2_fun,
    X = init2$centers,
    init_params = init2$parameterization,
    f_input = "vector"
  )

  baseline_msd <- calc_msd(
    x = data2,
    km = init2$km,
    f = f2_fun,
    t = init2$parameterization,
    D = D,
    d = d
  )

  thresholds <- c_values * baseline_msd

  eta_fit_msd <- rep(NA_real_, length(eta_vec))
  eta_lambda_star <- rep(NA_real_, length(eta_vec))
  eta_fit_list <- vector("list", length(eta_vec))
  eta_error_message <- rep(NA_character_, length(eta_vec))
  names(eta_fit_list) <- paste0("eta_", seq_along(eta_vec))

  selected_idx <- rep(NA_integer_, length(c_values))
  stop_idx <- rep(NA_integer_, length(c_values))
  c_error_message <- rep(NA_character_, length(c_values))
  names(selected_idx) <- names(stop_idx) <- as.character(c_values)

  eta_start <- 1L
  for (i in seq_along(c_values)) {
    c_now <- c_values[i]
    threshold_now <- thresholds[i]
    if (i > 1L) {
      selected_idx[i] <- selected_idx[i - 1L]
    }

    if (verbose) {
      message(sprintf(
        "Selecting eta for c = %g; threshold = %g * %g = %g",
        c_now, c_now, baseline_msd, threshold_now
      ))
    }

    if (eta_start > length(eta_vec)) {
      selected_idx[i] <- length(eta_vec)
      next
    }

    for (e in eta_start:length(eta_vec)) {
      eta_now <- eta_vec[e]

      if (is.null(eta_fit_list[[e]])) {
        fit_e <- tryCatch(
          SIME(
            f1_fun = f1_fun,
            init2 = init2,
            data2 = data2,
            eta = eta_now,
            lambda = lambda,
            d = d,
            epsilon = epsilon,
            max_iter = max_iter,
            SSD_ratio_threshold = SSD_ratio_threshold,
            verbose = FALSE
          ),
          error = function(err) err
        )

        if (inherits(fit_e, "error")) {
          eta_error_message[e] <- conditionMessage(fit_e)
          c_error_message[i] <- eta_error_message[e]
          stop_idx[i] <- e
          if (verbose) {
            message(sprintf(
              "eta = %g failed: %s. Return NA for this c and continue.",
              eta_now, eta_error_message[e]
            ))
          }
          eta_start <- e + 1L
          break
        }

        idx_star <- match(fit_e$tuning, fit_e$tuning_vec)
        if (is.na(idx_star)) idx_star <- which.min(fit_e$MSD)

        eta_fit_msd[e] <- fit_e$MSD[idx_star]
        eta_lambda_star[e] <- fit_e$tuning
        eta_fit_list[[e]] <- fit_e
      } else if (!is.na(eta_error_message[e])) {
        c_error_message[i] <- eta_error_message[e]
        stop_idx[i] <- e
        eta_start <- e + 1L
        break
      }

      if (eta_fit_msd[e] <= threshold_now) {
        selected_idx[i] <- e
      } else {
        stop_idx[i] <- e
        if (verbose) {
          message(sprintf(
            "eta = %g, lambda* = %g, MSD = %g > threshold = %g. Move to next c.",
            eta_now, eta_lambda_star[e], eta_fit_msd[e], threshold_now
          ))
        }
        break
      }

      if (verbose) {
        message(sprintf(
          "eta = %g, lambda* = %g, MSD = %g <= threshold = %g",
          eta_now, eta_lambda_star[e], eta_fit_msd[e], threshold_now
        ))
      }
    }

    if (!is.na(selected_idx[i])) {
      eta_start <- selected_idx[i] + 1L
    }
  }

  selected_eta <- ifelse(is.na(selected_idx), NA_real_, eta_vec[selected_idx])
  selected_lambda <- ifelse(is.na(selected_idx), NA_real_, eta_lambda_star[selected_idx])
  selected_msd <- ifelse(is.na(selected_idx), NA_real_, eta_fit_msd[selected_idx])

  fit_ref <- ifelse(is.na(selected_idx), NA_character_, names(eta_fit_list)[selected_idx])

  summary <- data.frame(
    c = c_values,
    baseline_msd = baseline_msd,
    threshold = thresholds,
    selected_eta = selected_eta,
    selected_lambda = selected_lambda,
    selected_msd = selected_msd,
    selected_fit_ref = as.character(fit_ref),
    error_message = c_error_message,
    stringsAsFactors = FALSE
  )

  evaluated <- !vapply(eta_fit_list, is.null, logical(1))
  eta_path <- data.frame(
    eta_index = seq_along(eta_vec),
    eta = eta_vec,
    fit_msd = eta_fit_msd,
    lambda_star = eta_lambda_star,
    evaluated = evaluated,
    error_message = eta_error_message,
    stringsAsFactors = FALSE
  )

  keep_idx <- which(evaluated)
  if (output_mode == "compact") {
    keep_idx <- sort(unique(selected_idx[!is.na(selected_idx)]))
  }
  fits <- eta_fit_list[keep_idx]
  saved_shrinkage <- compute_saved_shrinkage(fits)

  summary$sime_to_f1_msd <- NA_real_
  summary$sime_to_f2_msd <- NA_real_
  summary$f2_to_f1_msd <- NA_real_
  summary$shrinkage <- NA_real_
  if (nrow(saved_shrinkage) > 0L) {
    idx <- match(summary$selected_fit_ref, saved_shrinkage$fit_ref)
    keep <- !is.na(idx)
    summary$sime_to_f1_msd[keep] <- saved_shrinkage$sime_to_f1_msd[idx[keep]]
    summary$sime_to_f2_msd[keep] <- saved_shrinkage$sime_to_f2_msd[idx[keep]]
    summary$f2_to_f1_msd[keep] <- saved_shrinkage$f2_to_f1_msd[idx[keep]]
    summary$shrinkage[keep] <- saved_shrinkage$shrinkage[idx[keep]]
  }

  list(
    rule = "largest eta with SIME MSD <= c * baseline MSD; multiple c values share one increasing eta path",
    seed = seed,
    output_mode = output_mode,
    c = c_values,
    summary = summary,
    eta_path = eta_path,
    fits = fits,
    shrinkage_grid = shrinkage_grid
  )
}


#' Extract a selected SIME fit from SIME_select or SIME_select_path output
#'
#' @param result A \code{SIME_select()} or \code{SIME_select_path()} result.
#' @param c Optional c value. Required when \code{result} contains multiple
#'   selected c values.
#'
#' @return The selected SIME fit object, or \code{NULL} if that c has no
#'   selected fit.
#' @export
select_SIME_fit <- function(result, c = NULL) {
  selected_fit <- result[["selected_fit", exact = TRUE]]
  if (!is.null(selected_fit)) {
    return(selected_fit)
  }

  summary <- result[["summary", exact = TRUE]]
  fits <- result[["fits", exact = TRUE]]
  if (is.null(summary) || is.null(fits)) {
    stop("result must be a SIME_select() or SIME_select_path() object.")
  }

  if (is.null(c)) {
    ok <- which(!is.na(summary$selected_fit_ref))
    if (length(ok) != 1L) {
      stop("c must be supplied when the result contains multiple selected fits.")
    }
    row_idx <- ok
  } else {
    row_idx <- which(summary$c == as.numeric(c))
    if (length(row_idx) != 1L) {
      stop("Requested c was not found in result$summary$c.")
    }
  }

  fit_ref <- summary$selected_fit_ref[row_idx]
  if (is.na(fit_ref) || !nzchar(fit_ref)) {
    return(NULL)
  }
  fits[[fit_ref]]
}
