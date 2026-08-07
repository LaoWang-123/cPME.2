# SIME code

.make_sime_embedding <- function(coefs, params, Ua, d) {
  coefs <- as.matrix(coefs)
  params_all <- rbind(as.matrix(params), as.matrix(Ua))
  I_all <- nrow(params_all)
  d <- as.integer(d)

  env <- new.env(parent = parent.env(environment()))
  env$coefs <- coefs
  env$params_all <- params_all
  env$I_all <- I_all
  env$d <- d

  f <- function(parameters) {
    as.vector(
      (t(coefs[1:I_all, , drop = FALSE]) %*% etaFunc(parameters, params_all, 4 - d)) +
        (t(coefs[(I_all + 1):(I_all + d + 1), , drop = FALSE]) %*%
           matrix(c(1, parameters), ncol = 1))
    )
  }
  environment(f) <- env
  f
}

.calc_sime_correspondence_msd <- function(f, params, points) {
  params <- as.matrix(params)
  points <- as.matrix(points)
  pred <- t(vapply(seq_len(nrow(params)), function(i) {
    as.numeric(f(as.numeric(params[i, , drop = TRUE])))
  }, numeric(ncol(points))))
  mean(rowSums((pred - points)^2))
}

.nearest_parameter_init <- function(X, centers, center_params) {
  X <- as.matrix(X)
  centers <- as.matrix(centers)
  center_params <- as.matrix(center_params)

  nearest_idx <- vapply(seq_len(nrow(X)), function(i) {
    x_i <- matrix(
      X[i, ],
      nrow = nrow(centers),
      ncol = ncol(centers),
      byrow = TRUE
    )
    which.min(rowSums((centers - x_i)^2))
  }, integer(1))

  center_params[nearest_idx, , drop = FALSE]
}

#' SIME with anchor penalty + lambda screening
#' This code is run with given eta and chose lambda automatically
#'
#' @param f1_fun f1 function with vector input (u -> R^3)
#' @param init2 list: $centers (I x 3), $parameterization (I x 2),
#'   $theta_hat (optional), $km (required for calc_msd)
#' @param data2 raw data matrix (n x D) used in calc_msd (PME criterion)
#' @param eta anchor weight (penalty strength)
#' @param lambda tuning vector (can be length 1 or longer, like PME exp(-15:5))
#' @param lambda_selection Criterion used to select lambda. \code{"subject_only"}
#'   is the original SIME behavior and selects lambda by raw subject MSD only.
#'   \code{"anchor_weighted"} selects lambda by
#'   \code{(subject_MSD + eta * anchor_MSD) / (1 + eta)}, where
#'   \code{anchor_MSD} is computed on the fixed template anchors used in the
#'   SIME fit. \code{"template_all_weighted"} uses the same weighted criterion
#'   but replaces anchor MSD with direct-correspondence MSD over the raw
#'   subject points projected to \code{registered_fun}; their parameter values
#'   are then evaluated on \code{f1_fun} to get the corresponding template
#'   points.
#' @param registered_fun Registered subject surface used to build automatic
#'   all-template correspondences when
#'   \code{lambda_selection = "template_all_weighted"}.
#' @param template_points,template_params Optional advanced override for
#'   direct-correspondence template data. In the standard workflow these are
#'   built automatically from \code{registered_fun}.
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
    lambda_selection = c("subject_only", "anchor_weighted", "template_all_weighted"),
    registered_fun = NULL,
    template_points = NULL,
    template_params = NULL,
    d = 2,
    epsilon = 1,
    max_iter = 100,
    SSD_ratio_threshold = 5,
    verbose = FALSE
) {

  lambda_selection <- match.arg(lambda_selection)

  if (!exists("calc_msd", mode = "function"))
    stop("[SIMEpme] calc_msd() not found. Load PME utilities first.")
  if (is.null(init2$km))
    stop("[SIMEpme] init2$km is required for calc_msd().")
  if (lambda_selection == "template_all_weighted") {
    if ((is.null(template_points) || is.null(template_params)) && is.null(registered_fun)) {
      stop("registered_fun is required when lambda_selection = 'template_all_weighted' unless template_points and template_params are supplied.")
    }
    if (!is.null(template_points) || !is.null(template_params)) {
      if (is.null(template_points) || is.null(template_params)) {
        stop("template_points and template_params must be supplied together.")
      }
      template_points <- as.matrix(template_points)
      template_params <- as.matrix(template_params)
      if (nrow(template_points) != nrow(template_params)) {
        stop("template_points and template_params must have the same number of rows.")
      }
    }
  }

  data2 <- as.matrix(data2)
  D <- ncol(data2)

  # ----------------------------
  # Extract init2
  # ----------------------------
  X2 <- as.matrix(init2$centers)                # I x 3
  U2_init <- as.matrix(init2$parameterization)  # I x 2
  I <- nrow(X2)

  if (lambda_selection == "template_all_weighted" && is.null(template_points)) {
    template_params <- calc_params(
      f = registered_fun,
      X = data2,
      init_params = .nearest_parameter_init(data2, X2, U2_init),
      f_input = "vector"
    )
    template_points <- t(vapply(seq_len(nrow(template_params)), function(i) {
      as.numeric(f1_fun(as.numeric(template_params[i, , drop = TRUE])))
    }, numeric(D)))
  }

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
  subject_mse <- vector()
  anchor_mse <- rep(NA_real_, length(lambda))
  template_mse <- rep(NA_real_, length(lambda))
  selection_mse <- vector()
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

    # PME criterion: MSD on raw subject data.
    subject_mse[tuning_idx] <- calc_msd(data2, init2$km, f_embedding, params, D, d)
    anchor_mse[tuning_idx] <- .calc_sime_correspondence_msd(f_embedding, Ua, Xa)
    if (lambda_selection == "template_all_weighted") {
      template_mse[tuning_idx] <- .calc_sime_correspondence_msd(
        f = f_embedding,
        params = template_params,
        points = template_points
      )
    }

    selection_mse[tuning_idx] <- switch(
      lambda_selection,
      subject_only = subject_mse[tuning_idx],
      anchor_weighted = (subject_mse[tuning_idx] + eta * anchor_mse[tuning_idx]) / (1 + eta),
      template_all_weighted = (subject_mse[tuning_idx] + eta * template_mse[tuning_idx]) / (1 + eta)
    )

    if (verbose == TRUE) {
      message(sprintf("When lambda = %s, subject MSD = %s, selection MSD = %s.",
                      as.character(lambda[tuning_idx]),
                      as.character(subject_mse[tuning_idx]),
                      as.character(selection_mse[tuning_idx])))
    }

    # store path objects like PME
    coefs[[tuning_idx]] <- spline_coefs
    parameterization[[tuning_idx]] <- params

    embeddings[[tuning_idx]] <- .make_sime_embedding(
      coefs = coefs[[tuning_idx]],
      params = parameterization[[tuning_idx]],
      Ua = Ua,
      d = d
    )


    # PME early stop: last 4 MSD nondecreasing => break
    if (tuning_idx >= 4) {
      if (!is.unsorted(selection_mse[(tuning_idx - 3):tuning_idx])) {
        break
      }
    }
  }

  optimal_idx <- min(which(selection_mse == min(selection_mse, na.rm = TRUE)))

  coefs_opt <- coefs[[optimal_idx]]
  params_opt <- parameterization[[optimal_idx]]

  embedding_opt <- .make_sime_embedding(
    coefs = coefs_opt,
    params = params_opt,
    Ua = Ua,
    d = d
  )

  # return pme-like + sime-specific extras
  list(
    embedding_map = embedding_opt,
    params_opt = params_opt,
    centers = X2,
    theta_hat = theta,
    anchors = list(Ua = Ua, Xa = Xa, eta = eta),
    knots = init2$km,
    tuning = lambda[optimal_idx],
    MSD = subject_mse,
    subject_MSD = subject_mse,
    anchor_MSD = anchor_mse,
    template_MSD = template_mse,
    selection_MSD = selection_mse,
    lambda_selection = lambda_selection,
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
    lambda_selection = c("subject_only", "anchor_weighted", "template_all_weighted"),
    template_points = NULL,
    template_params = NULL,
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

  lambda_selection <- match.arg(lambda_selection)

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

  if (lambda_selection == "template_all_weighted" && is.null(template_points)) {
    template_params <- calc_params(
      f = f2_fun,
      X = data2,
      init_params = .nearest_parameter_init(data2, init2$centers, init2$parameterization),
      f_input = "vector"
    )
    template_points <- t(vapply(seq_len(nrow(template_params)), function(i) {
      as.numeric(f1_fun(as.numeric(template_params[i, , drop = TRUE])))
    }, numeric(D)))
  }

  threshold <- c * baseline_msd

  if (verbose) {
    message(sprintf("Baseline MSD = %g", baseline_msd))
    message(sprintf("Threshold = %g * %g = %g", c, baseline_msd, threshold))
  }

  # ----------------------------
  # storage
  # ----------------------------
  eta_fit_msd <- rep(NA_real_, length(eta_vec))
  eta_fit_anchor_msd <- rep(NA_real_, length(eta_vec))
  eta_fit_template_msd <- rep(NA_real_, length(eta_vec))
  eta_fit_selection_msd <- rep(NA_real_, length(eta_vec))
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
      lambda_selection = lambda_selection,
      registered_fun = f2_fun,
      template_points = template_points,
      template_params = template_params,
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
    eta_fit_anchor_msd[e] <- fit_e$anchor_MSD[idx_star]
    eta_fit_template_msd[e] <- fit_e$template_MSD[idx_star]
    eta_fit_selection_msd[e] <- fit_e$selection_MSD[idx_star]
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

  selected_eta <- if (is.na(selected_idx)) NA_real_ else eta_vec[selected_idx]
  selected_fit <- if (is.na(selected_idx)) NULL else eta_fit_list[[selected_idx]]

  if (verbose) {
    message("==================================================")
    if (is.na(selected_idx)) {
      message("No eta satisfied SIME MSD <= c * baseline MSD.")
    } else {
      message(sprintf("Selected eta = %g", selected_eta))
      message(sprintf("Selected lambda = %g", eta_lambda_star[selected_idx]))
      message(sprintf("Selected MSD = %g", eta_fit_msd[selected_idx]))
    }
    message("==================================================")
  }

  list(
    rule = "largest eta with SIME MSD <= c * baseline MSD; early stop at first exceedance",
    lambda_selection = lambda_selection,
    baseline_msd = baseline_msd,
    threshold = threshold,
    c = c,
    eta_fit_msd = eta_fit_msd,
    eta_fit_anchor_msd = eta_fit_anchor_msd,
    eta_fit_template_msd = eta_fit_template_msd,
    eta_fit_selection_msd = eta_fit_selection_msd,
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
#' @param lambda_selection Criterion used by \code{SIME()} to select lambda.
#' @param template_points,template_params Optional direct-correspondence
#'   template data passed to \code{SIME()} when
#'   \code{lambda_selection = "template_all_weighted"}.
#' @param d Intrinsic dimension.
#' @param epsilon,max_iter,SSD_ratio_threshold Same as in \code{SIME()}.
#' @param init_args_f2 Arguments passed to \code{initialize_pme()}.
#' @param seed Optional random seed.
#' @param output_mode Output size mode. \code{"compact"} keeps only selected
#'   unique eta fits; \code{"full"} keeps every evaluated eta fit.
#' @param shrinkage_grid Parameter grid used for shrinkage. If \code{NULL},
#'   a 40 by 40 disk grid from \code{make_uv_grid()} is used.
#' @param shrinkage_method Method used to evaluate shrinkage.
#'   \code{"data_projection"} is the default and
#'   projects all rows of \code{data2} onto \code{f2_fun} and compares
#'   manifolds at those projected parameter locations. \code{"grid"} compares
#'   manifolds on \code{shrinkage_grid}.
#' @param shrinkage_grid_name Optional label for a user-supplied
#'   \code{shrinkage_grid}. If omitted, grid shrinkage is labeled
#'   \code{"default_disk_40x40"} for the default grid or \code{"user_grid"} for
#'   a supplied grid.
#' @param verbose Logical.
#'
#' @return A list with \code{summary}, \code{eta_path}, and \code{fits}. In
#'   compact mode, \code{fits} stores only selected unique eta fits. In full
#'   mode, \code{fits} stores every evaluated eta fit. The top-level
#'   \code{shrinkage_method} records whether shrinkage used
#'   \code{"data_projection"} or a named grid.
#' @export
SIME_select_path <- function(
    f1_fun,
    f2_fun,
    data2,
    c = 1.3,
    eta_vec = c(exp(-15:-2), 0.1, exp(-1:5)),
    lambda = exp(-15:5),
    lambda_selection = c("subject_only", "anchor_weighted", "template_all_weighted"),
    template_points = NULL,
    template_params = NULL,
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
    shrinkage_method = c("data_projection", "grid"),
    shrinkage_grid_name = NULL,
    verbose = TRUE
) {

  output_mode <- match.arg(output_mode)
  shrinkage_method <- match.arg(shrinkage_method)
  lambda_selection <- match.arg(lambda_selection)

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

  shrinkage_grid_supplied <- !is.null(shrinkage_grid)
  if (shrinkage_method == "grid" && !shrinkage_grid_supplied) {
    shrinkage_grid <- make_uv_grid(n_u = 40, n_v = 40, grid_type = "disk")
  }
  shrinkage_label <- if (shrinkage_method == "data_projection") {
    "data_projection"
  } else if (!is.null(shrinkage_grid_name) &&
             length(shrinkage_grid_name) == 1L &&
             nzchar(shrinkage_grid_name)) {
    as.character(shrinkage_grid_name)
  } else if (!shrinkage_grid_supplied) {
    "default_disk_40x40"
  } else {
    "user_grid"
  }

  eval_on_params <- function(f, params) {
    t(vapply(seq_len(nrow(params)), function(i) {
      as.numeric(f(as.numeric(params[i, , drop = TRUE])))
    }, numeric(D)))
  }

  nearest_parameter_init <- function(X, centers, center_params) {
    X <- as.matrix(X)
    centers <- as.matrix(centers)
    center_params <- as.matrix(center_params)

    nearest_idx <- vapply(seq_len(nrow(X)), function(i) {
      x_i <- matrix(
        X[i, ],
        nrow = nrow(centers),
        ncol = ncol(centers),
        byrow = TRUE
      )
      which.min(rowSums((centers - x_i)^2))
    }, integer(1))

    center_params[nearest_idx, , drop = FALSE]
  }

  resolve_shrinkage_params <- function(init2) {
    if (shrinkage_method == "grid") {
      return(as.matrix(shrinkage_grid))
    }

    init_guess <- nearest_parameter_init(
      X = data2,
      centers = init2$centers,
      center_params = init2$parameterization
    )

    calc_params(
      f = f2_fun,
      X = data2,
      init_params = init_guess,
      f_input = "vector"
    )
  }

  compute_saved_shrinkage <- function(fits, shrinkage_params) {
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

    f1_data <- eval_on_params(f1_fun, shrinkage_params)
    f2_data <- eval_on_params(f2_fun, shrinkage_params)

    if (!all(dim(f1_data) == dim(f2_data))) {
      stop("Reference and baseline manifolds do not have matching output dimensions.")
    }

    f2_to_f1_msd <- mean(rowSums((f2_data - f1_data)^2))

    rows <- lapply(names(fits), function(fit_ref) {
      fit <- fits[[fit_ref]]
      sime_data <- eval_on_params(fit$embedding_map, shrinkage_params)
      if (!all(dim(sime_data) == dim(f1_data))) {
        stop("SIME and reference manifolds do not have matching output dimensions.")
      }

      sime_to_f1_msd <- mean(rowSums((sime_data - f1_data)^2))
      sime_to_f2_msd <- mean(rowSums((sime_data - f2_data)^2))
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

  if (lambda_selection == "template_all_weighted" && is.null(template_points)) {
    template_params <- calc_params(
      f = f2_fun,
      X = data2,
      init_params = nearest_parameter_init(
        X = data2,
        centers = init2$centers,
        center_params = init2$parameterization
      ),
      f_input = "vector"
    )
    template_points <- eval_on_params(f1_fun, template_params)
  }

  thresholds <- c_values * baseline_msd

  eta_fit_msd <- rep(NA_real_, length(eta_vec))
  eta_fit_anchor_msd <- rep(NA_real_, length(eta_vec))
  eta_fit_template_msd <- rep(NA_real_, length(eta_vec))
  eta_fit_selection_msd <- rep(NA_real_, length(eta_vec))
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
            lambda_selection = lambda_selection,
            registered_fun = f2_fun,
            template_points = template_points,
            template_params = template_params,
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
        eta_fit_anchor_msd[e] <- fit_e$anchor_MSD[idx_star]
        eta_fit_template_msd[e] <- fit_e$template_MSD[idx_star]
        eta_fit_selection_msd[e] <- fit_e$selection_MSD[idx_star]
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
  selected_anchor_msd <- ifelse(is.na(selected_idx), NA_real_, eta_fit_anchor_msd[selected_idx])
  selected_template_msd <- ifelse(is.na(selected_idx), NA_real_, eta_fit_template_msd[selected_idx])
  selected_selection_msd <- ifelse(is.na(selected_idx), NA_real_, eta_fit_selection_msd[selected_idx])

  fit_ref <- ifelse(is.na(selected_idx), NA_character_, names(eta_fit_list)[selected_idx])

  summary <- data.frame(
    c = c_values,
    baseline_msd = baseline_msd,
    threshold = thresholds,
    selected_eta = selected_eta,
    selected_lambda = selected_lambda,
    selected_msd = selected_msd,
    selected_anchor_msd = selected_anchor_msd,
    selected_template_msd = selected_template_msd,
    selected_selection_msd = selected_selection_msd,
    selected_fit_ref = as.character(fit_ref),
    error_message = c_error_message,
    stringsAsFactors = FALSE
  )

  evaluated <- !vapply(eta_fit_list, is.null, logical(1))
  eta_path <- data.frame(
    eta_index = seq_along(eta_vec),
    eta = eta_vec,
    fit_msd = eta_fit_msd,
    anchor_msd = eta_fit_anchor_msd,
    template_msd = eta_fit_template_msd,
    selection_msd = eta_fit_selection_msd,
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
  shrinkage_params <- resolve_shrinkage_params(init2)
  saved_shrinkage <- compute_saved_shrinkage(fits, shrinkage_params)

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
    lambda_selection = lambda_selection,
    seed = seed,
    output_mode = output_mode,
    shrinkage_method = shrinkage_label,
    c = c_values,
    summary = summary,
    eta_path = eta_path,
    fits = fits,
    shrinkage_grid = if (shrinkage_method == "grid") shrinkage_grid else NULL
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
