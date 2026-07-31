### Evaluation Utils
# This file include functions used to calculate MSD across manifolds
# theoretical variance estimator



#' Compute manifold approximation error (MSD) with respect to a true surface
#'
#' This function computes the mean squared distance (MSD) between a true
#' manifold \eqn{f_{\text{true}}} and an estimated manifold \eqn{\hat f}.
#' The metric is defined as
#' \deqn{
#' \mathrm{MSD}(\hat f) = \mathbb{E} \left\| X - \hat f(\pi_{\hat f}(X)) \right\|^2,
#' }
#' @param f_hat A function representing the estimated manifold \eqn{\hat f},
#'   mapping from parameter space \eqn{U} to \eqn{\mathbb{R}^3}.
#' @param f_true A function representing the true manifold \eqn{f_{\text{true}}},
#'   mapping from parameter space \eqn{U} to \eqn{\mathbb{R}^3}.
#' @param uv_grid A matrix or data frame of parameter values (u, v) used to
#'   generate samples from the true manifold.
#' @param seed Integer seed for reproducibility when generating data.
#' @param init_method Character string specifying the initialization method
#'   for projection (passed to \code{pme_initial_guess}). Default is \code{"pca"}.
#' @param return_xyz default FALSE
#'
#' @return A numeric scalar representing the mean squared distance (MSD)
#'
#' @export
calc_manifold_msd <- function(f_hat,
                              f_true,
                              uv_grid,
                              seed = 123,
                              init_method = "pca",
                              return_xyz = FALSE) {

  # 1. generate true surface (no noise)
  true_data <- generate_surface_data(
    f_true,
    noise_sd = 0,
    seed = seed,
    uv_grid = uv_grid
  )

  # 2. project true points onto estimated manifold
  proj <- calc_params(
    f = f_hat,
    X = true_data$XYZ,
    init_params = pme_initial_guess(true_data$XYZ, 2, init_method),
    f_input = "vector"
  )

  # 3. reconstruct points from projection
  fit <- generate_surface_data(
    f_hat,
    noise_sd = 0,
    seed = seed,
    uv_grid = proj
  )

  # 4. compute MSD
  msd <- mean(rowSums((fit$XYZ - true_data$XYZ)^2))

  if (return_xyz) {
    return(list(
      msd = msd,
      fit_xyz = fit$XYZ,
      true_xyz = true_data$XYZ
    ))
  }

  return(msd)
}


#' Compute Mean Squared Distance Between Two Manifolds Under Parameter Correspondence
#'
#' This function computes the mean squared distance (MSD) between two manifold
#' functions evaluated on a shared parameter grid \code{uv_grid}. It assumes that
#' both manifolds are defined on the same parameter domain and that pointwise
#' correspondence is given by identical \code{(u, v)} locations.
#'
#' @param f_ref A function representing the reference manifold. It should accept
#'   input in \code{uv_grid} format and return an \eqn{n \times 3} matrix of 3D coordinates.
#' @param f_cmp A function representing the comparison manifold (e.g., SIME or PME estimate).
#'   Same interface as \code{f_ref}.
#' @param uv_grid A data frame or matrix of parameter locations with columns \code{u} and \code{v}.
#' @param seed Integer random seed for reproducibility (passed to data generation).
#'
#' @return A numeric value representing the mean squared distance between the two manifolds.
#'
#' @export
calc_correspondence_msd <- function(f_ref,
                                    f_cmp,
                                    uv_grid,
                                    seed = 123) {
  ref_data <- generate_surface_data(
    f = f_ref,
    noise_sd = 0,
    seed = seed,
    uv_grid = uv_grid
  )

  cmp_data <- generate_surface_data(
    f = f_cmp,
    noise_sd = 0,
    seed = seed,
    uv_grid = uv_grid
  )

  if (!all(dim(ref_data$XYZ) == dim(cmp_data$XYZ))) {
    stop("Reference and comparison manifolds do not have matching output dimensions.")
  }

  msd <- mean(rowSums((ref_data$XYZ - cmp_data$XYZ)^2))
  return(msd)
}

########################################################################
#######################################################################
# theoretical variance estimator



#' Compute Truncated Prior Kernel Matrix
#'
#' Computes the spectral prior covariance kernel
#' \deqn{
#' K_0(s,t) = \sum_m
#' \frac{\psi_m(s)\psi_m(t)}{\lambda \alpha_m^2 + \eta}.
#' }
#'
#' @param basis_grid_1 Scalar basis grid evaluated at the first point set.
#' @param basis_grid_2 Scalar basis grid evaluated at the second point set.
#' @param lambda Smoothness tuning parameter.
#' @param eta Structural regularization tuning parameter.
#'
#' @return A kernel matrix with rows indexed by the first point set and
#'   columns indexed by the second point set.
#' @export
prior_kernel_matrix <- function(
    basis_grid_1,
    basis_grid_2,
    lambda,
    eta
) {
  keys1 <- sort(names(basis_grid_1))
  keys2 <- sort(names(basis_grid_2))

  if (!identical(keys1, keys2)) {
    stop("basis grids must use identical basis functions.")
  }

  Phi1 <- do.call(cbind, lapply(keys1, function(key) basis_grid_1[[key]]$psi))
  Phi2 <- do.call(cbind, lapply(keys2, function(key) basis_grid_2[[key]]$psi))

  alpha <- vapply(
    keys1,
    function(key) basis_grid_1[[key]]$lambda,
    numeric(1)
  )

  weights <- 1 / (lambda * alpha^2 + eta)

  sweep(Phi1, 2, weights, `*`) %*% t(Phi2)
}


#' Estimate Fixed-Projection Pointwise Posterior Variance
#'
#' Computes the fixed-projection posterior pointwise variance
#' \deqn{
#' Var\{f_\ell(u)\}
#' =
#' K_0(u,u)
#' -
#' k(u)^T (I + K_U)^{-1} k(u),
#' }
#' where \eqn{K_U = K_0(U_{\mathrm{obs}}, U_{\mathrm{obs}})}
#' and \eqn{k(u) = K_0(U_{\mathrm{obs}}, u)}.
#'
#' @param target_u A two-column matrix or data frame of target parameter
#'   locations where posterior variance is evaluated.
#' @param obs_u A two-column matrix or data frame of fixed projection
#'   locations \eqn{u_i = \pi_{\hat f_{\mathrm{SIME}}}(x_i)}. Required if
#'   \code{basis_obs} is not supplied.
#' @param basis_set A basis set returned by \code{build_basis_set()}.
#'   For UQ, this should usually be constructed with
#'   \code{include_constant = TRUE}.
#' @param lambda Selected smoothness tuning parameter.
#' @param eta Selected structural regularization parameter.
#' @param basis_obs Optional precomputed scalar basis grid for \code{obs_u},
#'   usually from \code{build_basis_grid(basis_set, obs_u, mode = "scalar")}.
#' @param KU Optional precomputed prior kernel matrix
#'   \eqn{K_0(U_{\mathrm{obs}}, U_{\mathrm{obs}})}. If supplied, it must be
#'   compatible with \code{basis_obs}.
#' @param total_3d Logical. If \code{FALSE}, returns the scalar coordinate-wise
#'   variance. If \code{TRUE}, returns the summed three-dimensional variance,
#'   equal to three times the scalar variance under independent identical
#'   coordinate priors.
#' @param jitter Small diagonal stabilization added to \eqn{I + K_U}.
#'
#' @return A numeric vector of posterior variances, one for each row of
#'   \code{target_u}.
#' @export
pointwise_variance_estimator <- function(
    target_u,
    obs_u = NULL,
    basis_set,
    lambda,
    eta,
    basis_obs = NULL,
    KU = NULL,
    total_3d = FALSE,
    jitter = 1e-8
) {
  target_u <- as.matrix(target_u)

  if (ncol(target_u) != 2) {
    stop("target_u must have two columns.")
  }
  if (lambda < 0 || eta <= 0) {
    stop("lambda must be nonnegative and eta must be positive.")
  }

  if (is.null(basis_obs)) {
    if (is.null(obs_u)) {
      stop("obs_u is required when basis_obs is not supplied.")
    }

    obs_u <- as.matrix(obs_u)

    if (ncol(obs_u) != 2) {
      stop("obs_u must have two columns.")
    }

    basis_obs <- build_basis_grid(basis_set, obs_u, mode = "scalar")
  }

  basis_target <- build_basis_grid(basis_set, target_u, mode = "scalar")

  if (is.null(KU)) {
    KU <- prior_kernel_matrix(basis_obs, basis_obs, lambda, eta)
  } else {
    KU <- as.matrix(KU)
  }

  n_obs <- length(basis_obs[[1]]$psi)

  if (!all(dim(KU) == c(n_obs, n_obs))) {
    stop("KU must be a square matrix with dimension equal to the number of observed parameter locations.")
  }

  KTU <- prior_kernel_matrix(basis_target, basis_obs, lambda, eta)
  KTT <- prior_kernel_matrix(basis_target, basis_target, lambda, eta)

  KTT_diag <- diag(KTT)
  M <- diag(nrow(KU)) + KU + jitter * diag(nrow(KU))

  solved <- t(solve(M, t(KTU)))
  var_scalar <- KTT_diag - rowSums(KTU * solved)
  var_scalar <- pmax(var_scalar, 0)

  if (total_3d) 3 * var_scalar else var_scalar
}

########################################################################
# Surface volume evaluation

# Volume estimation methods for fitted PME/SIME surfaces.
#
# These functions are the package version of the real-data workflow volume
# core. They intentionally exclude ADNI-specific surface extraction, plotting,
# simulation helpers, and the original LPME reference code.

#' Extract Surface Components for SIME/PME Volume Evaluation
#'
#' Extract the embedding map, coefficient matrix, and parameter locations needed
#' by volume and plotting utilities. The input can be either a PME object or a
#' `SIME_select()` result; for `SIME_select()`, the selected fit is used.
#'
#' @param fit A cPME.2 PME object or a `SIME_select()` result.
#'
#' @return A list containing `embedding_map`, `coefs`, `params`, the selected
#'   tuning index, and detected `surface_type` (`"PME"` or `"SIME"`).
#' @export
extract_sime_volume_components <- function(fit) {
  if (!is.null(fit$selected_fit)) {
    fit <- fit$selected_fit
  }

  if (is.null(fit$embedding_map)) {
    stop("fit must contain embedding_map, or be a SIME_select object.")
  }
  if (is.null(fit$coefs) || is.null(fit$parameterization)) {
    stop("fit must contain coefs and parameterization.")
  }

  idx <- NA_integer_
  if (!is.null(fit$tuning) && !is.null(fit$tuning_vec)) {
    idx <- match(fit$tuning, fit$tuning_vec)
  }
  if (is.na(idx) && !is.null(fit$MSD)) {
    idx <- which.min(fit$MSD)
  }
  if (is.na(idx)) {
    idx <- length(fit$coefs)
  }

  params <- as.matrix(fit$parameterization[[idx]])
  coefs <- as.matrix(fit$coefs[[idx]])
  surface_type <- "PME"

  if (!is.null(fit$anchors) && !is.null(fit$anchors$Ua)) {
    params <- rbind(params, fit$anchors$Ua)
    surface_type <- "SIME"
  }

  list(
    embedding_map = fit$embedding_map,
    coefs = coefs,
    params = params,
    selected_index = idx,
    surface_type = surface_type
  )
}

signed_distance_to_partition_plane <- function(points, partition_plane) {
  if (is.null(partition_plane)) {
    stop("partition_plane must not be NULL.")
  }
  if (is.null(partition_plane$center) || is.null(partition_plane$normal)) {
    stop("partition_plane must contain center and normal.")
  }

  points <- as.matrix(points)[, 1:3, drop = FALSE]
  center <- as.numeric(partition_plane$center)[1:3]
  normal <- as.numeric(partition_plane$normal)[1:3]
  normal <- normal / sqrt(sum(normal^2))

  as.numeric(sweep(points, 2, center, "-") %*% normal)
}

surface_point_orientation <- function(embedding_map,
                                      coefs,
                                      centers,
                                      params,
                                      x) {
  if (!requireNamespace("pracma", quietly = TRUE)) {
    stop("Package 'pracma' is required for surface volume estimation.")
  }

  d <- ncol(params)
  n <- nrow(params)
  x <- as.numeric(x)[1:3]

  nearest_center <- which.min(vapply(seq_len(nrow(centers)), function(i) {
    dist_euclidean(x, as.vector(centers[i, ]))
  }, numeric(1)))

  x_params <- projection_pme(x, embedding_map, params[nearest_center, ])
  if (is.null(x_params) || any(!is.finite(x_params))) {
    return(NA_real_)
  }
  x_proj <- embedding_map(x_params)

  param_diff <- -1 * sweep(params, 2, x_params)
  param_diff_norms <- sqrt(rowSums(param_diff^2))

  grad_factor <- rep(0, length(param_diff_norms))
  keep <- param_diff_norms > .Machine$double.eps
  grad_factor[keep] <- 2 * log(param_diff_norms[keep]) + 1

  jacobian <- t(coefs[seq_len(n), , drop = FALSE]) %*%
    (param_diff * grad_factor) +
    t(coefs[(n + 1):(n + d + 1), , drop = FALSE][-1, , drop = FALSE])

  normal_vector <- pracma::cross(
    as.vector(jacobian[1:3, 1]),
    as.vector(jacobian[1:3, 2])
  )

  as.numeric(sign(crossprod(as.numeric(x_proj)[1:3] - x, normal_vector)))
}

surface_interior_identification <- function(fit,
                                            x,
                                            ref = c(0, 0, 0)) {
  component <- extract_sime_volume_components(fit)
  x <- as.matrix(x)
  if (nrow(x) == 0) {
    return(logical(0))
  }

  param_centers <- t(vapply(seq_len(nrow(component$params)), function(i) {
    as.numeric(component$embedding_map(component$params[i, ]))[1:3]
  }, numeric(3)))

  orientation_ref <- surface_point_orientation(
    embedding_map = component$embedding_map,
    coefs = component$coefs,
    centers = param_centers,
    params = component$params,
    x = ref
  )[1]

  orientation_x <- vapply(seq_len(nrow(x)), function(i) {
    surface_point_orientation(
      embedding_map = component$embedding_map,
      coefs = component$coefs,
      centers = param_centers,
      params = component$params,
      x = x[i, ]
    )
  }, numeric(1))

  !is.na(orientation_x) & !is.na(orientation_ref) & orientation_x == orientation_ref
}

sample_volume_candidates <- function(data,
                                     n_points = 10000,
                                     limit_scaler = 0.05,
                                     data_max = NULL) {
  data <- as.matrix(data)[, 1:3, drop = FALSE]

  if (is.null(data_max)) {
    mins <- apply(data, 2, min, na.rm = TRUE)
    maxs <- apply(data, 2, max, na.rm = TRUE)
  } else {
    scales <- as.numeric(data_max)[1:3]
    mins <- -scales
    maxs <- scales
  }

  center <- (mins + maxs) / 2
  half_range <- (maxs - mins) / 2
  half_range <- half_range * (1 + limit_scaler)

  bounds <- cbind(
    min = center - half_range,
    max = center + half_range
  )
  rownames(bounds) <- c("x", "y", "z")

  candidates <- cbind(
    x = runif(n_points, bounds["x", "min"], bounds["x", "max"]),
    y = runif(n_points, bounds["y", "min"], bounds["y", "max"]),
    z = runif(n_points, bounds["z", "min"], bounds["z", "max"])
  )

  list(
    candidates = candidates,
    full_volume = prod(bounds[, "max"] - bounds[, "min"]),
    bounds = bounds
  )
}

empty_volume_result <- function(message,
                                stage,
                                n_points,
                                partition_index,
                                partition_plane,
                                ref,
                                classification_method,
                                volume_multiplier,
                                warning_enabled = TRUE) {
  if (isTRUE(warning_enabled)) {
    warning(message, call. = FALSE)
  }

  list(
    volume = NA_real_,
    interior = NA,
    interior_candidates = matrix(numeric(0), ncol = 3),
    candidates = matrix(numeric(0), ncol = 3),
    full_volume = NA_real_,
    bounds = matrix(NA_real_, nrow = 3, ncol = 2,
                    dimnames = list(c("x", "y", "z"), c("min", "max"))),
    n_points = n_points,
    partition_index = partition_index,
    partition_plane = partition_plane,
    ref = ref,
    classification_method = classification_method,
    fit_types = c(part1 = NA_character_, part2 = NA_character_),
    volume_multiplier = volume_multiplier,
    volume_success = FALSE,
    volume_log = data.frame(
      level = "warning",
      stage = stage,
      message = message,
      stringsAsFactors = FALSE
    )
  )
}

safe_extract_volume_components <- function(fit, label) {
  tryCatch(
    extract_sime_volume_components(fit),
    error = function(err) {
      structure(
        list(message = sprintf("%s: %s", label, conditionMessage(err))),
        class = "volume_component_error"
      )
    }
  )
}

#' Estimate Volume From Two Fitted Surface Pieces
#'
#' Estimate the volume enclosed by a shape represented by two fitted PME/SIME
#' surface pieces. Candidate 3D points are sampled from a padded bounding box
#' around `data`, classified as interior by surface normal orientation, and
#' converted to volume by Monte Carlo integration.
#'
#' The inputs `fit_part1` and `fit_part2` can be cPME.2 PME objects or
#' `SIME_select()` results. For a `SIME_select()` result, the selected SIME fit
#' is used automatically.
#'
#' @param fit_part1 Fitted PME object or `SIME_select()` result for the first
#'   surface piece.
#' @param fit_part2 Fitted PME object or `SIME_select()` result for the second
#'   surface piece.
#' @param data Numeric matrix/data frame with at least three columns giving the
#'   point cloud used to define the candidate sampling box.
#' @param n_points Number of Monte Carlo candidate points.
#' @param limit_scaler Proportional padding added to each side of the data
#'   bounding box.
#' @param partition_index Coordinate index used to split candidates when
#'   `partition_plane` is `NULL`.
#' @param partition_plane Optional list with `center` and `normal` fields. If
#'   supplied, candidates are assigned to part 1 or part 2 by signed distance to
#'   this plane.
#' @param ref A known interior reference point in the same coordinate system as
#'   `data`.
#' @param data_max Optional numeric length-3 scale. If supplied, candidates are
#'   sampled from `[-data_max, data_max]` rather than the data bounding box.
#' @param volume_multiplier Multiplier applied to the estimated volume, for
#'   example voxel volume when coordinates are voxel indices.
#' @param classification_method Either `"partition_side"` or `"both_surfaces"`.
#'   `"partition_side"` classifies each candidate using only the surface piece
#'   on its partition side. `"both_surfaces"` requires the candidate to be
#'   inside both fitted pieces.
#' @param seed Optional random seed for candidate sampling.
#' @param fail_on_error Logical. If \code{FALSE}, invalid or missing fits return
#'   an NA volume with a \code{volume_log} entry instead of stopping.
#' @param warn_on_failure Logical. If \code{TRUE}, non-fatal NA-volume failures
#'   also emit a warning.
#'
#' @return A list containing the estimated `volume`, candidate points, interior
#'   indicators, sampling bounds, metadata, and detected fit types.
#'
#' @export
estimate_partitioned_surface_volume <- function(
    fit_part1,
    fit_part2,
    data,
    n_points = 10000,
    limit_scaler = 0.05,
    partition_index = 3,
    partition_plane = NULL,
    ref = c(0, 0, 0),
    data_max = NULL,
    volume_multiplier = 1,
    classification_method = c("partition_side", "both_surfaces"),
    seed = NULL,
    fail_on_error = FALSE,
    warn_on_failure = TRUE) {
  classification_method <- match.arg(classification_method)
  if (!is.null(seed)) {
    set.seed(seed)
  }
  if (!partition_index %in% 1:3) {
    stop("partition_index must be 1, 2, or 3.")
  }

  component_part1 <- safe_extract_volume_components(fit_part1, "fit_part1")
  component_part2 <- safe_extract_volume_components(fit_part2, "fit_part2")
  component_errors <- c()
  if (inherits(component_part1, "volume_component_error")) {
    component_errors <- c(component_errors, component_part1$message)
  }
  if (inherits(component_part2, "volume_component_error")) {
    component_errors <- c(component_errors, component_part2$message)
  }

  if (length(component_errors) > 0) {
    msg <- paste(
      "Surface volume not estimated because at least one fitted surface is missing or invalid.",
      paste(component_errors, collapse = " | ")
    )
    if (isTRUE(fail_on_error)) {
      stop(msg, call. = FALSE)
    }
    return(empty_volume_result(
      message = msg,
      stage = "fit_validation",
      n_points = n_points,
      partition_index = partition_index,
      partition_plane = partition_plane,
      ref = ref,
      classification_method = classification_method,
      volume_multiplier = volume_multiplier,
      warning_enabled = warn_on_failure
    ))
  }

  candidate_info <- sample_volume_candidates(
    data = data,
    n_points = n_points,
    limit_scaler = limit_scaler,
    data_max = data_max
  )

  candidates <- candidate_info$candidates
  if (is.null(partition_plane)) {
    part1_index <- candidates[, partition_index] > ref[partition_index]
  } else {
    part1_index <- signed_distance_to_partition_plane(candidates, partition_plane) >= 0
  }
  part2_index <- !part1_index

  interior <- rep(FALSE, nrow(candidates))
  classification_result <- tryCatch({
    if (classification_method == "partition_side") {
      interior[part1_index] <- surface_interior_identification(
        fit = fit_part1,
        x = candidates[part1_index, , drop = FALSE],
        ref = ref
      )
      interior[part2_index] <- surface_interior_identification(
        fit = fit_part2,
        x = candidates[part2_index, , drop = FALSE],
        ref = ref
      )
    } else {
      interior_part1 <- surface_interior_identification(
        fit = fit_part1,
        x = candidates,
        ref = ref
      )
      interior_part2 <- surface_interior_identification(
        fit = fit_part2,
        x = candidates,
        ref = ref
      )
      interior <- interior_part1 & interior_part2
    }
    interior
  }, error = function(err) err)

  if (inherits(classification_result, "error")) {
    msg <- paste(
      "Surface volume not estimated because candidate classification failed.",
      conditionMessage(classification_result)
    )
    if (isTRUE(fail_on_error)) {
      stop(msg, call. = FALSE)
    }
    out <- empty_volume_result(
      message = msg,
      stage = "classification",
      n_points = n_points,
      partition_index = partition_index,
      partition_plane = partition_plane,
      ref = ref,
      classification_method = classification_method,
      volume_multiplier = volume_multiplier,
      warning_enabled = warn_on_failure
    )
    out$candidates <- candidates
    out$bounds <- candidate_info$bounds
    out$full_volume <- candidate_info$full_volume
    out$interior <- rep(NA, nrow(candidates))
    out$fit_types <- c(
      part1 = component_part1$surface_type,
      part2 = component_part2$surface_type
    )
    return(out)
  }
  interior <- classification_result

  list(
    volume = mean(interior) * candidate_info$full_volume * volume_multiplier,
    interior = interior,
    interior_candidates = candidates[interior, , drop = FALSE],
    candidates = candidates,
    full_volume = candidate_info$full_volume,
    bounds = candidate_info$bounds,
    n_points = n_points,
    partition_index = partition_index,
    partition_plane = partition_plane,
    ref = ref,
    classification_method = classification_method,
    volume_multiplier = volume_multiplier,
    volume_success = TRUE,
    volume_log = data.frame(
      level = "info",
      stage = "complete",
      message = "Surface volume estimated successfully.",
      stringsAsFactors = FALSE
    ),
    fit_types = c(
      part1 = component_part1$surface_type,
      part2 = component_part2$surface_type
    )
  )
}

