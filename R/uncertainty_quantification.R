

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
