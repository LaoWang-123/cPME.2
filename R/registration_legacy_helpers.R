# Legacy registration helper functions
#
# These functions support `Registration_legacy`. They are intentionally kept
# internal to the package; current workflows should use `Registration` and the
# grid-precomputed helpers in registration_helpers.R.
#### Utils.R
# This file contains all the public functions

# ==========================================================
# Vectorized computation of E = âˆ« ||f1(u,v) - f2(u,v)||^2 du dv
# ----------------------------------------------------------
# f1, f2 : functions (u,v) -> c(x,y,z)
# Ugrid  : data.frame with columns u,v  (regular grid)
# ==========================================================

# Ugrid <- expand.grid(
#   u = seq(0, 1, by = 0.01),
#   v = seq(0, 1, by = 0.01)
# )

compute_E_vectorized <- function(f1, f2, Ugrid) {

  du <- mean(diff(unique(Ugrid$u)))
  dv <- mean(diff(unique(Ugrid$v)))

  vals1 <- t(apply(Ugrid, 1, function(uv) f1(uv[1], uv[2])))
  vals2 <- t(apply(Ugrid, 1, function(uv) f2(uv[1], uv[2])))

  diff_sq <- rowSums((vals1 - vals2)^2)

  E_val <- sum(diff_sq) * du * dv
  return(E_val)
}

##################################################
# bi set function is stored in basis_functions.R
#######################################################
# Calculate dphi(b)
# ---------------------------
# Reference: Section 3.5, Eq. (11)-(12) in cPME paper
# For each basis field b_i (either âˆ‡Ïˆ_i or *âˆ‡Ïˆ_i), compute
#   dphi^k(b_i) = 1 / ||âˆ‡Ïˆ_i|| * (-0.5 * Î»_i * Ïˆ_i * f2^k + âˆ‡f2^k Â· âˆ‡Ïˆ_i)      # non-rotated
#   dphi^k(b_i) = 1 / ||âˆ‡Ïˆ_i|| * (âˆ‡f2^k Â· *âˆ‡Ïˆ_i)                              # rotated

# We get âˆ‡f2^k(u,v) from assemble_delta_f2k_from_state function defined previously

#
#' Factory that returns a callable dphi function
#'
#' @param basis bi
#' @param grad_f2k_fun grad_f2k_fun
#' @param f2_fun f2_fun
#' @param rotated True or False
#'
#' @returns a function (u,v) to 3x1 vector
#'
#' @noRd
make_dphi_fn <- function(basis, grad_f2k_fun, f2_fun, rotated = FALSE) {

  norm_pq <- basis$norm_pq  # ||âˆ‡Ïˆ_pq||_L2
  lambda_pq <- basis$lambda_pq
  if (!rotated) {
    # Gradient-type basis (curl-free)
    # Eq. (11): 1/norm * (-0.5 * Î» * Ïˆ * f2^k + âˆ‡f2^k Â· âˆ‡Ïˆ)
    psi_fun  <- basis$psi
    grad_fun <- basis$grad_psi

    # Return a function of (u,v)
    function(u, v) {
      (1 / norm_pq) * (
        -0.5 * lambda_pq * psi_fun(u, v) * as.numeric(f2_fun(u, v)) +
          as.numeric(grad_f2k_fun(u, v) %*% grad_fun(u, v))
      )
    }

  } else {
    # Rotated (divergence-free) basis
    # Eq. (12): 1/norm * (âˆ‡f2^k Â· *âˆ‡Ïˆ)
    rot_fun <- basis$rot_grad_psi

    function(u, v) {
      (1 / norm_pq) * (grad_f2k_fun(u,v) %*% rot_fun(u, v))
    }
  }
}

# ---------------------------
# (X) Build dphi^k(b_i) set
# ---------------------------
# Reference: Eq.(11)-(12)
# For each basis field b_i (gradient-type and rotated-type),
# generate a callable dphi^k(b_i)(u,v) function using make_dphi_fn().
# Naming follows the same convention as build_bi_set():
#   "p_q.grad" â†’ non-rotated
#   "p_q.rot"  â†’ rotated

build_dphi_set <- function(basis_set, grad_f2k_fun, f2_fun) {
  keys <- sort(names(basis_set))
  dphi_set <- list()

  for (key in keys) {
    bs <- basis_set[[key]]
    # Gradient-type (non-rotated)
    dphi_set[[paste0(key, ".grad")]] <-
      make_dphi_fn(basis = bs,
                   f2_fun = f2_fun,
                   grad_f2k_fun = grad_f2k_fun,
                   rotated = FALSE)

    # Rotated-type (rotated)
    dphi_set[[paste0(key, ".rot")]]  <-
      make_dphi_fn(basis = bs,
                   f2_fun = f2_fun,
                   grad_f2k_fun = grad_f2k_fun,
                   rotated = TRUE)
  }

  return(dphi_set)
}




# ---------------------------
# (10) Compute <f1 - f2^k, dphi(b_i)>_{L2} for all basis
# ---------------------------
# Inputs:
#   f1, f2_k   : functions mapping (u,v) -> numeric(3)
#   dphi_list  : list of functions (each maps (u,v) -> numeric(3))
#   uv_grid    : matrix [n, 2] or list of sampling coordinates
#   area_weight: scalar or vector weights for numerical integration (default uniform)
# Output:
#   named numeric vector of inner products (same order as dphi_list)

compute_inner_products <- function(f1, f2_k, dphi_list, uv_grid, area_weight = NULL) {

  uv_grid <- as.matrix(uv_grid)
  n <- nrow(uv_grid)

  # Uniform weights if none provided (L2 integral â‰ˆ mean)
  if (is.null(area_weight)) area_weight <- rep(1/n, n)

  # Evaluate f1 - f2^k on all points
  diff_vals <- t(apply(uv_grid, 1, function(x) f1(x[1], x[2]) - f2_k(x[1], x[2]))) # nÃ—3

  # Compute inner products for each dphi_i
  coef_vec <- vapply(names(dphi_list), function(nm) {
    dphi_i <- dphi_list[[nm]]

    # Evaluate dphi_i(u,v) over the grid
    dphi_vals <- t(apply(uv_grid, 1, function(x) dphi_i(x[1], x[2])))  # nÃ—3

    # L2 inner product: sum_i ( (f1 - f2)Â·dphi_i ) * Î”A
    inner_val <- sum(rowSums(diff_vals * dphi_vals) * area_weight)
    inner_val
  }, numeric(1))

  coef_vec
}



####################################################
# ---------------------------
# Single gamma update: Î³^{(k+1)} = Î³^{(k)} âˆ˜ (Î³_id + Îµ Î´Î³^{(k)})
# ---------------------------
# Inputs:
#   gamma_k      : current gamma function, takes matrix uv -> matrix [n,2]
#   delta_gamma_k: deformation field function, takes matrix uv -> matrix [n,2]
#   epsilon       : small step size
# Output:
#   A new function gamma_{k+1}(uv) that applies the composition

update_gamma_fn <- function(gamma_k, delta_gamma_k, epsilon) {

  # old_gamma <- gamma_k
  # gamma_id(u,v) = identity map on U
  gamma_id <- function(u, v) c(u,v)

  # Return a new callable function of (u,v) as matrix
  function(u, v) {
    # Compute incremental deformation: Î³_id + Îµ Î´Î³
    inc <- gamma_id(u, v) + epsilon * delta_gamma_k(u, v)

    # Apply right-composition: Î³_{k+1}(u,v) = Î³_k( inc )
    gamma_k(inc[1],inc[2])
  }
}




