#### Utils.R
# This file contains all the public functions

# ==========================================================
# Vectorized computation of E = ∫ ||f1(u,v) - f2(u,v)||^2 du dv
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
# For each basis field b_i (either ∇ψ_i or *∇ψ_i), compute
#   dphi^k(b_i) = 1 / ||∇ψ_i|| * (-0.5 * λ_i * ψ_i * f2^k + ∇f2^k · ∇ψ_i)      # non-rotated
#   dphi^k(b_i) = 1 / ||∇ψ_i|| * (∇f2^k · *∇ψ_i)                              # rotated

# We get ∇f2^k(u,v) from assemble_delta_f2k_from_state function defined previously

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
make_dphi_fn <- function(basis, grad_f2k_fun, f2_fun, rotated = FALSE) {

  norm_pq <- basis$norm_pq  # ||∇ψ_pq||_L2
  lambda_pq <- basis$lambda_pq
  if (!rotated) {
    # Gradient-type basis (curl-free)
    # Eq. (11): 1/norm * (-0.5 * λ * ψ * f2^k + ∇f2^k · ∇ψ)
    psi_fun  <- basis$psi
    grad_fun <- basis$grad_psi

    # Return a function of (u,v)
    function(u, v) {
      (1 / norm_pq) * (t(-0.5 * lambda_pq * psi_fun(u, v) * f2_fun(u,v)) + ( grad_f2k_fun(u,v) %*% grad_fun(u, v)))
    }

  } else {
    # Rotated (divergence-free) basis
    # Eq. (12): 1/norm * (∇f2^k · *∇ψ)
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
#   "p_q.grad" → non-rotated
#   "p_q.rot"  → rotated

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

  # Uniform weights if none provided (L2 integral ≈ mean)
  if (is.null(area_weight)) area_weight <- rep(1/n, n)

  # Evaluate f1 - f2^k on all points
  diff_vals <- t(apply(uv_grid, 1, function(x) f1(x[1], x[2]) - f2_k(x[1], x[2]))) # n×3

  # Compute inner products for each dphi_i
  coef_vec <- vapply(names(dphi_list), function(nm) {
    dphi_i <- dphi_list[[nm]]

    # Evaluate dphi_i(u,v) over the grid
    dphi_vals <- t(apply(uv_grid, 1, function(x) dphi_i(x[1], x[2])))  # n×3

    # L2 inner product: sum_i ( (f1 - f2)·dphi_i ) * ΔA
    inner_val <- sum(rowSums(diff_vals * dphi_vals) * area_weight)
    inner_val
  }, numeric(1))

  coef_vec
}


###################################################
# ---------------------------
# (C) Assemble δγ^k from coefficients and b_i
# ---------------------------

# dgamma_coefs : named numeric vector, names must match bi names
# bi_set       : list returned by build_bi_set()
# Returns      : function(u,v) -> vector of length 2
assemble_delta_gamma_fn <- function(dgamma_coefs, bi_set) {
  stopifnot(is.numeric(dgamma_coefs))
  bi_names <- names(bi_set)

  # Pre-extract coefficients as numeric vector in fixed order
  coefs <- as.numeric(dgamma_coefs[bi_names])

  # Return function(u,v)
  function(u, v) {
    # Evaluate all basis fields at (u,v) → 2D matrix per basis
    # vapply ensures output is [2 × n_basis] numeric matrix
    mat <- vapply(bi_names, function(nm) bi_set[[nm]](u, v), numeric(2))
    # Weighted sum across columns
    as.vector(mat %*% coefs)
  }
}

###################################################
# ---------------------------
# (C) Assemble Dδγ^k from coefficients and D_b_i
# ---------------------------
#
#' Assemble Dδγ^k from coefficients and D_b_i
#' Dδγ^k will be used in assembling grad_f2k and used in assemble_grad_f2k_from_state function
#' @param dgamma_coefs coefs
#' @param D_bi_set D_bi sets
#'
#' @returns a function (u,v) to 2x2 matrix
#'
#' @examples
#' assemble_D_delta_gamma_fn(dgamma_coefs,D_bi_set)
assemble_D_delta_gamma_fn <- function(dgamma_coefs, D_bi_set) {
  stopifnot(is.numeric(dgamma_coefs))
  bi_names <- names(D_bi_set)

  # Pre-extract coefficients as numeric vector in fixed order
  coefs <- as.numeric(dgamma_coefs[bi_names])

  function(u, v) {
    # Compute all c_i * D b_i(u,v) as a list of 2×2 matrices
    mats <- lapply(seq_along(D_bi_set), function(i) {
      coefs[i] * D_bi_set[[bi_names[i]]](u, v)
    })
    # Sum all matrices in the list elementwise
    Reduce(`+`, mats)
  }
}




####################################################
# ---------------------------
# Single gamma update: γ^{(k+1)} = γ^{(k)} ∘ (γ_id + ε δγ^{(k)})
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
    # Compute incremental deformation: γ_id + ε δγ
    inc <- gamma_id(u, v) + epsilon * delta_gamma_k(u, v)

    # Apply right-composition: γ_{k+1}(u,v) = γ_k( inc )
    gamma_k(inc[1],inc[2])
  }
}



#####################################################
# Calculate ∇f2^k
# In this step, we have state_0, state_1,... state_k-1 in state_list
# The state_k, we have gamma_k, but we don't have delta gamma^k and D delta gamma^k
#' Title
#'
#' @param state_list state_list
#' @param gamma_k gamma_k
#' @param f2_grad_fn f2_grad_fn
#' @param epsilon step to upgrade
#'
#' @returns a function (u,v) to 3x2 matrix, named grad_f2k
#'
assemble_grad_f2k_from_state <- function(state_list, gamma_k, f2_grad_fn, epsilon) {
  # Determine iteration depth
  k <- length(state_list) # state_0, state_1,... state_k-1

  # Extract all gamma and Ddelta_gamma functions
  gamma_list <- lapply(state_list, `[[`, "gamma_k")
  Ddelta_gamma_fns <- lapply(state_list, `[[`, "Ddelta_gamma_fn")

  # Return a function that maps (u,v) → 3×2 matrix
  function(u, v) {

    x_seq <- lapply(seq_len(k), function(t) gamma_list[[t]](u, v)) #x0=u, x1, x2,....x_k-1

    x_k <- gamma_k(u,v)
    grad_f <- f2_grad_fn(x_k[1], x_k[2]) # ∇f2(x_k) (3×2 matrix)

    # Step 3. chain multiply Jacobians (reversed order)
    # G <- diag(2)
    # for (t in seq_len(k)) {
    #   idx <- k - t + 1
    #   x_input <- x_seq[[idx]]
    #   Ddelta <- Ddelta_gamma_fns[[t]](x_input[1], x_input[2])
    #   G <- G %*% (diag(2) + epsilon * Ddelta)
    # }
    #
    # grad_f %*% G
    ##
    # ---- Begin optimized chain computation ----
    x_seq_rev <- rev(x_seq)
    D_list <- mapply(
      function(f, x) f(x[1], x[2]),
      Ddelta_gamma_fns,
      x_seq_rev,
      SIMPLIFY = FALSE
    )

    mat_list <- lapply(D_list, function(D) diag(2) + epsilon * D)
    #     Equivalent to: G = (I + ε·D_{k-1}) %*% ... %*% (I + ε·D_0)
    G <- Reduce(`%*%`, mat_list)

    grad_f2k <- grad_f %*% G
    return(grad_f2k)
  }
}



### We optimize some functions and workflow to accelerate the computation speed




#' Compute dphi based on grid points
#'
#' @param basis_grid a list of bi on grid points.
#' @param f2_grid surface function on grid points.
#' @param grad_f2_grid gradient of surface function on grid points.
#' @param mode mode c("full","div_free")
#'
#' @returns out, used to input in compute_inner_products_fast
#' @export
#'
compute_dphi_grid <- function(basis_grid, f2_grid, grad_f2_grid, mode = c("full", "div_free")) {
  # f2_grid: n × 3
  # grad_f2_grid: list(Gx = n×3, Gy = n×3), gradient of f2 for each coord
  mode <- match.arg(mode)
  out <- list()

  for (key in names(basis_grid)) {
    bg <- basis_grid[[key]]

    psi  <- bg$psi
    gx   <- bg$grad_psi_x
    gy   <- bg$grad_psi_y
    rx   <- bg$rot_x
    ry   <- bg$rot_y
    lam  <- bg$lambda
    norm <- bg$norm

    if (mode == "full") {
      ## ---- gradient-type (Eq. 11) ----
      term1 <- (-0.5 * lam) * psi * f2_grid                 # n×3
      term2 <- cbind(
        gx * grad_f2_grid$Gx[,1] + gy * grad_f2_grid$Gy[,1],
        gx * grad_f2_grid$Gx[,2] + gy * grad_f2_grid$Gy[,2],
        gx * grad_f2_grid$Gx[,3] + gy * grad_f2_grid$Gy[,3]
      )

      dphi_grad <- (term1 + term2) / norm   # n×3
      out[[paste0(key,".grad")]] <- dphi_grad
    }

    ## ---- rotated-type (Eq. 12) ----
    term_rot <- cbind(
      rx * grad_f2_grid$Gx[,1] + ry * grad_f2_grid$Gy[,1],
      rx * grad_f2_grid$Gx[,2] + ry * grad_f2_grid$Gy[,2],
      rx * grad_f2_grid$Gx[,3] + ry * grad_f2_grid$Gy[,3]
    )

    dphi_rot <- term_rot / norm          # n×3
    out[[paste0(key,".rot")]]  <- dphi_rot
  }
  out
}


#' Compute the inner products based on grid points to get coefficients
#'
#' @param diff_grid Two surface function difference on grid points.
#' @param dphi_grid_list The output of compute_dphi_grid.
#' @param weight default 1/n.
#'
#' @returns a named vector of coefficients
#' @export
#'
compute_inner_products_fast <- function(diff_grid, dphi_grid_list, weight) {

  vapply(dphi_grid_list, function(mat3) {
    sum(rowSums(diff_grid * mat3) * weight)
  }, numeric(1))

}


make_f2_grid <- function(f2_fun, Ugrid) {
  n <- nrow(Ugrid)

  mat <- t(vapply(
    1:n,
    function(i) f2_fun(Ugrid[i,1], Ugrid[i,2]),
    numeric(3)
  ))

  # n × 3
  colnames(mat) <- c("x","y","z")
  mat
}

make_grad_f2_grid <- function(grad_f2_fun, Ugrid) {
  n <- nrow(Ugrid)

  # 1. each element is 3x2 matrix
  G_list <- lapply(1:n, function(i) grad_f2_fun(Ugrid[i,1], Ugrid[i,2]))

  G0 <- G_list[[1]]
  if (!is.matrix(G0) || !all(dim(G0) == c(3,2))) {
    stop("grad_f2_fun must return a 3×2 matrix.")
  }


  G_array <- simplify2array(G_list)

  Gx <- t(G_array[,1,])
  Gy <- t(G_array[,2,])

  list(Gx = Gx, Gy = Gy)
}

# ==========================================================
# Vectorized computation of E = ∫ ||f1(u,v) - f2(u,v)||^2 du dv
# ----------------------------------------------------------
# f1, f2 : functions (u,v) -> c(x,y,z)
# f1_grid, f2_grid, nx3 matrix

compute_E_grid <- function(f1_grid, f2_grid) {

  diff_sq <- rowSums((f1_grid - f2_grid)^2)

  E_val <- mean(diff_sq)
  return(E_val)
}


################################################################
### Change the method to composite gamma_k
##############################################################
# delta_gamma_fns: list of all δγ^{(0)}, …, δγ^{(k-1)}
#' Composite gamma from delta_gamma history
#'
#' @param delta_gamma_fns a list from all delta gamma fns from state_list
#' @param eps default value=0.007
#'
#' @returns gamma^k function, (u,v) input (u,v) out
#' @export
#'
#' @examples
#' delta_gamma_fns <- lapply(state_list, function(s) s$delta_gamma_fn)
#' gamma_next <- make_gamma_from_history(delta_gamma_fns = delta_gamma_fns,eps = 0.007)
#' f2_next <- function(u, v) {
#' xy <- gamma_next.2(u, v)
#' self$f2(xy[1], xy[2])
#' }
make_gamma_from_history <- function(delta_gamma_fns, eps) {
  function(u, v) {
    x <- u
    y <- v
    for (g in rev(delta_gamma_fns)) {
      dxy <- g(x, y)
      x <- x + eps * dxy[1]
      y <- y + eps * dxy[2]
    }
    c(x, y)
  }
}




