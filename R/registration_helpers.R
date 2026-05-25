# Current registration helper functions
#
# These helpers support the grid-precomputed safe `Registration` class. They are
# package-internal unless explicitly exported elsewhere.
###################################################
# ---------------------------
# (C) Assemble Î´Î³^k from coefficients and b_i
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
    # Evaluate all basis fields at (u,v) â†’ 2D matrix per basis
    # vapply ensures output is [2 Ã— n_basis] numeric matrix
    mat <- vapply(bi_names, function(nm) bi_set[[nm]](u, v), numeric(2))
    # Weighted sum across columns
    as.vector(mat %*% coefs)
  }
}

###################################################
# ---------------------------
# (C) Assemble DÎ´Î³^k from coefficients and D_b_i
# ---------------------------
#
#' Assemble DÎ´Î³^k from coefficients and D_b_i
#' DÎ´Î³^k will be used in assembling grad_f2k and used in assemble_grad_f2k_from_state function
#' @param dgamma_coefs coefs
#' @param D_bi_set D_bi sets
#'
#' @returns a function (u,v) to 2x2 matrix
#'
#' @examples
#' assemble_D_delta_gamma_fn(dgamma_coefs,D_bi_set)
#' @noRd
assemble_D_delta_gamma_fn <- function(dgamma_coefs, D_bi_set) {
  stopifnot(is.numeric(dgamma_coefs))
  bi_names <- names(D_bi_set)

  # Pre-extract coefficients as numeric vector in fixed order
  coefs <- as.numeric(dgamma_coefs[bi_names])

  function(u, v) {
    # Compute all c_i * D b_i(u,v) as a list of 2Ã—2 matrices
    mats <- lapply(seq_along(D_bi_set), function(i) {
      coefs[i] * D_bi_set[[bi_names[i]]](u, v)
    })
    # Sum all matrices in the list elementwise
    Reduce(`+`, mats)
  }
}




#####################################################
# Calculate âˆ‡f2^k
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
#' @noRd
assemble_grad_f2k_from_state <- function(state_list, gamma_k, f2_grad_fn, epsilon) {
  # Determine iteration depth
  k <- length(state_list) # state_0, state_1,... state_k-1

  # Extract all gamma and Ddelta_gamma functions
  gamma_list <- lapply(state_list, `[[`, "gamma_k")
  Ddelta_gamma_fns <- lapply(state_list, `[[`, "Ddelta_gamma_fn")

  # Return a function that maps (u,v) â†’ 3Ã—2 matrix
  function(u, v) {

    x_seq <- lapply(seq_len(k), function(t) gamma_list[[t]](u, v)) #x0=u, x1, x2,....x_k-1

    x_k <- gamma_k(u,v)
    grad_f <- f2_grad_fn(x_k[1], x_k[2]) # âˆ‡f2(x_k) (3Ã—2 matrix)

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
    #     Equivalent to: G = (I + ÎµÂ·D_{k-1}) %*% ... %*% (I + ÎµÂ·D_0)
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
#' @noRd
#'
compute_dphi_grid <- function(basis_grid, f2_grid, grad_f2_grid, mode = c("full", "div_free")) {
  # f2_grid: n Ã— 3
  # grad_f2_grid: list(Gx = nÃ—3, Gy = nÃ—3), gradient of f2 for each coord
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
      term1 <- (-0.5 * lam) * psi * f2_grid                 # nÃ—3
      term2 <- cbind(
        gx * grad_f2_grid$Gx[,1] + gy * grad_f2_grid$Gy[,1],
        gx * grad_f2_grid$Gx[,2] + gy * grad_f2_grid$Gy[,2],
        gx * grad_f2_grid$Gx[,3] + gy * grad_f2_grid$Gy[,3]
      )

      dphi_grad <- (term1 + term2) / norm   # nÃ—3
      out[[paste0(key,".grad")]] <- dphi_grad
    }

    ## ---- rotated-type (Eq. 12) ----
    term_rot <- cbind(
      rx * grad_f2_grid$Gx[,1] + ry * grad_f2_grid$Gy[,1],
      rx * grad_f2_grid$Gx[,2] + ry * grad_f2_grid$Gy[,2],
      rx * grad_f2_grid$Gx[,3] + ry * grad_f2_grid$Gy[,3]
    )

    dphi_rot <- term_rot / norm          # nÃ—3
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
#' @noRd
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

  # n Ã— 3
  colnames(mat) <- c("x","y","z")
  mat
}

make_grad_f2_grid <- function(grad_f2_fun, Ugrid) {
  n <- nrow(Ugrid)

  # 1. each element is 3x2 matrix
  G_list <- lapply(1:n, function(i) grad_f2_fun(Ugrid[i,1], Ugrid[i,2]))

  G0 <- G_list[[1]]
  if (!is.matrix(G0) || !all(dim(G0) == c(3,2))) {
    stop("grad_f2_fun must return a 3Ã—2 matrix.")
  }


  G_array <- simplify2array(G_list)

  Gx <- t(G_array[,1,])
  Gy <- t(G_array[,2,])

  list(Gx = Gx, Gy = Gy)
}

# ==========================================================
# Vectorized computation of E = âˆ« ||f1(u,v) - f2(u,v)||^2 du dv
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
# delta_gamma_fns: list of all Î´Î³^{(0)}, â€¦, Î´Î³^{(k-1)}
#' Composite gamma from delta_gamma history
#'
#' @param delta_gamma_fns a list from all delta gamma fns from state_list
#' @param eps default value=0.007
#'
#' @returns gamma^k function, (u,v) input (u,v) out
#' @noRd
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





