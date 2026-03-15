### new Utils.R

### We optimize some functions and workflow to accelerate the computation speed

#' Build basis with grid points
#'
#' @param basis_set A list of basis functions.
#' @param Ugrid Data frame of (u,v) sample locations.
#' @param mode c("full","div_free")
#'
#' @returns a list of bi, calculated on grid points.
#' @export
#'
build_basis_grid <- function(basis_set, Ugrid, mode = c("full","div_free")) {
  mode <- match.arg(mode)
  out <- list()
  n  <- nrow(Ugrid)
  keys <- sort(names(basis_set))

  for (key in keys) {
    bs <- basis_set[[key]]

    psi_vals <- vapply(1:n, function(i) bs$psi(Ugrid[i,1], Ugrid[i,2]), numeric(1))

    # rot always needed if div_free
    rot_vals <- t(vapply(
      1:n,
      function(i) bs$rot_grad_psi(Ugrid[i,1], Ugrid[i,2]),
      numeric(2)
    ))
    rot_x <- rot_vals[,1]
    rot_y <- rot_vals[,2]

    if (mode == "full") {
      grad_vals <- t(vapply(
        1:n,
        function(i) bs$grad_psi(Ugrid[i,1], Ugrid[i,2]),
        numeric(2)
      ))
      grad_psi_x <- grad_vals[,1]
      grad_psi_y <- grad_vals[,2]
    } else {
      grad_psi_x <- NULL
      grad_psi_y <- NULL
    }

    out[[key]] <- list(
      psi        = psi_vals,
      grad_psi_x = grad_psi_x,
      grad_psi_y = grad_psi_y,
      rot_x      = rot_x,
      rot_y      = rot_y,
      lambda     = bs$lambda_pq,
      norm       = bs$norm_pq
    )
  }
  out
}

# build_basis_grid <- function(basis_set, Ugrid) {
#   out <- list()
#   n  <- nrow(Ugrid)
#
#   for (key in names(basis_set)) {
#     bs <- basis_set[[key]]
#
#     # psi(u,v)
#     psi_vals <- vapply(
#       1:n,
#       function(i) bs$psi(Ugrid[i,1], Ugrid[i,2]),
#       numeric(1)
#     )
#
#     # grad psi(u,v)
#     grad_vals <- t(vapply(
#       1:n,
#       function(i) bs$grad_psi(Ugrid[i,1], Ugrid[i,2]),
#       numeric(2)
#     ))
#     grad_psi_x <- grad_vals[,1]
#     grad_psi_y <- grad_vals[,2]
#
#     # rot grad psi(u,v)
#     rot_vals <- t(vapply(
#       1:n,
#       function(i) bs$rot_grad_psi(Ugrid[i,1], Ugrid[i,2]),
#       numeric(2)
#     ))
#     rot_x <- rot_vals[,1]
#     rot_y <- rot_vals[,2]
#
#     out[[key]] <- list(
#       psi        = psi_vals,     # n
#       grad_psi_x = grad_psi_x,   # n
#       grad_psi_y = grad_psi_y,   # n
#       rot_x      = rot_x,        # n
#       rot_y      = rot_y,        # n
#       lambda     = bs$lambda_pq,
#       norm       = bs$norm_pq
#     )
#   }
#
#   out
# }



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
