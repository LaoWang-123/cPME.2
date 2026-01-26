#### This file contains additional functions used for pme registration

############# PME surface function
#################################################################################################################
######################
#########################################
####################
# Here are some functions based on pme R package
# library(pme)

# ---------------------------------------------------------
# 1. PME INITIAL GUESS VIA ISOMAP or PCA
# ---------------------------------------------------------
#' PME INITIAL GUESS VIA ISOMAP or PCA
#'
#' @param X dataset
#' @param d default=2
#' @param method method = c("pca","isomap")
#'
#' @returns initialization of d embedding
#' @import vegan
#' @export
#'
pme_initial_guess <- function(X, d, method = c("pca","isomap")) {

  method <- match.arg(method)

  if (method == "isomap") {

    # ---- ISOMAP 初始值 ----
    dissimilarity <- as.matrix(stats::dist(X))

    iso_obj <- vegan::isomap(
      dissimilarity,
      ndim = d,
      k = floor(sqrt(nrow(dissimilarity)))
    )

    U0 <- as.matrix(iso_obj$points)[, 1:d]

  } else if (method == "pca") {


    pca_obj <- stats::prcomp(X, center = TRUE, scale. = FALSE)
    U0 <- pca_obj$x[, 1:d]
  }

  return(U0)
}

#####
# We copied this from pme.R,  projection_pme can be used directly, which is exported by pme package
#' Calculate a New Parameterization
#'
#' @param f Embedding map.
#' @param X Numeric matrix of high-dimensional data.
#' @param t Numeric matrix of initial low-dimensional parameterizations.
#'
#' @return A numeric matrix of parameterizations.
#'
calc_params <- function(f, X, init_params) {
  params <- purrr::map(1:nrow(X), ~ projection_pme(X[.x, ], f, init_params[.x, ])) %>%
    unlist() %>%
    matrix(nrow = nrow(X), byrow = TRUE)
  params
}

########################################
######### function to scale the 2d parameters of pme projection

scale_uniform_square_with_params <- function(U) {
  U <- as.matrix(U)
  x <- U[,1]; y <- U[,2]

  min_x <- min(x); max_x <- max(x)
  min_y <- min(y); max_y <- max(y)

  width  <- max_x - min_x
  height <- max_y - min_y

  k <- max(width, height)
  scale_factor <- 1 / k

  # forward transform:  s = A_s t + b_s
  A_s <- scale_factor * diag(2)
  b_s <- c(-min_x * scale_factor, -min_y * scale_factor)

  # inverse transform:  t = A s + b
  A   <- diag(2) * k
  b   <- c(min_x, min_y)

  s_x <- (x - min_x) * scale_factor
  s_y <- (y - min_y) * scale_factor

  list(
    U_scaled = cbind(s_x, s_y),
    A = A,
    b = b
  )
}

#############################################################################################
##############################################################################################
#### We need to define pme surface function and grad pme surface function

pme_embedding_factory <- function(pme_result,d=2,A=diag(2),b=0){
  # etaFunc export from pme package
  f_embedding <- function(u,v) {
    parameters=as.vector(A %*% c(u,v) +b)

    as.vector(
      (t(pme_result$kernel_coefs) %*% etaFunc(parameters,  pme_result$params_opt, 4 - d)) +
        (t(pme_result$polynomial_coefs) %*% matrix(c(1, parameters), ncol = 1))
    )
  }
  return(f_embedding)
}

pme_grad_factory <- function(pme_result,A=diag(2),b=0) {

  parameterization <- pme_result$params_opt        # n × 2
  kernel_coefs <- pme_result$kernel_coefs            # n × 3
  poly_coefs <- pme_result$polynomial_coefs

  # gradient of η₂(r)
  grad_eta <- function(r) {
    rnorm <- sqrt(sum(r^2))
    if (rnorm == 0) return(c(0,0))
    return(2 * r * (1 + log(rnorm)))
  }

  grad_f <- function(u, v) {
    s    <- c(u, v)
    tvec <- as.vector(A %*% s + b)

    n <- nrow(parameterization)

    # G is n×2 matrix: every row = ∇η₂(t - parameterization[j,])
    G <- t(apply(parameterization, 1, function(uj) grad_eta(tvec - uj)))  # (n × 2)

    # kernel part: kernel_coefs^T G   → (3 × 2)
    grad_kernel <- t(kernel_coefs) %*% G

    # polynomial gradient part: rows = l=1..3, cols = (∂/∂u, ∂/∂v)
    grad_poly <- t(poly_coefs[2:3, ])   # (3 × 2)

    grad_t <- grad_kernel + grad_poly                  # 3×2  (∇ₜf)

    # chain rule: ∇ₛf = ∇ₜf ⋅ A
    grad_s <- grad_t %*% A
    return(grad_s)
  }

  return(grad_f)
}


evaluate_embedding <- function(grid, f) {
  pts <- t(apply(grid, 1, f))
  pts <- matrix(pts, ncol = 3)
  colnames(pts) <- c("x","y","z")
  cbind(grid, pts)
}
