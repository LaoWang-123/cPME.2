# ==========================================================
# γ function factory —  "smooth" / "rotate" / "identity"
# ==========================================================
#' Ganna function factory
#'
#' @param mode mode choice
#' @param au default 0.3
#' @param av default 0.25
#' @param angle
#'
#' @returns a chosen gamma function
#' @export
#'
gamma_function_factory <- function(mode = c("smooth", "rotate", "identity"),
                                   au = 0.3, av = 0.25, angle = 0.5*pi) {
  mode <- match.arg(mode)
  .au <- au; .av <- av; .angle <- angle
  force(.au); force(.av); force(.angle)

  if (mode == "smooth") {
    fun <- function(uv, au = .au, av = .av) {
      u <- uv[,1]; v <- uv[,2]
      bu <- u * (1 - u); bv <- v * (1 - v)
      du <- au * bu * sin(2*pi*v)
      dv <- av * bv * sin(2*pi*u)
      u2 <- pmin(pmax(u + du, 0), 1)
      v2 <- pmin(pmax(v + dv, 0), 1)
      cbind(u2, v2)
    }
  } else if (mode == "rotate") { # counterclockwise rotation
    fun <- function(uv, angle = .angle) {
      u <- uv[,1]; v <- uv[,2]
      uc <- u - 0.5; vc <- v - 0.5
      cosA <- cos(angle); sinA <- sin(angle)
      u2 <- pmin(pmax(cosA*uc - sinA*vc + 0.5, 0), 1)
      v2 <- pmin(pmax(sinA*uc + cosA*vc + 0.5, 0), 1)
      cbind(u2, v2)
    }
  } else { # identity
    fun <- function(uv) {
      cbind(uv[,1], uv[,2])
    }
  }

  attr(fun, "params") <- list(au = .au, av = .av, angle = .angle, mode = mode)
  return(fun)
}

# ==========================================================
# bowl function
# ==========================================================
f_bowl <- function(uv, R = 1.0, c = 1.5) {
  u <- uv[,1]; v <- uv[,2]
  x <- (2*u - 1) * R
  y <- (2*v - 1) * R
  z <- c * ((2*u - 1)^2 + (2*v - 1)^2)
  cbind(x, y, z)
}

###################
# surface function factor
##################
make_surface_function <- function(R=1,
                                  c = 1.5,
                                  gamma_mode = "identity",
                                  au = 0.3, av = 0.25,
                                  angle = 0.5*pi) {

  gamma_fun <- gamma_function_factory(mode = gamma_mode,
                                      au = au, av = av, angle = angle)

  function(u,v){
    uv <- cbind(u, v)
    f_bowl(gamma_fun(uv), R = R, c = c)
  }

}

###### Example
# # f_m (rotate reparam)
# f1 <- make_surface_function(R = 1.2,gamma_mode = "rotate")
#
# # f_tau (identity)
# f2 <- make_surface_function(R = 1, gamma_mode = "identity")


##### grad_surface_function
grad_f_bowl <- function(u, v, R = 1.0, c = 1.5) {
  # ensure vectorized operation
  du <- 2 * R
  dv <- 2 * R
  dz_du <- 4 * c * (2 * u - 1)
  dz_dv <- 4 * c * (2 * v - 1)

  matrix(c(
    du, 0,
    0,  dv,
    dz_du, dz_dv
  ), nrow = 3, byrow = TRUE)
}

# f2_grad_fn <- grad_f_bowl



# We define the bowl surface function and gamma function from surface_functions.R
# ==========================================================
# Generate simulated surface data
# ==========================================================
generate_surface_data <- function(f,
                                  n_u = 50,
                                  n_v = 50,
                                  noise_sd = 0,
                                  seed = NULL) {
  # Optional reproducibility
  if (!is.null(seed)) set.seed(seed)

  # Generate UV grid
  u_seq <- seq(0, 1, length.out = n_u)
  v_seq <- seq(0, 1, length.out = n_v)
  UV_grid <- expand.grid(u = u_seq, v = v_seq)

  # Evaluate surface function f(u, v)
  xyz <- t(apply(UV_grid, 1, function(uv) {
    f(uv[1], uv[2])
  }))
  colnames(xyz) <- c("x", "y", "z")

  # Add Gaussian noise if needed
  if (noise_sd > 0) {
    xyz <- xyz + matrix(rnorm(length(xyz), sd = noise_sd), ncol = 3)
  }

  # Combine results into a list
  data_list <- list(
    UV_grid = UV_grid,
    XYZ = xyz,
    n_points = nrow(UV_grid),
    noise_sd = noise_sd
  )

  class(data_list) <- "surface_data"
  return(data_list)
}

# Example usage:
# ==========================================================
# generate one dataset for f1 and f2
# ==========================================================
# data_f1 <- generate_surface_data(f1, n_u = 60, n_v = 60, noise_sd = 0, seed = 123) # R=1.2
# data_f2 <- generate_surface_data(f2, n_u = 60, n_v = 60, noise_sd = 0, seed = 123) # R=1
#
# # Inspect results
# str(data_f1)
# head(data_f1$XYZ)


################################################################################################################
##################################################################################################################

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


