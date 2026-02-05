#' Smoothing Kernel for Density Estimation
#'
#' Implements Gaussian kernel smoothing.
#'
#' @param x A vector of numeric values.
#' @param mu The mean of the Gaussian density.
#' @param sigma The standard deviation of the Gaussian density.
#'
#' @return A vector of numeric values.
#' @export
#'
smoothing_kernel <- function(x, mu, sigma) {
  yseq <- stats::dnorm((x - mu) / sigma)
  return((sigma^(-length(x))) * prod(yseq))
}


#' Smoothing Kernel for Density Estimation
#'
#' Implements Gaussian kernel smoothing on log-scale
#'
#' @param x A vector of numeric values.
#' @param mu The mean of the Gaussian density
#' @param sigma The standard deviation of the Gaussian density
#'
#' @return A numeric value
#' @export
#'
log_smoothing_kernel_r <- function(x, mu, sigma) {
  yseq <- stats::dnorm((x - mu) / sigma)
  output <- (-length(x) * log(sigma)) + sum(log(yseq))
  return(output)
}

#' Find the Coefficients of a Weighted Spline Function
#'
#' @param E A numeric matrix.
#' @param W A numeric matrix.
#' @param t_val A numeric matrix.
#' @param X A numeric matrix.
#' @param w The smoothing parameter.
#' @param d The intrinsic dimension.
#' @param D The dimension of the higher dimensional space.
#'
#' @return A numeric matrix.
#' @export
solve_weighted_spline_r <- function(E, W, t_val, X, w, d, D) {
  M1 <- cbind(
    2 * E %*% W %*% E + 2 * w * E,
    2 * E %*% W %*% t_val,
    t_val
  )
  M2 <- cbind(
    2 * t(t_val) %*% W %*% E,
    2 * t(t_val) %*% W %*% t_val,
    matrix(0, ncol = d + 1, nrow = d + 1)
  )
  M3 <- cbind(
    t(t_val),
    matrix(0, ncol = d + 1, nrow = d + 1),
    matrix(0, ncol = d + 1, nrow = d + 1)
  )
  M <- rbind(M1, M2, M3)

  b <- rbind(
    2 * E %*% W %*% X,
    2 * t(t_val) %*% W %*% X,
    matrix(0, nrow = d + 1, ncol = D)
  )
  sol <- MASS::ginv(M) %*% b
  sol
}

#' Solve for Smoothing Spline Coefficients
#'
#' @param E A matrix of values.
#' @param t_val A numeric matrix of input values.
#' @param X A numeric matrix of outputs.
#' @param w The smoothing parameter.
#' @param d The dimension of the input.
#' @param D The dimension of the output.
#'
#' @return A numeric matrix of coefficients.
#' @export
solve_spline <- function(E, t_val, X, w, d, D) {
  M1 <- cbind(E + (w * diag(rep(1, nrow(t_val)))), t_val)
  M2 <- cbind(t(t_val), matrix(0, ncol = d + 1, nrow = d + 1))
  M <- rbind(M1, M2)
  b <- rbind(X, matrix(0, nrow = d + 1, ncol = D))
  sol <- MASS::ginv(M) %*% b
  sol
}


#' Project onto Low-Dimensional Manifold
#'
#' @param x A data point in high-dimensional space.
#' @param f An embedding function.
#' @param initial_guess Guess of the parameterization.
#'
#' @return A vector describing the data point in low-dimensional space.
#' @export
projection_pme <- function(x, f, initial_guess) {
  nlm_est <- try(
    stats::nlm(
      function(t) dist_euclidean(x = x, f(t)),
      p = initial_guess
    ),
    silent = TRUE
  )

  if (inherits(nlm_est, "try-error")) {
    opts <- list("algorithm" = "NLOPT_LN_COBYLA", "xtol_rel" = 1e-10)
    nlopt_est <- try(
      nloptr::nloptr(
        x0 <- initial_guess,
        function(t) dist_euclidean(x = x, f(t)),
        opts = opts
      ),
      silent = TRUE
    )
    if (inherits(nlopt_est, "try-error")) {
      return(NULL)
    } else {
      return(nlopt_est$solution)
    }
  } else {
    return(nlm_est$estimate)
  }
}

#' Project onto Low-Dimensional Manifold
#'
#' @param x A value
#' @param f A value
#' @param initial_guess A value
#' @param n_knots A value
#' @param d_new A value
#' @param gamma A value
#'
#' @return A value
#' @export
projection_lpme <- function(x, f, initial_guess, n_knots, d_new, gamma) {
  nlm_est <- try(
    stats::nlm(
      function(t) dist_euclidean(x = x, f(matrix(c(initial_guess[1], t), nrow = 1))),
      p = initial_guess[-1]
    ),
    silent = TRUE
  )
  if (inherits(nlm_est, "try-error")) {
    opts <- list("algorithm" = "NLOPT_LN_COBYLA", "xtol_rel" = 1e-07)
    nlopt_est <- try(
      nloptr::nloptr(
        x0 = initial_guess[-1],
        function(t) dist_euclidean(x = x, f(c(initial_guess[1], t))),
        opts = opts
      ),
      silent = TRUE
    )
    if (inherits(nlopt_est, "try-error")) {
      return(NULL)
    } else {
      return(c(initial_guess[1], nlopt_est$solution))
    }
  } else {
    return(c(initial_guess[1], nlm_est$estimate))
  }
}

#' Compute Spline Coefficients
#'
#' @param X Numeric matrix of high-dimensional data.
#' @param t Numeric matrix of low-dimensional paramterizations.
#' @param weights Numeric matrix of cluster weights.
#' @param w smoothing parameter.
#'
#' @return A matrix of spline coefficients.
#'
#' @noRd
calc_coefficients <- function(X, t, weights, w) {
  t_val <- cbind(rep(1, nrow(t)), t)
  E <- calcE(t, 4 - ncol(t))
  solve_weighted_spline(E, weights, t_val, X, w, ncol(t), ncol(X))
}

#' Calculate a New Parameterization
#'
#' @param f Embedding map.
#' @param X Numeric matrix of high-dimensional data.
#' @param init_params Numeric matrix of initial low-dimensional parameterizations.
#' @param f_input default c("vector", "uv"), when chose "uv", it will divide the input into u and v.
#' @return A numeric matrix of parameterizations.
#'
#' @export
calc_params <- function(f, X, init_params, f_input = c("vector", "uv")) {
  f_input <- match.arg(f_input)

  f_use <- if (f_input == "uv") {
    function(p) {
      p <- as.numeric(p)
      f(p[1], p[2])
    }
  } else {
    f
  }

  params <- purrr::map(1:nrow(X), ~ projection_pme(X[.x, ], f_use, init_params[.x, ])) %>%
    unlist() %>%
    matrix(nrow = nrow(X), byrow = TRUE)
  params
}

#' Calculate Sum of Squared Distances
#'
#' @param f Embedding map.
#' @param X Numeric matrix of high-dimensional data.
#' @param t Numeric matrix of low-dimensional paramterizations.
#'
#' @return A numeric value.
#'
#' @noRd
calc_SSD <- function(f, X, t) {
  SSD_val <- purrr::map(1:nrow(X), ~ dist_euclidean(X[.x, ], f(t[.x, ]))^2) %>%
    unlist() %>%
    sum()
  SSD_val
}

#' Print the Results of a PME Iteration
#'
#' @param w The smoothing parameter.
#' @param SSD_new The most recent SSD calculation.
#' @param SSD_ratio The ratio of the new and old SSD values.
#' @param count The iteration number.
#'
#' @noRd
print_SSD <- function(tuning_val, SSD_new, SSD_ratio, count) {
  print(
    paste0(
      "For tuning parameter ",
      as.character(round(tuning_val, 4)),
      ", iteration #",
      as.character(count),
      " gives SSD = ",
      as.character(round(SSD_new, 4)),
      " and SSD_ratio = ",
      as.character(round(SSD_ratio, 4)),
      "."
    )
  )
}

#' Plot a PME Object
#'
#' @param f A function of the embedding map.
#' @param x A numeric matrix of the unreduced data.
#' @param centers A numeric matrix of the mixture component centers.
#' @param sol A numeric matrix of the embedding map coefficients.
#' @param t A numeric matrix of the parameterization of the component centers.
#' @param d The intrinsic dimension.
#'
#' @noRd
plot_pme <- function(f, x, centers, sol, t, d) {
  # pred_grid <- calc_tnew(centers, t, sol, I, d, lambda)
  pred_grid <- calc_params(f, centers, t)
  r_bounds <- Rfast::colMinsMaxs(pred_grid)
  r_list <- list()
  for (idx in 1:dim(r_bounds)[2]) {
    r_list[[idx]] <- seq(
      r_bounds[1, idx],
      r_bounds[2, idx],
      length.out = nrow(centers)
    )
  }
  r_mat <- as.matrix(expand.grid(r_list))

  pred_grid <- r_mat
  f_pred <- purrr::map(
    1:nrow(pred_grid),
    ~ f(unlist(as.vector(pred_grid[.x, ]))) %>%
      as.vector()
  ) %>%
    unlist() %>%
    matrix(nrow = nrow(pred_grid), byrow = TRUE)

  f_pred_full <- cbind(pred_grid, f_pred)

  if (dim(x)[2] == 2) {
    plt <- ggplot2::ggplot() +
      ggplot2::geom_point(
        ggplot2::aes(
          x = x[, 1],
          y = x[, 2]
        ),
        alpha = 0.5
      ) +
      ggplot2::geom_point(
        ggplot2::aes(
          x = f_pred_full[, d + 1],
          y = f_pred_full[, d + 2]
        ),
        color = "red"
      )
    print(plt)
  } else if (dim(x)[2] >= 3) {
    plt <- plotly::plot_ly(
      x = f_pred_full[, d + 1],
      y = f_pred_full[, d + 2],
      z = f_pred_full[, d + 3],
      type = "scatter3d",
      mode = "markers",
      opacity = 0.5
    ) %>%
      plotly::add_markers(
        x = x[, 1],
        y = x[, 2],
        z = x[, 3],
        opacity = 0.15,
        marker = list(size = 3)
      )
    print(plt)
  }
}

#' Calculate Mean Squared Distance
#'
#' @param x A numeric matrix of data
#' @param km A kmeans object describing reduced data.
#' @param f An embedding map.
#' @param t A numeric matrix of parameterizations for the mixture components
#' described in `km`.
#' @param D The dimension of the original dataset.
#' @param d The intrinsic dimension.
#'
#' @return A numeric value.
#'
#' @noRd
calc_msd <- function(x, km, f, t, D, d) {
  data_initial <- matrix(0, nrow = 1, ncol = D + d)
  center_order <- order(km$centers[, 1])
  for (i in 1:max(km$cluster)) {
    index_temp <- which(km$cluster == center_order[i])
    length_temp <- length(index_temp)
    temp_x <- matrix(x[index_temp, ], nrow = length_temp)
    t_temp <- matrix(rep(t[i, 1], length_temp))
    for (j in 1:d) {
      t_temp <- cbind(t_temp, rep(t[i, j], length_temp))
    }
    t_temp <- matrix(t_temp[, -1], nrow = length_temp)
    data_initial <- rbind(data_initial, cbind(temp_x, t_temp))
  }
  data_initial <- data_initial[-1, ]
  proj_para <- purrr::map(
    1:nrow(data_initial),
    ~ projection_pme(data_initial[.x, 1:D], f, data_initial[.x, (D + 1):(D + d)]) %>%
      t()
  ) %>%
    purrr::reduce(rbind)
  proj_points <- purrr::map(
    1:nrow(proj_para),
    ~ f(proj_para[.x, ]) %>%
      t()
  ) %>%
    purrr::reduce(rbind)

  mse <- purrr::map(
    1:nrow(data_initial),
    ~ dist_euclidean(data_initial[.x, 1:D], proj_points[.x, ])^2
  ) %>%
    unlist() %>%
    mean()

  mse
}
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


########################################

#' Generate surface points with embedding function
#'
#' @param grid a data.frame generated by expand.grid
#' @param f embedding function
#'
#' @returns a data.frame with the surface points
#' @export
#'
evaluate_embedding <- function(grid, f) {
  pts <- t(apply(grid, 1, f))
  pts <- matrix(pts, ncol = 3)
  colnames(pts) <- c("x","y","z")
  cbind(grid, pts)
}

########################################
######### function to scale the 2d parameters of pme projection

#' Scale the 2d parameters to [0,1]^2 by a linear transformation
#'
#' @param U original 2d parameters
#'
#' @returns a new 2d parameters in domain
#' @export
#'
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

#### Thesse two functions are used for post-pme fitting and if we want to rescale its domain to [0,1]^2,
#### then we need to use them.

### Latest note on 2/4/2026, I decide to rescale the parameterization domain in the initialization step instead.

#' Make a pme function with pme result and rescale matrix
#'
#' @param pme_result output from pme
#' @param d default 2
#' @param A a 2x2 matrix for scale
#' @param b a 2x1 vector
#'
#' @returns a f(u,v) function
#' @export
#'
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

#' Get a grad function with pme and scale matrix
#'
#' @param pme_result output from pme
#' @param A a 2x2 matrix for scale
#' @param b a 2x1 vector
#'
#' @returns a g(u,v) function
#' @export
#'
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


