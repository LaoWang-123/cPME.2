### Simulation helper functions

#' Make a UV Grid
#'
#' @param n_u Number of grid values along the u direction.
#' @param n_v Number of grid values along the v direction.
#' @param grid_type Either `"disk"` or `"square"`.
#' @param center Optional disk center.
#' @param radius Optional disk radius.
#' @param xmin Minimum u value.
#' @param xmax Maximum u value.
#' @param ymin Minimum v value.
#' @param ymax Maximum v value.
#'
#' @return A data frame with columns `u` and `v`.
#' @export
make_uv_grid <- function(n_u = 60,
                         n_v = 60,
                         grid_type = c("disk", "square"),
                         center = NULL,
                         radius = NULL,
                         xmin = 0,
                         xmax = 1,
                         ymin = 0,
                         ymax = 1) {
  grid_type <- match.arg(grid_type)

  uv <- expand.grid(
    u = seq(xmin, xmax, length.out = n_u),
    v = seq(ymin, ymax, length.out = n_v)
  )

  if (is.null(center)) {
    center <- c((xmin + xmax) / 2, (ymin + ymax) / 2)
  }
  if (is.null(radius)) {
    radius <- min(xmax - xmin, ymax - ymin) / 2
  }

  if (grid_type == "disk") {
    cu <- center[1]
    cv <- center[2]
    uv <- subset(uv, (u - cu)^2 + (v - cv)^2 <= radius^2)
  }

  uv
}

#' Generate Simulated Surface Data
#'
#' @param f Surface function.
#' @param uv_grid Optional data frame with `u` and `v` columns.
#' @param n_u Number of grid values along the u direction when `uv_grid` is not
#' supplied.
#' @param n_v Number of grid values along the v direction when `uv_grid` is not
#' supplied.
#' @param n_points Optional number of points sampled from `uv_grid`.
#' @param noise_sd Gaussian noise standard deviation added to xyz coordinates.
#' @param seed Optional random seed.
#' @param f_input Function input convention: `"auto"`, `"uv"`, or `"vector"`.
#'
#' @return A `surface_data` list with `UV_grid`, `XYZ`, `n_points`, and
#' `noise_sd`.
#' @export
generate_surface_data <- function(f,
                                  uv_grid = NULL,
                                  n_u = 50,
                                  n_v = 50,
                                  n_points = NULL,
                                  noise_sd = 0,
                                  seed = NULL,
                                  f_input = c("auto", "uv", "vector")) {
  f_input <- match.arg(f_input)

  if (!is.null(seed)) {
    set.seed(seed)
  }

  if (is.null(uv_grid)) {
    uv_grid <- expand.grid(
      u = seq(0, 1, length.out = n_u),
      v = seq(0, 1, length.out = n_v)
    )
  }

  if (!is.null(n_points)) {
    idx <- sample(seq_len(nrow(uv_grid)), size = n_points, replace = TRUE)
    uv_grid <- uv_grid[idx, , drop = FALSE]
  }

  f_use <- switch(
    f_input,
    uv = function(uv) {
      uv <- as.numeric(uv)
      f(uv[1], uv[2])
    },
    vector = function(uv) {
      uv <- as.numeric(uv)
      f(uv)
    },
    auto = function(uv) {
      uv <- as.numeric(uv)
      tryCatch(
        f(uv[1], uv[2]),
        error = function(e) f(uv)
      )
    }
  )

  xyz <- t(apply(uv_grid, 1, f_use))
  colnames(xyz) <- c("x", "y", "z")

  if (noise_sd > 0) {
    xyz <- xyz + matrix(rnorm(length(xyz), sd = noise_sd), ncol = 3)
  }

  data_list <- list(
    UV_grid = uv_grid,
    XYZ = xyz,
    n_points = nrow(uv_grid),
    noise_sd = noise_sd
  )

  class(data_list) <- "surface_data"
  data_list
}

#' Contamination Field From Gaussian Bumps
#'
#' @param A Numeric vector of bump amplitudes.
#' @param u0 Numeric vector of bump centers in the u direction.
#' @param v0 Numeric vector of bump centers in the v direction.
#' @param sigma Numeric vector of positive bump bandwidths.
#'
#' @return A function mapping `(u, v)` to an xyz deformation matrix.
#' @export
make_bump_field <- function(A, u0, v0, sigma) {
  K <- length(A)
  stopifnot(length(u0) == K, length(v0) == K, length(sigma) == K)

  function(u, v) {
    z <- sum(
      A * exp(-((u - u0)^2 + (v - v0)^2) / (2 * sigma^2))
    )
    cbind(0, 0, z)
  }
}

#' Combine Two Surface Functions
#'
#' @param f1 First function mapping `(u, v)` to an n by 3 matrix.
#' @param f2 Second function mapping `(u, v)` to an n by 3 matrix.
#'
#' @return A function mapping `(u, v)` to the pointwise sum of `f1` and `f2`.
#' @export
combine_surface_functions <- function(f1, f2) {
  function(u, v) {
    out1 <- f1(u, v)
    out2 <- f2(u, v)

    stopifnot(
      is.matrix(out1),
      is.matrix(out2),
      ncol(out1) == 3,
      ncol(out2) == 3,
      nrow(out1) == nrow(out2)
    )

    out1 + out2
  }
}
