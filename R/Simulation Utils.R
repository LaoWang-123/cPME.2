### Simulation Datasets

#' Translate surface with a vector
#'
#' @param f u,v input function
#' @param translate translate vector
#'
#' @returns a translated function
#' @export
#'
translate_surface_function <- function(f, translate = c(0,0,0)) {
  stopifnot(is.numeric(translate), length(translate) == 3)

  function(u, v) {
    xyz <- f(u, v)
    xyz <- sweep(xyz, 2, translate, "+")
    xyz
  }
}


#' Make UV grid
#'
#' @param n_u 60
#' @param n_v 60
#' @param grid_type "disk"
#' @param center c(0.5, 0.5)
#' @param radius 0.5
#' @param xmin 0
#' @param xmax 1
#' @param ymin 0
#' @param ymax 1
#'
#' @returns a data.frame contains grid points
#' @export
make_uv_grid <- function(n_u = 60, n_v = 60,
                         grid_type = c("disk", "square"),
                         center = c(0.5, 0.5),
                         radius = 0.5,xmin=0,xmax=1,ymin=0,ymax=1) {
  grid_type <- match.arg(grid_type)

  uv <- expand.grid(
    u = seq(xmin, xmax, length.out = n_u),
    v = seq(ymin, ymax, length.out = n_v)
  )

  if (grid_type == "disk") {
    cu <- center[1]; cv <- center[2]
    uv <- subset(uv, (u - cu)^2 + (v - cv)^2 <= radius^2)
  }

  uv
}



# We define the bowl surface function and gamma function from surface_functions.R
# ==========================================================
# Generate simulated surface data
# ==========================================================
#' Generate simulated surface data
#'
#' @param f surface function, (u,v) input and (x,y,z) output
#' @param uv_grid a data.frame with u,v columns
#' @param n_u default = 50
#' @param n_v default = 50
#' @param n_points = NULL
#' @param noise_sd default = 0
#' @param seed default = NULL
#' @param f_input  c("auto", "uv", "vector")
#'
#' @returns a list of data to generate plots
#' @export
#'
#' @examples
#' data_f1 <- generate_surface_data(f1, n_u = 60, n_v = 60, noise_sd = 0, seed = 123) # R=1.2
generate_surface_data <- function(f,
                                  uv_grid = NULL,
                                  n_u = 50,
                                  n_v = 50,
                                  n_points = NULL,
                                  noise_sd = 0,
                                  seed = NULL,
                                  f_input = c("auto", "uv", "vector")) {
  f_input <- match.arg(f_input)

  # Optional reproducibility
  if (!is.null(seed)) set.seed(seed)

  if(is.null(uv_grid)){
    # Generate UV grid
    u_seq <- seq(0, 1, length.out = n_u)
    v_seq <- seq(0, 1, length.out = n_v)
    uv_grid <- expand.grid(u = u_seq, v = v_seq)
  }

  # ---- NEW: sample from uv_grid ----
  if (!is.null(n_points)) {
    idx <- sample(seq_len(nrow(uv_grid)), size = n_points, replace = TRUE)
    uv_grid <- uv_grid[idx, , drop = FALSE]
  }
  # ---------------------------------

  f_use <- switch(
    f_input,
    "uv" = function(uv) {
      uv <- as.numeric(uv)
      f(uv[1], uv[2])
    },
    "vector" = function(uv) {
      uv <- as.numeric(uv)
      f(uv)
    },
    "auto" = function(uv) {
      uv <- as.numeric(uv)
      tryCatch(
        f(uv[1], uv[2]),
        error = function(e) f(uv)
      )
    }
  )

  # Evaluate surface function f(u, v)
  xyz <- t(apply(uv_grid, 1, f_use))
  colnames(xyz) <- c("x", "y", "z")

  # Add Gaussian noise if needed
  if (noise_sd > 0) {
    xyz <- xyz + matrix(rnorm(length(xyz), sd = noise_sd), ncol = 3)
  }

  # Combine results into a list
  data_list <- list(
    UV_grid = uv_grid,
    XYZ = xyz,
    n_points = nrow(uv_grid),
    noise_sd = noise_sd
  )

  class(data_list) <- "surface_data"
  return(data_list)
}


##### Contamination in simulations


#' Contamination field defined by a sum of Gaussian bumps (z-axis deformation)
#'
#' This function constructs a smooth contamination vector field
#' \eqn{c : [0,1]^2 \to \mathbb{R}^3} acting along the z-axis,
#' defined as a linear combination of Gaussian radial basis functions:
#' \deqn{
#' h(u,v) = \sum_{k=1}^K A_k \exp\left(
#' -\frac{(u - u_k)^2 + (v - v_k)^2}{2\sigma_k^2}
#' \right).
#' }
#'
#' This construction is commonly used to model structured contamination.
#'
#' @param A Numeric vector of amplitudes \eqn{A_k}. Can contain both positive
#'   and negative values to represent upward or downward deformation.
#' @param u0 Numeric vector of bump centers in the \eqn{u}-direction.
#' @param v0 Numeric vector of bump centers in the \eqn{v}-direction.
#' @param sigma Numeric vector of positive bandwidth parameters \eqn{\sigma_k},
#'   controlling the spatial spread of each bump.
#'
#' @returns A function \code{f(u, v)} that evaluates the contamination field.
#' @export
#'
#' @examples
#' # Two bumps: one positive, one negative
#' A     <- c(0.15, -0.12)
#' u0    <- c(0.3, 0.7)
#' v0    <- c(0.4, 0.6)
#' sigma <- c(0.15, 0.2)
#'
#' contam_fn <- make_bump_field(A, u0, v0, sigma)
#'
make_bump_field <- function(A, u0, v0, sigma) {
  K <- length(A)
  stopifnot(length(u0) == K,
            length(v0) == K,
            length(sigma) == K)

  function(u, v) {
    z <- sum(
      A * exp(-((u - u0)^2 + (v - v0)^2) / (2 * sigma^2))
    )
    cbind(0, 0, z)
  }
}

#' Combine two surface/field functions by pointwise addition
#'
#' Given two functions f1 and f2 mapping (u, v) -> R^3,
#' construct a new function
#'
#'   f(u,v) = f1(u,v) + f2(u,v)
#'
#' @param f1 function(u, v) returning n x 3 matrix
#' @param f2 function(u, v) returning n x 3 matrix
#'
#' @returns A function (u, v) -> n x 3 matrix
#' @export
combine_surface_functions <- function(f1, f2) {

  function(u, v) {
    out1 <- f1(u, v)
    out2 <- f2(u, v)

    stopifnot(
      is.matrix(out1), is.matrix(out2),
      ncol(out1) == 3, ncol(out2) == 3,
      nrow(out1) == nrow(out2)
    )

    out1 + out2
  }
}
