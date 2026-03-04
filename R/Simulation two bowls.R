### Simulation Datasets

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

# ==========================================================
# bowl function
# ==========================================================

#' Make bowl surface function
#'
#' @param R 1
#' @param c 1.5
#' @param p 1
#' @param flip 1
#'
#' @returns a function with (u,v) input
#' @export
#'
make_bowl_surface_function <- function(R = 1.0, c = 1.5, p = 1.0, flip = 1) {

  f_bowl <- function(uv, R = 1.0, c = 1.5, p = 1.0, flip = 1) {
    u <- uv[,1]; v <- uv[,2]
    x <- (2*u - 1) * R
    y <- (2*v - 1) * R
    r2 <- ((2*u - 1)^2 + (2*v - 1)^2)
    z  <- flip * c * (r2^p)
    cbind(x, y, z)
  }

  function(u, v) {
    uv <- cbind(u, v)
    f_bowl(uv, R = R, c = c, p = p, flip = flip)
  }
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

translate_surface_function <- function(f, translate = c(0,0,0)) {
  stopifnot(is.numeric(translate), length(translate) == 3)
  function(u, v) {
    xyz <- f(u, v)
    sweep(xyz, 2, translate, "+")
  }
}

generate_two_surface_simulation <- function(
    f1,
    f2,
    translate2 = c(0, 0, 0),   # ONLY translate surface2
    uv_grid = NULL,
    uv_common = list(n_u = 60, n_v = 60, grid_type = "disk", center = c(0.5,0.5), radius = 0.5),
    noise1_sd = 0,
    noise2_sd = 0,
    seed = 123
) {

  if (is.null(uv_grid)) uv_grid <- do.call(make_uv_grid, uv_common)

  f2_use <- translate_surface_function(f2, translate = translate2)

  data1 <- generate_surface_data(f1,     uv_grid = uv_grid, noise_sd = noise1_sd, seed = seed)
  data2 <- generate_surface_data(f2_use, uv_grid = uv_grid, noise_sd = noise2_sd, seed = seed)

  list(surface1 = data1, surface2 = data2, seed = seed, translate2=translate2)
}
