# ==========================================================
# γ function factory —  "smooth" / "rotate" / "identity"
# ==========================================================
#' Gamma function factory nx2 matrix input
#'
#' @param mode mode choice
#' @param au default 0.3
#' @param av default 0.25
#' @param angle default 0.5*pi
#'
#' @returns a chosen gamma function with input with nx2 matrix
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

#' Gamma function factory u,v input
#'
#' @param mode mode choice
#' @param au default 0.3
#' @param av default 0.25
#' @param angle default 0.5*pi
#'
#' @returns a chosen gamma function with input with u,v
#' @export
#'
gamma_function_factory_scalar <- function(
    mode = c("smooth", "rotate", "identity"),
    au = 0.3, av = 0.25, angle = 0.5*pi
) {
  mode <- match.arg(mode)
  .au <- au; .av <- av; .angle <- angle
  force(.au); force(.av); force(.angle)

  if (mode == "smooth") {

    fun <- function(u, v, au = .au, av = .av) {
      bu <- u*(1-u)
      bv <- v*(1-v)
      du <- au * bu * sin(2*pi*v)
      dv <- av * bv * sin(2*pi*u)
      u2 <- pmin(pmax(u + du, 0), 1)
      v2 <- pmin(pmax(v + dv, 0), 1)
      c(u2, v2)
    }

  } else if (mode == "rotate") {

    fun <- function(u, v, angle = .angle) {
      uc <- u - 0.5
      vc <- v - 0.5
      cosA <- cos(angle); sinA <- sin(angle)
      u2 <- cosA*uc - sinA*vc + 0.5
      v2 <- sinA*uc + cosA*vc + 0.5
      # clamp to [0,1]
      u2 <- pmin(pmax(u2, 0), 1)
      v2 <- pmin(pmax(v2, 0), 1)
      c(u2, v2)
    }

  } else if (mode == "identity") {

    fun <- function(u, v) {
      c(u, v)
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
#' Generate simulated surface data
#'
#' @param f surface function, (u,v) input and (x,y,z) output
#' @param uv_grid a data.frame with u,v columns
#' @param n_u default = 50
#' @param n_v default = 50
#' @param noise_sd default = 0
#' @param seed default = NULL
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
                                  noise_sd = 0,
                                  seed = NULL) {
  # Optional reproducibility
  if (!is.null(seed)) set.seed(seed)

  if(is.null(uv_grid)){
    # Generate UV grid
    u_seq <- seq(0, 1, length.out = n_u)
    v_seq <- seq(0, 1, length.out = n_v)
    uv_grid <- expand.grid(u = u_seq, v = v_seq)
  }

  # Evaluate surface function f(u, v)
  xyz <- t(apply(uv_grid, 1, function(uv) {
    f(uv[1], uv[2])
  }))
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

