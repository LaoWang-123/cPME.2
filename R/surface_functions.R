######## This file is mostly used for the surface consturction for simulations in registration

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





################################################################################################################
##################################################################################################################

