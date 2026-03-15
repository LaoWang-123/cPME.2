cpq_cos <- function(p, q) {
  if (p == 0 && q == 0) return(1)
  if (p == 0 || q == 0) return(sqrt(2))
  return(2)
}

cpq_sin <- function(p, q) {
  return(2)  # Dirichlet sin–sin normalization
}


#' Fourier Laplacian Basis
#'
#' A collection of functions defining the periodic Fourier eigenbasis,
#' including scalar, gradient, rotation, and associated Jacobians.
#'
#' @details
#' This object is a named list containing function generators:
#' \describe{
#'   \item{scalar}{Generate the normalized scalar basis \eqn{\psi_{pq}}.}
#'   \item{grad}{Gradient field \eqn{\nabla\psi_{pq}}.}
#'   \item{rot}{Divergence-free rotated field.}
#'   \item{jac_g}{Jacobian of the gradient field.}
#'   \item{jac_r}{Jacobian of the rotated field.}
#'   \item{lambda}{Eigenvalue \eqn{\lambda_{pq}}.}
#'   \item{norm}{Gradient L2 norm.}
#' }
#'
#' @format A named list of 7 functions.
#' @examples
#' basis <- fourier_basis
#' psi23 <- basis$scalar(2,3)
#' psi23(0.5,0.2)
#'
#' @export
fourier_basis <- list(

  scalar = function(p, q) {
    cpq <- cpq_cos(p, q)
    A <- 2*pi*p; B <- 2*pi*q
    function(u, v) cpq * cos(A*u) * cos(B*v)
  },

  grad = function(p, q) {
    cpq <- cpq_cos(p, q)
    A <- 2*pi*p; B <- 2*pi*q
    function(u, v) c(
      -cpq*A*sin(A*u)*cos(B*v),
      -cpq*B*cos(A*u)*sin(B*v)
    )
  },

  rot = function(p, q) {
    cpq <- cpq_cos(p, q)
    A <- 2*pi*p; B <- 2*pi*q
    function(u, v) c(
      cpq*B*cos(A*u)*sin(B*v),
      -cpq*A*sin(A*u)*cos(B*v)
    )
  },

  jac_g = function(p, q) {
    cpq <- cpq_cos(p, q)
    A <- 2*pi*p; B <- 2*pi*q
    function(u, v) matrix(c(
      -A^2*cpq*cos(A*u)*cos(B*v),
      A*B*cpq*sin(A*u)*sin(B*v),
      A*B*cpq*sin(A*u)*sin(B*v),
      -B^2*cpq*cos(A*u)*cos(B*v)
    ), nrow=2, byrow=TRUE)
  },

  jac_r = function(p, q) {
    cpq <- cpq_cos(p, q)
    A <- 2*pi*p; B <- 2*pi*q
    function(u, v) matrix(c(
      -A*B*cpq*sin(A*u)*sin(B*v),
      B^2*cpq*cos(A*u)*cos(B*v),
      -A^2*cpq*cos(A*u)*cos(B*v),
      A*B*cpq*sin(A*u)*sin(B*v)
    ), nrow=2, byrow=TRUE)
  },

  lambda = function(p, q) 4*pi^2*(p^2 + q^2),
  norm   = function(p, q) 2*pi*sqrt(p^2 + q^2),
  p_min=0,
  q_min=0
)

#' Neumann Laplacian Basis
#'
#' A collection of functions implementing the Neumann Laplacian eigenbasis
#' on the unit square \eqn{U=[0,1]^2}. Each basis element is constructed from
#' the normalized cosine modes
#' \deqn{\psi_{pq}(u,v)=c_{pq}\cos(p\pi u)\cos(q\pi v)}
#' which satisfy \eqn{\partial_n\psi_{pq}=0} on the boundary.
#'
#' The basis supports generation of scalar fields, gradient fields,
#' divergence-free rotated fields, and their Jacobians, along with the
#' eigenvalue \eqn{\lambda_{pq}} and gradient norm \eqn{\|\nabla\psi_{pq}\|}.
#'
#' @details
#' This object is a named list containing the following components:
#' \describe{
#'   \item{scalar}{Function generating the normalized scalar eigenfunction \eqn{\psi_{pq}}.}
#'   \item{grad}{Function generating the gradient field \eqn{\nabla\psi_{pq}}.}
#'   \item{rot}{Function generating the rotated (divergence-free) field \eqn{\ast\nabla\psi_{pq}}.}
#'   \item{jac_g}{Jacobian of the gradient field.}
#'   \item{jac_r}{Jacobian of the rotated field.}
#'   \item{lambda}{Laplacian eigenvalue \eqn{\lambda_{pq}=\pi^2(p^2+q^2)}.}
#'   \item{norm}{L2 norm of the gradient field, equal to \eqn{\pi\sqrt{p^2+q^2}}.}
#'   \item{p_min}{Minimum valid index for \eqn{p}. For Neumann, \eqn{p_{\min}=0}.}
#'   \item{q_min}{Minimum valid index for \eqn{q}. For Neumann, \eqn{q_{\min}=0}.}
#' }
#'
#' @format A named list of 9 functions/values.
#'
#' @examples
#' basis <- neumann_basis
#' psi <- basis$scalar(1,2)
#' psi(0.4, 0.7)
#'
#' grad12 <- basis$grad(1,2)
#' grad12(0.3, 0.6)
#'
#' @export
neumann_basis <- list(

  scalar = function(p, q) {
    cpq <- cpq_cos(p, q)
    A <- pi*p; B <- pi*q
    function(u, v) cpq * cos(A*u) * cos(B*v)
  },

  grad = function(p, q) {
    cpq <- cpq_cos(p, q)
    A <- pi*p; B <- pi*q
    function(u, v) c(
      -cpq*A*sin(A*u)*cos(B*v),
      -cpq*B*cos(A*u)*sin(B*v)
    )
  },

  rot = function(p, q) {
    cpq <- cpq_cos(p, q)
    A <- pi*p; B <- pi*q
    function(u, v) c(
      cpq*B*cos(A*u)*sin(B*v),
      -cpq*A*sin(A*u)*cos(B*v)
    )
  },

  jac_g = function(p, q) {
    cpq <- cpq_cos(p, q)
    A <- pi*p; B <- pi*q
    function(u, v) matrix(c(
      -A^2*cpq*cos(A*u)*cos(B*v),
      A*B*cpq*sin(A*u)*sin(B*v),
      A*B*cpq*sin(A*u)*sin(B*v),
      -B^2*cpq*cos(A*u)*cos(B*v)
    ), nrow=2, byrow=TRUE)
  },

  jac_r = function(p, q) {
    cpq <- cpq_cos(p, q)
    A <- pi*p; B <- pi*q
    function(u, v) matrix(c(
      -A*B*cpq*sin(A*u)*sin(B*v),
      B^2*cpq*cos(A*u)*cos(B*v),
      -A^2*cpq*cos(A*u)*cos(B*v),
      A*B*cpq*sin(A*u)*sin(B*v)
    ), nrow=2, byrow=TRUE)
  },

  lambda = function(p, q) pi^2*(p^2 + q^2),
  norm   = function(p, q) pi*sqrt(p^2 + q^2),
  p_min=0,
  q_min=0
)

#' Dirichlet Laplacian Basis
#'
#' A collection of functions implementing the Dirichlet Laplacian eigenbasis
#' on the unit square \eqn{U=[0,1]^2}. Each basis element is constructed from
#' the normalized sine modes
#' \deqn{\psi_{pq}(u,v)=2\sin(p\pi u)\sin(q\pi v)}
#' which satisfy \eqn{\psi_{pq}=0} on the boundary.
#'
#' Unlike Neumann and Fourier bases, Dirichlet modes require
#' \eqn{p\ge 1,\ q\ge 1}. No constant or axis-aligned modes exist.
#'
#' The basis supports generation of scalar fields, gradient fields,
#' divergence-free rotated fields, their Jacobians, the eigenvalue
#' \eqn{\lambda_{pq}}, and gradient norm \eqn{\|\nabla\psi_{pq}\|}.
#'
#' @details
#' This object is a named list containing:
#' \describe{
#'   \item{scalar}{Function generating the normalized scalar eigenfunction \eqn{\psi_{pq}}.}
#'   \item{grad}{Function generating the gradient field \eqn{\nabla\psi_{pq}}.}
#'   \item{rot}{Function generating the rotated (divergence-free) field.}
#'   \item{jac_g}{Jacobian of the gradient field.}
#'   \item{jac_r}{Jacobian of the rotated field.}
#'   \item{lambda}{Laplacian eigenvalue \eqn{\lambda_{pq}=\pi^2(p^2+q^2)}.}
#'   \item{norm}{L2 norm of the gradient field \eqn{\pi\sqrt{p^2+q^2}}.}
#'   \item{p_min}{Minimum valid index \eqn{p_{\min}=1}.}
#'   \item{q_min}{Minimum valid index \eqn{q_{\min}=1}.}
#' }
#'
#' @format A named list of 9 functions/values.
#'
#' @examples
#' basis <- dirichlet_basis
#' psi <- basis$scalar(2,3)
#' psi(0.2, 0.5)
#'
#' grad23 <- basis$grad(2,3)
#' grad23(0.6, 0.1)
#'
#' @export
dirichlet_basis <- list(

  scalar = function(p, q) {
    cpq <- cpq_sin(p, q)  # always = 2
    A <- pi*p; B <- pi*q
    function(u, v) cpq * sin(A*u) * sin(B*v)
  },

  grad = function(p, q) {
    cpq <- cpq_sin(p, q)
    A <- pi*p; B <- pi*q
    function(u, v) c(
      cpq*A*cos(A*u)*sin(B*v),
      cpq*B*sin(A*u)*cos(B*v)
    )
  },

  rot = function(p, q) {
    cpq <- cpq_sin(p, q)
    A <- pi*p; B <- pi*q
    function(u, v) c(
      -cpq*B*sin(A*u)*cos(B*v),
      cpq*A*cos(A*u)*sin(B*v)
    )
  },

  jac_g = function(p, q) {
    cpq <- cpq_sin(p, q)
    A <- pi*p; B <- pi*q
    function(u, v) matrix(c(
      -A^2*cpq*sin(A*u)*sin(B*v),
      A*B*cpq*cos(A*u)*cos(B*v),
      A*B*cpq*cos(A*u)*cos(B*v),
      -B^2*cpq*sin(A*u)*sin(B*v)
    ), nrow=2, byrow=TRUE)
  },

  jac_r = function(p, q) {
    cpq <- cpq_sin(p, q)
    A <- pi*p; B <- pi*q
    function(u, v) matrix(c(
      -A*B*cpq*cos(A*u)*cos(B*v),
      B^2*cpq*sin(A*u)*sin(B*v),
      -A^2*cpq*sin(A*u)*sin(B*v),
      A*B*cpq*cos(A*u)*cos(B*v)
    ), nrow=2, byrow=TRUE)
  },

  lambda = function(p, q) pi^2*(p^2 + q^2),
  norm   = function(p, q) pi*sqrt(p^2 + q^2),
  p_min=1,
  q_min=1
)

