### basis_construction.R

#' Build a Set of Basis Functions
#'
#' Constructs scalar, gradient, rotated gradient, and Jacobian basis
#' functions on the tensor grid \eqn{0 \le p \le Pmax}, \eqn{0 \le q \le Qmax},
#' respecting the index constraints defined in the basis object
#' (e.g. Dirichlet requires \eqn{p,q \ge 1}).
#'
#' @param Pmax Maximum frequency index in the first dimension.
#' @param Qmax Maximum frequency index in the second dimension.
#' @param include_constant Logical; include the constant basis term when
#'   supported by the selected basis family.
#' @param basis A basis object such as `fourier_basis`, `neumann_basis`,
#'   or `dirichlet_basis`, containing scalar/grad/rot/Jacobian generators,
#'   eigenvalue and norm functions, and minimally valid indices `p_min`, `q_min`.
#'
#' @return A named list of basis entries indexed by "p\_q".
#'
#' @examples
#' Bf <- build_basis_set(5, 5, basis = fourier_basis)
#' Bd <- build_basis_set(5, 5, basis = dirichlet_basis)
#'
#' @export
build_basis_set <- function(Pmax, Qmax,include_constant = FALSE, basis=c(fourier_basis,neumann_basis,dirichlet_basis))
{
  basis_set <- list()

  # Default to 0 if basis does not explicitly define p_min or q_min
  p_min <- if (!is.null(basis$p_min)) basis$p_min else 0
  q_min <- if (!is.null(basis$q_min)) basis$q_min else 0

  for (p in p_min:Pmax) {
    for (q in q_min:Qmax) {

      # Skip constant mode for cosine (Fourier/Neumann)
      # Dirichlet will never hit (0,0) because p_min=q_min=1
      # default: skip (0,0), same as old behavior
      # for UQ kernel: set include_constant = TRUE
      if (!include_constant && p == 0 && q == 0) next

      key <- paste(p, q, sep = "_")

      basis_set[[key]] <- list(
        psi   = basis$scalar(p, q),
        grad_psi  = basis$grad(p, q),
        rot_grad_psi   = basis$rot(p, q),
        Dgrad_psi = basis$jac_g(p, q),
        Drot_grad_psi  = basis$jac_r(p, q),
        lambda_pq = basis$lambda(p, q),
        norm_pq   = basis$norm(p, q) # L2 norm of the gradient field ||∇ψ_pq||_L2, used for normalized tangent vector construction.

      )
    }
  }

  return(basis_set)
}



#' Evaluate Basis Functions on a Grid
#'
#' Evaluates a pre-constructed basis system on a set of parameter-domain
#' locations. The function supports three modes. The `"full"` mode evaluates
#' scalar, gradient, and rotated-gradient basis functions. The `"div_free"`
#' mode evaluates scalar and rotated-gradient basis functions for
#' divergence-free reparameterization updates. The `"scalar"` mode evaluates
#' only scalar basis functions and is intended for uncertainty quantification,
#' spectral kernel construction, and pointwise posterior variance estimation.
#'
#' @param basis_set A named list of basis functions returned by
#'   \code{build_basis_set()}.
#' @param Ugrid A two-column matrix or data frame of parameter locations.
#'   Each row is one point \eqn{(u,v)} in the parameter domain.
#' @param mode Character string specifying which quantities to evaluate.
#'   Options are \code{"full"}, \code{"div_free"}, and \code{"scalar"}.
#'   The default is \code{"full"}.
#'
#' @return A named list indexed by basis labels. Each entry contains the
#'   evaluated scalar basis values \code{psi}, the Laplacian eigenvalue
#'   \code{lambda}, and the normalization constant \code{norm}. In
#'   \code{"full"} mode, gradient and rotated-gradient components are also
#'   returned. In \code{"div_free"} mode, only rotated-gradient components
#'   are returned in addition to scalar quantities. In \code{"scalar"} mode,
#'   no vector-field quantities are evaluated.
#'
#' @export
build_basis_grid <- function(
    basis_set,
    Ugrid,
    mode = c("full", "div_free", "scalar")
) {
  mode <- match.arg(mode)
  Ugrid <- as.matrix(Ugrid)

  out <- list()
  n <- nrow(Ugrid)
  keys <- sort(names(basis_set))

  for (key in keys) {
    bs <- basis_set[[key]]

    # Scalar basis values are needed in all modes.
    psi_vals <- vapply(
      seq_len(n),
      function(i) bs$psi(Ugrid[i, 1], Ugrid[i, 2]),
      numeric(1)
    )

    # UQ only needs scalar basis values and eigenvalues.
    if (mode == "scalar") {
      out[[key]] <- list(
        psi = psi_vals,
        lambda = bs$lambda_pq,
        norm = bs$norm_pq
      )
      next
    }

    # Rotated-gradient basis is needed for divergence-free updates.
    rot_vals <- t(vapply(
      seq_len(n),
      function(i) bs$rot_grad_psi(Ugrid[i, 1], Ugrid[i, 2]),
      numeric(2)
    ))

    if (mode == "full") {
      grad_vals <- t(vapply(
        seq_len(n),
        function(i) bs$grad_psi(Ugrid[i, 1], Ugrid[i, 2]),
        numeric(2)
      ))

      grad_psi_x <- grad_vals[, 1]
      grad_psi_y <- grad_vals[, 2]
    } else {
      grad_psi_x <- NULL
      grad_psi_y <- NULL
    }

    out[[key]] <- list(
      psi = psi_vals,
      grad_psi_x = grad_psi_x,
      grad_psi_y = grad_psi_y,
      rot_x = rot_vals[, 1],
      rot_y = rot_vals[, 2],
      lambda = bs$lambda_pq,
      norm = bs$norm_pq  # Gradient-field normalization constant.
      # For normalized scalar bases, this is not ||ψ_pq||, but ||∇ψ_pq||_L2.

    )
  }

  out
}


# Returns a named list of functions;
# Names are "p_q.grad" and "p_q.rot" (both included).

#' Build bi set from basis_set
#'
#' @param basis_set A list of basis functions.
#' @param mode c("full","div_free")
#'
#' @returns a list of bi, normailized of grad and grad rot basis functions
#' @export
#'
#' @examples
#' basis_set <- build_basis_set(2, 2, basis = neumann_basis)
#' bi_set <- build_bi_set(basis_set)
build_bi_set <- function(basis_set, mode = c("full", "div_free")) {
  mode <- match.arg(mode)
  keys <- sort(names(basis_set))

  rot_list <- lapply(keys, function(key) {
    bs <- basis_set[[key]]
    norm_pq <- bs$norm_pq
    function(u, v) bs$rot_grad_psi(u, v) / norm_pq
  })
  names(rot_list) <- paste0(keys, ".rot")

  if (mode == "div_free") return(rot_list)

  grad_list <- lapply(keys, function(key) {
    bs <- basis_set[[key]]
    norm_pq <- bs$norm_pq
    function(u, v) bs$grad_psi(u, v) / norm_pq
  })
  names(grad_list) <- paste0(keys, ".grad")

  bi_list <- c(grad_list, rot_list)
  return(bi_list)
}



# Returns a named list of functions Dbi, (u,v) to 2x2 matrix
# Names are "p_q.grad" and "p_q.rot" (both included).

#' Build D_bi set from basis_set
#'
#' @param basis_set A list of basis functions.
#' @param mode c("full","div_free")
#'
#' @returns a list of D bi sets
#' @export
#'
#' @examples
#' basis_set <- build_basis_set(2, 2, basis = neumann_basis)
#' D_bi_set <- build_D_bi_set(basis_set)
build_D_bi_set <- function(basis_set, mode = c("full", "div_free")) {
  mode <- match.arg(mode)
  keys <- sort(names(basis_set))

  D_rot_list <- lapply(keys, function(key) {
    bs <- basis_set[[key]]
    norm_pq <- bs$norm_pq
    function(u, v) bs$Drot_grad_psi(u, v) / norm_pq
  })
  names(D_rot_list) <- paste0(keys, ".rot")

  if (mode == "div_free") return(D_rot_list)

  D_grad_list <- lapply(keys, function(key) {
    bs <- basis_set[[key]]
    norm_pq <- bs$norm_pq
    function(u, v) bs$Dgrad_psi(u, v) / norm_pq
  })
  names(D_grad_list) <- paste0(keys, ".grad")

  D_bi_set <- c(D_grad_list, D_rot_list)
  return(D_bi_set)
}


######## Basis Formula

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


