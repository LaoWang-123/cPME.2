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
#' @param basis A basis object such as `fourier_basis`, `neumann_basis`,
#'   or `dirichlet_basis`, containing scalar/grad/rot/Jacobian generators,
#'   eigenvalue and norm functions, and minimally valid indices `p_min`, `q_min`.
#'
#' @return A named list of basis entries indexed by "p\_q".
#'
#' @examples
#' Bf <- build_basis_set(5, 5, fourier_basis)
#' Bd <- build_basis_set(5, 5, dirichlet_basis)
#'
#' @export
build_basis_set <- function(Pmax, Qmax, basis=c(fourier_basis,neumann_basis,dirichlet_basis))
{
  basis_set <- list()

  # Default to 0 if basis does not explicitly define p_min or q_min
  p_min <- if (!is.null(basis$p_min)) basis$p_min else 0
  q_min <- if (!is.null(basis$q_min)) basis$q_min else 0

  for (p in p_min:Pmax) {
    for (q in q_min:Qmax) {

      # Skip constant mode for cosine (Fourier/Neumann)
      # Dirichlet will never hit (0,0) because p_min=q_min=1
      if (p == 0 && q == 0) next

      key <- paste(p, q, sep = "_")

      basis_set[[key]] <- list(
        psi   = basis$scalar(p, q),
        grad_psi  = basis$grad(p, q),
        rot_grad_psi   = basis$rot(p, q),
        Dgrad_psi = basis$jac_g(p, q),
        Drot_grad_psi  = basis$jac_r(p, q),
        lambda_pq = basis$lambda(p, q),
        norm_pq   = basis$norm(p, q)
      )
    }
  }

  return(basis_set)
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
  bi_list
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
  D_bi_set
}

