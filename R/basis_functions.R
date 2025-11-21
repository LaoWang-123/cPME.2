#### Basis Construction
# ---------------- basis.R ----------------

# ---------------------------
# (1) Scalar basis ψ_pq(u,v)
# ---------------------------
scalar_basis_fn <- function(p, q) {
  cpq <- if (p == 0 && q == 0) 1 else if (p == 0 || q == 0) sqrt(2) else 2
  A <- 2 * pi * p
  B <- 2 * pi * q
  function(u, v) {
    cpq * cos(A * u) * cos(B * v)
  }
}

# ---------------------------
# (2) Gradient basis ∇ψ_pq(u,v)
# ---------------------------
grad_basis_fn <- function(p, q) {
  cpq <- if (p == 0 && q == 0) 1 else if (p == 0 || q == 0) sqrt(2) else 2
  A <- 2 * pi * p
  B <- 2 * pi * q
  function(u, v) {
    c(
      -cpq * A * sin(A * u) * cos(B * v),
      -cpq * B * cos(A * u) * sin(B * v)
    )
  }
}

# ---------------------------
# (3) Rotated gradient basis *∇ψ_pq(u,v)
# ---------------------------
rot_basis_fn <- function(p, q) {
  cpq <- if (p == 0 && q == 0) 1 else if (p == 0 || q == 0) sqrt(2) else 2
  A <- 2 * pi * p
  B <- 2 * pi * q
  function(u, v) {
    c(
      cpq * B * cos(A * u) * sin(B * v),
      -cpq * A * sin(A * u) * cos(B * v)
    )
  }
}

# ---------------------------
# (4) Jacobian of gradient basis D(∇ψ_pq)
# ---------------------------
jacobian_grad_fn <- function(p, q) {
  cpq <- if (p == 0 && q == 0) 1 else if (p == 0 || q == 0) sqrt(2) else 2
  A <- 2 * pi * p
  B <- 2 * pi * q
  function(u, v) {
    matrix(c(
      -A^2 * cpq * cos(A * u) * cos(B * v),
      A * B * cpq * sin(A * u) * sin(B * v),
      A * B * cpq * sin(A * u) * sin(B * v),
      -B^2 * cpq * cos(A * u) * cos(B * v)
    ), nrow = 2, byrow = TRUE)
  }
}

# ---------------------------
# (5) Jacobian of rotated basis D(*∇ψ_pq)
# ---------------------------
jacobian_rot_fn <- function(p, q) {
  cpq <- if (p == 0 && q == 0) 1 else if (p == 0 || q == 0) sqrt(2) else 2
  A <- 2 * pi * p
  B <- 2 * pi * q
  function(u, v) {
    matrix(c(
      -A * B * cpq * sin(A * u) * sin(B * v),
      B^2 * cpq * cos(A * u) * cos(B * v),
      -A^2 * cpq * cos(A * u) * cos(B * v),
      A * B * cpq * sin(A * u) * sin(B * v)
    ), nrow = 2, byrow = TRUE)
  }
}

# ---------------------------
# (6) Build all basis sets
# ---------------------------
#' Build all basis sets
#'
#' @param Pmax default 3
#' @param Qmax default 3
#'
#' @returns a list of basis functions
#' @export
#'
#' @examples
#' basis_set <- build_basis_set(3,3)
build_basis_set <- function(Pmax, Qmax) {
  basis_set <- list()
  for (p in 0:Pmax) {
    for (q in 0:Qmax) {
      # skip the constant mode (p=0, q=0)
      if (p == 0 && q == 0) next

      ### Calculate the norm ||\nambla psi_pq||_L2 = pi * cpq * \sqrt(p2+q2)
      cpq <- if (p == 0 && q == 0) 1 else if (p == 0 || q == 0) sqrt(2) else 2
      norm_pq <- pi * abs(cpq) * sqrt(p^2+q^2)

      key <- paste(p, q, sep = "_")
      basis_set[[key]] <- list(
        psi   = scalar_basis_fn(p, q),
        grad_psi  = grad_basis_fn(p, q),
        rot_grad_psi   = rot_basis_fn(p, q),
        Dgrad_psi = jacobian_grad_fn(p, q),
        Drot_grad_psi  = jacobian_rot_fn(p, q),
        norm_pq = norm_pq,
        lambda_pq = 4*pi^2*(p^2+q^2)
      )
    }
  }
  return(basis_set)
}





# Returns a named list of functions;
# Names are "p_q.grad" and "p_q.rot" (both included).
#' Title
#'
#' @param basis_set
#'
#' @returns a list of bi set (grad and normalized of scalar basis)
#' @export
#'
#' @examples
#' bi_set <- build_bi_set(basis_set)
build_bi_set <- function(basis_set) {
  keys <- sort(names(basis_set))

  grad_list <- lapply(keys, function(key) {
    bs <- basis_set[[key]]
    norm_pq <- bs$norm_pq
    function(u, v) bs$grad_psi(u, v) / norm_pq
  })

  rot_list <- lapply(keys, function(key) {
    bs <- basis_set[[key]]
    norm_pq <- bs$norm_pq
    function(u, v) bs$rot_grad_psi(u, v) / norm_pq
  })

  bi_list <- c(grad_list, rot_list)
  names(bi_list) <- c(paste0(keys, ".grad"), paste0(keys, ".rot"))
  bi_list
}


# Returns a named list of functions Dbi, (u,v) to 2x2 matrix
# Names are "p_q.grad" and "p_q.rot" (both included).
#' Title
#'
#' @param basis_set
#'
#' @returns a list of D bi sets
#' @export
#'
#' @examples
#' D_bi_set <- build_D_bi_set(basis_set)
build_D_bi_set <- function(basis_set) {
  keys <- sort(names(basis_set))

  D_grad_list <- lapply(keys, function(key) {
    bs <- basis_set[[key]]
    norm_pq <- bs$norm_pq
    function(u, v) bs$Dgrad_psi(u, v) / norm_pq
  })

  D_rot_list <- lapply(keys, function(key) {
    bs <- basis_set[[key]]
    norm_pq <- bs$norm_pq
    function(u, v) bs$Drot_grad_psi(u, v) / norm_pq
  })


  D_bi_set <- c(D_grad_list, D_rot_list)
  names(D_bi_set) <- c(paste0(keys, ".grad"), paste0(keys, ".rot"))

  D_bi_set
}

