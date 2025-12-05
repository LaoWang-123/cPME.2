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


