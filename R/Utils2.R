### new Utils.R

### We optimize some functions and workflow to accelerate the computation speed

build_basis_grid <- function(basis_set, Ugrid) {
  out <- list()
  n  <- nrow(Ugrid)

  for (key in names(basis_set)) {
    bs <- basis_set[[key]]

    # psi(u,v)
    psi_vals <- vapply(
      1:n,
      function(i) bs$psi(Ugrid[i,1], Ugrid[i,2]),
      numeric(1)
    )

    # grad psi(u,v)
    grad_vals <- t(vapply(
      1:n,
      function(i) bs$grad_psi(Ugrid[i,1], Ugrid[i,2]),
      numeric(2)
    ))
    grad_psi_x <- grad_vals[,1]
    grad_psi_y <- grad_vals[,2]

    # rot grad psi(u,v)
    rot_vals <- t(vapply(
      1:n,
      function(i) bs$rot_grad_psi(Ugrid[i,1], Ugrid[i,2]),
      numeric(2)
    ))
    rot_x <- rot_vals[,1]
    rot_y <- rot_vals[,2]

    out[[key]] <- list(
      psi        = psi_vals,     # n
      grad_psi_x = grad_psi_x,   # n
      grad_psi_y = grad_psi_y,   # n
      rot_x      = rot_x,        # n
      rot_y      = rot_y,        # n
      lambda     = bs$lambda_pq,
      norm       = bs$norm_pq
    )
  }

  out
}



compute_dphi_grid <- function(basis_grid, f2_grid, grad_f2_grid) {
  # f2_grid: n × 3
  # grad_f2_grid: list(Gx = n×3, Gy = n×3), gradient of f2 for each coord

  out <- list()

  for (key in names(basis_grid)) {
    bg <- basis_grid[[key]]

    psi  <- bg$psi
    gx   <- bg$grad_psi_x
    gy   <- bg$grad_psi_y
    rx   <- bg$rot_x
    ry   <- bg$rot_y
    lam  <- bg$lambda
    norm <- bg$norm

    ## ---- gradient-type (Eq. 11) ----
    term1 <- (-0.5 * lam) * psi * f2_grid                 # n×3
    term2 <- cbind(
      gx * grad_f2_grid$Gx[,1] + gy * grad_f2_grid$Gy[,1],
      gx * grad_f2_grid$Gx[,2] + gy * grad_f2_grid$Gy[,2],
      gx * grad_f2_grid$Gx[,3] + gy * grad_f2_grid$Gy[,3]
    )

    dphi_grad <- (term1 + term2) / norm   # n×3

    ## ---- rotated-type (Eq. 12) ----
    term_rot <- cbind(
      rx * grad_f2_grid$Gx[,1] + ry * grad_f2_grid$Gy[,1],
      rx * grad_f2_grid$Gx[,2] + ry * grad_f2_grid$Gy[,2],
      rx * grad_f2_grid$Gx[,3] + ry * grad_f2_grid$Gy[,3]
    )

    dphi_rot <- term_rot / norm          # n×3

    out[[paste0(key,".grad")]] <- dphi_grad
    out[[paste0(key,".rot")]]  <- dphi_rot
  }
  out
}


compute_inner_products_fast <- function(diff_grid, dphi_grid_list, weight) {

  vapply(dphi_grid_list, function(mat3) {
    sum(rowSums(diff_grid * mat3) * weight)
  }, numeric(1))
}


make_f2_grid <- function(f2_fun, Ugrid) {
  n <- nrow(Ugrid)

  # 每个点调用 f2_fun(u,v)，返回长度 3 vector
  mat <- t(vapply(
    1:n,
    function(i) f2_fun(Ugrid[i,1], Ugrid[i,2]),
    numeric(3)
  ))

  # n × 3
  colnames(mat) <- c("x","y","z")
  mat
}

make_grad_f2_grid <- function(grad_f2_fun, Ugrid) {
  n <- nrow(Ugrid)

  # 1. 对所有 grid 点生成一个 list，每个元素是一个 3×2 matrix
  G_list <- lapply(1:n, function(i) grad_f2_fun(Ugrid[i,1], Ugrid[i,2]))

  # 2. 安全检查（只检查第一个）
  G0 <- G_list[[1]]
  if (!is.matrix(G0) || !all(dim(G0) == c(3,2))) {
    stop("grad_f2_fun must return a 3×2 matrix.")
  }

  # 3. 将 list of 3×2 转成 3×2×n 数组
  G_array <- simplify2array(G_list)
  # 现在 G_array 的维度是：3 × 2 × n

  # 4. 提取两列偏导
  # G_array[,1,] 是 3 × n，转置后 n × 3
  Gx <- t(G_array[,1,])
  Gy <- t(G_array[,2,])

  list(Gx = Gx, Gy = Gy)
}

