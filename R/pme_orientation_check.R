### This file is used to check PME function orientation
# =========================
# Initialization diagnostics for PME fitting
# =========================

# =========================
# Joint-PCA corner order check for two initializations
# init structure: $centers (n x 3), $parameterization (n x 2)
# =========================

# --- min-max normalize to [0,1] ---
.range01 <- function(x) {
  r <- range(x, na.rm = TRUE)
  if (!is.finite(r[1]) || !is.finite(r[2]) || r[1] == r[2]) return(rep(0.5, length(x)))
  (x - r[1]) / (r[2] - r[1])
}

.normalize_uv <- function(UV) {
  UV <- as.matrix(UV)
  if (ncol(UV) != 2) stop("UV must be n x 2")
  cbind(u = .range01(UV[,1]), v = .range01(UV[,2]))
}

# --- pick corners by nearest to canonical corners ---
.pick_uv_corners_nearest <- function(UV01, XYZ) {
  corners <- rbind(c(0,0), c(1,0), c(1,1), c(0,1))
  labels  <- c("00","10","11","01")
  idx <- integer(4)
  for (k in 1:4) {
    d2 <- (UV01[,1] - corners[k,1])^2 + (UV01[,2] - corners[k,2])^2
    idx[k] <- which.min(d2)
  }
  XYZ <- as.matrix(XYZ)
  list(
    labels = labels,
    idx = idx,
    UV4 = UV01[idx, , drop=FALSE],
    XYZ4 = XYZ[idx, , drop=FALSE]
  )
}

# --- same up to cyclic rotation ---
.same_up_to_rotation <- function(a, b) {
  if (length(a) != length(b)) return(FALSE)
  n <- length(a)
  for (s in 0:(n-1)) {
    if (all(a == c(b[(s+1):n], b[1:s]))) return(TRUE)
  }
  FALSE
}

# --- compute cyclic order by angles around centroid in 2D ---
.cyclic_order_2d <- function(P2) {
  ctr <- colMeans(P2)
  ang <- atan2(P2[,2] - ctr[2], P2[,1] - ctr[1])  # (-pi, pi]
  order(ang)  # CCW
}

# --- build a shared 2D PCA basis from pooled XYZ ---
.joint_pca_basis <- function(XYZ_pool) {
  X <- as.matrix(XYZ_pool)
  Xc <- scale(X, center = TRUE, scale = FALSE)
  pc <- stats::prcomp(Xc)
  # shared basis in original 3D coords:
  B <- pc$rotation[,1:2, drop=FALSE]  # 3x2
  center <- attr(Xc, "scaled:center")
  list(B = B, center = center)
}

# --- project XYZ to shared PCA plane ---
.project_shared <- function(XYZ, basis) {
  X <- sweep(as.matrix(XYZ), 2, basis$center, FUN = "-")
  X %*% basis$B
}

# =========================
# Main function: joint PCA
# =========================
#' Use pca to check pme initialization orientation
#'
#' @param init1 initialization for pme1
#' @param init2 initialization for pme2
#' @param pca_source c("all_centers", "eight_corners")
#' @param verbose default TRUE
#'
#' @returns a list of whether rot true
#' @export
#'
check_pme_orientation <- function(init1, init2,
                                    pca_source = c("all_centers", "eight_corners"),
                                    verbose = TRUE) {
  pca_source <- match.arg(pca_source)

  # Extract
  XYZ1 <- as.matrix(init1$centers); UV1 <- as.matrix(init1$parameterization)
  XYZ2 <- as.matrix(init2$centers); UV2 <- as.matrix(init2$parameterization)
  if (ncol(XYZ1) != 3 || ncol(XYZ2) != 3) stop("$centers must be n x 3")
  if (ncol(UV1)  != 2 || ncol(UV2)  != 2) stop("$parameterization must be n x 2")

  # Normalize UVs separately (isomap -> [0,1]^2)
  UV101 <- .normalize_uv(UV1)
  UV201 <- .normalize_uv(UV2)

  # Pick 4 labeled corners for each
  c1 <- .pick_uv_corners_nearest(UV101, XYZ1)
  c2 <- .pick_uv_corners_nearest(UV201, XYZ2)

  # Build joint PCA basis
  XYZ_pool <- if (pca_source == "all_centers") {
    rbind(XYZ1, XYZ2)
  } else {
    rbind(c1$XYZ4, c2$XYZ4)
  }
  basis <- .joint_pca_basis(XYZ_pool)

  # Project corners using shared basis
  P1 <- .project_shared(c1$XYZ4, basis)  # 4x2
  P2 <- .project_shared(c2$XYZ4, basis)

  # Compute cyclic order (permutation of 1:4) in fixed label order (00,10,11,01)
  o1 <- .cyclic_order_2d(P1)
  o2 <- .cyclic_order_2d(P2)

  # Comparisons
  strict_same <- identical(o1, o2)             # no rotation allowed (label-locked order)
  rot_same    <- .same_up_to_rotation(o1, o2)  # allow cyclic shift (diagnose rotation mismatch)
  mirror_same <- .same_up_to_rotation(o1, rev(o2))

  # A simple "final" label
  final <- if (strict_same) {
    "same_strict"
  } else if (rot_same) {
    "same_up_to_rotation"
  } else if (mirror_same) {
    "mirror_reversed"
  } else {
    "different"
  }

  out <- list(
    final = final,
    order1 = o1,
    order2 = o2,
    corners1_UV01 = c1$UV4,
    corners2_UV01 = c2$UV4,
    corners1_XYZ = c1$XYZ4,
    corners2_XYZ = c2$XYZ4,
    corners1_2D = P1,
    corners2_2D = P2,
    pca_source = pca_source
  )

  if (isTRUE(verbose)) {
    message(sprintf(
      "Joint-PCA corner check: final=%s",out$final)
    )
  }

  out
}
