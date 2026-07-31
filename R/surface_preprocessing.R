# Surface preprocessing helpers for labeled 3D arrays.
#
# These functions are intentionally generic: they operate on arrays or point
# tables and do not read NIfTI files directly. Cluster workflows can read image
# data however they prefer, then use these helpers before PME registration.

#' Extract Surface Voxels From a 3D Label Array
#'
#' A voxel is treated as a surface voxel when it belongs to the target ROI and
#' at least one selected neighboring voxel is background. Use 6-connectivity for
#' face neighbors or 26-connectivity for face, edge, and corner neighbors.
#'
#' @param label_space A 3D array containing ROI labels or a binary mask.
#' @param label_values Optional vector of labels to keep. If `NULL`, all values
#'   not equal to `background` are treated as ROI.
#' @param background Background value. Defaults to `0`.
#' @param connectivity Neighbor definition: `6` for face neighbors or `26` for
#'   all face/edge/corner neighbors.
#' @param output Either `"coords"` for a data frame of surface coordinates or
#'   `"mask"` for a logical 3D surface mask.
#' @param voxel_spacing Optional numeric length-3 vector. When provided with
#'   `output = "coords"`, `x_mm`, `y_mm`, and `z_mm` columns are added.
#' @param separate_labels If `TRUE`, surface voxels are extracted for each label
#'   separately and then combined. This prevents adjacent labels from being
#'   treated as one merged object.
#'
#' @return A data frame of surface voxel coordinates or a logical 3D mask.
#' @export
extract_surface_voxels <- function(label_space,
                                   label_values = NULL,
                                   background = 0,
                                   connectivity = 6,
                                   output = c("coords", "mask"),
                                   voxel_spacing = NULL,
                                   separate_labels = FALSE) {
  output <- match.arg(output)
  connectivity <- as.integer(connectivity)

  if (length(dim(label_space)) != 3) {
    stop("label_space must be a 3D array.")
  }
  if (!connectivity %in% c(6L, 26L)) {
    stop("connectivity must be either 6 or 26.")
  }
  if (!is.null(voxel_spacing) && length(voxel_spacing) != 3) {
    stop("voxel_spacing must be NULL or a numeric vector of length 3.")
  }
  if (separate_labels && is.null(label_values)) {
    stop("separate_labels = TRUE requires label_values.")
  }

  if (separate_labels) {
    surface_list <- lapply(label_values, function(label_value) {
      extract_surface_voxels(
        label_space = label_space,
        label_values = label_value,
        background = background,
        connectivity = connectivity,
        output = "coords",
        voxel_spacing = voxel_spacing,
        separate_labels = FALSE
      )
    })

    surface_df <- do.call(rbind, surface_list)
    rownames(surface_df) <- NULL

    if (output == "mask") {
      surface_mask <- array(FALSE, dim = dim(label_space))
      surface_mask[as.matrix(surface_df[, c("x", "y", "z"), drop = FALSE])] <- TRUE
      return(surface_mask)
    }

    return(surface_df)
  }

  label_space <- as.array(label_space)
  roi_mask <- if (is.null(label_values)) {
    label_space != background
  } else {
    label_space %in% label_values
  }
  dim(roi_mask) <- dim(label_space)

  if (!any(roi_mask)) {
    if (output == "mask") {
      return(array(FALSE, dim = dim(label_space)))
    }
    empty <- data.frame(
      x = integer(),
      y = integer(),
      z = integer(),
      label_value = numeric()
    )
    if (!is.null(voxel_spacing)) {
      empty$x_mm <- numeric()
      empty$y_mm <- numeric()
      empty$z_mm <- numeric()
    }
    return(empty)
  }

  offsets <- as.matrix(expand.grid(dx = -1:1, dy = -1:1, dz = -1:1))
  offsets <- offsets[rowSums(abs(offsets)) > 0, , drop = FALSE]
  if (connectivity == 6L) {
    offsets <- offsets[rowSums(abs(offsets)) == 1, , drop = FALSE]
  }

  dims <- dim(roi_mask)
  padded_mask <- array(FALSE, dim = dims + 2L)
  padded_mask[
    2:(dims[1] + 1L),
    2:(dims[2] + 1L),
    2:(dims[3] + 1L)
  ] <- roi_mask

  has_background_neighbor <- array(FALSE, dim = dims)
  for (idx in seq_len(nrow(offsets))) {
    dx <- offsets[idx, 1]
    dy <- offsets[idx, 2]
    dz <- offsets[idx, 3]

    neighbor_mask <- padded_mask[
      (2 + dx):(dims[1] + 1L + dx),
      (2 + dy):(dims[2] + 1L + dy),
      (2 + dz):(dims[3] + 1L + dz)
    ]
    has_background_neighbor <- has_background_neighbor | !neighbor_mask
  }

  surface_mask <- roi_mask & has_background_neighbor

  if (output == "mask") {
    return(surface_mask)
  }

  coords <- which(surface_mask, arr.ind = TRUE)
  surface_df <- as.data.frame(coords)
  colnames(surface_df) <- c("x", "y", "z")
  surface_df$label_value <- as.numeric(label_space[coords])

  if (!is.null(voxel_spacing)) {
    xyz_mm <- sweep(
      as.matrix(surface_df[, c("x", "y", "z")]),
      2,
      voxel_spacing,
      `*`
    )
    surface_df$x_mm <- xyz_mm[, 1]
    surface_df$y_mm <- xyz_mm[, 2]
    surface_df$z_mm <- xyz_mm[, 3]
  }

  surface_df
}

#' Estimate a Midplane From Surface Coordinates
#'
#' The plane is defined as `normal' * (point - center) = 0`. The normal is
#' chosen from a PCA/SVD direction of the point cloud, making the rule reusable
#' for both template and subject shapes.
#'
#' @param surface_points Data frame or matrix containing 3D points.
#' @param coord_cols Coordinate columns to use when `surface_points` is a data
#'   frame.
#' @param normal_pc Which PCA direction is used as the plane normal. PC1 cuts
#'   across the longest axis, PC2 across the second axis, and PC3 is the
#'   least-variance direction.
#' @param center Either `"mean"` or `"median"`.
#' @param reference_normal Optional numeric length-3 vector used to make the
#'   normal direction sign consistent across template/subjects.
#'
#' @return A list with center, unit normal, plane coefficients, PCA rotation,
#'   standard deviations along PCs, and signed distances for the input points.
#' @export
estimate_midplane <- function(surface_points,
                              coord_cols = c("x", "y", "z"),
                              normal_pc = 3,
                              center = c("mean", "median"),
                              reference_normal = NULL) {
  center <- match.arg(center)
  normal_pc <- as.integer(normal_pc)

  if (!normal_pc %in% 1:3) {
    stop("normal_pc must be 1, 2, or 3.")
  }

  coords <- if (is.matrix(surface_points)) {
    surface_points[, seq_len(3), drop = FALSE]
  } else {
    missing_cols <- setdiff(coord_cols, names(surface_points))
    if (length(missing_cols) > 0) {
      stop(
        "surface_points is missing coordinate columns: ",
        paste(missing_cols, collapse = ", ")
      )
    }
    as.matrix(surface_points[, coord_cols, drop = FALSE])
  }
  storage.mode(coords) <- "double"

  keep <- stats::complete.cases(coords)
  coords <- coords[keep, , drop = FALSE]
  if (nrow(coords) < 3) {
    stop("At least three complete 3D points are required.")
  }

  plane_center <- if (center == "mean") {
    colMeans(coords)
  } else {
    apply(coords, 2, stats::median)
  }

  centered_coords <- sweep(coords, 2, plane_center, "-")
  svd_fit <- svd(centered_coords, nu = 0, nv = 3)
  normal <- svd_fit$v[, normal_pc]
  normal <- as.numeric(normal / sqrt(sum(normal^2)))

  if (!is.null(reference_normal)) {
    if (length(reference_normal) != 3) {
      stop("reference_normal must be NULL or a numeric vector of length 3.")
    }
    reference_normal <- as.numeric(reference_normal)
    reference_normal <- reference_normal / sqrt(sum(reference_normal^2))
    if (sum(normal * reference_normal) < 0) {
      normal <- -normal
    }
  }

  signed_distance <- as.numeric(centered_coords %*% normal)
  names(plane_center) <- coord_cols
  names(normal) <- coord_cols

  list(
    center = plane_center,
    normal = normal,
    coefficients = c(normal, d = -sum(normal * plane_center)),
    normal_pc = normal_pc,
    pc_rotation = svd_fit$v,
    pc_sd = svd_fit$d / sqrt(nrow(coords) - 1),
    signed_distance = signed_distance
  )
}

#' Partition Surface Points by a Midplane
#'
#' Add partition indicators to surface points using signed distance to an
#' estimated midplane. Points inside the overlap band are assigned to both
#' partitions.
#'
#' @param surface_points Data frame or matrix containing 3D points.
#' @param midplane Output from `estimate_midplane()`.
#' @param coord_cols Coordinate columns to use when `surface_points` is a data
#'   frame.
#' @param overlap_prop Non-negative proportion of the largest absolute signed
#'   distance from points to the plane. For example, `0.10` uses an overlap band
#'   with half-width 10 percent of `max(abs(signed_distance))`.
#'
#' @return A data frame with `signed_distance`, `core_partition`, `in_overlap`,
#'   `partition1`, and `partition2` columns added.
#' @export
partition_by_midplane <- function(surface_points,
                                  midplane,
                                  coord_cols = c("x", "y", "z"),
                                  overlap_prop = 0.15) {
  if (overlap_prop < 0) {
    stop("overlap_prop must be non-negative.")
  }
  if (is.null(midplane$center) || is.null(midplane$normal)) {
    stop("midplane must be the output of estimate_midplane().")
  }

  coords <- if (is.matrix(surface_points)) {
    surface_points[, seq_len(3), drop = FALSE]
  } else {
    missing_cols <- setdiff(coord_cols, names(surface_points))
    if (length(missing_cols) > 0) {
      stop(
        "surface_points is missing coordinate columns: ",
        paste(missing_cols, collapse = ", ")
      )
    }
    as.matrix(surface_points[, coord_cols, drop = FALSE])
  }
  storage.mode(coords) <- "double"

  signed_distance <- as.numeric(
    sweep(coords, 2, as.numeric(midplane$center), "-") %*%
      as.numeric(midplane$normal)
  )
  overlap_width <- max(abs(signed_distance), na.rm = TRUE) * overlap_prop

  out <- as.data.frame(surface_points)
  out$signed_distance <- signed_distance
  out$core_partition <- ifelse(signed_distance >= 0, 1L, 2L)
  out$in_overlap <- abs(signed_distance) <= overlap_width
  out$partition1 <- signed_distance >= -overlap_width
  out$partition2 <- signed_distance <= overlap_width
  attr(out, "overlap_width") <- overlap_width
  attr(out, "overlap_prop") <- overlap_prop
  out
}
