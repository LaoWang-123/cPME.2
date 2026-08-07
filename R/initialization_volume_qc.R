#' Compute a proxy energy between template and subject initial PME fits
#'
#' This is the lightweight initialization-ranking criterion used before full
#' registration. It compares the two fitted PME embeddings on coarse grids over
#' their current parameter domains.
#'
#' @param pme1 Template PME fit.
#' @param pme2 Subject PME fit.
#' @param init1 Template initialization.
#' @param init2 Subject initialization.
#' @param n_grid Number of grid points per parameter direction.
#'
#' @return Mean squared discrepancy between coarse template and subject grids.
#' @export
compute_initialization_proxy_energy <- function(pme1, pme2, init1, init2, n_grid = 10L) {
  U1 <- as.matrix(init1$parameterization)
  U2 <- as.matrix(init2$parameterization)

  r1_min <- apply(U1, 2L, min)
  r1_max <- apply(U1, 2L, max)
  r2_min <- apply(U2, 2L, min)
  r2_max <- apply(U2, 2L, max)

  U_eval1 <- as.matrix(expand.grid(
    seq(r1_min[1], r1_max[1], length.out = n_grid),
    seq(r1_min[2], r1_max[2], length.out = n_grid)
  ))
  U_eval2 <- as.matrix(expand.grid(
    seq(r2_min[1], r2_max[1], length.out = n_grid),
    seq(r2_min[2], r2_max[2], length.out = n_grid)
  ))

  vals1 <- t(apply(U_eval1, 1L, function(u) pme1$embedding_map(as.numeric(u))))
  vals2 <- t(apply(U_eval2, 1L, function(u) pme2$embedding_map(as.numeric(u))))

  mean(rowSums((vals1 - vals2)^2))
}

#' Compute voxel target volumes for two partition sides
#'
#' @param voxel_points ROI voxel coordinates as a data frame or matrix. If this
#'   is \code{NULL}, coordinates are extracted from \code{label_space}.
#' @param label_space Optional 3D label array.
#' @param label_values Label value or values identifying the ROI in
#'   \code{label_space}.
#' @param partition_plane Partition plane with \code{center} and \code{normal}.
#' @param coord_cols Coordinate-column names when \code{voxel_points} is a data
#'   frame.
#' @param voxel_volume Volume represented by one voxel.
#'
#' @return A two-row data frame with \code{part}, \code{n_voxels}, and
#'   \code{volume}.
#' @export
compute_partition_voxel_targets <- function(
    voxel_points = NULL,
    label_space = NULL,
    label_values = NULL,
    partition_plane,
    coord_cols = c("x", "y", "z"),
    voxel_volume = 1) {
  if (is.null(voxel_points)) {
    if (is.null(label_space) || is.null(label_values)) {
      stop("Provide either voxel_points or both label_space and label_values.", call. = FALSE)
    }
    mask <- label_space %in% label_values
    dim(mask) <- dim(label_space)
    voxel_points <- as.data.frame(which(mask, arr.ind = TRUE))
    colnames(voxel_points) <- coord_cols
  } else {
    voxel_points <- as.data.frame(voxel_points)
    if (!all(coord_cols %in% colnames(voxel_points))) {
      if (ncol(voxel_points) < length(coord_cols)) {
        stop("voxel_points must have at least three coordinate columns.", call. = FALSE)
      }
      colnames(voxel_points)[seq_along(coord_cols)] <- coord_cols
    }
  }

  points <- as.matrix(voxel_points[, coord_cols, drop = FALSE])
  if (nrow(points) == 0L) {
    stop("No ROI voxels were found for the requested partition target.", call. = FALSE)
  }
  part1_index <- signed_distance_to_partition_plane(points, partition_plane) >= 0
  n_part1 <- sum(part1_index)
  n_part2 <- sum(!part1_index)

  data.frame(
    part = c("part1", "part2"),
    n_voxels = c(n_part1, n_part2),
    volume = c(n_part1, n_part2) * voxel_volume,
    stringsAsFactors = FALSE
  )
}

#' Estimate one partition surface volume against a partition side
#'
#' @param fit PME/SIME-like fitted surface object.
#' @param data Full ROI point cloud used to define the sampling box.
#' @param partition_plane Optional partition plane with \code{center} and
#'   \code{normal}. If omitted, \code{partition_index} and the resolved
#'   reference point are used.
#' @param part_label Either \code{"part1"} or \code{"part2"}.
#' @param n_points Number of Monte Carlo candidate points.
#' @param limit_scaler Sampling-box padding passed to the volume sampler.
#' @param partition_index Coordinate index used when \code{partition_plane} is
#'   \code{NULL}.
#' @param ref Optional reference point for interior classification.
#' @param data_max Optional sampling-box scale.
#' @param volume_multiplier Unit multiplier for the returned volume.
#' @param seed Optional random seed.
#' @param fail_on_error If \code{TRUE}, throw classification errors.
#'
#' @return A volume-result list with \code{volume}, \code{volume_success}, and
#'   diagnostic fields.
#' @export
estimate_single_partition_surface_volume <- function(
    fit,
    data,
    partition_plane = NULL,
    part_label = c("part1", "part2"),
    n_points = 5000,
    limit_scaler = 0.05,
    partition_index = 3,
    ref = NULL,
    data_max = NULL,
    volume_multiplier = 1,
    seed = NULL,
    fail_on_error = FALSE) {
  part_label <- match.arg(part_label)
  if (!is.null(seed)) {
    set.seed(seed)
  }

  ref_info <- resolve_volume_ref_point(
    ref = ref,
    data = data,
    partition_plane = partition_plane
  )
  candidate_info <- sample_volume_candidates(
    data = data,
    n_points = n_points,
    limit_scaler = limit_scaler,
    data_max = data_max
  )

  candidates <- candidate_info$candidates
  if (is.null(partition_plane)) {
    part1_index <- candidates[, partition_index] > ref_info$ref[partition_index]
  } else {
    part1_index <- signed_distance_to_partition_plane(candidates, partition_plane) >= 0
  }
  side_index <- if (part_label == "part1") part1_index else !part1_index

  interior <- rep(FALSE, nrow(candidates))
  classification_result <- tryCatch(
    surface_interior_identification(
      fit = fit,
      x = candidates[side_index, , drop = FALSE],
      ref = ref_info$ref
    ),
    error = function(err) err
  )

  if (inherits(classification_result, "error")) {
    msg <- paste(
      "Single partition surface volume not estimated because candidate classification failed.",
      conditionMessage(classification_result)
    )
    if (isTRUE(fail_on_error)) {
      stop(msg, call. = FALSE)
    }
    return(list(
      volume = NA_real_,
      volume_success = FALSE,
      message = msg,
      n_side_candidates = sum(side_index),
      candidates = candidates,
      full_volume = candidate_info$full_volume,
      bounds = candidate_info$bounds,
      ref = ref_info$ref,
      ref_method = ref_info$ref_method
    ))
  }

  interior[side_index] <- classification_result
  list(
    volume = mean(interior) * candidate_info$full_volume * volume_multiplier,
    volume_success = TRUE,
    message = NA_character_,
    n_side_candidates = sum(side_index),
    interior = interior,
    candidates = candidates,
    full_volume = candidate_info$full_volume,
    bounds = candidate_info$bounds,
    ref = ref_info$ref,
    ref_method = ref_info$ref_method
  )
}

#' Select an initial subject PME using orientation, proxy energy, and volume QC
#'
#' This function runs repeated subject PME initializations, filters mirror-
#' reversed parameterizations relative to a template initialization, ranks the
#' surviving candidates by a registration proxy energy, and selects the first
#' proxy-energy-ranked candidate whose single-part surface-volume relative
#' error is no larger than \code{volume_rel_error_threshold}. If no candidate
#' passes the threshold, the candidate with the smallest volume relative error
#' is selected and flagged by \code{selection_rule}.
#'
#' @param template_pme Template PME fit.
#' @param template_initialization Template initialization object.
#' @param subject_data Subject partition point cloud used to fit PME.
#' @param volume_data Full ROI point cloud used for single-part volume
#'   estimation.
#' @param partition_plane Partition plane with \code{center} and \code{normal}.
#' @param part_label Either \code{"part1"} or \code{"part2"}.
#' @param target_volume Target single-part volume in output units.
#' @param n_trials Number of initialization trials.
#' @param volume_rel_error_threshold Relative volume error threshold.
#' @param init_args Arguments passed to \code{initialize_pme()}.
#' @param pme_args Arguments passed to \code{pme()}.
#' @param volume_n_points Number of Monte Carlo points for volume QC.
#' @param volume_seed Optional base seed for volume estimates.
#' @param init_seed Optional seed before initialization generation.
#' @param proxy_energy_grid Number of grid points per parameter direction for
#'   proxy energy.
#' @param orientation_pca_source Passed to \code{check_pme_orientation()}.
#' @param volume_multiplier Unit multiplier for volume estimates.
#' @param limit_scaler Sampling-box padding for volume estimation.
#' @param ref Optional reference point for volume classification.
#' @param data_max Optional sampling-box scale.
#' @param verbose Logical.
#'
#' @return A list containing the selected initialization, selected PME,
#'   selected summary row, all candidate rows, and selection metadata.
#' @export
select_initial_pme_volume_qc <- function(
    template_pme,
    template_initialization,
    subject_data,
    volume_data,
    partition_plane,
    part_label = c("part1", "part2"),
    target_volume,
    n_trials = 100L,
    volume_rel_error_threshold = 0.5,
    init_args = list(d = 2),
    pme_args = list(d = 2, lambda = exp(-15:5), max_iter = 100, verbose = FALSE, print_plots = FALSE),
    volume_n_points = 5000L,
    volume_seed = NULL,
    init_seed = NULL,
    proxy_energy_grid = 10L,
    orientation_pca_source = "all_centers",
    volume_multiplier = 1,
    limit_scaler = 0.05,
    ref = NULL,
    data_max = NULL,
    verbose = FALSE) {
  part_label <- match.arg(part_label)
  subject_data <- as.matrix(subject_data)
  volume_data <- as.matrix(volume_data)

  if (!is.null(init_seed)) {
    set.seed(init_seed)
  }

  candidate_inits_all <- vector("list", n_trials)
  candidate_checks_all <- vector("list", n_trials)

  for (i in seq_len(n_trials)) {
    candidate_inits_all[[i]] <- do.call(
      initialize_pme,
      c(list(x = subject_data), init_args)
    )
    candidate_checks_all[[i]] <- check_pme_orientation(
      init1 = template_initialization,
      init2 = candidate_inits_all[[i]],
      pca_source = orientation_pca_source,
      verbose = verbose
    )
  }

  keep <- which(vapply(
    candidate_checks_all,
    function(chk) !identical(chk$final, "mirror_reversed"),
    logical(1)
  ))
  if (length(keep) == 0L) {
    stop("All ", part_label, " initialization candidates were mirror_reversed.")
  }

  candidate_inits <- candidate_inits_all[keep]
  candidate_checks <- candidate_checks_all[keep]
  n_candidates <- length(candidate_inits)
  candidate_pmes <- vector("list", n_candidates)
  candidate_energy <- rep(Inf, n_candidates)

  candidate_pmes[[1]] <- do.call(
    pme,
    c(list(data = subject_data, initialization = candidate_inits[[1]]), pme_args)
  )
  lambda_reuse <- candidate_pmes[[1]]$tuning
  candidate_energy[1] <- compute_initialization_proxy_energy(
    template_pme,
    candidate_pmes[[1]],
    template_initialization,
    candidate_inits[[1]],
    n_grid = proxy_energy_grid
  )

  if (n_candidates >= 2L) {
    for (i in 2:n_candidates) {
      pme_args_i <- pme_args
      pme_args_i$lambda <- lambda_reuse
      candidate_pmes[[i]] <- do.call(
        pme,
        c(list(data = subject_data, initialization = candidate_inits[[i]]), pme_args_i)
      )
      candidate_energy[i] <- compute_initialization_proxy_energy(
        template_pme,
        candidate_pmes[[i]],
        template_initialization,
        candidate_inits[[i]],
        n_grid = proxy_energy_grid
      )
    }
  }

  candidate_rows <- lapply(seq_len(n_candidates), function(i) {
    seed_i <- if (is.null(volume_seed)) NULL else volume_seed + i
    volume_result <- estimate_single_partition_surface_volume(
      fit = candidate_pmes[[i]],
      data = volume_data,
      partition_plane = partition_plane,
      part_label = part_label,
      n_points = volume_n_points,
      limit_scaler = limit_scaler,
      ref = ref,
      data_max = data_max,
      volume_multiplier = volume_multiplier,
      seed = seed_i
    )

    volume_rel_error <- abs(volume_result$volume - target_volume) / target_volume

    data.frame(
      part = part_label,
      candidate = i,
      original_trial = keep[i],
      orientation_status = candidate_checks[[i]]$final,
      candidate_E = candidate_energy[i],
      lambda = candidate_pmes[[i]]$tuning,
      pme_msd = min(candidate_pmes[[i]]$MSD, na.rm = TRUE),
      target_volume = target_volume,
      volume = volume_result$volume,
      volume_success = isTRUE(volume_result$volume_success),
      volume_rel_error = volume_rel_error,
      volume_error_le_threshold = isTRUE(volume_result$volume_success) &&
        is.finite(volume_rel_error) &&
        volume_rel_error <= volume_rel_error_threshold,
      error_message = volume_result$message,
      stringsAsFactors = FALSE
    )
  })

  candidates <- do.call(rbind, candidate_rows)
  candidates <- candidates[order(candidates$candidate_E), , drop = FALSE]
  candidates$proxy_energy_rank <- seq_len(nrow(candidates))

  passing <- candidates[
    candidates$volume_success & candidates$volume_rel_error <= volume_rel_error_threshold,
    ,
    drop = FALSE
  ]

  if (nrow(passing) > 0L) {
    selected_candidate <- passing$candidate[1]
    selection_rule <- paste0(
      "first proxy-energy-ranked candidate with volume_rel_error <= ",
      volume_rel_error_threshold
    )
  } else {
    ord <- order(candidates$volume_rel_error, candidates$candidate_E)
    selected_candidate <- candidates$candidate[ord[1]]
    selection_rule <- paste0(
      "fallback: no candidate passed volume_rel_error <= ",
      volume_rel_error_threshold,
      "; selected smallest volume_rel_error"
    )
  }

  candidates$selected <- candidates$candidate == selected_candidate
  selected_summary <- candidates[candidates$selected, , drop = FALSE]

  list(
    part = part_label,
    threshold = volume_rel_error_threshold,
    target_volume = target_volume,
    n_trials_requested = n_trials,
    n_candidates_after_orientation = n_candidates,
    lambda_reuse = lambda_reuse,
    selection_rule = selection_rule,
    selected_candidate = selected_candidate,
    selected_initialization = candidate_inits[[selected_candidate]],
    selected_pme = candidate_pmes[[selected_candidate]],
    selected_summary = selected_summary,
    candidates = candidates
  )
}
