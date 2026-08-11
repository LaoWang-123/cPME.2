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
#' surviving candidates by a registration proxy energy, and evaluates volume as
#' a guardrail. By default, candidates are checked in increasing proxy-energy
#' order and the first one passing the volume-error threshold is selected. This
#' keeps proxy energy as the primary criterion while excluding implausible
#' initial PME surfaces.
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
#' @param energy_window_multiplier Keep candidates whose proxy energy is no
#'   larger than this multiplier times the minimum proxy energy for volume QC
#'   and selection. Use \code{Inf} to evaluate all orientation-valid candidates.
#' @param selection_strategy Either \code{"first_energy_passing_volume"} or
#'   \code{"smallest_volume_error_in_energy_window"}. The first strategy selects
#'   the lowest-energy candidate whose volume relative error is below
#'   \code{volume_rel_error_threshold}; the second is the previous behavior.
#' @param save_candidate_objects If \code{TRUE}, return all orientation-valid
#'   initializations, PME fits, and orientation checks in the output object.
#' @param evaluate_volume Logical. If \code{TRUE}, evaluate the legacy
#'   single-part volume guardrail inside this function. If \code{FALSE}, only
#'   generate orientation-valid candidate PME fits and proxy-energy rankings;
#'   the minimum-proxy-energy candidate is stored as a provisional selection.
#'   Use \code{FALSE} when candidates will be passed to
#'   \code{select_initial_pme_pair_volume_qc()} or
#'   \code{refine_initial_pme_pair_normal_ray_qc()}.
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
    energy_window_multiplier = 1.2,
    selection_strategy = c("first_energy_passing_volume", "smallest_volume_error_in_energy_window"),
    save_candidate_objects = TRUE,
    evaluate_volume = TRUE,
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
  selection_strategy <- match.arg(selection_strategy)
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
  if (!is.numeric(energy_window_multiplier) ||
      length(energy_window_multiplier) != 1L ||
      is.na(energy_window_multiplier) ||
      energy_window_multiplier < 1) {
    stop("energy_window_multiplier must be a numeric scalar >= 1.", call. = FALSE)
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
    data.frame(
      part = part_label,
      candidate = i,
      original_trial = keep[i],
      orientation_status = candidate_checks[[i]]$final,
      candidate_E = candidate_energy[i],
      lambda = candidate_pmes[[i]]$tuning,
      pme_msd = min(candidate_pmes[[i]]$MSD, na.rm = TRUE),
      target_volume = target_volume,
      volume = NA_real_,
      volume_success = FALSE,
      volume_rel_error = NA_real_,
      volume_error_le_threshold = FALSE,
      in_energy_window = FALSE,
      error_message = NA_character_,
      stringsAsFactors = FALSE
    )
  })

  candidates <- do.call(rbind, candidate_rows)
  candidates <- candidates[order(candidates$candidate_E), , drop = FALSE]
  candidates$proxy_energy_rank <- seq_len(nrow(candidates))
  min_energy <- min(candidates$candidate_E, na.rm = TRUE)
  energy_cutoff <- min_energy * energy_window_multiplier
  candidates$in_energy_window <- FALSE

  if (selection_strategy == "first_energy_passing_volume") {
    volume_candidate_ids <- candidates$candidate
  } else {
    candidates$in_energy_window <- candidates$candidate_E <= energy_cutoff
    volume_candidate_ids <- candidates$candidate[candidates$in_energy_window]
    if (length(volume_candidate_ids) == 0L) {
      volume_candidate_ids <- candidates$candidate[1L]
      candidates$in_energy_window[candidates$candidate == volume_candidate_ids] <- TRUE
    }
  }

  volume_rows <- list()
  selected_candidate <- NA_integer_
  if (isTRUE(evaluate_volume)) {
    for (i in volume_candidate_ids) {
      candidates$in_energy_window[candidates$candidate == i] <- TRUE
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

      row <- list(
        candidate = i,
        volume = volume_result$volume,
        volume_success = isTRUE(volume_result$volume_success),
        volume_rel_error = volume_rel_error,
        volume_error_le_threshold = isTRUE(volume_result$volume_success) &&
          is.finite(volume_rel_error) &&
          volume_rel_error <= volume_rel_error_threshold,
        error_message = volume_result$message
      )
      volume_rows[[length(volume_rows) + 1L]] <- row

      if (selection_strategy == "first_energy_passing_volume" &&
          isTRUE(row$volume_error_le_threshold)) {
        selected_candidate <- i
        break
      }
    }
  }

  for (row in volume_rows) {
    idx <- candidates$candidate == row$candidate
    candidates$volume[idx] <- row$volume
    candidates$volume_success[idx] <- row$volume_success
    candidates$volume_rel_error[idx] <- row$volume_rel_error
    candidates$volume_error_le_threshold[idx] <- row$volume_error_le_threshold
    candidates$error_message[idx] <- row$error_message
  }

  if (!isTRUE(evaluate_volume)) {
    selected_candidate <- candidates$candidate[1]
    candidates$in_energy_window[candidates$candidate == selected_candidate] <- TRUE
    selection_rule <- "provisional minimum proxy energy; volume not evaluated"
  } else if (selection_strategy == "first_energy_passing_volume" &&
      is.finite(selected_candidate)) {
    selection_rule <- paste0(
      "first candidate by increasing proxy energy with volume_rel_error <= ",
      volume_rel_error_threshold
    )
  } else {
    window_success <- candidates[
      candidates$in_energy_window & candidates$volume_success & is.finite(candidates$volume_rel_error),
      ,
      drop = FALSE
    ]

    if (selection_strategy == "smallest_volume_error_in_energy_window" &&
        nrow(window_success) > 0L) {
      ord <- order(window_success$volume_rel_error, window_success$candidate_E)
      selected_candidate <- window_success$candidate[ord[1]]
      selection_rule <- paste0(
        "smallest volume_rel_error among candidates with candidate_E <= ",
        energy_window_multiplier,
        " * min(candidate_E)"
      )
    } else {
      selected_candidate <- candidates$candidate[1]
      selection_rule <- paste0(
        "fallback: no candidate passed volume_rel_error <= ",
        volume_rel_error_threshold,
        "; selected minimum proxy energy"
      )
    }
  }

  candidates$selected <- candidates$candidate == selected_candidate
  selected_summary <- candidates[candidates$selected, , drop = FALSE]

  out <- list(
    part = part_label,
    threshold = volume_rel_error_threshold,
    energy_window_multiplier = energy_window_multiplier,
    min_candidate_E = min_energy,
    energy_cutoff = energy_cutoff,
    target_volume = target_volume,
    n_trials_requested = n_trials,
    n_candidates_after_orientation = n_candidates,
    n_candidates_in_energy_window = sum(candidates$in_energy_window),
    selection_strategy = selection_strategy,
    lambda_reuse = lambda_reuse,
    selection_rule = selection_rule,
    selected_candidate = selected_candidate,
    selected_initialization = candidate_inits[[selected_candidate]],
    selected_pme = candidate_pmes[[selected_candidate]],
    selected_summary = selected_summary,
    candidates = candidates
  )
  if (isTRUE(save_candidate_objects)) {
    out$candidate_initializations <- candidate_inits
    out$candidate_pmes <- candidate_pmes
    out$candidate_orientation_checks <- candidate_checks
    out$kept_original_trials <- keep
  }
  out
}

#' Select paired initial PME fits with normal-ray volume QC
#'
#' This is the preferred initialization selector for two-part real-data
#' workflows. It generates orientation-valid PME candidates and PME-based proxy
#' energy rankings for part 1 and part 2, then selects the final pair using
#' paired normal-ray volume QC. Proxy energy remains the primary ranking
#' criterion; volume is used as a guardrail to reject implausible candidates.
#'
#' @param template_part1,template_part2 Template partition objects. Each object
#'   must contain \code{pme1} and \code{initialization_f1}.
#' @param subject_part1_data,subject_part2_data Subject partition point clouds.
#' @param volume_data Full ROI point cloud used for volume estimation.
#' @param partition_plane Partition plane with \code{center} and \code{normal}.
#' @param target_volumes Named numeric vector or list with \code{part1} and
#'   \code{part2} target volumes.
#' @param n_trials Number of initialization trials per part.
#' @param volume_rel_error_threshold Relative volume error threshold.
#' @param init_seed_part1,init_seed_part2 Optional initialization seeds.
#' @param volume_seed Optional base seed for paired volume estimates.
#' @param volume_seed_offsets Named or length-2 numeric vector giving part1 and
#'   part2 seed offsets.
#' @param init_args Arguments passed to \code{initialize_pme()}.
#' @param pme_args Arguments passed to \code{pme()}.
#' @param volume_n_points Number of Monte Carlo candidate points for volume QC.
#' @param volume_multiplier Unit multiplier for volume estimates.
#' @param limit_scaler Sampling-box padding for volume estimation.
#' @param proxy_energy_grid Number of grid points per parameter direction for
#'   PME-based proxy energy.
#' @param orientation_pca_source Passed to \code{check_pme_orientation()}.
#' @param grid_n Number of parameter-grid points per direction used to seed
#'   normal-ray intersection search.
#' @param ray_tol Maximum perpendicular distance for accepting a ray hit.
#' @param ray_min_t Minimum distance along the ray before a hit is treated as a
#'   crossing rather than a partition-plane boundary contact.
#' @param max_passes Number of alternating paired refinement passes.
#' @param verbose Logical.
#'
#' @return A list with refined \code{part1}, refined \code{part2}, and combined
#'   \code{candidates}. Each part contains all candidate objects and the selected
#'   initialization/PME.
#' @export
select_initial_pme_pair_volume_qc <- function(
    template_part1,
    template_part2,
    subject_part1_data,
    subject_part2_data,
    volume_data,
    partition_plane,
    target_volumes,
    n_trials = 100L,
    volume_rel_error_threshold = 0.5,
    init_seed_part1 = NULL,
    init_seed_part2 = NULL,
    volume_seed = NULL,
    volume_seed_offsets = c(part1 = 10000, part2 = 20000),
    init_args = list(d = 2),
    pme_args = list(d = 2, lambda = exp(-15:5), max_iter = 100, verbose = FALSE, print_plots = FALSE),
    volume_n_points = 5000L,
    volume_multiplier = 1,
    limit_scaler = 0.05,
    proxy_energy_grid = 10L,
    orientation_pca_source = "all_centers",
    grid_n = 75L,
    ray_tol = 0.45,
    ray_min_t = 0.05,
    max_passes = 2L,
    verbose = FALSE) {
  target_volumes <- unlist(target_volumes)
  if (!all(c("part1", "part2") %in% names(target_volumes))) {
    stop("target_volumes must contain named entries 'part1' and 'part2'.", call. = FALSE)
  }
  required_template_fields <- c("pme1", "initialization_f1")
  for (field in required_template_fields) {
    if (is.null(template_part1[[field]]) || is.null(template_part2[[field]])) {
      stop("template_part1 and template_part2 must each contain '", field, "'.", call. = FALSE)
    }
  }

  init_qc_part1 <- select_initial_pme_volume_qc(
    template_pme = template_part1$pme1,
    template_initialization = template_part1$initialization_f1,
    subject_data = subject_part1_data,
    volume_data = volume_data,
    partition_plane = partition_plane,
    part_label = "part1",
    target_volume = target_volumes[["part1"]],
    n_trials = n_trials,
    volume_rel_error_threshold = volume_rel_error_threshold,
    energy_window_multiplier = Inf,
    selection_strategy = "first_energy_passing_volume",
    save_candidate_objects = TRUE,
    evaluate_volume = FALSE,
    init_args = init_args,
    pme_args = pme_args,
    volume_n_points = volume_n_points,
    volume_seed = volume_seed,
    init_seed = init_seed_part1,
    proxy_energy_grid = proxy_energy_grid,
    orientation_pca_source = orientation_pca_source,
    volume_multiplier = volume_multiplier,
    limit_scaler = limit_scaler,
    verbose = verbose
  )

  init_qc_part2 <- select_initial_pme_volume_qc(
    template_pme = template_part2$pme1,
    template_initialization = template_part2$initialization_f1,
    subject_data = subject_part2_data,
    volume_data = volume_data,
    partition_plane = partition_plane,
    part_label = "part2",
    target_volume = target_volumes[["part2"]],
    n_trials = n_trials,
    volume_rel_error_threshold = volume_rel_error_threshold,
    energy_window_multiplier = Inf,
    selection_strategy = "first_energy_passing_volume",
    save_candidate_objects = TRUE,
    evaluate_volume = FALSE,
    init_args = init_args,
    pme_args = pme_args,
    volume_n_points = volume_n_points,
    volume_seed = volume_seed,
    init_seed = init_seed_part2,
    proxy_energy_grid = proxy_energy_grid,
    orientation_pca_source = orientation_pca_source,
    volume_multiplier = volume_multiplier,
    limit_scaler = limit_scaler,
    verbose = verbose
  )

  refined <- refine_initial_pme_pair_normal_ray_qc(
    init_qc_part1 = init_qc_part1,
    init_qc_part2 = init_qc_part2,
    volume_data = volume_data,
    partition_plane = partition_plane,
    volume_n_points = volume_n_points,
    volume_seed = volume_seed,
    volume_seed_offsets = volume_seed_offsets,
    volume_multiplier = volume_multiplier,
    limit_scaler = limit_scaler,
    grid_n = grid_n,
    ray_tol = ray_tol,
    ray_min_t = ray_min_t,
    max_passes = max_passes
  )
  refined$selection_strategy <- "paired_proxy_energy_normal_ray_volume_qc"
  refined$n_trials_requested <- n_trials
  refined$volume_rel_error_threshold <- volume_rel_error_threshold
  refined$target_volumes <- target_volumes[c("part1", "part2")]
  refined
}

#' Refine paired initial PME choices with normal-ray volume QC
#'
#' Re-select two partition initial PME fits after both parts have already been
#' initialized and proxy-energy ranked. Each part is scanned in increasing proxy
#' energy order, and the first candidate passing single-part normal-ray volume
#' QC is selected. The other part's current selected PME is used as the
#' counterpart surface for normal-ray intersection-order diagnostics.
#'
#' @param init_qc_part1 Output from \code{select_initial_pme_volume_qc()} for
#'   part 1, saved with \code{save_candidate_objects = TRUE}.
#' @param init_qc_part2 Output from \code{select_initial_pme_volume_qc()} for
#'   part 2, saved with \code{save_candidate_objects = TRUE}.
#' @param volume_data Full ROI point cloud used for volume estimation.
#' @param partition_plane Partition plane with \code{center} and \code{normal}.
#' @param volume_n_points Number of Monte Carlo candidate points.
#' @param volume_seed Optional base seed for volume estimates.
#' @param volume_seed_offsets Named or length-2 numeric vector giving part1 and
#'   part2 seed offsets.
#' @param volume_multiplier Unit multiplier for volume estimates.
#' @param limit_scaler Sampling-box padding for volume estimation.
#' @param grid_n Number of parameter-grid points per direction used to seed
#'   normal-ray intersection search.
#' @param ray_tol Maximum perpendicular distance for accepting a ray hit.
#' @param ray_min_t Minimum distance along the ray before a hit is treated as a
#'   crossing rather than a partition-plane boundary contact.
#' @param max_passes Number of alternating refinement passes. Two is usually
#'   enough because part 2 can use the updated part 1 choice in the first pass.
#'
#' @return A list with refined \code{part1}, refined \code{part2}, and a
#'   combined \code{candidates} data frame.
refine_initial_pme_pair_normal_ray_qc <- function(
    init_qc_part1,
    init_qc_part2,
    volume_data,
    partition_plane,
    volume_n_points = 5000L,
    volume_seed = NULL,
    volume_seed_offsets = c(part1 = 10000, part2 = 20000),
    volume_multiplier = 1,
    limit_scaler = 0.05,
    grid_n = 75L,
    ray_tol = 0.45,
    ray_min_t = 0.05,
    max_passes = 2L) {
  volume_data <- as.matrix(volume_data)
  if (length(volume_seed_offsets) < 2L) {
    stop("volume_seed_offsets must contain at least two values.", call. = FALSE)
  }
  if (is.null(names(volume_seed_offsets))) {
    names(volume_seed_offsets) <- c("part1", "part2")
  }
  if (!all(c("part1", "part2") %in% names(volume_seed_offsets))) {
    names(volume_seed_offsets)[seq_len(2L)] <- c("part1", "part2")
  }

  required <- c("candidate_pmes", "candidate_initializations", "candidates")
  for (nm in required) {
    if (is.null(init_qc_part1[[nm]]) || is.null(init_qc_part2[[nm]])) {
      stop(
        "Both init_qc_part1 and init_qc_part2 must contain '", nm,
        "'. Run select_initial_pme_volume_qc(save_candidate_objects = TRUE).",
        call. = FALSE
      )
    }
  }

  refine_one <- function(init_qc, counterpart_qc, part_label, seed_offset) {
    candidates <- as.data.frame(init_qc$candidates)
    add_col <- function(name, value) {
      if (!name %in% names(candidates)) {
        candidates[[name]] <<- value
      }
    }
    add_col("normal_ray_volume", NA_real_)
    add_col("normal_ray_volume_rel_error", NA_real_)
    add_col("normal_ray_volume_success", FALSE)
    add_col("normal_ray_flip_this_part", NA)
    add_col("normal_ray_flip_part1", NA)
    add_col("normal_ray_flip_part2", NA)
    add_col("normal_ray_n_positive", NA_integer_)
    add_col("normal_ray_n_negative", NA_integer_)

    candidate_order <- candidates$candidate[order(candidates$proxy_energy_rank)]
    selected_candidate <- NA_integer_

    for (candidate_id in candidate_order) {
      volume_result <- estimate_single_partition_surface_volume_normal_ray(
        fit = init_qc$candidate_pmes[[candidate_id]],
        counterpart_fit = counterpart_qc$selected_pme,
        data = volume_data,
        partition_plane = partition_plane,
        part_label = part_label,
        n_points = volume_n_points,
        limit_scaler = limit_scaler,
        volume_multiplier = volume_multiplier,
        seed = if (is.null(volume_seed)) NULL else volume_seed + seed_offset + candidate_id,
        grid_n = grid_n,
        ray_tol = ray_tol,
        ray_min_t = ray_min_t
      )

      idx <- candidates$candidate == candidate_id
      rel_error <- abs(volume_result$volume - init_qc$target_volume) / init_qc$target_volume
      candidates$volume[idx] <- volume_result$volume
      candidates$volume_success[idx] <- isTRUE(volume_result$volume_success)
      candidates$volume_rel_error[idx] <- rel_error
      candidates$volume_error_le_threshold[idx] <- isTRUE(volume_result$volume_success) &&
        is.finite(rel_error) &&
        rel_error <= init_qc$threshold
      candidates$error_message[idx] <- volume_result$message
      candidates$in_energy_window[idx] <- TRUE
      candidates$normal_ray_volume[idx] <- volume_result$volume
      candidates$normal_ray_volume_rel_error[idx] <- rel_error
      candidates$normal_ray_volume_success[idx] <- isTRUE(volume_result$volume_success)
      candidates$normal_ray_flip_this_part[idx] <- isTRUE(volume_result$flip_this_part)
      candidates$normal_ray_flip_part1[idx] <- isTRUE(volume_result$flip_part1)
      candidates$normal_ray_flip_part2[idx] <- isTRUE(volume_result$flip_part2)
      candidates$normal_ray_n_positive[idx] <- volume_result$normal_ray$direction_counts[1]
      candidates$normal_ray_n_negative[idx] <- volume_result$normal_ray$direction_counts[2]

      if (isTRUE(candidates$volume_error_le_threshold[idx])) {
        selected_candidate <- candidate_id
        break
      }
    }

    if (!is.finite(selected_candidate)) {
      selected_candidate <- candidates$candidate[order(candidates$proxy_energy_rank)][1]
      selection_rule <- paste0(
        "paired normal-ray fallback: no candidate passed volume_rel_error <= ",
        init_qc$threshold,
        "; selected minimum proxy energy"
      )
    } else {
      selection_rule <- paste0(
        "paired normal-ray first candidate by increasing proxy energy with volume_rel_error <= ",
        init_qc$threshold
      )
    }

    candidates$selected <- candidates$candidate == selected_candidate
    selected_summary <- candidates[candidates$selected, , drop = FALSE]

    init_qc$selection_rule <- selection_rule
    init_qc$selection_strategy <- "first_energy_passing_volume+paired_normal_ray"
    init_qc$selected_candidate <- selected_candidate
    init_qc$selected_initialization <- init_qc$candidate_initializations[[selected_candidate]]
    init_qc$selected_pme <- init_qc$candidate_pmes[[selected_candidate]]
    init_qc$selected_summary <- selected_summary
    init_qc$candidates <- candidates
    init_qc
  }

  part1 <- init_qc_part1
  part2 <- init_qc_part2
  for (pass in seq_len(max_passes)) {
    previous_selected <- c(part1$selected_candidate, part2$selected_candidate)
    part1 <- refine_one(part1, part2, "part1", volume_seed_offsets[["part1"]])
    part2 <- refine_one(part2, part1, "part2", volume_seed_offsets[["part2"]])
    if (identical(previous_selected, c(part1$selected_candidate, part2$selected_candidate))) {
      break
    }
  }

  candidates <- rbind(part1$candidates, part2$candidates)
  rownames(candidates) <- NULL
  list(part1 = part1, part2 = part2, candidates = candidates)
}
