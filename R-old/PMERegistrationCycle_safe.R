#' Safe PME + Registration Iterative Cycle
#'
#' Experimental variant of `PMERegistrationCycle` that calls
#' `Registration_safe`. Only accepted registration states are used for updating
#' the subject PME. If the first registration step is rejected, the cycle keeps
#' the initial subject PME and does not update/refit PME2.
#'
#' @docType class
#' @format An \code{R6Class} generator object
#' @export
PMERegistrationCycle_safe <- R6::R6Class(
  classname = "PMERegistrationCycle_safe",
  inherit = PMERegistrationCycle,

  public = list(
    #' @field init_strategy_f2 Subject initialization strategy.
    init_strategy_f2 = "isomap",

    #' @field refit_after_domain_stop Continue outer cycles after a domain stop.
    refit_after_domain_stop = TRUE,

    #' @field refit_after_energy_stop Continue outer cycles after an energy stop.
    refit_after_energy_stop = TRUE,

    #' @field best_state Best accepted outer-cycle state.
    best_state = NULL,

    #' @field energy_grid_mode Which grid is used to evaluate registration energy.
    energy_grid_mode = "same",

    #' @field optimization_grid_mode Which grid is used for basis/update steps.
    optimization_grid_mode = "provided",

    #' @field scale_mode How projected PME parameters are scaled before registration.
    scale_mode = "square",

    #' @field grid_n_u Resolution for generated disk/square optimization grids.
    grid_n_u = 60L,

    #' @field grid_n_v Resolution for generated disk/square optimization grids.
    grid_n_v = 60L,

    #' @field disk_compression Linear shrink factor used by `scale_mode = "disk"`.
    disk_compression = 1 / sqrt(2),

    #' @field f2_registered Registered subject function f2(gamma(u,v)), before PME refit.
    f2_registered = NULL,

    #' @field pme2_refit Subject PME object after refitting from f2_registered.
    pme2_refit = NULL,

    #' @field f2_refit Subject PME function after refitting. `f2_fun` is a legacy alias.
    f2_refit = NULL,

    #' @field f2_refit_grad Gradient of `f2_refit`.
    f2_refit_grad = NULL,

    initialize = function(...,
                          init_strategy_f2 = c("isomap",
                                               "template_projection_original"),
                          energy_grid_mode = c("same",
                                               "template_projected"),
                          optimization_grid_mode = c("provided",
                                                     "disk",
                                                     "square",
                                                     "template_projected"),
                          scale_mode = c("square", "disk"),
                          grid_n_u = 60L,
                          grid_n_v = 60L,
                          disk_compression = 1 / sqrt(2),
                          refit_after_domain_stop = TRUE,
                          refit_after_energy_stop = TRUE) {
      self$init_strategy_f2 <- match.arg(init_strategy_f2)
      self$energy_grid_mode <- match.arg(energy_grid_mode)
      self$optimization_grid_mode <- match.arg(optimization_grid_mode)
      self$scale_mode <- match.arg(scale_mode)
      self$grid_n_u <- as.integer(grid_n_u)
      self$grid_n_v <- as.integer(grid_n_v)
      self$disk_compression <- disk_compression
      self$refit_after_domain_stop <- refit_after_domain_stop
      self$refit_after_energy_stop <- refit_after_energy_stop
      self$best_state <- NULL
      super$initialize(...)
      invisible(self)
    },

    #' Restore active fields to the best accepted outer-cycle state.
    restore_best_state = function() {
      private$.restore_best_state()
      invisible(self)
    },

    #' Fit initial PME models, with optional template-projection f2 initialization.
    fit_initial = function() {
      if (self$init_strategy_f2 == "isomap") {
        return(super$fit_initial())
      }

      if (self$init_strategy_f2 != "template_projection_original") {
        stop("Unsupported init_strategy_f2: ", self$init_strategy_f2)
      }

      if (is.null(self$initialization_f1)) {
        self$initialization_f1 <- private$.init_pme(
          dataX = self$data1,
          init_args = self$init_args_f1,
          rescale = FALSE
        )
      }

      if (is.null(self$pme1)) {
        self$pme1 <- private$.fit_pme(
          dataX = self$data1,
          pme_args = self$pme_args_f1,
          initialization = self$initialization_f1
        )
      }

      if (is.null(self$initialization_f2)) {
        self$initialization_f2 <- private$.init_pme(
          dataX = self$data2,
          init_args = self$init_args_f2,
          rescale = FALSE
        )

        template_centers <- as.matrix(self$initialization_f1$centers)
        template_params <- as.matrix(self$pme1$params_opt)
        subject_centers <- as.matrix(self$initialization_f2$centers)

        if (nrow(template_centers) != nrow(template_params)) {
          stop(
            "Template center count does not match template PME parameter count."
          )
        }

        nearest_idx <- vapply(
          seq_len(nrow(subject_centers)),
          function(i) {
            center_i <- matrix(
              subject_centers[i, ],
              nrow = nrow(template_centers),
              ncol = ncol(template_centers),
              byrow = TRUE
            )
            which.min(rowSums((template_centers - center_i)^2))
          },
          integer(1)
        )

        init_guess <- template_params[nearest_idx, , drop = FALSE]
        projected_params <- calc_params(
          f = self$pme1$embedding_map,
          X = subject_centers,
          init_params = init_guess,
          f_input = "vector"
        )

        self$initialization_f2$parameterization <- projected_params
      }

      if (is.null(self$pme2)) {
        self$pme2 <- private$.fit_pme(
          dataX = self$data2,
          pme_args = self$pme_args_f2,
          initialization = self$initialization_f2
        )
      }

      self$scale_f1 <- private$.compute_projection_and_scale(self$pme1, self$data1)
      self$f1_fun <- private$.make_scaled_embedding(
        self$pme1,
        self$scale_f1,
        d = self$d
      )

      self$scale_f2 <- private$.compute_projection_and_scale(self$pme2, self$data2)
      self$f2_fun <- private$.make_scaled_embedding(
        self$pme2,
        self$scale_f2,
        d = self$d
      )
      self$f2_grad <- private$.make_scaled_grad(self$pme2, self$scale_f2)
      self$pme2_refit <- self$pme2
      self$f2_refit <- self$f2_fun
      self$f2_refit_grad <- self$f2_grad

      self$history$initial <- list(
        k = 0,
        stage = "init",
        init_strategy_f2 = self$init_strategy_f2,
        pme1 = self$pme1,
        pme2 = self$pme2,
        initialization_f1 = self$initialization_f1,
        initialization_f2 = self$initialization_f2,
        f1_fun = self$f1_fun,
        f2_fun = self$f2_fun,
        scale_f1 = self$scale_f1,
        scale_f2 = self$scale_f2,
        energy_grid_mode = self$energy_grid_mode,
        optimization_grid_mode = self$optimization_grid_mode,
        scale_mode = self$scale_mode
      )

      invisible(self)
    },

    #' Run one safe alternating registration cycle.
    run_cycle = function() {
      if (is.null(self$f1_fun) || is.null(self$f2_fun)) {
        self$fit_initial()
      }

      k <- self$cycle_idx + 1L

      self$reg <- private$.register_once(
        f1_fun = self$f1_fun,
        f2_fun = self$f2_fun,
        grad_f2_fun = self$f2_grad
      )

      n_accepted_states <- length(self$reg$state_list)
      no_accepted_update <- n_accepted_states <= 1L

      if (no_accepted_update) {
        last_state <- self$reg$state_list[[1]]
        # `f2_registered` is the registration output before PME refit:
        # f2_registered(u,v) = f2_previous(gamma(u,v)). With no accepted update
        # this is the identity-warped current f2.
        self$f2_warped <- last_state$f2_k
        self$f2_registered <- self$f2_warped
        self$pme2_refit <- self$pme2
        self$f2_refit <- self$f2_fun
        self$f2_refit_grad <- self$f2_grad
        lambda_to_use <- private$.resolve_lambda_for_f2()

      private$.record_history(
        k = k,
        reg = self$reg,
        new_pme2 = self$pme2,
        scale_f2 = self$scale_f2,
        f2_fun = self$f2_fun,
        lambda = lambda_to_use
      )
      private$.annotate_latest_cycle_outputs()
      self$history$cycles[[length(self$history$cycles)]]$safe_domain_refit_triggered <- FALSE
      self$history$cycles[[length(self$history$cycles)]]$safe_energy_refit_triggered <- FALSE
      self$history$cycles[[length(self$history$cycles)]]$safe_refit_triggered <- FALSE

      self$cycle_idx <- k

        if (self$verbose) {
          message(
            "[PMERegistrationCycle_safe] No accepted registration update; ",
            "kept current PME2 without refit."
          )
        }

        return(invisible(self))
      }

      last_state <- private$.extract_last_state(self$reg)
      # This is the real registration output for the current cycle, before PME
      # compression/refit. It is often the best SIME baseline when refit PME
      # loses the registered geometry.
      self$f2_warped <- last_state$f2_k
      self$f2_registered <- self$f2_warped

      private$.update_f2_initialization_from_f2k_inplace(
        f2k_fun = self$f2_warped,
        scale_f2 = self$scale_f2,
        center_points = self$initialization_f2$centers,
        init_params_for_centers = self$pme2$params_opt
      )
      updated_init <- self$initialization_f2

      lambda_to_use <- private$.resolve_lambda_for_f2()

      self$pme2 <- private$.fit_pme(
        dataX = self$data2,
        pme_args = self$pme_args_f2,
        initialization = updated_init,
        lambda = lambda_to_use
      )

      self$scale_f2 <- private$.compute_projection_and_scale(self$pme2, self$data2)
      self$f2_fun <- private$.make_scaled_embedding(
        self$pme2,
        self$scale_f2,
        d = self$d
      )
      self$f2_grad <- private$.make_scaled_grad(self$pme2, self$scale_f2)
      # `f2_fun` is kept for backwards compatibility with the original
      # PMERegistrationCycle interface. In the safe code, prefer the explicit
      # names below:
      #   f2_registered = f2(gamma(u,v)), before PME refit
      #   f2_refit      = PME refit after registration
      self$pme2_refit <- self$pme2
      self$f2_refit <- self$f2_fun
      self$f2_refit_grad <- self$f2_grad

      private$.record_history(
        k = k,
        reg = self$reg,
        new_pme2 = self$pme2,
        scale_f2 = self$scale_f2,
        f2_fun = self$f2_fun,
        lambda = lambda_to_use
      )
      private$.annotate_latest_cycle_outputs()
      self$history$cycles[[length(self$history$cycles)]]$safe_domain_refit_triggered <-
        isTRUE(self$refit_after_domain_stop) &&
        identical(self$reg$stop_reason, "rejected domain violation")
      self$history$cycles[[length(self$history$cycles)]]$safe_energy_refit_triggered <-
        isTRUE(self$refit_after_energy_stop) &&
        identical(self$reg$stop_reason, "rejected energy increase")
      self$history$cycles[[length(self$history$cycles)]]$safe_refit_triggered <-
        isTRUE(self$history$cycles[[length(self$history$cycles)]]$safe_domain_refit_triggered) ||
        isTRUE(self$history$cycles[[length(self$history$cycles)]]$safe_energy_refit_triggered)

      self$cycle_idx <- k
      invisible(self)
    },

    #' Run multiple safe alternating cycles.
    run = function(n_cycles = 5,
                   save_dir = NULL,
                   filename = NULL,
                   reinit = FALSE) {
      completed <- length(self$history$cycles)
      start_cycle <- completed + 1L
      end_cycle <- completed + as.integer(n_cycles)

      if (self$verbose) {
        ra <- private$.null_to_list(self$reg_args)
        cat("========================================\n")
        cat("PMERegistrationCycle_safe START\n")
        cat(sprintf("already_done : %d\n", as.integer(completed)))
        cat(sprintf("n_cycles     : %d\n", as.integer(n_cycles)))
        cat(sprintf("will_run     : %d -> %d (overall)\n",
                    as.integer(start_cycle), as.integer(end_cycle)))
        cat(sprintf("d            : %d\n", as.integer(self$d)))
        cat(sprintf("eps_step     : %s\n", as.character(ra$eps_step)))
        cat(sprintf("eps_energy   : %s\n", as.character(ra$eps_energy)))
        cat(sprintf("max_iter     : %s\n", as.character(ra$max_iter)))
        cat("========================================\n")
        flush.console()
      }

      if (!is.null(save_dir)) self$save_dir <- save_dir
      if (!is.null(filename)) self$filename <- filename

      if (reinit) {
        private$.reset_state()
        self$best_state <- NULL
        self$f2_registered <- NULL
        self$pme2_refit <- NULL
        self$f2_refit <- NULL
        self$f2_refit_grad <- NULL
      }

      if (is.null(self$pme1) || is.null(self$pme2)) {
        self$fit_initial()
      }

      self$converged <- FALSE
      self$stop_reason <- NULL

      for (i in seq_len(n_cycles)) {
        if (self$verbose) {
          overall_i <- start_cycle + i - 1L
          cat(sprintf(
            "\n[PMEReg safe] cycle %d / %d (this run %d / %d)\n",
            as.integer(overall_i), as.integer(end_cycle),
            as.integer(i), as.integer(n_cycles)
          ))
        }

        self$run_cycle()

        no_accepted_update <- length(self$reg$state_list) <= 1L
        if (!no_accepted_update) {
          current_is_best <- private$.update_best_state_from_latest_cycle()

          if (!isTRUE(current_is_best)) {
            private$.restore_best_state()
            self$converged <- TRUE
            self$stop_reason <- "cycle energy increased; restored best state"

            if (!is.null(self$save_dir) && !is.null(self$filename)) {
              private$.save_all_overwrite(
                save_dir = self$save_dir,
                filename = self$filename
              )
            }

            break
          }
        }

        if (!is.null(self$save_dir) && !is.null(self$filename)) {
          private$.save_all_overwrite(
            save_dir = self$save_dir,
            filename = self$filename
          )
        }

        if (no_accepted_update) {
          private$.restore_best_state()
          self$converged <- TRUE
          self$stop_reason <- "no accepted registration update"
          break
        }

        reg_stopped_early <- length(self$reg$E_history) < self$reg$max_iter + 1L
        domain_stop_with_update <- identical(
          self$reg$stop_reason,
          "rejected domain violation"
        ) && length(self$reg$state_list) > 1L
        energy_stop_with_update <- identical(
          self$reg$stop_reason,
          "rejected energy increase"
        ) && length(self$reg$state_list) > 1L

        if (reg_stopped_early &&
            isTRUE(self$refit_after_domain_stop) &&
            domain_stop_with_update) {
          if (self$verbose) {
            message(
              "[PMERegistrationCycle_safe] Registration reached the domain ",
              "boundary after accepted updates; refit PME2 and continue."
            )
          }
          next
        }

        if (reg_stopped_early &&
            isTRUE(self$refit_after_energy_stop) &&
            energy_stop_with_update) {
          if (self$verbose) {
            message(
              "[PMERegistrationCycle_safe] Registration could not find a ",
              "further energy-decreasing update after accepted updates; ",
              "refit PME2 and continue."
            )
          }
          next
        }

        if (reg_stopped_early) {
          self$converged <- TRUE
          self$stop_reason <- self$reg$stop_reason
          break
        }
      }

      invisible(self)
    }
  ),

  private = list(
    .register_once = function(f1_fun, f2_fun, grad_f2_fun) {
      args <- modifyList(
        private$.null_to_list(self$reg_args),
        list(
          f1 = f1_fun,
          f2 = f2_fun,
          f2_grad_fn = grad_f2_fun,
          verbose = self$verbose
        )
      )
      args <- private$.apply_registration_grid_modes(args)

      reg <- do.call(Registration_safe$new, args)
      reg$run()
      reg
    },

    .apply_registration_grid_modes = function(args) {
      # The safe workflow can now decouple two roles that used to share Ugrid:
      # 1. optimization/update grid: where basis inner products define gamma.
      # 2. energy grid: where the accept/reject objective is evaluated.
      #
      # Defaults preserve the original behavior: use the provided Ugrid for both.
      if (self$optimization_grid_mode == "disk") {
        args$Ugrid <- make_uv_grid(
          self$grid_n_u,
          self$grid_n_v,
          grid_type = "disk"
        )
      } else if (self$optimization_grid_mode == "square") {
        args$Ugrid <- make_uv_grid(
          self$grid_n_u,
          self$grid_n_v,
          grid_type = "square"
        )
      } else if (self$optimization_grid_mode == "template_projected") {
        args$Ugrid <- private$.as_uv_grid(self$scale_f1$U_scaled)
      }

      if (self$energy_grid_mode == "template_projected") {
        args$energy_Ugrid <- private$.as_uv_grid(self$scale_f1$U_scaled)
        args$energy_grid_source <- "template_projected"
      }

      args
    },

    .as_uv_grid = function(U) {
      U <- as.data.frame(U)
      if (ncol(U) < 2L) {
        stop("Grid must have at least two columns.")
      }
      U <- U[, 1:2, drop = FALSE]
      names(U) <- c("u", "v")
      U
    },

    .compute_projection_and_scale = function(pme_result,
                                             dataX,
                                             init_guess_method = "pca") {
      init_param <- pme_initial_guess(
        X = dataX,
        d = self$d,
        method = init_guess_method
      )
      U_proj <- calc_params(
        f = pme_result$embedding_map,
        X = dataX,
        init_params = init_param,
        f_input = "vector"
      )

      scaled_square <- scale_uniform_square_with_params(U_proj)
      scaled <- if (self$scale_mode == "disk") {
        private$.compress_square_scale_to_disk(
          scaled_square,
          alpha = self$disk_compression
        )
      } else {
        scaled_square
      }

      list(
        U_proj = U_proj,
        U_scaled = scaled$U_scaled,
        A = scaled$A,
        b = scaled$b,
        scale_mode = self$scale_mode,
        disk_compression = if (self$scale_mode == "disk") {
          self$disk_compression
        } else {
          NULL
        }
      )
    },

    .compress_square_scale_to_disk = function(square_scale,
                                              alpha = 1 / sqrt(2)) {
      if (!is.finite(alpha) || alpha <= 0 || alpha > 1) {
        stop("disk compression alpha must be in (0, 1].")
      }

      center <- c(0.5, 0.5)
      U_square <- as.matrix(square_scale$U_scaled)
      U_scaled <- sweep(U_square, 2L, center, "-") * alpha + 0.5

      # If s_disk = 0.5 + alpha * (s_square - 0.5), then
      # s_square = (s_disk - 0.5) / alpha + 0.5. Compose this inverse with the
      # original square inverse so identity-energy comparisons stay unchanged
      # after applying the same compression to template and subject functions.
      A <- square_scale$A / alpha
      b <- as.numeric(
        square_scale$A %*% (center - center / alpha) + square_scale$b
      )

      list(
        U_scaled = U_scaled,
        A = A,
        b = b,
        alpha = alpha
      )
    },

    .capture_active_state = function(cycle_entry) {
      list(
        source = "cycle",
        cycle = cycle_entry$k,
        final_E = cycle_entry$final_E,
        reg_E_history = cycle_entry$reg_E_history,
        reg_stop_reason = cycle_entry$reg$stop_reason,
        safe_domain_refit_triggered = isTRUE(
          cycle_entry$safe_domain_refit_triggered
        ),
        safe_energy_refit_triggered = isTRUE(
          cycle_entry$safe_energy_refit_triggered
        ),
        safe_refit_triggered = isTRUE(cycle_entry$safe_refit_triggered),
        pme2 = self$pme2,
        scale_f2 = self$scale_f2,
        pme2_refit = self$pme2_refit,
        f2_refit = self$f2_refit,
        f2_refit_grad = self$f2_refit_grad,
        f2_registered = self$f2_registered,
        f2_fun = self$f2_fun,
        f2_grad = self$f2_grad,
        f2_warped = self$f2_warped,
        initialization_f2 = self$initialization_f2,
        reg = self$reg
      )
    },

    .restore_best_state = function() {
      if (is.null(self$best_state)) {
        return(invisible(FALSE))
      }

      self$pme2 <- self$best_state$pme2
      self$scale_f2 <- self$best_state$scale_f2
      self$pme2_refit <- if (!is.null(self$best_state$pme2_refit)) {
        self$best_state$pme2_refit
      } else {
        self$best_state$pme2
      }
      self$f2_refit <- if (!is.null(self$best_state$f2_refit)) {
        self$best_state$f2_refit
      } else {
        self$best_state$f2_fun
      }
      self$f2_refit_grad <- if (!is.null(self$best_state$f2_refit_grad)) {
        self$best_state$f2_refit_grad
      } else {
        self$best_state$f2_grad
      }
      self$f2_registered <- if (!is.null(self$best_state$f2_registered)) {
        self$best_state$f2_registered
      } else {
        self$best_state$f2_warped
      }
      self$f2_fun <- self$best_state$f2_fun
      self$f2_grad <- self$best_state$f2_grad
      self$f2_warped <- self$best_state$f2_warped
      self$initialization_f2 <- self$best_state$initialization_f2
      self$reg <- self$best_state$reg

      invisible(TRUE)
    },

    .annotate_latest_cycle_outputs = function() {
      cycle_idx <- length(self$history$cycles)
      if (cycle_idx == 0L) {
        return(invisible(FALSE))
      }

      cycle_entry <- self$history$cycles[[cycle_idx]]
      cycle_entry$pme2_refit <- self$pme2_refit
      cycle_entry$f2_refit <- self$f2_refit
      cycle_entry$f2_refit_grad <- self$f2_refit_grad
      cycle_entry$f2_registered <- self$f2_registered
      cycle_entry$f2_warped <- self$f2_registered
      cycle_entry$f2_fun_legacy_alias <- self$f2_fun
      cycle_entry$output_name_note <- paste(
        "f2_registered is f2(gamma(u,v)) before PME refit;",
        "f2_refit is the PME refit after registration;",
        "f2_fun is kept as a legacy alias for f2_refit."
      )
      self$history$cycles[[cycle_idx]] <- cycle_entry

      invisible(TRUE)
    },

    .update_best_state_from_latest_cycle = function() {
      cycle_idx <- length(self$history$cycles)
      if (cycle_idx == 0L) {
        return(FALSE)
      }

      current_cycle <- self$history$cycles[[cycle_idx]]
      current_E <- current_cycle$final_E
      best_E <- if (is.null(self$best_state)) Inf else self$best_state$final_E
      current_is_best <- is.finite(current_E) && current_E <= best_E

      current_cycle$is_best_after_cycle <- current_is_best
      current_cycle$is_rejected_by_outer_energy <- !current_is_best
      self$history$cycles[[cycle_idx]] <- current_cycle

      if (current_is_best) {
        self$best_state <- private$.capture_active_state(current_cycle)
      }

      current_is_best
    }
  )
)
