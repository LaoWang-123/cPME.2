#' PME + Registration iterative cycle
#'
#' Current cPME alternating workflow. This class fits initial PME models,
#' registers the moving surface to the fixed surface with safe accepted-state
#' registration, refits the moving PME after accepted registration updates, and
#' keeps the best outer-cycle state as the active output.
#'
#' @docType class
#' @format An \code{R6Class} generator object
#' @export
PMERegistrationCycle <- R6::R6Class(
  classname = "PMERegistrationCycle",

  public = list(
    #' @field data1 Fixed/template data matrix.
    data1 = NULL,

    #' @field data2 Moving/subject data matrix.
    data2 = NULL,

    #' @field d Intrinsic dimension of the PME parameter domain.
    d = NULL,

    #' @field pme_args_f1 Arguments passed to `pme()` for `data1`.
    pme_args_f1 = NULL,

    #' @field pme_args_f2 Arguments passed to `pme()` for `data2`.
    pme_args_f2 = NULL,

    #' @field init_args_f1 Arguments passed to `initialize_pme()` for `data1`.
    init_args_f1 = NULL,

    #' @field init_args_f2 Arguments passed to `initialize_pme()` for `data2`.
    init_args_f2 = NULL,

    #' @field init_trials Number of subject initialization trials.
    init_trials = NULL,

    #' @field reg_args Arguments passed to `Registration$new()`.
    reg_args = NULL,

    #' @field lambda_policy_f2 Lambda policy for subject PME refits.
    lambda_policy_f2 = NULL,

    #' @field fixed_lambda_f2 Fixed lambda used when `lambda_policy_f2 = "fixed"`.
    fixed_lambda_f2 = NULL,

    #' @field cycle_idx Current outer-cycle index.
    cycle_idx = 0L,

    #' @field pme1_initial Initial fixed/template PME object.
    pme1_initial = NULL,

    #' @field pme2_initial Initial moving/subject PME object before registration.
    pme2_initial = NULL,

    #' @field f1_initial Wrapped/scaled fixed/template PME function.
    f1_initial = NULL,

    #' @field f2_initial Wrapped/scaled initial moving/subject PME function.
    f2_initial = NULL,

    #' @field f2_initial_grad Gradient of `f2_initial`.
    f2_initial_grad = NULL,

    #' @field initialization_f1 Cached PME initialization object for `data1`.
    initialization_f1 = NULL,

    #' @field initialization_f2 Cached PME initialization object for `data2`.
    initialization_f2 = NULL,

    #' @field scale_f1 Parameter scaling information for `pme1_initial`.
    scale_f1 = NULL,

    #' @field scale_f2 Parameter scaling information for the active subject PME.
    scale_f2 = NULL,

    #' @field reg Most recent `Registration` object.
    reg = NULL,

    #' @field converged Logical convergence flag.
    converged = FALSE,

    #' @field stop_reason Text description of the stopping condition.
    stop_reason = NULL,

    #' @field history Stored initial state and cycle states.
    history = NULL,

    #' @field init_history Optional initialization history.
    init_history = NULL,

    #' @field save_dir Directory for autosaving full objects.
    save_dir = NULL,

    #' @field filename File name for autosaving full objects.
    filename = NULL,

    #' @field verbose Logical; print progress messages.
    verbose = TRUE,

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

    #' @field pme2_refit Active moving PME object after the latest accepted registration/refit.
    pme2_refit = NULL,

    #' @field f2_refit Wrapped/scaled function for `pme2_refit`.
    f2_refit = NULL,

    #' @field f2_refit_grad Gradient of `f2_refit`.
    f2_refit_grad = NULL,

    #' @description
    #' Create a PME registration cycle object.
    #' @param data1 Fixed/template data matrix.
    #' @param data2 Moving/subject data matrix.
    #' @param pme1 Optional pre-fitted fixed/template PME object. Stored as
    #' `pme1_initial`.
    #' @param pme2 Optional pre-fitted moving/subject PME object. Stored as
    #' `pme2_initial`.
    #' @param initialization_f1 Optional cached initialization for `data1`.
    #' @param initialization_f2 Optional cached initialization for `data2`.
    #' @param d Intrinsic dimension of the PME parameter domain.
    #' @param pme_args_f1 Arguments passed to `pme()` for `data1`.
    #' @param pme_args_f2 Arguments passed to `pme()` for `data2`.
    #' @param init_args_f1 Arguments passed to `initialize_pme()` for `data1`.
    #' @param init_args_f2 Arguments passed to `initialize_pme()` for `data2`.
    #' @param init_trials Number of subject initialization trials.
    #' @param reg_args Arguments passed to `Registration$new()`.
    #' @param lambda_policy_f2 Lambda policy for subject PME refits:
    #' `"reuse_prev"`, `"fixed"`, or `"retune"`.
    #' @param fixed_lambda_f2 Fixed lambda used when
    #' `lambda_policy_f2 = "fixed"`.
    #' @param default_args Optional list of default PME, initialization, and
    #' registration arguments.
    #' @param save_dir Optional directory for autosaving the full object.
    #' @param filename Optional filename for autosaving the full object.
    #' @param verbose Logical; print progress messages.
    #' @param init_strategy_f2 Initialization strategy for the moving/subject
    #' PME. `"isomap"` uses the ordinary subject initialization;
    #' `"template_projection_original"` initializes subject center parameters
    #' from nearest template projected parameters.
    #' @param energy_grid_mode Grid used to evaluate registration energy.
    #' `"template_projected"` uses data-projected template parameters;
    #' `"same"` uses the registration optimization grid.
    #' @param optimization_grid_mode Grid used for the registration basis/update
    #' calculation. The default `"template_projected"` is the current recommended
    #' real-data workflow.
    #' @param scale_mode Parameter scaling mode before registration. `"square"`
    #' maps projected parameters to `[0,1]^2`; `"disk"` linearly compresses that
    #' square into the inscribed disk.
    #' @param grid_n_u Number of grid values used when generating square/disk
    #' grids.
    #' @param grid_n_v Number of grid values used when generating square/disk
    #' grids.
    #' @param disk_compression Linear compression factor used by
    #' `scale_mode = "disk"`.
    #' @param refit_after_domain_stop Logical; if accepted registration steps are
    #' followed by a domain-violation proposal, refit PME2 from the last accepted
    #' registered function and continue outer cycles.
    #' @param refit_after_energy_stop Logical; if accepted registration steps are
    #' followed by an energy-increasing proposal, refit PME2 from the last
    #' accepted registered function and continue outer cycles.
    #' @return A new `PMERegistrationCycle` object.
    initialize = function(data1,
                          data2,
                          pme1 = NULL,
                          pme2 = NULL,
                          initialization_f1 = NULL,
                          initialization_f2 = NULL,
                          d = 2,
                          pme_args_f1 = list(),
                          pme_args_f2 = list(),
                          init_args_f1 = list(),
                          init_args_f2 = list(),
                          init_trials = 5,
                          reg_args = list(),
                          lambda_policy_f2 = c("reuse_prev",
                                               "fixed",
                                               "retune"),
                          fixed_lambda_f2 = NULL,
                          default_args = NULL,
                          save_dir = NULL,
                          filename = NULL,
                          verbose = TRUE,
                          init_strategy_f2 = c("isomap",
                                               "template_projection_original"),
                          energy_grid_mode = c("template_projected",
                                               "same"),
                          optimization_grid_mode = c("template_projected",
                                                     "provided",
                                                     "disk",
                                                     "square"),
                          scale_mode = c("square", "disk"),
                          grid_n_u = 60L,
                          grid_n_v = 60L,
                          disk_compression = 1 / sqrt(2),
                          refit_after_domain_stop = TRUE,
                          refit_after_energy_stop = TRUE) {
      self$data1 <- data1
      self$data2 <- data2
      self$pme1_initial <- pme1
      self$pme2_initial <- pme2
      self$initialization_f1 <- initialization_f1
      self$initialization_f2 <- initialization_f2
      self$d <- as.integer(d)
      self$init_trials <- init_trials

      defaults <- private$.get_default_args(default_args)
      self$pme_args_f1 <- modifyList(
        private$.null_to_list(defaults$pme_args_f1),
        private$.null_to_list(pme_args_f1),
        keep.null = TRUE
      )
      self$pme_args_f2 <- modifyList(
        private$.null_to_list(defaults$pme_args_f2),
        private$.null_to_list(pme_args_f2),
        keep.null = TRUE
      )
      self$init_args_f1 <- modifyList(
        private$.null_to_list(defaults$init_args_f1),
        private$.null_to_list(init_args_f1),
        keep.null = TRUE
      )
      self$init_args_f2 <- modifyList(
        private$.null_to_list(defaults$init_args_f2),
        private$.null_to_list(init_args_f2),
        keep.null = TRUE
      )
      self$reg_args <- modifyList(
        private$.null_to_list(defaults$reg_args),
        private$.null_to_list(reg_args),
        keep.null = TRUE
      )
      self$lambda_policy_f2 <- match.arg(lambda_policy_f2)
      self$fixed_lambda_f2 <- fixed_lambda_f2
      self$save_dir <- save_dir
      self$filename <- filename
      self$verbose <- isTRUE(verbose)

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
      self$history <- list(initial = NULL, cycles = list())
      self$init_history <- list()
      invisible(self)
    },

    #' @description
    #' Restore active fields to the best accepted outer-cycle state.
    #' @return The object, invisibly.
    restore_best_state = function() {
      private$.restore_best_state()
      invisible(self)
    },

    #' @description
    #' Fit initial PME models, with optional template-projection f2 initialization.
    #' @return The object, invisibly.
    fit_initial = function() {
      compute_E_from_inits <- function(pme1, pme2, init1, init2, n_grid = 10L) {
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

        vals1 <- t(apply(U_eval1, 1L, function(u) {
          pme1$embedding_map(as.numeric(u))
        }))
        vals2 <- t(apply(U_eval2, 1L, function(u) {
          pme2$embedding_map(as.numeric(u))
        }))

        mean(rowSums((vals1 - vals2)^2))
      }

      if (is.null(self$initialization_f1)) {
        self$initialization_f1 <- private$.init_pme(
          dataX = self$data1,
          init_args = self$init_args_f1,
          rescale = FALSE
        )
      }

      if (is.null(self$pme1_initial)) {
        self$pme1_initial <- private$.fit_pme(
          dataX = self$data1,
          pme_args = self$pme_args_f1,
          initialization = self$initialization_f1
        )
      }

      if (self$init_strategy_f2 == "template_projection_original") {
        if (is.null(self$initialization_f2)) {
          self$initialization_f2 <- private$.init_pme(
            dataX = self$data2,
            init_args = self$init_args_f2,
            rescale = FALSE
          )
        }

        template_centers <- as.matrix(self$initialization_f1$centers)
        template_params <- as.matrix(self$pme1_initial$params_opt)
        subject_centers <- as.matrix(self$initialization_f2$centers)

        if (nrow(template_centers) != nrow(template_params)) {
          stop("Template center count does not match template PME parameter count.")
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
        self$initialization_f2$parameterization <- calc_params(
          f = self$pme1_initial$embedding_map,
          X = subject_centers,
          init_params = init_guess,
          f_input = "vector"
        )
      } else if (self$init_strategy_f2 == "isomap" &&
                 is.null(self$pme2_initial)) {
        if (is.null(self$initialization_f2)) {
          max_tries <- max(1L, as.integer(self$init_trials))
          cand_inits <- vector("list", max_tries)
          cand_checks <- vector("list", max_tries)

          for (i in seq_len(max_tries)) {
            init2_i <- private$.init_pme(
              dataX = self$data2,
              init_args = self$init_args_f2,
              rescale = FALSE
            )
            chk_i <- check_pme_orientation(
              init1 = self$initialization_f1,
              init2 = init2_i,
              pca_source = "all_centers",
              verbose = self$verbose
            )

            cand_inits[[i]] <- init2_i
            cand_checks[[i]] <- chk_i
          }

          keep <- which(vapply(
            cand_checks,
            function(chk) !identical(chk$final, "mirror_reversed"),
            logical(1)
          ))

          if (length(keep) == 0L) {
            stop(sprintf(
              "[fit_initial] all f2 initializations are mirror_reversed across %d trials.",
              max_tries
            ))
          }

          cand_inits <- cand_inits[keep]
          cand_checks <- cand_checks[keep]

          if (self$verbose) {
            message(sprintf(
              "[fit_initial] init pool %d -> %d after orientation filter.",
              max_tries, length(cand_inits)
            ))
          }
        } else {
          cand_inits <- list(self$initialization_f2)
          cand_checks <- list(list(final = "user_provided"))
        }

        n_cand <- length(cand_inits)
        cand_pme2 <- vector("list", n_cand)
        cand_E <- rep(Inf, n_cand)

        cand_pme2[[1]] <- private$.fit_pme(
          dataX = self$data2,
          pme_args = self$pme_args_f2,
          initialization = cand_inits[[1]]
        )
        cand_E[1] <- compute_E_from_inits(
          self$pme1_initial,
          cand_pme2[[1]],
          self$initialization_f1,
          cand_inits[[1]]
        )

        lambda_reuse <- cand_pme2[[1]]$tuning
        if (n_cand >= 2L) {
          for (i in 2:n_cand) {
            cand_pme2[[i]] <- private$.fit_pme(
              dataX = self$data2,
              pme_args = self$pme_args_f2,
              initialization = cand_inits[[i]],
              lambda = lambda_reuse
            )
            cand_E[i] <- compute_E_from_inits(
              self$pme1_initial,
              cand_pme2[[i]],
              self$initialization_f1,
              cand_inits[[i]]
            )
          }
        }

        best_idx <- which.min(cand_E)
        self$initialization_f2 <- cand_inits[[best_idx]]
        self$pme2_initial <- cand_pme2[[best_idx]]
        self$init_history$fit_initial <- list(
          init_strategy_f2 = self$init_strategy_f2,
          n_trials_requested = if (exists("max_tries")) max_tries else n_cand,
          n_candidates = n_cand,
          best_idx = best_idx,
          candidate_E = cand_E,
          candidate_checks = cand_checks,
          lambda_reuse = lambda_reuse
        )

        if (self$verbose) {
          message(sprintf(
            "[fit_initial] selected candidate %d/%d with E=%.6g",
            best_idx, n_cand, cand_E[best_idx]
          ))
        }
      }

      if (is.null(self$pme2_initial)) {
        self$pme2_initial <- private$.fit_pme(
          dataX = self$data2,
          pme_args = self$pme_args_f2,
          initialization = self$initialization_f2
        )
      }

      self$scale_f1 <- private$.compute_projection_and_scale(
        self$pme1_initial,
        self$data1
      )
      self$f1_initial <- private$.make_scaled_embedding(
        self$pme1_initial,
        self$scale_f1,
        d = self$d
      )

      self$scale_f2 <- private$.compute_projection_and_scale(
        self$pme2_initial,
        self$data2
      )
      self$f2_initial <- private$.make_scaled_embedding(
        self$pme2_initial,
        self$scale_f2,
        d = self$d
      )
      self$f2_initial_grad <- private$.make_scaled_grad(
        self$pme2_initial,
        self$scale_f2
      )

      self$pme2_refit <- self$pme2_initial
      self$f2_refit <- self$f2_initial
      self$f2_refit_grad <- self$f2_initial_grad
      self$f2_registered <- NULL

      self$history$initial <- list(
        k = 0,
        stage = "init",
        init_strategy_f2 = self$init_strategy_f2,
        pme1_initial = self$pme1_initial,
        pme2_initial = self$pme2_initial,
        initialization_f1 = self$initialization_f1,
        initialization_f2 = self$initialization_f2,
        f1_initial = self$f1_initial,
        f2_initial = self$f2_initial,
        scale_f1 = self$scale_f1,
        scale_f2 = self$scale_f2,
        energy_grid_mode = self$energy_grid_mode,
        optimization_grid_mode = self$optimization_grid_mode,
        scale_mode = self$scale_mode
      )

      invisible(self)
    },
    #' @description
    #' Run one safe alternating registration cycle.
    #' @return The object, invisibly.
    run_cycle = function() {
      if (is.null(self$f1_initial) || is.null(self$f2_refit)) {
        self$fit_initial()
      }

      k <- self$cycle_idx + 1L

      self$reg <- private$.register_once(
        f1_fun = self$f1_initial,
        f2_fun = self$f2_refit,
        grad_f2_fun = self$f2_refit_grad
      )

      n_accepted_states <- length(self$reg$state_list)
      no_accepted_update <- n_accepted_states <= 1L
      last_state <- self$reg$state_list[[length(self$reg$state_list)]]
      self$f2_registered <- last_state$f2_k
      lambda_to_use <- private$.resolve_lambda_for_f2()

      if (no_accepted_update) {
        private$.record_history(
          k = k,
          reg = self$reg,
          pme2_refit = self$pme2_refit,
          scale_f2 = self$scale_f2,
          f2_refit = self$f2_refit,
          lambda = lambda_to_use
        )
        private$.annotate_latest_cycle_outputs()
        self$history$cycles[[length(self$history$cycles)]]$safe_domain_refit_triggered <- FALSE
        self$history$cycles[[length(self$history$cycles)]]$safe_energy_refit_triggered <- FALSE
        self$history$cycles[[length(self$history$cycles)]]$safe_refit_triggered <- FALSE

        self$cycle_idx <- k

        if (self$verbose) {
          message(
            "[PMERegistrationCycle] No accepted registration update; ",
            "kept current pme2_refit without another refit."
          )
        }

        return(invisible(self))
      }

      private$.update_f2_initialization_from_f2k_inplace(
        f2k_fun = self$f2_registered,
        scale_f2 = self$scale_f2,
        center_points = self$initialization_f2$centers,
        init_params_for_centers = self$pme2_refit$params_opt
      )
      updated_init <- self$initialization_f2

      self$pme2_refit <- private$.fit_pme(
        dataX = self$data2,
        pme_args = self$pme_args_f2,
        initialization = updated_init,
        lambda = lambda_to_use
      )

      self$scale_f2 <- private$.compute_projection_and_scale(
        self$pme2_refit,
        self$data2
      )
      self$f2_refit <- private$.make_scaled_embedding(
        self$pme2_refit,
        self$scale_f2,
        d = self$d
      )
      self$f2_refit_grad <- private$.make_scaled_grad(
        self$pme2_refit,
        self$scale_f2
      )

      private$.record_history(
        k = k,
        reg = self$reg,
        pme2_refit = self$pme2_refit,
        scale_f2 = self$scale_f2,
        f2_refit = self$f2_refit,
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
    #' @description
    #' Run multiple safe alternating cycles.
    #' @param n_cycles Number of outer PME/registration cycles to run.
    #' @param save_dir Optional directory for autosaving the full object.
    #' @param filename Optional RDS filename for autosaving the full object.
    #' @param reinit Logical; when `TRUE`, clear fitted state and rerun from
    #' initialization.
    #' @return The object, invisibly.
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
        cat("PMERegistrationCycle START\n")
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

      if (is.null(self$pme1_initial) || is.null(self$pme2_initial) ||
          is.null(self$pme2_refit)) {
        self$fit_initial()
      }

      self$converged <- FALSE
      self$stop_reason <- NULL

      for (i in seq_len(n_cycles)) {
        if (self$verbose) {
          overall_i <- start_cycle + i - 1L
          cat(sprintf(
            "\n[PMEReg] cycle %d / %d (this run %d / %d)\n",
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
              "[PMERegistrationCycle] Registration reached the domain ",
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
              "[PMERegistrationCycle] Registration could not find a ",
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
    .null_to_list = function(x) {
      if (is.null(x)) list() else x
    },

    .get_default_args = function(default_args = NULL) {
      if (!is.null(default_args)) {
        return(default_args)
      }

      list(
        pme_args_f1 = list(
          initialization_rescale = FALSE
        ),
        pme_args_f2 = list(
          initialization_rescale = FALSE
        ),
        init_args_f1 = list(
          min_clusters = 10,
          alpha = 0.01,
          max_clusters = 100,
          algorithm = "isomap",
          rescale = FALSE
        ),
        init_args_f2 = list(
          min_clusters = 10,
          alpha = 0.01,
          max_clusters = 100,
          algorithm = "isomap",
          rescale = FALSE
        ),
        reg_args = list(
          eps_step = 0.005,
          eps_energy = 0.005,
          max_iter = 10,
          basis_mode = "div_free",
          basis_set = build_basis_set(5, 5, basis = neumann_basis),
          Ugrid = subset(
            expand.grid(
              u = seq(0, 1, length.out = 60),
              v = seq(0, 1, length.out = 60)
            ),
            (u - 0.5)^2 + (v - 0.5)^2 <= 0.5^2
          )
        )
      )
    },

    .reset_state = function() {
      self$cycle_idx <- 0L
      self$pme1_initial <- NULL
      self$pme2_initial <- NULL
      self$f1_initial <- NULL
      self$f2_initial <- NULL
      self$f2_initial_grad <- NULL
      self$pme2_refit <- NULL
      self$f2_refit <- NULL
      self$f2_refit_grad <- NULL
      self$f2_registered <- NULL
      self$scale_f1 <- NULL
      self$scale_f2 <- NULL
      self$reg <- NULL
      self$converged <- FALSE
      self$stop_reason <- NULL
      self$best_state <- NULL
      self$history <- list(initial = NULL, cycles = list())
      invisible(TRUE)
    },

    .init_pme = function(dataX, init_args = NULL, rescale = FALSE) {
      args <- modifyList(
        list(x = dataX, d = self$d, rescale = rescale),
        private$.null_to_list(init_args)
      )

      do.call(initialize_pme, args)
    },

    .fit_pme = function(dataX,
                        pme_args = NULL,
                        initialization = NULL,
                        lambda = NULL) {
      args <- modifyList(
        list(data = dataX, d = self$d),
        private$.null_to_list(pme_args)
      )

      if (!is.null(initialization)) {
        args$initialization <- initialization
      }
      if (!is.null(lambda)) {
        args$lambda <- lambda
      }

      do.call(pme, args)
    },

    .make_scaled_embedding = function(pme_result, scale_info, d = 2) {
      pme_embedding_factory(
        pme_result = pme_result,
        d = d,
        A = scale_info$A,
        b = scale_info$b
      )
    },

    .make_scaled_grad = function(pme_result, scale_info) {
      pme_grad_factory(
        pme_result = pme_result,
        A = scale_info$A,
        b = scale_info$b
      )
    },

    .update_f2_initialization_from_f2k_inplace = function(f2k_fun,
                                                          scale_f2,
                                                          center_points,
                                                          init_params_for_centers) {
      if (is.null(self$initialization_f2)) {
        stop("initialization_f2 is NULL. Call $fit_initial() first.")
      }
      if (is.null(init_params_for_centers)) {
        stop("init_params_for_centers is NULL. Expected self$pme2_refit$params_opt.")
      }

      new_params_centers <- calc_params(
        f = f2k_fun,
        X = center_points,
        init_params = init_params_for_centers,
        f_input = "uv"
      )

      A <- scale_f2$A
      b <- scale_f2$b
      new_params_back_toPME <- t(A %*% t(new_params_centers) + b)

      self$initialization_f2$parameterization <- new_params_back_toPME
      invisible(TRUE)
    },

    .resolve_lambda_for_f2 = function() {
      if (self$lambda_policy_f2 == "reuse_prev") {
        return(self$pme2_refit$tuning)
      }
      if (self$lambda_policy_f2 == "fixed") {
        if (is.null(self$fixed_lambda_f2)) {
          stop("fixed_lambda_f2 must be provided for lambda_policy_f2='fixed'.")
        }
        return(self$fixed_lambda_f2)
      }
      if (self$lambda_policy_f2 == "retune") {
        return(exp(-15:5))
      }
    },

    .save_all_overwrite = function(save_dir, filename) {
      if (!dir.exists(save_dir)) {
        dir.create(save_dir, recursive = TRUE)
      }

      saveRDS(self, file = file.path(save_dir, filename))
      invisible(TRUE)
    },

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

      reg <- do.call(Registration$new, args)
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

    .record_history = function(k,
                               reg,
                               pme2_refit,
                               scale_f2,
                               f2_refit,
                               lambda) {
      E_hist <- reg$E_history
      final_E <- if (!is.null(E_hist) && length(E_hist) > 0L) {
        tail(E_hist, 1)
      } else {
        NA_real_
      }

      self$history$cycles[[length(self$history$cycles) + 1L]] <- list(
        k = k,
        reg = reg,
        lambda_f2 = lambda,
        n_iter = length(E_hist),
        final_E = final_E,
        scale_f2 = scale_f2,
        pme2_refit = pme2_refit,
        f2_refit = f2_refit,
        f2_registered = self$f2_registered,
        reg_E_history = E_hist
      )

      invisible(TRUE)
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
        scale_f2 = self$scale_f2,
        pme2_refit = self$pme2_refit,
        f2_refit = self$f2_refit,
        f2_refit_grad = self$f2_refit_grad,
        f2_registered = self$f2_registered,
        initialization_f2 = self$initialization_f2,
        reg = self$reg
      )
    },

    .restore_best_state = function() {
      if (is.null(self$best_state)) {
        return(invisible(FALSE))
      }

      self$scale_f2 <- self$best_state$scale_f2
      self$pme2_refit <- self$best_state$pme2_refit
      self$f2_refit <- self$best_state$f2_refit
      self$f2_refit_grad <- self$best_state$f2_refit_grad
      self$f2_registered <- self$best_state$f2_registered
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
      cycle_entry$output_name_note <- paste(
        "f2_registered is f2(gamma(u,v)) before PME refit;",
        "f2_refit is the PME refit after registration."
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


