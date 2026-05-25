#' Safe Registration Class for cPME Algorithm
#'
#' This experimental registration class inherits from `Registration_new` but
#' rejects any registration step that increases the energy. If the first
#' proposed update is worse, the registration remains at the identity state.
#' It can also backtrack the step size within each registration update.
#'
#' @docType class
#' @format An \code{R6Class} generator object
#' @export
Registration_safe <- R6::R6Class(
  classname = "Registration_safe",
  inherit = Registration_new,

  public = list(
    #' @field rejected_steps List of proposed but rejected registration steps.
    rejected_steps = NULL,

    #' @field last_step_accepted Logical flag for the most recent step.
    last_step_accepted = TRUE,

    #' @field stop_reason Text description of the safe stopping condition.
    stop_reason = NULL,

    #' @field eps_shrink Multiplicative backtracking shrink factor.
    eps_shrink = 0.5,

    #' @field eps_min Smallest step size considered in backtracking.
    eps_min = 0.0005,

    #' @field max_backtracking Maximum number of step-size attempts per update.
    max_backtracking = 12L,

    #' @field domain_lower Lower allowed coordinate for gamma(Ugrid).
    domain_lower = -0.15,

    #' @field domain_upper Upper allowed coordinate for gamma(Ugrid).
    domain_upper = 1.15,

    #' @field max_displacement Optional cap on max gamma displacement per step.
    max_displacement = 0.05,

    #' @field backtracking_history List of attempted step sizes and energies.
    backtracking_history = NULL,

    #' @field energy_Ugrid Grid used only for energy evaluation.
    energy_Ugrid = NULL,

    #' @field f1_energy_grid Template surface values on `energy_Ugrid`.
    f1_energy_grid = NULL,

    #' @field energy_grid_source Description of the energy grid choice.
    energy_grid_source = NULL,

    initialize = function(...,
                          energy_Ugrid = NULL,
                          energy_grid_source = NULL,
                          eps_shrink = 0.5,
                          eps_min = 0.0005,
                          max_backtracking = 12L,
                          domain_lower = -0.15,
                          domain_upper = 1.15,
                          max_displacement = 0.05) {
      super$initialize(...)
      self$rejected_steps <- list()
      self$last_step_accepted <- TRUE
      self$stop_reason <- NULL
      self$eps_shrink <- eps_shrink
      self$eps_min <- eps_min
      self$max_backtracking <- as.integer(max_backtracking)
      self$domain_lower <- domain_lower
      self$domain_upper <- domain_upper
      self$max_displacement <- max_displacement
      self$backtracking_history <- list()

      # `Ugrid` is the optimization/update grid inherited from Registration_new.
      # `energy_Ugrid` lets us evaluate the registration objective on a
      # different set of parameter locations without changing the basis update
      # locations. When omitted, safe registration is identical to the original
      # single-grid behavior.
      if (is.null(energy_Ugrid)) {
        self$energy_Ugrid <- self$Ugrid
        self$energy_grid_source <- "optimization_Ugrid"
      } else {
        self$energy_Ugrid <- private$.standardize_uv_grid(energy_Ugrid)
        self$energy_grid_source <- if (is.null(energy_grid_source)) {
          "custom_energy_Ugrid"
        } else {
          energy_grid_source
        }

        self$f1_energy_grid <- make_f2_grid(self$f1, self$energy_Ugrid)
        f2_energy_grid <- make_f2_grid(self$f2, self$energy_Ugrid)
        E0 <- compute_E_grid(self$f1_energy_grid, f2_energy_grid)

        self$E_history <- E0
        if (length(self$state_list) >= 1L) {
          self$state_list[[1L]]$E <- E0
        }
      }

      if (is.null(self$f1_energy_grid)) {
        self$f1_energy_grid <- self$f1_grid
      }

      invisible(self)
    },

    #' Perform one safe registration update.
    #'
    #' A proposed update is accepted only when its energy is not larger than the
    #' current accepted energy and gamma(Ugrid) remains in the allowed domain.
    #' If the first proposed step fails, the method shrinks eps and retries.
    #' Rejected updates are logged but not appended to `state_list` or
    #' `E_history`, so downstream code only sees accepted states.
    step = function() {
      E_prev <- self$E_history[length(self$E_history)]
      delta_gamma_fns <- lapply(self$state_list, function(s) s$delta_gamma_fn)
      proposed_iter <- self$iter + 1L

      current_delta <- self$state_list[[length(self$state_list)]]$delta_gamma_fn
      delta_grid <- t(vapply(
        seq_len(nrow(self$Ugrid)),
        function(i) current_delta(self$Ugrid[i, 1], self$Ugrid[i, 2]),
        numeric(2)
      ))
      max_delta <- max(sqrt(rowSums(delta_grid^2)), na.rm = TRUE)

      eps_start <- self$eps_step
      if (!is.null(self$max_displacement) &&
          is.finite(max_delta) &&
          max_delta > 0) {
        eps_start <- min(eps_start, self$max_displacement / max_delta)
      }
      eps_start <- max(eps_start, self$eps_min)

      attempts <- list()
      accepted <- NULL
      eps_try <- eps_start

      for (attempt_idx in seq_len(self$max_backtracking)) {
        gamma_next <- make_gamma_from_history(
          delta_gamma_fns = delta_gamma_fns,
          eps = eps_try
        )

        gamma_uv_update <- private$.eval_gamma_on_grid(gamma_next, self$Ugrid)
        gamma_uv_energy <- private$.eval_gamma_on_grid(
          gamma_next,
          self$energy_Ugrid
        )
        gamma_uv <- rbind(gamma_uv_update, gamma_uv_energy)

        domain_ok <- all(is.finite(gamma_uv)) &&
          min(gamma_uv[, 1]) >= self$domain_lower &&
          max(gamma_uv[, 1]) <= self$domain_upper &&
          min(gamma_uv[, 2]) >= self$domain_lower &&
          max(gamma_uv[, 2]) <= self$domain_upper

        f2_next <- function(u, v) {
          xy <- gamma_next(u, v)
          self$f2(xy[1], xy[2])
        }

        f2_next_update_grid <- make_f2_grid(f2_next, self$Ugrid)
        f2_next_energy_grid <- make_f2_grid(f2_next, self$energy_Ugrid)
        E_curr <- compute_E_grid(self$f1_energy_grid, f2_next_energy_grid)
        dE <- E_curr - E_prev
        energy_ok <- is.finite(E_curr) && E_curr <= E_prev

        attempt <- list(
          iter = proposed_iter,
          attempt = attempt_idx,
          eps_try = eps_try,
          E_prev = E_prev,
          E_proposed = E_curr,
          dE = dE,
          domain_ok = domain_ok,
          energy_ok = energy_ok,
          gamma_u_min = min(gamma_uv[, 1]),
          gamma_u_max = max(gamma_uv[, 1]),
          gamma_v_min = min(gamma_uv[, 2]),
          gamma_v_max = max(gamma_uv[, 2])
        )
        attempts[[length(attempts) + 1L]] <- attempt

        if (isTRUE(self$verbose)) {
          message(sprintf(
            paste0(
              "Iter %d attempt %d: eps = %.6g, E = %.6f, ",
              "dE = %.3e, domain_ok = %s"
            ),
            proposed_iter, attempt_idx, eps_try, E_curr, dE, domain_ok
          ))
        }

        if (energy_ok && domain_ok) {
          accepted <- list(
            gamma_next = gamma_next,
            f2_next = f2_next,
            f2_next_update_grid = f2_next_update_grid,
            f2_next_energy_grid = f2_next_energy_grid,
            E_curr = E_curr,
            dE = dE,
            eps_used = eps_try,
            attempts = attempts
          )
          break
        }

        next_eps <- eps_try * self$eps_shrink
        if (next_eps < self$eps_min) {
          break
        }
        eps_try <- next_eps
      }

      self$backtracking_history[[length(self$backtracking_history) + 1L]] <- list(
        iter = proposed_iter,
        max_delta = max_delta,
        eps_start = eps_start,
        attempts = attempts,
        accepted = !is.null(accepted)
      )

      if (is.null(accepted)) {
        last_attempt <- attempts[[length(attempts)]]
        self$last_step_accepted <- FALSE
        self$stop_reason <- if (!isTRUE(last_attempt$domain_ok)) {
          "rejected domain violation"
        } else if (!is.finite(last_attempt$E_proposed)) {
          "rejected non-finite energy"
        } else {
          "rejected energy increase"
        }
        self$rejected_steps[[length(self$rejected_steps) + 1L]] <- list(
          iter = proposed_iter,
          E_prev = E_prev,
          attempts = attempts,
          stop_reason = self$stop_reason
        )

        if (isTRUE(self$verbose)) {
          message(sprintf(
            "Rejected iter %d after %d attempts: %s.",
            proposed_iter, length(attempts), self$stop_reason
          ))
        }

        if (!is.null(self$folder)) {
          self$save_state()
        }

        return(invisible(E_prev))
      }

      grad_f2k_fun <- assemble_grad_f2k_from_state(
        state_list = self$state_list,
        gamma_k = accepted$gamma_next,
        f2_grad_fn = self$f2_grad_fn,
        epsilon = accepted$eps_used
      )

      grad_f2_grid <- make_grad_f2_grid(grad_f2k_fun, self$Ugrid)

      dphi_grid_list <- compute_dphi_grid(
        basis_grid = self$basis_grid,
        f2_grid = accepted$f2_next_update_grid,
        grad_f2_grid = grad_f2_grid,
        mode = self$basis_mode
      )

      dgamma_coefs <- compute_inner_products_fast(
        diff_grid = self$f1_grid - accepted$f2_next_update_grid,
        dphi_grid_list = dphi_grid_list,
        weight = 1 / self$n
      )

      delta_gamma_fn <- assemble_delta_gamma_fn(dgamma_coefs, self$bi_set)
      Ddelta_gamma_fn <- assemble_D_delta_gamma_fn(dgamma_coefs, self$D_bi_set)

      self$iter <- proposed_iter
      self$E_history <- c(self$E_history, accepted$E_curr)
      self$gamma_k <- accepted$gamma_next
      self$f2_k <- accepted$f2_next
      self$grad_f2k_fun <- grad_f2k_fun
      self$dgamma_coefs <- dgamma_coefs
      self$delta_gamma_fn <- delta_gamma_fn
      self$Ddelta_gamma_fn <- Ddelta_gamma_fn
      self$last_step_accepted <- TRUE
      self$stop_reason <- NULL

      self$state_list[[length(self$state_list) + 1L]] <- list(
        iter = self$iter,
        E = accepted$E_curr,
        gamma_k = accepted$gamma_next,
        f2_k = accepted$f2_next,
        dgamma_coefs = dgamma_coefs,
        delta_gamma_fn = delta_gamma_fn,
        Ddelta_gamma_fn = Ddelta_gamma_fn,
        eps_used = accepted$eps_used,
        energy_grid_source = self$energy_grid_source,
        n_energy_grid = nrow(self$energy_Ugrid)
      )

      if (isTRUE(self$verbose)) {
        message(sprintf(
          "Accepted iter %d: E = %.6f, eps = %.6g",
          self$iter, accepted$E_curr, accepted$eps_used
        ))
      }

      if (!is.null(self$folder)) {
        self$save_state()
      }

      invisible(accepted$E_curr)
    },

    #' Run safe registration up to `max_iter`.
    run = function() {
      for (i in seq_len(self$max_iter)) {
        E_before <- self$E_history[length(self$E_history)]
        self$step()

        if (!isTRUE(self$last_step_accepted)) {
          if (isTRUE(self$verbose)) {
            message(sprintf(
              "Stopped at iteration %d: %s.",
              i, self$stop_reason
            ))
          }
          break
        }

        E_after <- self$E_history[length(self$E_history)]
        dE_abs <- abs(E_after - E_before)
        if (dE_abs < self$eps_energy) {
          self$stop_reason <- "converged by energy tolerance"
          if (isTRUE(self$verbose)) {
            message(sprintf(
              "Stopped at iteration %d: |dE| < %.1e.",
              i, self$eps_energy
            ))
          }
          break
        }
      }

      invisible(self)
    },

    #' Continue safe registration from the current accepted state.
    continue = function(n_steps = NULL,
                        max_iter_total = NULL,
                        eps_step = NULL,
                        eps_energy = NULL) {
      if (!is.null(eps_step)) self$eps_step <- eps_step
      if (!is.null(eps_energy)) self$eps_energy <- eps_energy

      if (!is.null(max_iter_total)) {
        target_iter <- max_iter_total
      } else if (!is.null(n_steps)) {
        target_iter <- self$iter + n_steps
      } else {
        target_iter <- self$max_iter
      }

      while (self$iter < target_iter) {
        E_before <- self$E_history[length(self$E_history)]
        self$step()

        if (!isTRUE(self$last_step_accepted)) {
          break
        }

        E_after <- self$E_history[length(self$E_history)]
        if (abs(E_after - E_before) < self$eps_energy) {
          self$stop_reason <- "converged by energy tolerance"
          break
        }
      }

      invisible(self)
    }
  ),

  private = list(
    .standardize_uv_grid = function(Ugrid) {
      Ugrid <- as.data.frame(Ugrid)
      if (ncol(Ugrid) < 2L) {
        stop("energy_Ugrid must have at least two columns.")
      }
      Ugrid <- Ugrid[, 1:2, drop = FALSE]
      names(Ugrid) <- c("u", "v")
      Ugrid
    },

    .eval_gamma_on_grid = function(gamma_fn, Ugrid) {
      gamma_uv <- t(vapply(
        seq_len(nrow(Ugrid)),
        function(i) gamma_fn(Ugrid[i, 1], Ugrid[i, 2]),
        numeric(2)
      ))
      colnames(gamma_uv) <- c("u", "v")
      gamma_uv
    }
  )
)
