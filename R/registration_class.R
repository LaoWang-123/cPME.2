

#' Registration Class for cPME Algorithm
#'
#' This is the current registration implementation used by cPME. It combines
#' the fast grid-precomputed registration path with safe step acceptance:
#' each proposed update is accepted only if it does not increase the selected
#' registration energy and remains inside the allowed parameter-domain
#' tolerance. Backtracking is used to shrink the step size when needed.
#' It supports:
#' - building basis functions
#' - computing tangent updates
#' - separate optimization and energy grids
#' - safe backtracking and accepted-state history
#' - iterative optimization
#' - saving state to RDS
#'
#' @docType class
#' @format An \code{R6Class} generator object
#' @export
Registration <- R6::R6Class("Registration",

                            public = list(
                              #' @field f1 Function f1(u, v). The target surface embedding or image.
                              #' @field f2 Function f2(u, v). The source surface to be warped.
                              #' @field f1_grid The output of f1 function with uv-grid.
                              #' @field f2_grid The output of f2 function with uv-grid.
                              #' @field f2_grad_fn Gradient of f2, a function returning (df/du, df/dv).
                              #' @field grad_f2_grid The output of f2_grad_fn function with uv-grid.
                              #' @field basis_set A list of basis functions used to build the tangent basis.
                              #' @field basis_grid A list of basis functions' output with grid input.
                              #' @field Ugrid A data frame of (u, v) evaluation points on [0,1]^2.
                              #' @field n the length of grid.

                              #' @field bi_set Precomputed basis functions b_i(u,v).
                              #' @field D_bi_set Precomputed gradients Db_i(u,v).

                              #' @field gamma_k Current reparameterization map Î³_k(u,v).
                              #' @field f2_k Current warped surface f2 âˆ˜ Î³_k.
                              #' @field grad_f2k_fun Current gradient âˆ‡(f2 âˆ˜ Î³_k).
                              #' @field dgamma_coefs Coefficients Î±_i in Î´Î³ = Î£ Î±_i b_i.
                              #' @field delta_gamma_fn Function Î´Î³_k(u,v).
                              #' @field Ddelta_gamma_fn Jacobian of Î´Î³_k(u,v).

                              #' @field state_list List storing all iteration states (Î³_k, Î±_k, E_k).
                              #' @field E_history Numeric vector storing energy values E_k.
                              #' @field iter Current iteration index.

                              #' @field eps_step Step size for Î³ update.
                              #' @field eps_energy Convergence threshold on |E_k - E_{k-1}|.
                              #' @field max_iter Maximum number of iterations.

                              #' @field basis_mode "full" or "div_free"

                              #' @field folder Optional folder path for autosaving state.
                              #' @field filename default "reg_state.rds"
                              #' @field verbose If TRUE, print per-iteration energy updates as the optimizer runs.
                              #' @field energy_Ugrid Grid used to evaluate the registration objective.
                              #' @field f1_energy_grid Template surface values on `energy_Ugrid`.
                              #' @field energy_grid_source Text label describing the energy grid choice.
                              #' @field rejected_steps Proposed steps rejected by energy or domain checks.
                              #' @field last_step_accepted Logical flag for the most recent registration step.
                              #' @field stop_reason Text description of the current stopping condition.
                              #' @field eps_shrink Multiplicative backtracking shrink factor.
                              #' @field eps_min Smallest step size tried during backtracking.
                              #' @field max_backtracking Maximum number of step-size attempts per update.
                              #' @field domain_lower Lower allowed coordinate for gamma grids.
                              #' @field domain_upper Upper allowed coordinate for gamma grids.
                              #' @field max_displacement Optional cap on max gamma displacement per step.
                              #' @field backtracking_history Attempted eps values and energies for each step.

                              # ---------------------------------------------------------
                              # Fields (User-provided functions and algorithm settings)
                              # ---------------------------------------------------------
                              f1 = NULL,
                              f2 = NULL,
                              f1_grid = NULL,
                              f2_grid = NULL,
                              f2_grad_fn = NULL,
                              grad_f2_grid = NULL,
                              basis_set = NULL,
                              basis_grid = NULL,
                              Ugrid = NULL,
                              n = NULL,

                              bi_set = NULL,
                              D_bi_set = NULL,

                              # Current state
                              gamma_k = NULL,
                              f2_k = NULL,
                              grad_f2k_fun = NULL,
                              dgamma_coefs = NULL,
                              delta_gamma_fn = NULL,
                              Ddelta_gamma_fn = NULL,

                              # Logs
                              state_list = NULL,
                              E_history = NULL,
                              iter = 0,

                              # Hyperparameters
                              eps_step = NULL,
                              eps_energy = NULL,
                              max_iter = NULL,

                              # basis selection mode (divergence or not)
                              basis_mode = NULL,

                              # Folder for autosave
                              folder = NULL,
                              filename = "reg_state.rds",
                              verbose = TRUE,

                              # Safe registration controls
                              rejected_steps = NULL,
                              last_step_accepted = TRUE,
                              stop_reason = NULL,
                              eps_shrink = 0.5,
                              eps_min = 0.0005,
                              max_backtracking = 12L,
                              domain_lower = -0.15,
                              domain_upper = 1.15,
                              max_displacement = 0.05,
                              backtracking_history = NULL,
                              energy_Ugrid = NULL,
                              f1_energy_grid = NULL,
                              energy_grid_source = NULL,


                              # ---------------------------------------------------------
                              # Initialization
                              # ---------------------------------------------------------
                              #' @description Initialization of the registration object.
                              #' Sets user-provided surfaces, gradients, basis functions, and algorithm
                              #' hyperparameters. Precomputes basis fields and initializes Î³_0 and f2_0.
                              #'
                              #' @param f1 Function f1(u,v), target surface.
                              #' @param f2 Function f2(u,v), source surface.
                              #' @param f2_grad_fn Gradient function of f2.
                              #' @param basis_set A list of basis functions.
                              #' @param Ugrid Data frame of (u,v) sample locations.
                              #' @param eps_step Numeric. Step size for updates.
                              #' @param eps_energy Numeric. Convergence threshold.
                              #' @param max_iter Maximum number of iterations.
                              #' @param basis_mode c("full", "div_free".
                              #' @param folder Folder for autosave (optional).
                              #' @param filename filename for autosave, default "reg_state.rds"
                              #' @param verbose If TRUE, print per-iteration energy updates as the optimizer runs.
                              #' @param energy_Ugrid Optional grid used only for objective evaluation.
                              #' @param energy_grid_source Optional label describing `energy_Ugrid`.
                              #' @param eps_shrink Multiplicative shrink factor for backtracking.
                              #' @param eps_min Smallest eps tried during backtracking.
                              #' @param max_backtracking Maximum number of eps attempts per step.
                              #' @param domain_lower Lower allowed coordinate for gamma grids.
                              #' @param domain_upper Upper allowed coordinate for gamma grids.
                              #' @param max_displacement Optional cap on the maximum update displacement per step.
                              #' @return The initialized Registration object (invisibly).
                              #'
                              initialize = function(
    f1, f2, f2_grad_fn,
    basis_set, Ugrid,
    eps_step = 0.05,
    eps_energy = 1e-5,
    max_iter = 100,
    basis_mode = c("full", "div_free"),
    folder = NULL,
    filename = "reg_state.rds",
    verbose = TRUE,
    energy_Ugrid = NULL,
    energy_grid_source = NULL,
    eps_shrink = 0.5,
    eps_min = 0.0005,
    max_backtracking = 12L,
    domain_lower = -0.15,
    domain_upper = 1.15,
    max_displacement = 0.05) {
                                # =======================
                                # INITIALIZATION MODE
                                # =======================
                                # User-supplied functions
                                self$f1 = f1
                                self$f2 = f2
                                self$f2_grad_fn = f2_grad_fn
                                self$basis_set = basis_set
                                self$Ugrid = Ugrid

                                # Settings
                                self$eps_step = eps_step
                                self$eps_energy = eps_energy
                                self$max_iter = max_iter
                                self$folder = folder
                                self$filename = filename
                                self$verbose <- isTRUE(verbose)

                                # NEW: save mode
                                self$basis_mode <- match.arg(basis_mode)

                                # Precompute basis fields
                                self$bi_set   <- build_bi_set(basis_set = self$basis_set, mode = self$basis_mode) ## need to add mode parameter
                                self$D_bi_set <- build_D_bi_set(basis_set = self$basis_set, mode = self$basis_mode) ## need to add mode parameter

                                # Precompute grid elements
                                self$n=nrow(self$Ugrid)
                                self$basis_grid <- build_basis_grid(basis_set = self$basis_set, Ugrid = self$Ugrid, mode = self$basis_mode) ## need to add mode parameter
                                self$f1_grid <- make_f2_grid(self$f1, self$Ugrid)

                                # Safe registration settings
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
                                }
                                self$f1_energy_grid <- make_f2_grid(self$f1, self$energy_Ugrid)

                                # Initialize the algorithm state (k = 0)
                                self$initialize_state()
                              },


    # ---------------------------------------------------------
    # Initialize iteration state (k = 0)
    # ---------------------------------------------------------
    #' @description Initialize iteration state for k = 0.
    #' Builds Î³_0, f2_0, initial basis derivatives, Î± coefficients,
    #' and stores the first state.
    #' @return Nothing. Updates internal fields.
    initialize_state = function() {
      gamma_id <- function(u, v) c(u, v)
      self$gamma_k <- gamma_id
      self$f2_k    <- self$f2

      # gradient of f2^0 = Df2
      self$grad_f2k_fun <- self$f2_grad_fn

      ####################################
      #### New edited part
      self$f2_grid <- make_f2_grid(self$f2_k, self$Ugrid)
      f2_energy_grid <- make_f2_grid(self$f2_k, self$energy_Ugrid)

      # Compute initial energy E_0
      E0 <- compute_E_grid(self$f1_energy_grid, f2_energy_grid)
      if (isTRUE(self$verbose)) {
        message(sprintf("Iter %d: E = %.6f", 0L, E0))
      }

      self$grad_f2_grid <- make_grad_f2_grid(self$grad_f2k_fun, self$Ugrid)

      dphi_grid_list <- compute_dphi_grid(  ## need to add mode parameter
        basis_grid = self$basis_grid,
        f2_grid = self$f2_grid,
        grad_f2_grid = self$grad_f2_grid,
        mode = self$basis_mode
      )

      self$dgamma_coefs <- compute_inner_products_fast(
        diff_grid = self$f1_grid-self$f2_grid,
        dphi_grid_list = dphi_grid_list,
        weight = 1/self$n
      )
      ####################################
      ####################################

      # Assemble Î´Î³ and DÎ´Î³
      self$delta_gamma_fn  <- assemble_delta_gamma_fn(self$dgamma_coefs, self$bi_set)
      self$Ddelta_gamma_fn <- assemble_D_delta_gamma_fn(self$dgamma_coefs, self$D_bi_set)

      # Initialize storage
      self$state_list <- list()
      self$state_list[[1]] <- list(
        iter             = 0,
        E                = E0,
        gamma_k          = self$gamma_k,
        f2_k             = self$f2_k,
        dgamma_coefs     = self$dgamma_coefs,
        delta_gamma_fn   = self$delta_gamma_fn,
        Ddelta_gamma_fn  = self$Ddelta_gamma_fn
      )

      self$E_history <- E0
      self$iter <- 0
    },


    # ---------------------------------------------------------
    # Perform ONE gradient-descent iteration
    # ---------------------------------------------------------
    #' @description Perform ONE gradient-descent update step.
    #' Updates Î³_{k+1}, f2_{k+1}, computes energy, convergence check,
    #' recomputes basis derivatives and coefficients.
    #' @return The current energy E_k (invisibly).
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
        gamma_uv_energy <- private$.eval_gamma_on_grid(gamma_next, self$energy_Ugrid)
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


    # ---------------------------------------------------------
    # Run all iterations up to max_iter
    # ---------------------------------------------------------
    #' @description Run safe registration up to `max_iter`.
    #' Stops at convergence, domain violation, non-finite energy, or rejected
    #' energy increase. Rejected states are logged but not accepted.
    #' @return The current registration object, invisibly.
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



    # ---------------------------------------------------------
    # Continue optimization from a saved state
    # ---------------------------------------------------------
    #' @description Continue optimization from a previously saved state.
    #' Allows overriding step size and energy threshold.
    #'
    #' @param n_steps Number of additional steps to run (optional).
    #' @param max_iter_total Absolute iteration target (optional).
    #' @param eps_step Optional override of step size.
    #' @param eps_energy Optional new convergence tolerance.
    #' @return The current registration object, invisibly.
    continue = function(
    n_steps = NULL,
    max_iter_total = NULL,
    eps_step = NULL,
    eps_energy = NULL
    ) {

      if (!is.null(eps_step))      self$eps_step <- eps_step
      if (!is.null(eps_energy))    self$eps_energy <- eps_energy

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
    },

    # ---------------------------------------------------------
    # Save state to folder as RDS file
    # ---------------------------------------------------------
    #' @description Save the full registration object into an RDS file.
    #' @return Nothing. Writes `reg_state.rds` to `self$folder`.
    save_state = function() {
      if (!dir.exists(self$folder)) {
        dir.create(self$folder, recursive = TRUE)
      }
      saveRDS(self, file.path(self$folder, self$filename))
    }

                            ), # end public list

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

# Compatibility alias for older scripts. New code should call Registration$new().
Registration_new <- Registration

# Compatibility alias for older scripts. New code should call Registration$new().
Registration_safe <- Registration


