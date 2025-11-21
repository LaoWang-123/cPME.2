

#' Registration Class for cPME Algorithm
#'
#' This class implements the surface registration algorithm used in cPME.
#' It supports:
#' - building basis functions
#' - computing tangent updates
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
                              #' @field f2_grad_fn Gradient of f2, a function returning (df/du, df/dv).
                              #' @field basis_set A list of basis functions used to build the tangent basis.
                              #' @field Ugrid A data frame of (u, v) evaluation points on [0,1]^2.

                              #' @field bi_set Precomputed basis functions b_i(u,v).
                              #' @field D_bi_set Precomputed gradients Db_i(u,v).

                              #' @field gamma_k Current reparameterization map γ_k(u,v).
                              #' @field f2_k Current warped surface f2 ∘ γ_k.
                              #' @field grad_f2k_fun Current gradient ∇(f2 ∘ γ_k).
                              #' @field dgamma_coefs Coefficients α_i in δγ = Σ α_i b_i.
                              #' @field delta_gamma_fn Function δγ_k(u,v).
                              #' @field Ddelta_gamma_fn Jacobian of δγ_k(u,v).

                              #' @field state_list List storing all iteration states (γ_k, α_k, E_k).
                              #' @field E_history Numeric vector storing energy values E_k.
                              #' @field iter Current iteration index.

                              #' @field eps_step Step size for γ update.
                              #' @field eps_energy Convergence threshold on |E_k - E_{k-1}|.
                              #' @field max_iter Maximum number of iterations.

                              #' @field folder Optional folder path for autosaving state.

                              # ---------------------------------------------------------
                              # Fields (User-provided functions and algorithm settings)
                              # ---------------------------------------------------------
                              f1 = NULL,
                              f2 = NULL,
                              f2_grad_fn = NULL,
                              basis_set = NULL,
                              Ugrid = NULL,

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

                              # Folder for autosave
                              folder = NULL,


                              # ---------------------------------------------------------
                              # Initialization
                              # ---------------------------------------------------------
                              #' @description Initialization of the registration object.
                              #' Sets user-provided surfaces, gradients, basis functions, and algorithm
                              #' hyperparameters. Precomputes basis fields and initializes γ_0 and f2_0.
                              #'
                              #' @param f1 Function f1(u,v), target surface.
                              #' @param f2 Function f2(u,v), source surface.
                              #' @param f2_grad_fn Gradient function of f2.
                              #' @param basis_set A list of basis functions.
                              #' @param Ugrid Data frame of (u,v) sample locations.
                              #' @param eps_step Numeric. Step size for updates.
                              #' @param eps_energy Numeric. Convergence threshold.
                              #' @param max_iter Maximum number of iterations.
                              #' @param folder Folder for autosave (optional).
                              #' @return The initialized Registration object (invisibly).
                              #'
                              initialize = function(
    f1, f2, f2_grad_fn,
    basis_set, Ugrid,
    eps_step = 0.05,
    eps_energy = 1e-5,
    max_iter = 100,
    folder = NULL) {
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

                                # Precompute basis fields
                                self$bi_set   <- build_bi_set(basis_set)
                                self$D_bi_set <- build_D_bi_set(basis_set)

                                # Initialize the algorithm state (k = 0)
                                self$initialize_state()
                              },


    # ---------------------------------------------------------
    # Initialize iteration state (k = 0)
    # ---------------------------------------------------------
    #' @description Initialize iteration state for k = 0.
    #' Builds γ_0, f2_0, initial basis derivatives, α coefficients,
    #' and stores the first state.
    #' @return Nothing. Updates internal fields.
    initialize_state = function() {
      gamma_id <- function(u, v) c(u, v)
      self$gamma_k <- gamma_id
      self$f2_k    <- self$f2

      # Compute initial energy E_0
      E0 <- compute_E_vectorized(self$f1, self$f2_k, self$Ugrid)

      # gradient of f2^0 = Df2
      self$grad_f2k_fun <- self$f2_grad_fn

      # Build dphi basis for iteration 0
      dphi_set <- build_dphi_set(
        basis_set     = self$basis_set,
        grad_f2k_fun  = self$grad_f2k_fun,
        f2_fun        = self$f2_k
      )

      # Compute coefficients α_i
      self$dgamma_coefs <- compute_inner_products(
        f1        = self$f1,
        f2_k      = self$f2_k,
        dphi_list = dphi_set,
        uv_grid   = self$Ugrid
      )

      # Assemble δγ and Dδγ
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
    #' Updates γ_{k+1}, f2_{k+1}, computes energy, convergence check,
    #' recomputes basis derivatives and coefficients.
    #' @return The current energy E_k (invisibly).
    step = function() {

      # --- 1) Determine step size (either scheduled or eps_step) ---
      eps <- self$eps_step

      # --- 2) Update γ^{k+1} = γ_k ∘ (id + eps * δγ_k) ---
      gamma_next <- update_gamma_fn(
        gamma_k       = self$gamma_k,
        delta_gamma_k = self$delta_gamma_fn,
        epsilon       = eps
      )

      # --- 3) Update f2^{k+1} = f2 ∘ γ^{k+1} ---
      f2_next <- function(u, v) {
        xy <- gamma_next(u, v)
        self$f2(xy[1], xy[2])
      }

      # --- 4) Compute new energy ---
      E_curr <- compute_E_vectorized(self$f1, f2_next, self$Ugrid)

      self$iter <- self$iter + 1
      self$E_history <- c(self$E_history, E_curr)

      # --- 5) Convergence check: |E_k - E_{k-1}| ---
      if (length(self$E_history) >= 2) {
        dE <- abs(self$E_history[length(self$E_history)] -
                    self$E_history[length(self$E_history) - 1])
        message(sprintf("Iter %d: E = %.6f, ΔE = %.3e", self$iter, E_curr, dE))

        if (dE < self$eps_energy) {
          message(sprintf("Converged at iter %d (|ΔE| < %.1e)",
                          self$iter, self$eps_energy))
        }
      }

      # --- 6) Compute ∇(f2 ∘ γ^{k+1}) using chain rule ---
      grad_f2k_fun <- assemble_grad_f2k_from_state(
        state_list  = self$state_list,
        gamma_k     = gamma_next,
        f2_grad_fn  = self$f2_grad_fn,
        epsilon     = eps
      )

      # --- Build new dϕ(b_i) ---
      dphi_set <- build_dphi_set(
        basis_set     = self$basis_set,
        grad_f2k_fun  = grad_f2k_fun,
        f2_fun        = f2_next
      )

      # --- Compute new coefficients α_i ---
      dgamma_coefs <- compute_inner_products(
        f1        = self$f1,
        f2_k      = f2_next,
        dphi_list = dphi_set,
        uv_grid   = self$Ugrid
      )

      # --- Assemble δγ and Dδγ ---
      delta_gamma_fn  <- assemble_delta_gamma_fn(dgamma_coefs, self$bi_set)
      Ddelta_gamma_fn <- assemble_D_delta_gamma_fn(dgamma_coefs, self$D_bi_set)

      # --- Update internal state ---
      self$gamma_k        <- gamma_next
      self$f2_k           <- f2_next
      self$grad_f2k_fun   <- grad_f2k_fun
      self$dgamma_coefs   <- dgamma_coefs
      self$delta_gamma_fn <- delta_gamma_fn
      self$Ddelta_gamma_fn<- Ddelta_gamma_fn

      # --- Save iteration record (state_list) ---
      self$state_list[[self$iter + 1]] <- list(
        iter             = self$iter,
        E                = E_curr,
        gamma_k          = gamma_next,
        f2_k             = f2_next,
        dgamma_coefs     = dgamma_coefs,
        delta_gamma_fn   = delta_gamma_fn,
        Ddelta_gamma_fn  = Ddelta_gamma_fn
      )

      # --- Autosave to folder ---
      if (!is.null(self$folder)) {
        self$save_state()
      }

      invisible(E_curr)
    },


    # ---------------------------------------------------------
    # Run all iterations up to max_iter
    # ---------------------------------------------------------
    #' @description Run the full optimization loop up to `max_iter`.
    #' Automatically stops if the energy change falls below `eps_energy`.
    #' @return Nothing. Updates the entire optimization history.
    run = function() {
      for (i in seq_len(self$max_iter)) {

        self$step()

        # Auto-stop when converged
        if (length(self$E_history) >= 2) {
          dE <- abs(
            self$E_history[length(self$E_history)] -
              self$E_history[length(self$E_history) - 1]
          )
          if (dE < self$eps_energy) {
            message(sprintf(
              "Stopped at iteration %d: |ΔE| < %.1e",
              self$iter, self$eps_energy
            ))
            break
          }
        }
      }
    },



    # ---------------------------------------------------------
    # Continue optimization from a saved state
    # ---------------------------------------------------------
    #' @description Continue optimization from a previously saved state.
    #' Allows overriding step size, schedule, energy threshold, etc.
    #'
    #' @param n_steps Number of additional steps to run (optional).
    #' @param max_iter_total Absolute iteration target (optional).
    #' @param eps_step Optional override of step size.
    #' @param step_schedule Optional new step schedule.
    #' @param eps_energy Optional new convergence tolerance.
    #' @return Nothing. Continues modifying the registration state.
    continue = function(
    n_steps = NULL,
    max_iter_total = NULL,
    eps_step = NULL,
    eps_energy = NULL
    ) {

      # Parameter overrides
      if (!is.null(eps_step))      self$eps_step <- eps_step
      if (!is.null(eps_energy))    self$eps_energy <- eps_energy

      start_iter <- self$iter

      # Determine target iteration
      if (!is.null(max_iter_total)) {
        target_iter <- max_iter_total
      } else if (!is.null(n_steps)) {
        target_iter <- self$iter + n_steps
      } else {
        target_iter <- self$max_iter
      }

      while (self$iter < target_iter) {

        self$step()

        # Auto-stop when converged
        if (length(self$E_history) >= 2) {
          dE <- abs(
            self$E_history[length(self$E_history)] -
              self$E_history[length(self$E_history) - 1]
          )
          if (dE < self$eps_energy) {
            message(sprintf(
              "Early stop at iter %d: |ΔE| < %.1e",
              self$iter, self$eps_energy
            ))
            break
          }
        }
      }
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
      saveRDS(self, file.path(self$folder, "reg_state.rds"))
    }

    ) # end public list

)
