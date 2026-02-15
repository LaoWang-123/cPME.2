#' PME + Registration iterative cycle (fixed f1, update f2)
#'
#' This R6 class orchestrates:
#' 1) Fit PME for data1 and data2
#' 2) Rescale both parameterizations to [0,1]^d
#' 3) Run registration to warp f2 toward f1
#' 4) Update f2 initialization (parameterization of centers) using final warped embedding f2_k
#' 5) Refit PME for data2, repeat
#'
#' Requirements (already in your package):
#' - pme(), initialize_pme(), pme_initial_guess()
#' - calc_params(), scale_uniform_square_with_params()
#' - pme_embedding_factory(), pme_grad_factory()
#' - Registration_new R6 class
#'
#' @docType class
#' @format An \code{R6Class} generator object
#' @export
PMERegistrationCycle <- R6::R6Class(
  classname = "PMERegistrationCycle",
  public = list(

    # -------------------------
    # Data / config
    # -------------------------
    data1 = NULL,
    data2 = NULL,
    d = NULL,

    pme_args_f1 = NULL,
    pme_args_f2 = NULL,
    init_args_f2 = NULL,

    reg_args = NULL,

    lambda_policy_f2 = NULL,  # "reuse_prev" (default) / "fixed" / "retune"
    fixed_lambda_f2 = NULL,   # used if lambda_policy_f2 == "fixed"

    # -------------------------
    # Current state
    # -------------------------
    cycle_idx = 0L,

    pme1 = NULL,
    pme2 = NULL,

    scale_f1 = NULL,          # cached scaling info for current pme1
    scale_f2 = NULL,          # scaling info for current pme2 (changes each cycle)

    reg = NULL,
    gamma = NULL,
    f2_warped = NULL,         # last state's f2_k (warped embedding)

    converged = FALSE,
    stop_reason = NULL,

    # -------------------------
    # History
    # -------------------------
    history = NULL,

    # -------------------------
    # Constructor
    # -------------------------
    initialize = function(
    data1,
    data2,
    d = 2,
    pme_args_f1 = list(),
    pme_args_f2 = list(),
    init_args_f2 = list(),
    reg_args = list(),
    lambda_policy_f2 = c("reuse_prev", "fixed", "retune"),
    fixed_lambda_f2 = NULL
    ) {
      self$data1 <- data1
      self$data2 <- data2
      self$d <- as.integer(d)

      self$pme_args_f1 <- pme_args_f1
      self$pme_args_f2 <- pme_args_f2
      self$init_args_f2 <- init_args_f2

      self$reg_args <- reg_args

      self$lambda_policy_f2 <- match.arg(lambda_policy_f2)
      self$fixed_lambda_f2 <- fixed_lambda_f2

      self$cycle_idx <- 0L
      self$history <- list()
      invisible(self)
    },

    # -------------------------
    # Public API
    # -------------------------

    #' Fit initial PME for both datasets
    fit_initial = function() {
      private$.fit_initial_pmes()
      # cache scaling for f1 now (f1 fixed across cycles)
      self$scale_f1 <- private$.compute_projection_and_scale(self$pme1, self$data1)
      invisible(self)
    },

    #' Run a single cycle: register + update f2 init + refit f2
    run_cycle = function(record_full = TRUE) {
      if (is.null(self$pme1) || is.null(self$pme2)) {
        stop("Call $fit_initial() before running cycles.")
      }

      k <- self$cycle_idx + 1L

      # 1) compute scaling for current f2 (must be recomputed each cycle)
      self$scale_f2 <- private$.compute_projection_and_scale(self$pme2, self$data2)

      # 2) make scaled embedding functions for registration
      f1_fun <- private$.make_scaled_embedding(self$pme1, self$scale_f1)
      f2_pack <- private$.make_scaled_embedding_and_grad(self$pme2, self$scale_f2)

      # 3) registration
      self$reg <- private$.register_once(
        f1_fun = f1_fun,
        f2_fun = f2_pack$f,
        grad_f2_fun = f2_pack$grad
      )

      last_state <- private$.extract_last_state(self$reg)
      self$gamma <- last_state$gamma_k
      self$f2_warped <- last_state$f2_k

      # 4) update f2 initialization using final f2_k
      updated_init <- private$.update_f2_initialization_from_f2k(
        f2k_fun = self$f2_warped,
        scale_f2 = self$scale_f2
      )

      # 5) refit f2
      lambda_to_use <- private$.resolve_lambda_for_f2()
      old_pme2 <- self$pme2

      self$pme2 <- private$.refit_f2(
        updated_init = updated_init,
        lambda = lambda_to_use
      )

      # 6) record history
      private$.record_history(
        k = k,
        old_pme2 = old_pme2,
        new_pme2 = self$pme2,
        scale_f1 = self$scale_f1,
        scale_f2 = self$scale_f2,
        reg = self$reg,
        last_state = last_state,
        lambda = lambda_to_use,
        record_full = record_full
      )

      self$cycle_idx <- k
      invisible(self)
    },

    #' Run multiple cycles with (optional) early stop
    run = function(
    n_cycles = 5,
    tol_E = NULL,
    stop_rule = c("none", "delta_E"),
    record_full = TRUE
    ) {
      stop_rule <- match.arg(stop_rule)

      self$converged <- FALSE
      self$stop_reason <- NULL

      for (i in seq_len(n_cycles)) {
        self$run_cycle(record_full = record_full)

        if (stop_rule == "delta_E" && !is.null(tol_E)) {
          if (private$.check_convergence_delta_E(tol_E = tol_E)) {
            self$converged <- TRUE
            self$stop_reason <- sprintf("delta_E < %.3g", tol_E)
            break
          }
        }
      }

      invisible(self)
    },

    #' Get current key objects
    get_current = function() {
      list(
        cycle = self$cycle_idx,
        pme1 = self$pme1,
        pme2 = self$pme2,
        reg = self$reg,
        gamma = self$gamma,
        f2_warped = self$f2_warped,
        converged = self$converged,
        stop_reason = self$stop_reason
      )
    },

    #' Get history (compact or full)
    get_history = function(compact = TRUE) {
      if (!compact) return(self$history)

      lapply(self$history, function(h) {
        list(
          k = h$k,
          lambda_f2 = h$lambda_f2,
          final_E = h$final_E,
          n_iter = h$n_iter
        )
      })
    }
  ),

  private = list(

    # -------------------------
    # Internal steps
    # -------------------------

    .fit_initial_pmes = function() {
      # Fit f1
      self$pme1 <- do.call(
        pme,
        c(list(data = self$data1, d = self$d), self$pme_args_f1)
      )
      # Fit f2
      self$pme2 <- do.call(
        pme,
        c(list(data = self$data2, d = self$d), self$pme_args_f2)
      )
      invisible(TRUE)
    },

    .compute_projection_and_scale = function(pme_result, dataX, init_guess_method = "pca") {
      init_param <- pme_initial_guess(X = dataX, d = self$d, method = init_guess_method)
      U_proj <- calc_params(
        f = pme_result$embedding_map,
        X = dataX,
        init_params = init_param
      )
      scaled <- scale_uniform_square_with_params(U = U_proj)

      # Standardize return fields
      list(
        init_param = init_param,
        U_proj = U_proj,
        A = scaled$A,
        b = scaled$b,
        U_scaled = scaled$U_scaled
      )
    },

    .make_scaled_embedding = function(pme_result, scale_info) {
      # embedding factory only
      pme_embedding_factory(
        pme_result = pme_result,
        d = self$d,
        A = scale_info$A,
        b = scale_info$b
      )
    },

    .make_scaled_embedding_and_grad = function(pme_result, scale_info) {
      f <- pme_embedding_factory(
        pme_result = pme_result,
        d = self$d,
        A = scale_info$A,
        b = scale_info$b
      )
      grad <- pme_grad_factory(
        pme_result = pme_result,
        A = scale_info$A,
        b = scale_info$b
      )
      list(f = f, grad = grad)
    },

    .register_once = function(f1_fun, f2_fun, grad_f2_fun) {
      reg <- do.call(
        Registration_new$new,
        c(list(
          f1 = f1_fun,
          f2 = f2_fun,
          f2_grad_fn = grad_f2_fun
        ), self$reg_args)
      )
      reg$run()
      reg
    },

    .extract_last_state = function(reg) {
      if (is.null(reg$state_list) || length(reg$state_list) == 0L) {
        stop("Registration object has empty state_list.")
      }
      reg$state_list[[length(reg$state_list)]]
    },

    .update_f2_initialization_from_f2k = function(f2k_fun, scale_f2) {
      # Build initialization for data2 in ORIGINAL parameter domain (rescale=FALSE)
      init <- do.call(
        initialize_pme,
        c(list(x = self$data2, d = self$d, rescale = FALSE), self$init_args_f2)
      )

      # Find new params of centers in SCALED domain ([0,1]^d) under warped embedding f2_k
      new_params_scaled <- calc_params(
        f = f2k_fun,
        X = init$centers,
        init_params = init$parameterization,
        f_input = "uv"
      )

      # Map scaled params back to ORIGINAL domain using (A, b) from current pme2 scaling
      # new_params_back = A %*% u + b  (row-wise)
      A <- scale_f2$A
      b <- scale_f2$b
      new_params_back <- t(A %*% t(new_params_scaled) + b)

      init$parameterization <- new_params_back
      init
    },

    .resolve_lambda_for_f2 = function() {
      if (self$lambda_policy_f2 == "reuse_prev") {
        # assumes pme2 has $tuning
        return(self$pme2$tuning)
      }
      if (self$lambda_policy_f2 == "fixed") {
        if (is.null(self$fixed_lambda_f2)) stop("fixed_lambda_f2 must be provided for lambda_policy_f2='fixed'.")
        return(self$fixed_lambda_f2)
      }
      # "retune": let pme pick its own (set lambda NULL)
      NULL
    },

    .refit_f2 = function(updated_init, lambda) {
      args <- c(list(
        data = self$data2,
        d = self$d,
        initialization = updated_init
      ), self$pme_args_f2)

      # Override lambda depending on policy
      if (!is.null(lambda)) args$lambda <- lambda

      do.call(pme, args)
    },

    .record_history = function(
    k,
    old_pme2,
    new_pme2,
    scale_f1,
    scale_f2,
    reg,
    last_state,
    lambda,
    record_full = TRUE
    ) {
      E_hist <- reg$E_history
      final_E <- if (!is.null(E_hist) && length(E_hist) > 0L) tail(E_hist, 1) else NA_real_

      entry <- list(
        k = k,
        lambda_f2 = lambda,
        n_iter = if (!is.null(E_hist)) length(E_hist) else NA_integer_,
        final_E = final_E,
        scale_f1 = if (record_full) scale_f1 else NULL,
        scale_f2 = if (record_full) scale_f2 else NULL,
        gamma = last_state$gamma_k,
        reg_E_history = if (record_full) E_hist else NULL,
        reg_last_state = if (record_full) last_state else NULL,
        pme2_before = if (record_full) old_pme2 else NULL,
        pme2_after = if (record_full) new_pme2 else NULL
      )

      self$history[[length(self$history) + 1L]] <- entry
      invisible(TRUE)
    },

    # -------------------------
    # Convergence rules
    # -------------------------
    .check_convergence_delta_E = function(tol_E) {
      if (length(self$history) < 2L) return(FALSE)
      E1 <- self$history[[length(self$history)]]$final_E
      E0 <- self$history[[length(self$history) - 1L]]$final_E
      if (is.na(E1) || is.na(E0)) return(FALSE)
      abs(E1 - E0) < tol_E
    }
  )
)
