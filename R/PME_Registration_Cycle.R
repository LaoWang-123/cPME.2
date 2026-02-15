#' PME + Registration iterative cycle (fixed f1, update f2)
#'
#' This R6 class orchestrates:
#' 1) Fit PME for data1 and data2 (f1 fixed; f2 updated)
#' 2) Rescale both parameterizations to [0,1]^d
#' 3) Run registration to warp f2 toward f1
#' 4) Update cached f2 initialization (parameterization of centers) using final warped embedding f2_k
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

    # init args for initialize_pme (run ONCE)
    init_args_f1 = NULL,
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

    f1_fun = NULL,
    f2_fun = NULL,
    f2_grad = NULL,

    initialization_f1 = NULL,
    initialization_f2 = NULL,  # cached initialization for f2 (centers fixed; parameterization updated)

    scale_f1 = NULL,           # cached scaling info for current pme1
    scale_f2 = NULL,           # scaling info for current pme2 (changes each cycle)

    reg = NULL,
    gamma = NULL,
    f2_warped = NULL,          # last state's f2_k (warped embedding)

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
    init_args_f1 = NULL,
    init_args_f2 = NULL,
    reg_args = list(),
    lambda_policy_f2 = c("reuse_prev", "fixed", "retune"),
    fixed_lambda_f2 = NULL,
    default_args = NULL
    ) {
      self$data1 <- data1
      self$data2 <- data2
      self$d <- as.integer(d)

      # ---- defaults (optional) ----
      defaults <- private$.get_default_args(default_args)

      # merge lists: user args override defaults
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

      # init args: if init_args_f1 is NULL, keep the original behavior (fallback to init_args_f2)
      init_args_f2_final <- if (is.null(init_args_f2)) defaults$init_args_f2 else init_args_f2
      init_args_f1_final <- if (is.null(init_args_f1)) defaults$init_args_f1 else init_args_f1

      self$init_args_f1 <- init_args_f1_final
      self$init_args_f2 <- init_args_f2_final

      self$reg_args <- modifyList(
        private$.null_to_list(defaults$reg_args),
        private$.null_to_list(reg_args),
        keep.null = TRUE
      )

      self$lambda_policy_f2 <- match.arg(lambda_policy_f2)
      self$fixed_lambda_f2 <- fixed_lambda_f2

      # allow setting policy defaults via default_args (only when user didn't pass these args)
      if (!is.null(defaults$lambda_policy_f2) && missing(lambda_policy_f2)) {
        self$lambda_policy_f2 <- defaults$lambda_policy_f2
      }
      if (!is.null(defaults$fixed_lambda_f2) && missing(fixed_lambda_f2)) {
        self$fixed_lambda_f2 <- defaults$fixed_lambda_f2
      }

      self$cycle_idx <- 0L
      self$history <- list()

      invisible(self)
    },

    # -------------------------
    # Public API
    # -------------------------

    #' Fit initial PME for both datasets (and cache initialization_f2)
    fit_initial = function() {

      # f1 init + fit
      self$initialization_f1 <- private$.init_pme(
        dataX = self$data1,
        init_args = self$init_args_f1,
        rescale = FALSE
      )
      self$pme1 <- private$.fit_pme(
        dataX = self$data1,
        pme_args = self$init_args_f1,
        initialization = self$initialization_f1
      )

      # f2 init + fit
      self$initialization_f2 <- private$.init_pme(
        dataX = self$data2,
        init_args = self$init_args_f2,
        rescale = FALSE
      )
      self$pme2 <- private$.fit_pme(
        dataX = self$data2,
        pme_args = self$pme_args_f2,
        initialization = self$initialization_f2
      )

      self$scale_f1 <- private$.compute_projection_and_scale(self$pme1, self$data1)
      self$f1_fun <- private$.make_scaled_embedding(self$pme1, self$scale_f1)

      invisible(self)
    },


    #' Run a single cycle: register + update cached f2 init + refit f2
    run_cycle = function(record_full = TRUE) {

      # auto-initialize if user didn't call fit_initial()
      if (is.null(self$pme1) || is.null(self$pme2)) {
        self$fit_initial()
      }
      if (is.null(self$initialization_f2)) {
        stop("initialization_f2 is NULL after fit_initial(). Check initialize_pme() call.")
      }

      k <- self$cycle_idx + 1L

      # 1) compute scaling for current f2 (must be recomputed each cycle)
      self$scale_f2 <- private$.compute_projection_and_scale(self$pme2, self$data2)

      # 2) make scaled embedding functions for registration

      self$f2_fun <- private$.make_scaled_embedding(self$pme2, self$scale_f2)
      self$f2_grad <- private$.make_scaled_grad(self$pme2, self$scale_f2)

      # 3) registration
      self$reg <- private$.register_once(
        f1_fun = self$f1_fun,
        f2_fun = self$f2_fun,
        grad_f2_fun = self$f2_grad
      )

      last_state <- private$.extract_last_state(self$reg)
      self$gamma <- last_state$gamma_k
      self$f2_warped <- last_state$f2_k

      # 4) update cached f2 initialization using final f2_k
      #    IMPORTANT: use pme2$params_opt as init_params (your updated design)
      private$.update_f2_initialization_from_f2k_inplace(
        f2k_fun = self$f2_warped,
        scale_f2 = self$scale_f2,
        init_params_for_centers = self$pme2$params_opt
      )
      updated_init <- self$initialization_f2

      # 5) refit f2
      lambda_to_use <- private$.resolve_lambda_for_f2()
      old_pme2 <- self$pme2

      self$pme2 <- private$.fit_pme(
        dataX = self$data2,
        pme_args = self$pme_args_f2,
        initialization = updated_init, # Use the new initialization (old centers but new params)
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

    #' Run multiple cycles (auto fit_initial if needed)
    run = function(
    n_cycles = 5,
    tol_E = NULL,
    stop_rule = c("none", "delta_E"),
    record_full = TRUE,
    reinit = FALSE
    ) {
      stop_rule <- match.arg(stop_rule)

      # if reinit==TRUE, force redo initial fit (useful if user changed args after construction)
      if (reinit) {
        private$.reset_state()
      }

      # auto fit if needed
      if (is.null(self$pme1) || is.null(self$pme2)) {
        self$fit_initial()
      }

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
        initialization_f2 = self$initialization_f2,
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
    # Helpers
    # -------------------------
    .null_to_list = function(x) {
      if (is.null(x)) list() else x
    },

    # Default argument set used when `default_args` is not provided.
    # User-provided args always override these defaults (via modifyList()).
    .get_default_args = function(default_args = NULL) {
      if (!is.null(default_args)) return(default_args)

      list(
        # PME defaults
        pme_args_f1 = list(
          initialization_rescale = FALSE
        ),
        pme_args_f2 = list(
          initialization_rescale = FALSE
        ),

        # initialize_pme() defaults (typical choices used in your qmd)
        init_args_f1 = NULL,
        init_args_f2 = list(
          min_clusters = 10,
          alpha = 0.01,
          max_clusters = 100,
          algorithm = "isomap"
          # rescale is controlled by .init_pme(rescale = FALSE) in this class
        ),

        # registration defaults (leave basis_set / Ugrid to user)
        reg_args = list(
          eps_step = 0.005,
          eps_energy = 0.005,
          max_iter = 15,
          basis_mode = "div_free"
        ),

        # optional: allow setting these through defaults as well
        lambda_policy_f2 = NULL,
        fixed_lambda_f2 = NULL
      )
    },

    .reset_state = function() {
      self$cycle_idx <- 0L
      self$pme1 <- NULL
      self$pme2 <- NULL
      self$initialization_f2 <- NULL
      self$scale_f1 <- NULL
      self$scale_f2 <- NULL
      self$reg <- NULL
      self$gamma <- NULL
      self$f2_warped <- NULL
      self$converged <- FALSE
      self$stop_reason <- NULL
      self$history <- list()
      invisible(TRUE)
    },

    # -------------------------
    # Internal steps
    # -------------------------

    # Generic initializer (returns an initialization object)
    .init_pme = function(dataX, init_args = NULL, rescale = FALSE) {
      args <- c(
        list(x = dataX, d = self$d, rescale = rescale),
        private$.null_to_list(init_args)
      )
      # avoid user overwriting
      args$x <- dataX
      args$d <- self$d
      args$rescale <- rescale

      do.call(initialize_pme, args)
    },

    # Generic pme fitter (returns a pme result)
    .fit_pme = function(dataX, pme_args = NULL, initialization = NULL, lambda = NULL) {
      args <- c(
        list(data = dataX, d = self$d),
        private$.null_to_list(pme_args)
      )
      # avoid user overwriting
      args$data <- dataX
      args$d <- self$d

      if (!is.null(initialization)) args$initialization <- initialization
      if (!is.null(lambda)) args$lambda <- lambda

      do.call(pme, args)
    },

    .compute_projection_and_scale = function(pme_result, dataX, init_guess_method = "pca") {
      init_param <- pme_initial_guess(X = dataX, d = self$d, method = init_guess_method)
      U_proj <- calc_params(
        f = pme_result$embedding_map,
        X = dataX,
        init_params = init_param,
        f_input = "vector"
      )  # when dealing with pme_result, we used f_input = "vector"
      scaled <- scale_uniform_square_with_params(U = U_proj) # to [0,1]^2

      list(
        A = scaled$A,
        b = scaled$b
      )
    },

    .make_scaled_embedding = function(pme_result, scale_info) {
      pme_embedding_factory(
        pme_result = pme_result,
        d = self$d,
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

    .register_once = function(f1_fun, f2_fun, grad_f2_fun) {
      reg <- do.call(
        Registration_new$new,
        c(private$.null_to_list(self$reg_args),
          list(
          f1 = f1_fun,
          f2 = f2_fun,
          f2_grad_fn = grad_f2_fun)
          )
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

    # In-place update: initialization_f2$parameterization <- new_params_back
    # using X = initialization_f2$centers and init_params = pme2$params_opt (your new design).
    .update_f2_initialization_from_f2k_inplace = function(f2k_fun, scale_f2, init_params_for_centers) {
      if (is.null(self$initialization_f2)) {
        stop("initialization_f2 is NULL. Call $fit_initial() first.")
      }
      if (is.null(init_params_for_centers)) {
        stop("init_params_for_centers is NULL. Expected self$pme2$params_opt.")
      }

      new_params_centers <- calc_params(
        f = f2k_fun,
        X = self$initialization_f2$centers,
        init_params = init_params_for_centers,
        f_input = "uv"
      ) # when dealing with reg output, we used f_input = "uv"

      A <- scale_f2$A
      b <- scale_f2$b
      new_params_back <- t(A %*% t(new_params_centers) + b)

      self$initialization_f2$parameterization <- new_params_back
      invisible(TRUE)
    },

    .resolve_lambda_for_f2 = function() {
      if (self$lambda_policy_f2 == "reuse_prev") {
        return(self$pme2$tuning)
      }
      if (self$lambda_policy_f2 == "fixed") {
        if (is.null(self$fixed_lambda_f2)) {
          stop("fixed_lambda_f2 must be provided for lambda_policy_f2='fixed'.")
        }
        return(self$fixed_lambda_f2)
      }
      NULL
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
