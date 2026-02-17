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
PMERegistrationCycle <- R6::R6Class(classname = "PMERegistrationCycle",

    public = list(
      #' @field data1 Numeric matrix or array for the first dataset (fixed).
      #' @field data2 Numeric matrix or array for the second dataset (moving).
      #' @field d Integer. Intrinsic dimension of the parameterization domain.
      #'
      #' @field pme_args_f1 List of arguments passed to \code{pme()} when fitting the first dataset.
      #' @field pme_args_f2 List of arguments passed to \code{pme()} when fitting the second dataset.
      #'
      #' @field init_args_f1 List of arguments passed to \code{initialize_pme()} for the first dataset (executed once).
      #' @field init_args_f2 List of arguments passed to \code{initialize_pme()} for the second dataset (executed once).
      #'
      #' @field reg_args List of arguments controlling the registration step.
      #'
      #' @field lambda_policy_f2 Character. Strategy for selecting lambda in PME fitting of data2.
      #'        One of \code{"reuse_prev"}, \code{"fixed"}, or \code{"retune"}.
      #' @field fixed_lambda_f2 Numeric. Fixed lambda value used when
      #'        \code{lambda_policy_f2 == "fixed"}.
      #'
      #' @field cycle_idx Integer. Current registration cycle index.
      #'
      #' @field pme1 Fitted PME object for \code{data1}.
      #' @field pme2 Fitted PME object for \code{data2}.
      #'
      #' @field f1_fun Function. Current embedding function associated with \code{pme1}.
      #' @field f2_fun Function. Current embedding function associated with \code{pme2}.
      #' @field f2_grad Function. Gradient function corresponding to \code{f2_fun}.
      #'
      #' @field initialization_f1 List. Cached initialization object for PME fitting of \code{data1}.
      #' @field initialization_f2 List. Cached initialization object for PME fitting of \code{data2}.
      #'        The centers remain fixed while parameterization may be updated across cycles.
      #'
      #' @field scale_f1 List. Affine scaling information mapping the parameterization
      #'        of \code{pme1} to the unit domain.
      #' @field scale_f2 List. Affine scaling information mapping the parameterization
      #'        of \code{pme2} to the unit domain (updated each cycle).
      #'
      #' @field reg Registration object storing the current registration model.
      #' @field f2_warped Function. Warped embedding of \code{f2} under the current \code{gamma}.
      #'
      #' @field converged Logical. Indicates whether the alternating optimization has converged.
      #' @field stop_reason Character. Text description of the stopping condition.
      #'
      #' @field history List. Stores per-cycle state information.
      #' @field init_history List. Stores initialization-stage information.
      #' @field save_dir Character. Directory path for saving intermediate states.
      #' @field filename Character. File name for autosaving registration state.
      #' @field verbose printing log.

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
      f2_warped = NULL,          # last state's f2_k (warped embedding)

      converged = FALSE,
      stop_reason = NULL,

      # -------------------------
      # History
      # -------------------------
      history = NULL,
      init_history = NULL,
      save_dir = NULL,
      filename = NULL,
      verbose = TRUE,

      # -------------------------
      # Constructor/Initialize
      # -------------------------
      #' Initialize a PME Registration Cycle Object
      #'
      #' Constructs a new \code{PME_Registration_Cycle} object and sets up
      #' configuration, initialization parameters, and registration controls
      #' for the alternating PME–registration procedure.
      #'
      #' @param data1 Numeric matrix or array representing the first dataset (fixed).
      #' @param data2 Numeric matrix or array representing the second dataset (moving).
      #' @param d Integer. Intrinsic dimension of the parameterization domain.
      #' @param pme_args_f1 List of arguments passed to \code{pme()} when fitting \code{data1}.
      #' @param pme_args_f2 List of arguments passed to \code{pme()} when fitting \code{data2}.
      #' @param init_args_f1 List of arguments passed to \code{initialize_pme()} for \code{data1}.
      #' @param init_args_f2 List of arguments passed to \code{initialize_pme()} for \code{data2}.
      #' @param reg_args List of arguments controlling the registration step.
      #' @param lambda_policy_f2 Character string specifying the lambda selection
      #' strategy for PME fitting of \code{data2}. One of
      #' \code{"reuse_prev"}, \code{"fixed"}, or \code{"retune"}.
      #' @param fixed_lambda_f2 Numeric value used when
      #' \code{lambda_policy_f2 = "fixed"}.
      #' @param default_args Optional list of default arguments applied to PME fitting.
      #' @param save_dir Optional character string specifying a directory for
      #' autosaving intermediate states.
      #' @param filename Optional character string specifying the filename used
      #' @param verbose printing log
      #' when saving registration state.
      #'
      #' @return A new \code{PME_Registration_Cycle} R6 object.
      #'
      #' @examples
      #' \dontrun{
      #' obj <- PME_Registration_Cycle$new(
      #'   data1 = X1,
      #'   data2 = X2,
      #'   d = 2
      #' )
      #' }
      initialize = function(
      data1,
      data2,
      d = 2,
      pme_args_f1 = list(),
      pme_args_f2 = list(),
      init_args_f1 = list(),
      init_args_f2 = list(),
      reg_args = list(),
      lambda_policy_f2 = c("reuse_prev", "fixed", "retune"),
      fixed_lambda_f2 = NULL,
      default_args = NULL,
      save_dir = NULL,
      filename = NULL,
      verbose = TRUE
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

        # init args:
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

        # reg args:
        self$reg_args <- modifyList(
          private$.null_to_list(defaults$reg_args),
          private$.null_to_list(reg_args),
          keep.null = TRUE
        )

        self$lambda_policy_f2 <- match.arg(lambda_policy_f2)
        self$fixed_lambda_f2 <- fixed_lambda_f2

        self$cycle_idx <- 0L
        self$history <- list(
          initial = NULL,
          cycles  = list()
        )

        self$save_dir <- save_dir
        self$filename <- filename
        self$verbose <- isTRUE(verbose)

        invisible(self)
      },

      # -------------------------
      # Public API
      # -------------------------

      #' Fit Initial PME Models for Both Datasets
      #'
      #' Performs the initialization and fitting steps of PME for
      #' \code{data1} (fixed) and \code{data2} (moving).
      #'
      #' For each dataset, this method:
      #' \enumerate{
      #'   \item Calls the internal PME initialization routine,
      #'   \item Fits the PME model using the provided arguments,
      #'   \item Stores the fitted PME objects in \code{pme1} and \code{pme2}.
      #' }
      #'
      #' The initialization object for \code{data2} is cached in
      #' \code{initialization_f2} for reuse across registration cycles.
      #'
      #' After fitting, affine scaling information for \code{pme1}
      #' is computed and the corresponding embedding function
      #' \code{f1_fun} is constructed.
      #'
      #' The initial state is recorded in \code{init_history}.
      #'
      #' @return The updated \code{PME_Registration_Cycle} object (invisibly).
      fit_initial = function() {

        # f1 init + fit
        self$initialization_f1 <- private$.init_pme(
          dataX = self$data1,
          init_args = self$init_args_f1,
          rescale = FALSE
        )
        self$pme1 <- private$.fit_pme(
          dataX = self$data1,
          pme_args = self$pme_args_f1,
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
        self$f1_fun <- private$.make_scaled_embedding(self$pme1, self$scale_f1, d = self$d)

        self$scale_f2 <- private$.compute_projection_and_scale(self$pme2, self$data2)
        self$f2_fun <- private$.make_scaled_embedding(self$pme2, self$scale_f2, d = self$d)

        # initial history
        self$history$initial <- list(
          k = 0,
          stage = "init",
          pme1 = self$pme1,
          pme2 = self$pme2,
          f1_fun = self$f1_fun,
          f2_fun = self$f2_fun,
          scale_f1 = self$scale_f1,
          scale_f2 = self$scale_f2
        )

        invisible(self)
      },


      #' Run One Alternating Registration Cycle
      #' register + update cached f2 init + refit f2
      #' Executes a single iteration of the alternating
      #' PME–registration procedure.
      #'
      #' If initial PME models have not yet been fitted,
      #' this method automatically calls \code{fit_initial()}.
      #'
      #' A cycle consists of:
      #' \enumerate{
      #'   \item Performing registration between the current embeddings,
      #'   \item Updating the cached initialization of \code{data2}
      #'         according to the estimated reparameterization,
      #'   \item Refitting PME for \code{data2} using the updated initialization.
      #' }
      #'
      #' The cycle index \code{cycle_idx} is incremented and
      #' internal state variables (including \code{pme2}, \code{gamma},
      #' scaling information, and embedding functions) are updated.
      #'
      #'
      #' @return The updated \code{PME_Registration_Cycle} object (invisibly).
      run_cycle = function() {

        # auto-initialize if user didn't call fit_initial()
        if (is.null(self$pme1) || is.null(self$pme2)) {
          self$fit_initial()
        }
        if (is.null(self$initialization_f2)) {
          stop("initialization_f2 is NULL after fit_initial(). Check initialize_pme() call.")
        }

        k <- self$cycle_idx + 1L

        # 1) compute scaling for current f2 (must be recomputed each cycle)
        # 2) make scaled embedding functions for registration
        # These two steps have completed in initial_fit once.

        # 3) registration
        self$reg <- private$.register_once(
          f1_fun = self$f1_fun,
          f2_fun = self$f2_fun,
          grad_f2_fun = self$f2_grad  # self$reg_args are included in the private function
        )

        last_state <- private$.extract_last_state(self$reg)
        self$f2_warped <- last_state$f2_k

        # 4) update cached f2 initialization using final f2_k
        #    IMPORTANT: use pme2$params_opt as init_params (your updated design)
        private$.update_f2_initialization_from_f2k_inplace(
          f2k_fun = self$f2_warped,
          scale_f2 = self$scale_f2,
          center_points = self$initialization_f2$centers,
          init_params_for_centers = self$pme2$params_opt
        )
        updated_init <- self$initialization_f2

        # 5) refit f2
        lambda_to_use <- private$.resolve_lambda_for_f2()

        self$pme2 <- private$.fit_pme(
          dataX = self$data2,
          pme_args = self$pme_args_f2,
          initialization = updated_init, # Use the new initialization (old centers but new params)
          lambda = lambda_to_use # Here input the initialization and lambda will override the pme_args_f2
        )

        # 1) compute scaling for current f2 (must be recomputed each cycle)
        self$scale_f2 <- private$.compute_projection_and_scale(self$pme2, self$data2)

        # 2) make scaled embedding functions for registration

        self$f2_fun <- private$.make_scaled_embedding(self$pme2, self$scale_f2,d = self$d)
        self$f2_grad <- private$.make_scaled_grad(self$pme2, self$scale_f2)

        # 6) record history
        private$.record_history(
          k = k,
          reg = self$reg,
          new_pme2 = self$pme2,
          scale_f2 = self$scale_f2,
          f2_fun = self$f2_fun,
          lambda = lambda_to_use
        )


        self$cycle_idx <- k
        invisible(self)
      },


      #' Executes multiple alternating PME–registration cycles.
      #'
      #' This method serves as the main driver of the registration
      #' algorithm. It optionally reinitializes the object state,
      #' performs initial PME fitting if necessary, and then iteratively
      #' calls \code{run_cycle()}.
      #'
      #' During each cycle, intermediate results may be saved to disk,
      #' and convergence can be checked based on a user-specified
      #' stopping rule.
      #'
      #' @param n_cycles Integer. Maximum number of alternating cycles.
      #' @param save_dir Optional character string specifying a directory
      #'        for saving intermediate states.
      #' @param filename Optional character string specifying the filename
      #'        used when saving registration state.
      #' @param tol_E Optional numeric tolerance used when
      #'        \code{stop_rule = "delta_E"}.
      #' @param stop_rule Character string specifying the stopping criterion.
      #'        One of \code{"none"} (no early stopping) or
      #'        \code{"delta_E"} (stop when energy change is below \code{tol_E}).
      #' @param reinit Logical. If \code{TRUE}, resets the internal state
      #'        and refits the initial PME models before running cycles.
      #'
      #' @details
      #' If PME models have not yet been fitted, \code{fit_initial()}
      #' is called automatically. The object fields \code{converged}
      #' and \code{stop_reason} are updated if early stopping occurs.
      #'
      #' @return The updated \code{PME_Registration_Cycle} object (invisibly).
      run = function(
      n_cycles = 5,
      save_dir = NULL,
      filename = NULL,
      tol_E = NULL,
      stop_rule = c("none", "delta_E"),
      reinit = FALSE
      ) {
        stop_rule <- match.arg(stop_rule)
        # ---- how many cycles already completed? ----
        completed <- length(self$history$cycles)
        start_cycle <- completed + 1L
        end_cycle   <- completed + as.integer(n_cycles)

        # print settings
        if (self$verbose) {
          ra <- private$.null_to_list(self$reg_args)

          cat("========================================\n")
          cat("PMERegistrationCycle START\n")
          cat(sprintf("already_done : %d\n", as.integer(completed)))
          cat(sprintf("n_cycles     : %d\n", as.integer(n_cycles)))
          cat(sprintf("will_run     : %d -> %d (overall)\n", as.integer(start_cycle), as.integer(end_cycle)))
          cat(sprintf("d            : %d\n", as.integer(self$d)))
          cat(sprintf("eps_step     : %s\n", as.character(ra$eps_step)))
          cat(sprintf("eps_energy   : %s\n", as.character(ra$eps_energy)))
          cat(sprintf("max_iter     : %s\n", as.character(ra$max_iter)))
          cat("========================================\n")
          flush.console()
        }

        # new save_dir and new filename
        if (!is.null(save_dir)) {
          self$save_dir <- save_dir
        }

        if (!is.null(filename)) {
          self$filename <- filename
        }

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

          # print cycle number
          if (self$verbose) {
            overall_i <- start_cycle + i - 1L
            cat(sprintf("\n[PMEReg] cycle %d / %d (this run %d / %d)\n",
                        as.integer(overall_i), as.integer(end_cycle),
                        as.integer(i), as.integer(n_cycles)))
          }

          self$run_cycle()

          # save the result of our PME registration cycle
          if (!is.null(self$save_dir) && !is.null(self$filename)) {
            private$.save_history_overwrite(save_dir=self$save_dir,filename=self$filename)
          }

          if (stop_rule == "delta_E" && !is.null(tol_E)) {
            if (private$.check_convergence_delta_E(tol_E = tol_E)) {
              self$converged <- TRUE
              self$stop_reason <- sprintf("delta_E < %.3g", tol_E)
              break
            }
          }
        }

        invisible(self)
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

          # registration defaults (leave basis_set / Ugrid to user)
          reg_args = list(
            eps_step = 0.005,
            eps_energy = 0.005,
            max_iter = 10,
            basis_mode = "div_free",
            basis_set = build_basis_set(5,5,basis = neumann_basis),
            Ugrid = subset(expand.grid(
              u = seq(0, 1, length.out = 60),
              v = seq(0, 1, length.out = 60)),
              (u - 0.5)^2 + (v - 0.5)^2 <= 0.5^2))# We need to define basis and Ugrid manually
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
        self$f2_warped <- NULL
        self$converged <- FALSE
        self$stop_reason <- NULL
        self$history <- list(
          initial = NULL,
          cycles  = list()
        )
        invisible(TRUE)
      },

      # -------------------------
      # Internal steps
      # -------------------------

      # Generic initializer (returns an initialization object)
      .init_pme = function(dataX, init_args = NULL, rescale = FALSE) {

        args <- modifyList(
          list(x = dataX, d = self$d, rescale = rescale), # The parameter names must be corresponded to the function parameters correctly.
          private$.null_to_list(init_args)
        )

        do.call(initialize_pme, args)
      },

      # Generic pme fitter (returns a pme result)
      .fit_pme = function(dataX, pme_args = NULL, initialization = NULL, lambda = NULL) {

        args <- modifyList(
          list(data = dataX, d = self$d),
          private$.null_to_list(pme_args) # The parameter names must be corresponded to the function parameters correctly.
        )

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

      .make_scaled_embedding = function(pme_result, scale_info,d=2) {
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

        reg <- do.call(Registration_new$new, args)
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
      .update_f2_initialization_from_f2k_inplace = function(f2k_fun, scale_f2, center_points, init_params_for_centers) {
        if (is.null(self$initialization_f2)) {
          stop("initialization_f2 is NULL. Call $fit_initial() first.")
        }
        if (is.null(init_params_for_centers)) {
          stop("init_params_for_centers is NULL. Expected self$pme2$params_opt.")
        }

        new_params_centers <- calc_params(
          f = f2k_fun,
          X = center_points,
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
        if (self$lambda_policy_f2 == "retune") {
          return(exp(-15:5))
        }
      },

      .record_history = function(
      k,
      reg,
      new_pme2,
      scale_f2,
      f2_fun,
      lambda
      ) {
        E_hist <- reg$E_history
        final_E <- if (!is.null(E_hist) && length(E_hist) > 0L) tail(E_hist, 1) else NA_real_

        entry <- list(
          k = k,
          reg = reg,
          lambda_f2 = lambda,
          n_iter = length(E_hist),
          final_E = final_E,
          scale_f2 = scale_f2,
          f2_fun = f2_fun,
          reg_E_history = E_hist,
          new_pme2 = new_pme2
        )

        # cycle history
        self$history$cycles[[length(self$history$cycles) + 1L]] <- entry
        invisible(TRUE)
      },

      ### save function
      .save_history_overwrite = function(save_dir,filename) {
        if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)

        saveRDS(self$history, file = file.path(save_dir, filename))
        invisible(TRUE)
      },

      # -------------------------
      # Convergence rules
      # -------------------------
      .check_convergence_delta_E = function(tol_E) {
        if (length(self$history$cycles) < 2L) return(FALSE)
        E1 <- self$history$cycles[[length(self$history$cycles)]]$final_E
        E0 <- self$history$cycles[[length(self$history$cycles) - 1L]]$final_E
        if (is.na(E1) || is.na(E0)) return(FALSE)
        abs(E1 - E0) < tol_E
      }
    )
  )
