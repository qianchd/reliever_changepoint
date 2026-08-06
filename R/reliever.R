#' Reliever changepoint detection with a custom interval-loss function
#'
#' Run a changepoint search with a user-supplied interval loss. This is the
#' primary fitting interface for custom models. The supplied \code{reg_fun}
#' fits a model using observations \code{l:r} and returns the loss of every
#' observation in \code{l_end:r_end}. For a built-in loss family, use
#' \code{reliever()} instead. Model-specific \code{reliever_*()} functions
#' expose the built-in implementations with complete formal argument lists.
#'
#' @param data A vector, matrix, data frame, or model-specific data object with
#'   observations in rows. The object must be
#'   accepted by `reg_fun`.
#' @param reg_fun A function following the interval-loss contract
#'   below. It fits on \code{l:r} and returns observation-level losses on
#'   \code{l_end:r_end}.
#' @param cpn_max,dm,cov_rate,method,cpn_crit Common path controls; see
#'   \code{\link{reliever}()}.
#' @param pen_val,prune_value,M,wbs_seed,wbs_stop_crit Common penalty and
#'   search controls; see \code{\link{reliever}()}.
#' @param detail,cache_backend,owner_key,echo Common output and cache controls;
#'   see \code{\link{reliever}()}.
#' @param dc_grid_size,dc_grid Common candidate-grid controls; see
#'   \code{\link{reliever}()}.
#' @param cache_profile \code{NULL} or a reusable cache-profile object returned
#'   by a previous call
#'   made with \code{detail = TRUE}. Reuse is valid only when the data,
#'   \code{reg_fun}, its arguments, search grid, \code{cov_rate}, and
#'   \code{cache_backend} are unchanged.
#' @param run_cpd_ids \code{NULL} or an integer vector of identifiers from the
#'   \code{loss_output_id} column of the virtual metadata table, selecting
#'   which output losses to search for changepoints. The default \code{NULL}
#'   fits changepoints for all outputs.
#' @param ... Additional arguments passed unchanged to \code{reg_fun}.
#'
#' @section Writing a custom reg_fun:
#' A single \code{reg_fun} definition supplies its loss-output count, optional
#' metadata, and all interval losses. \code{reliever_generic()} uses it in two
#' ways:
#' \itemize{
#'   \item Before the search, it calls
#'   \code{reg_fun(data, l = 1L, r = 1L, is_virtual_run = TRUE, ...)}.
#'   This call must return either a positive integer giving the number of loss
#'   outputs, or a list with \code{n_loss_outputs} and optional
#'   \code{loss_output_meta}.
#'   \item During the search, \code{reg_fun} receives fitting endpoints,
#'   evaluation endpoints, and \code{save_model = FALSE}.
#'   Here \code{l:r} is the fitting interval and \code{l_end:r_end} is the
#'   interval whose observations must be scored. The result must be a list with
#'   a numeric \code{loss} matrix having one row for every observation in
#'   \code{l_end:r_end} and one column for every declared loss output. For one
#'   loss output, \code{loss} may instead be a numeric vector of that length.
#' }
#' This is the complete minimum contract: a positive output count for the
#' virtual call and correctly sized losses for each fitted interval. All
#' metadata are optional.
#'
#' A virtual call may additionally return \code{loss_output_meta}, with one row
#' per output. Its standard columns are:
#' \itemize{
#'   \item \code{loss_output_id} and \code{row_type};
#'   \item \code{hyper_value} and its display label \code{hyper_name};
#'   \item \code{default_selection}; and
#'   \item \code{loss_kind}. Built-in labels are \code{"rss"},
#'   \code{"negative_log_likelihood"}, and \code{"crossfit_loss"}.
#'   Custom losses may use another non-empty label.
#' }
#' These fields support output filtering, readable summaries, default plotting,
#' and compatible model selection. Omitted output identifiers are generated
#' automatically, unrecognized columns are retained, and missing
#' \code{loss_kind} leaves scale compatibility to the author.
#' This model-independent contract gives custom \code{reg_fun} implementations
#' the common result format and public post-fit selectors. The fitted result also
#' retains the function and its resolved arguments for external segment
#' evaluation; \code{select_holdout()} therefore does not require them to be
#' entered a second time.
#'
#' @return A list of class \code{reliever_result}. Its main components are:
#'   \describe{
#'     \item{\code{summary}}{A data frame containing any explicitly requested fit-time selection;
#'     empty by default for a K-indexed path.}
#'     \item{\code{cpd_path}}{A list containing every distinct candidate segmentation, including
#'     its K, changepoints, and stored loss.}
#'     \item{\code{run_meta}}{A data frame containing the identifiers and optional metadata for each
#'     loss output returned by \code{reg_fun}.}
#'     \item{\code{settings}}{A named list containing the resolved search, grid, cache, and loss
#'     arguments.}
#'     \item{\code{timing}}{A data frame containing model-fit counts and elapsed times by run.}
#'   }
#'   Requested diagnostics and reusable cache objects are included when
#'   available. Use \code{summary()}, \code{plot()}, and the post-fit selectors
#'   to work with the result.
#' @references
#' Qian, C., Wang, G., and Zou, C. (2025). Reliever: Relieving the burden of
#' costly model fits for changepoint detection. \emph{Journal of Machine
#' Learning Research}, 26(203), 1--57.
#'
#' Qian, C., Wang, G., Wang, Z., and Zou, C. (2024). Changepoint detection in
#' complex models: Cross-fitting is needed. arXiv:2411.07874.
#' @seealso \code{\link{reliever}()},
#'   \code{\link{cv.reliever_generic}()}, and
#'   \code{\link{reg_fun_crossfit_template}()}.
#' @export
#'
#' @examples
#' set.seed(2026)
#' n_seg <- 150
#' x <- rbind(
#'   matrix(rnorm(n_seg * 5, mean = 0, sd = 0.5), n_seg, 5),
#'   matrix(rnorm(n_seg * 5, mean = 4, sd = 0.5), n_seg, 5),
#'   matrix(rnorm(n_seg * 5, mean = -4, sd = 0.5), n_seg, 5)
#' )
#' custom_mean_loss <- function(data, l, r, l_end = l, r_end = r,
#'                              save_model = FALSE,
#'                              is_virtual_run = FALSE) {
#'   if (is_virtual_run) {
#'     return(1L)
#'   }
#'   data <- as.matrix(data)
#'   mu <- colMeans(data[l:r, , drop = FALSE])
#'   residual <- sweep(data[l_end:r_end, , drop = FALSE], 2L, mu)
#'   list(loss = rowMeans(residual^2))
#' }
#' res <- reliever_generic(
#'   data = x, reg_fun = custom_mean_loss,
#'   cpn_max = 7, dm = 30, cov_rate = 0.8, method = "SN"
#' )
#' summary(res) # empty: the generic fit does not infer a selection rule
#'
#' # Explicit RSS-SIC is available when the custom loss is known to be RSS.
#' selected_sic <- select_by_run(result = res, cpn_crit = "rss_sic")
#' selected_sic
#' stopifnot(identical(selected_sic$K_hat, 2L))
#' res$cpd_path$candidates
reliever_generic <- function(data, reg_fun, cpn_max = 3, dm = 50, cov_rate = 0.8,
                             method = "SN", cpn_crit = "none",
                             pen_val = 1, prune_value = 0,
                             M = 100, wbs_seed = NULL,
                             wbs_stop_crit = NULL,
                             detail = FALSE,
                             cache_backend = "by_loss_block",
                             owner_key = TRUE,
                             echo = FALSE,
                             dc_grid_size = NULL,
                             dc_grid = NULL,
                             cache_profile = NULL, run_cpd_ids = NULL,
                             ...) {
  if (missing(reg_fun)) {
    stop("reg_fun must be supplied to reliever_generic().", call. = FALSE)
  }
  if (!is.logical(owner_key) || length(owner_key) != 1L ||
      is.na(owner_key)) {
    stop("owner_key must be TRUE or FALSE.", call. = FALSE)
  }
  dc_grid <- .reliever_resolve_dc_grid(
    data = data,
    dc_grid_size = dc_grid_size,
    dc_grid = dc_grid
  )
  if (!is.null(dc_grid_size)) {
    dc_grid_size <- as.integer(dc_grid_size)
  }
  para_list <- list(...)
  .reliever_reject_renamed_cpn_max(para_list, "reliever_generic")
  setup <- .reliever_prepare(
    data = data,
    cpn_max = cpn_max,
    dm = dm,
    pen_val = pen_val,
    cov_rate = cov_rate,
    reg_fun = reg_fun,
    method = method,
    M = M,
    prune_value = prune_value,
    run_cpd_ids = run_cpd_ids,
    dc_grid = dc_grid,
    cache_backend = cache_backend,
    wbs_stop_crit = wbs_stop_crit,
    para_list = para_list
  )
  cpn_penalty_info <- .cpn_penalty(cpn_crit, setup$n_original)
  detail_mode <- .reliever_detail_mode(detail)
  return_cache_profile <- detail_mode != "none"
  individual_loss_fun <- .reliever_make_individual_loss_fun(
    data = data,
    reg_fun = reg_fun,
    dc_grid = setup$dc_grid,
    para_list = para_list,
    n_loss_outputs = setup$n_loss_outputs
  )
  loss_context <- list(
    data = data,
    reg_fun = reg_fun,
    para_list = para_list,
    dc_grid = setup$dc_grid
  )

  if (echo) timer_start <- proc.time()[["elapsed"]]
  use_owner_key <- setup$cache_backend == "by_loss_block" && owner_key
  request <- list(
    setup = setup,
    individual_loss_fun = individual_loss_fun,
    pen_val = pen_val,
    wbs_seed = wbs_seed,
    backend_cache = .reliever_backend_cache_from_profile(
      cache_profile, setup, loss_context
    ),
    loss_context = loss_context,
    return_cache_profile = return_cache_profile,
    owner_key = use_owner_key
  )
  backend_res <- if (setup$cache_backend == "by_loss_block") {
    .reliever_run_by_loss_block_backend(request)
  } else {
    .reliever_run_by_cost_mat_backend(request)
  }

  if (echo) {
    message(sprintf(
      "reliever-%s: %.3f sec elapsed",
      setup$method, proc.time()[["elapsed"]] - timer_start
    ))
  }

  out <- list(
    cpd_path = backend_res$cpd_path,
    run_meta = setup$run_meta,
    n_model_fit = backend_res$n_model_fit,
    model_fit_time = backend_res$model_fit_time,
    cpd_time = backend_res$cpd_time,
    total_time = backend_res$total_time,
    diagnostics = if (detail_mode == "none") NULL else backend_res$diagnostics,
    cache_profile = backend_res$cache_profile,
    settings = .reliever_result_settings(
      setup = setup,
      cpn_crit = cpn_crit,
      cpn_penalty_info = cpn_penalty_info,
      pen_val = pen_val,
      wbs_seed = wbs_seed,
      owner_key = use_owner_key,
      dc_grid_size = dc_grid_size
    ),
    loss_spec = list(reg_fun = reg_fun, args = para_list)
  )
  .reliever_finalize_result(out)
}

.reliever_response_input_spec <- function(X, y) {
  X <- as.matrix(X)
  list(
    type = "response_predictor",
    n_predictors = if (is.null(y)) ncol(X) - 1L else ncol(X),
    original_form = if (is.null(y)) "response_first" else "separate_xy"
  )
}

.reliever_set_input_spec <- function(result, input_spec) {
  if (!is.list(result) || is.null(result$settings)) {
    return(result)
  }
  result$settings$input_spec <- input_spec
  result
}

.reliever_match_cpd_family <- function(cpd_family, choices, caller) {
  tryCatch(
    match.arg(cpd_family, choices),
    error = function(e) {
      stop(
        "cpd_family must be one of ",
        paste0('"', choices, '"', collapse = ", "),
        " in ", caller, "().",
        call. = FALSE
      )
    }
  )
}

.reliever_prepare_input <- function(X, y, cpd_family, response_families) {
  if (cpd_family %in% response_families) {
    if (is.null(y)) {
      warning(
        "cpd_family = \"", cpd_family, "\" expects a response y; ",
        "using the first column of X as the response. Pass y explicitly ",
        "to make the response/predictor split unambiguous.",
        call. = FALSE
      )
      return(as.matrix(X))
    }
    X <- as.matrix(X)
    y <- as.numeric(y)
    if (NROW(X) != length(y)) {
      stop("X and y must have the same number of observations.", call. = FALSE)
    }
    return(cbind(y, X))
  }

  if (!is.null(y)) {
    warning(
      "y is ignored for cpd_family = \"", cpd_family,
      "\". Omit y, or choose a response-based cpd_family if y should ",
      "define the response.",
      call. = FALSE
    )
  }
  X
}

#' Reliever changepoint detection with built-in loss families
#'
#' Primary fitting interface for built-in Reliever methods.
#' \code{cpd_family = "mean"} fits mean changes by default; the alternatives
#' cover regression, density, nonparametric, and classifier losses.
#' \code{cpn_crit = "none"} retains the candidate path for later selection.
#'
#' Use \code{cv.reliever()} for CPSS outer-CV selection with its built-in
#' single-path parametric losses or lasso path.
#'
#' Cross-fitted changepoint families support two ReCV workflows:
#' \itemize{
#'   \item For interval-adaptive ReCV, use \code{select_by_run()} with run type
#'   \code{"recv"} and criterion \code{"loss"}.
#'   \item Homogeneous-hyperparameter ReCV jointly selects one hyperparameter,
#'   K, and the changepoints from stored fixed-hyperparameter rows.
#' }
#' Apply \code{select_across_runs()} for the homogeneous workflow.
#'
#' Its run type is \code{"crossfit_homo_hyper"}; request both crossfit outputs
#' during fitting.
#'
#' If these unpenalized-loss rules select too many changepoints, replace
#' \code{cpn_crit = "loss"} with \code{cpn_crit = "sic"} to add a
#' \eqn{\log(n)K} penalty.
#'
#' Information criteria are available through
#' \code{select_by_run()}. Independent hold-out selection supports ordinary
#' losses plus the built-in lasso- and KDE-NLL-crossfit scorers; other crossfit
#' templates require a model-specific external scorer.
#'
#' Focused \code{reliever_*()} functions expose the complete model-specific
#' arguments for the same built-in methods. Use \code{reliever_generic()} for
#' a custom interval-loss function.
#'
#' @param X A numeric vector, matrix, data frame, distance/kernel matrix, or
#'   family-specific data object required by the selected \code{cpd_family}.
#'   For regression families
#'   with \code{y = NULL}, the first column of \code{X} is treated as the
#'   response; lasso families require at least two predictor columns because
#'   they use \code{glmnet}. KDE-NLL families ordinarily accept the original
#'   observations and compute one reusable distance matrix; a precomputed
#'   squared-distance matrix is also accepted. KDE-L2 accepts one precomputed
#'   kernel-evaluation matrix, or raw observations together with one fixed
#'   \code{kernel}; it evaluates an empirical-grid kernel-density L2 loss.
#'   Mean, covariance, mean-and-covariance, exponential-family, NMCD, and
#'   classifier families accept observations in rows.
#' @param y \code{NULL}, a numeric response vector, or for a grouped binomial
#'   GLM a two-column numeric matrix of successes and failures. Used by
#'   \code{"lm"}, \code{"glm"},
#'   \code{"lasso"}, and \code{"lasso_crossfit"}. A grouped binomial GLM may
#'   use a two-column matrix of successes and failures. It is ignored with a
#'   warning for non-regression families.
#' @param cpd_family A character scalar selecting the built-in loss family:
#'   \itemize{
#'     \item \code{"mean"}: ordinary mean changes.
#'     \item \code{"mean_crossfit"}: cross-fitted mean changes.
#'     \item \code{"var"} and \code{"meanvar"}: covariance changes around a
#'     fixed common mean, or simultaneous mean and covariance changes, using
#'     multivariate Gaussian likelihood; see \code{\link{reliever_var}()} and
#'     \code{\link{reliever_meanvar}()}.
#'     \item \code{"lm"} and \code{"glm"}: changes in linear or generalized
#'     linear regression relationships; see \code{\link{reliever_lm}()} and
#'     \code{\link{reliever_glm}()}.
#'     \item \code{"em"}: changes in one selected exponential-family
#'     distribution.
#'     \item \code{"lasso"} and \code{"lasso_crossfit"}: regression-coefficient
#'     changes using an in-sample lasso path or cross-fitted lasso losses; see
#'     \code{\link{reliever_lasso}()} and
#'     \code{\link{reliever_lasso_crossfit}()}. Built-in lasso fits use no
#'     intercept and no automatic predictor standardization.
#'     \item \code{"kde_nll"} and \code{"kde_nll_crossfit"}: distribution
#'     changes using Gaussian, radial-Laplace, or multivariate-Student KDE
#'     negative log-likelihood; see
#'     \code{\link{reliever_kde_nll}()} and
#'     \code{\link{reliever_kde_nll_crossfit}()}.
#'     \item \code{"kde_l2"}: distributional changes under an empirical-grid
#'     kernel-density L2 loss with one fixed kernel. See
#'     \code{\link{reliever_kde_l2}()}.
#'     \item \code{"nmcd"}: arbitrary univariate distribution changes using an
#'     empirical-CDF loss; see \code{\link{reliever_nmcd}()}.
#'     \item \code{"ranger_crossfit"}: cross-fitted random-forest changes.
#'     \item \code{"mlp_crossfit"}: cross-fitted MLP classifier changes.
#'   }
#' @param cpn_max A non-negative integer scalar giving the maximum number of
#'   changepoints considered by fixed-K
#'   searches. Set it above the largest plausible changepoint count.
#' @param dm A positive integer scalar giving the minimum segment length on the
#'   original sample scale. For a gridded
#'   search, Reliever uses a conservative grid-cell minimum that guarantees
#'   every returned segment contains at least \code{dm} observations.
#' @param cov_rate A numeric scalar in \code{(0, 1]} giving the fraction of
#'   eligible intervals covered by the deterministic
#'   representative-interval collection, in \code{(0, 1]}. Use 1 for a full
#'   search.
#' @param method A character scalar selecting the changepoint search algorithm.
#'   Choices are segment neighbourhood
#'   (\code{"SN"}), Reliever's globally greedy wild binary segmentation
#'   (\code{"WBS"}), recursive wild binary segmentation
#'   (\code{"WBS_recursive"}), seeded binary segmentation
#'   (\code{"SeedBS"}), binary segmentation (\code{"BS"}), optimal
#'   partitioning (\code{"OP"}), and PELT (\code{"PELT"}).
#'   PELT applies pruning; OP runs the same penalized dynamic program without
#'   pruning.
#' @param cpn_crit A character scalar or non-negative numeric scalar specifying
#'   the optional rule used to select K while constructing the fit.
#'   Every candidate has a stored \code{loss} and changepoint count K:
#'   \itemize{
#'     \item \code{"loss"} chooses the smallest stored loss without a K
#'     penalty. For recycled cross-validation, select the
#'     appropriate loss rows with \code{run_type} in a post-fit selector.
#'     \item \code{"aic"}, \code{"hqc"}, and \code{"sic"} minimize
#'     \eqn{loss + \gamma K}, with \eqn{\gamma=2},
#'     \eqn{2\log\log(n)}, and \eqn{\log(n)}, respectively. They apply directly
#'     to the stored loss, including RSS, and therefore retain its scale.
#'     \item \code{"rss_aic"}, \code{"rss_hqc"}, and \code{"rss_sic"}
#'     minimize \eqn{n\log(RSS/n)/2 + \gamma K}. This alternative log-RSS
#'     formulation may be preferable when residual variance is unknown or
#'     potentially heterogeneous.
#'     \item A non-negative number supplies \eqn{\gamma} directly in
#'     \eqn{loss + \gamma K}.
#'     \item \code{"none"}, the default, applies no loss-based K selection.
#'     For paths indexed directly by K, the compact summary is empty. For WBS
#'     thresholds or PELT/OP penalties, the summary still reports the
#'     segmentation associated with each supplied search value.
#'   }
#'   A criterion is applied separately within every selected loss path.
#'   Comparable paths can instead be selected jointly with
#'   \code{\link{select_across_runs}()}. In every formula, \eqn{n} is the original number of
#'   observations, including when a candidate grid compresses the search grid.
#'   The named presets penalize K only; they do not estimate
#'   model-family degrees of freedom. If \code{wbs_stop_crit} is supplied,
#'   selection compares only its mapped segmentations and the no-change
#'   baseline.
#' @param pen_val Non-negative numeric penalties for \code{"PELT"} or
#'   \code{"OP"}. Each value produces one candidate segmentation. Other
#'   methods ignore this argument.
#' @param prune_value A numeric scalar giving the constant in the PELT pruning
#'   inequality. The default 0
#'   gives the standard PELT rule; \code{method = "OP"} disables pruning.
#' @param M A positive integer scalar giving the maximum number of random
#'   intervals retained by \code{"WBS"} and
#'   \code{"WBS_recursive"}.
#' @param wbs_seed \code{NULL} or an integer scalar used as the seed for random
#'   WBS interval generation.
#' @param wbs_stop_crit Optional numeric stopping thresholds for WBS-family
#'   methods.
#'   Larger thresholds tend to retain fewer changepoints. When supplied,
#'   selection compares their mapped segmentations and the no-change baseline.
#' @param detail A logical or character scalar. Use \code{TRUE} or
#'   \code{"cache"} to retain cache state and
#'   search-method diagnostics; use \code{FALSE} or \code{"none"} for a compact
#'   result. For non-native families the returned cache can also be supplied as
#'   \code{cache_profile} in a later compatible fit; the native mean cache is
#'   retained for inspection only.
#' @param cache_backend A character scalar selecting the cost-cache
#'   implementation used by every built-in
#'   family. \code{"by_loss_block"} stores losses for fitted representative
#'   intervals, while \code{"by_cost_mat"} stores reconstructed interval
#'   costs in an expanded matrix. A full search automatically uses
#'   \code{"by_cost_mat"}.
#' @param owner_key A logical scalar. With
#'   \code{cache_backend = "by_loss_block"}, cache the
#'   mapping from each exact interval to its representative interval. This
#'   usually improves speed at a modest memory cost.
#' @param echo A logical scalar indicating whether to print timing messages.
#' @param dc_grid_size Regular candidate-grid spacing, or \code{NULL} to search
#'   every boundary. This coarse-grid option is motivated by the
#'   divide step of DCDP (Li, Wang, and Rinaldo, 2023); it does not apply the
#'   subsequent DCDP local-refinement step.
#' @param dc_grid \code{NULL} or an increasing integer vector providing an
#'   advanced explicit interface for an irregular grid
#'   ending at the number of observations. Supply only one of
#'   \code{dc_grid_size} and \code{dc_grid}.
#'
#' @param ... Model-specific arguments passed to the selected \code{reg_fun};
#'   see that function's help page.
#'
#' @return A list of class \code{reliever_result} with:
#'   \describe{
#'     \item{\code{summary}}{A data frame containing the compact selection shown when the object is
#'     printed and returned by \code{summary()}. With the ordinary
#'     \code{cpn_crit = "none"} default, this is empty for paths indexed
#'     directly by K. An explicit criterion selects K separately within each
#'     eligible loss output. A \code{hyper_value} identifies that row's model
#'     setting; it does not imply that the hyperparameter was compared across
#'     runs. For eligible loss outputs, supplied WBS thresholds or PELT/OP
#'     penalties are reported with their corresponding segmentations even when
#'     \code{cpn_crit = "none"}.}
#'     \item{\code{cpd_path}}{A list containing all distinct candidate segmentations. Its
#'     \code{candidates} table contains K, loss, and changepoint locations. For
#'     PELT/OP penalties or WBS stopping thresholds, \code{selector_map}
#'     records which supplied value produced each candidate.}
#'     \item{\code{run_meta}}{A data frame containing metadata describing each searched loss output.
#'     It always links the exact numeric \code{run_id} to
#'     \code{loss_output_id}; built-in and custom losses may also declare
#'     \code{row_type}, a hyperparameter, and \code{loss_kind}. Post-fit
#'     selectors can match declared \code{row_type} labels through
#'     \code{run_type}, so ordinary ReCV workflows need not inspect numeric
#'     identifiers. Built-in \code{loss_kind} metadata describes the stored
#'     loss and can flag unlike kinds in across-run selection.}
#'     \item{\code{settings}}{A named list containing the effective fitting configuration, including
#'     sample and search settings, cache controls, the resolved
#'     changepoint-number criterion, and model-specific loss arguments.
#'     Built-in response-based and kernel fits also record an
#'     \code{input_spec}. It validates later hold-out inputs and reconstructs
#'     pairwise matrices when the original kernel fit used raw data.}
#'     \item{\code{timing}}{A data frame containing model-fitting counts and elapsed times by
#'     \code{run_id}.}
#'     \item{\code{diagnostics}}{A list of method-specific details, such as WBS split
#'     gains or PELT pruning counts, when requested.}
#'     \item{\code{cache_profile}}{A list of reusable cache objects when requested by
#'     \code{detail}.}
#'   }
#'   Use \code{\link{select_by_run}()} for K selection within a stored loss
#'   path, \code{\link{select_across_runs}()} for a joint choice across
#'   comparable paths, \code{\link{cv.reliever}()} or
#'   \code{\link{cv.reliever_generic}()} for CPSS-style outer CV, and
#'   \code{\link{select_holdout}()} for compatible independent hold-out
#'   selection. The hold-out help page lists the supported crossfit families.
#' @references
#' Qian, C., Wang, G., and Zou, C. (2025). Reliever: Relieving the burden of
#' costly model fits for changepoint detection. \emph{Journal of Machine
#' Learning Research}, 26(203), 1--57.
#'
#' Qian, C., Wang, G., Wang, Z., and Zou, C. (2024). Changepoint detection in
#' complex models: Cross-fitting is needed. arXiv:2411.07874.
#'
#' Li, W., Wang, D., and Rinaldo, A. (2023). Divide and conquer dynamic
#' programming: An almost linear time change point detection methodology in
#' high dimensions. \emph{Proceedings of Machine Learning Research}, 202,
#' 20065--20148.
#' @seealso \code{\link{reliever_generic}()}, \code{\link{cv.reliever}()},
#'   and \code{\link{select_by_run}()}.
#' @export
#'
#' @examples
#' set.seed(2026)
#' n_seg <- 300
#' x <- rbind(
#'   matrix(rnorm(n_seg * 5, mean = 0, sd = 0.5), n_seg, 5),
#'   matrix(rnorm(n_seg * 5, mean = 4, sd = 0.5), n_seg, 5),
#'   matrix(rnorm(n_seg * 5, mean = -4, sd = 0.5), n_seg, 5)
#' )
#' fit <- reliever(
#'   X = x, cpn_max = 5, dm = 15, cov_rate = 0.6, method = "SN"
#' )
#' summary(fit) # empty: fitting did not choose K
#'
#' # CPSS-style outer sample splitting selects K by held-out loss.
#' fit_cv <- cv.reliever(
#'   X = x, cpn_max = 5, dm = 15, cov_rate = 0.6, nfolds = 3
#' )
#' selected <- summary(fit_cv)
#' selected
#' stopifnot(identical(selected$K_hat, 2L))
#'
#' \donttest{
#' # The model-based examples below take more than five seconds together.
#' # Regression-coefficient changes use the same entry point, with X and y
#' # supplied separately.
#' p <- 20
#' b0 <- c(3, -2.5, 2, -1.5, 1.5, rep(0, p - 5))
#' delta <- cbind(-2 * b0, 1.8 * b0)
#' set.seed(2026)
#' reg_data <- dgp_linear_regression(
#'   n = 450, p = p, tau = c(150, 300),
#'   b0 = b0, delta = delta, sig = 1
#' )$data
#' reg_y <- reg_data[, 1]
#' reg_x <- reg_data[, -1, drop = FALSE]
#'
#' # Cross-fitted lasso paths support both interval-adaptive ReCV and
#' # homogeneous-lambda ReCV without an additional outer-CV fit.
#' fit_lasso_cf <- reliever(
#'   X = reg_x, y = reg_y, cpd_family = "lasso_crossfit",
#'   cpn_max = 5, dm = 15, cov_rate = 0.7, method = "SN",
#'   nfolds = 2, nlambda = 25,
#'   loss_output_types = c("recv", "crossfit_homo_hyper")
#' )
#' selected_lasso_cf <- select_by_run(
#'   result = fit_lasso_cf, run_type = "recv", cpn_crit = "loss"
#' )
#' selected_lasso_cf
#' stopifnot(identical(selected_lasso_cf$K_hat, 2L))
#' selected_lasso_fixed <- select_across_runs(
#'   result = fit_lasso_cf,
#'   run_type = "crossfit_homo_hyper", cpn_crit = "loss"
#' )
#' selected_lasso_fixed
#' stopifnot(identical(selected_lasso_fixed$K_hat, 2L))
#' plot(fit_lasso_cf, run_type = "recv")
#' plot(
#'   fit_lasso_cf, x_axis = "hyperparameter", K = 2,
#'   run_type = "crossfit_homo_hyper", cpn_crit = "loss", log = "x"
#' )
#'
#' # KDE-NLL also receives the original observations. One distance matrix is
#' # reused across the bandwidth path; no bandwidth-specific Gram path is kept.
#' x_kernel <- x[seq(1, nrow(x), by = 3), , drop = FALSE]
#' fit_kde_cf <- reliever(
#'   X = x_kernel, cpd_family = "kde_nll_crossfit",
#'   kernel = "laplace",
#'   cpn_max = 3, dm = 10, cov_rate = 0.6, method = "SN", nfolds = 2
#' )
#' selected_kde_cf <- select_by_run(
#'   result = fit_kde_cf, run_type = "recv", cpn_crit = "sic"
#' )
#' selected_kde_cf
#' stopifnot(identical(selected_kde_cf$K_hat, 2L))
#' # Multivariate Student KDE is selected with, for example,
#' # kernel = "student", kernel_args = list(df = 5).
#'
#' # KDE-L2 evaluates a fixed density-kernel L2 loss on the common sample grid.
#' # One Gram matrix is constructed and reused; bandwidths are not compared.
#' fit_kernel_l2 <- reliever(
#'   X = x_kernel, cpd_family = "kde_l2",
#'   kernel = "gaussian", bandwidth = 1,
#'   cpn_max = 3, dm = 10, cov_rate = 0.6, method = "SN"
#' )
#' selected_kde_l2 <- select_by_run(
#'   result = fit_kernel_l2, cpn_crit = "rss_sic"
#' )
#' selected_kde_l2
#' stopifnot(identical(selected_kde_l2$K_hat, 2L))
#' }
reliever <- function(X, y = NULL, cpd_family = "mean",
                     cpn_max = 3, dm = 50, cov_rate = 0.8,
                     method = "SN", cpn_crit = "none",
                     pen_val = 1, prune_value = 0,
                     M = 100, wbs_seed = NULL, wbs_stop_crit = NULL,
                     detail = FALSE, cache_backend = "by_loss_block",
                     owner_key = TRUE, echo = FALSE,
                     dc_grid_size = NULL, dc_grid = NULL, ...) {
  cpd_family <- .reliever_match_cpd_family(
    cpd_family,
    c(
      "mean", "var", "meanvar", "lm", "glm", "em",
      "mean_crossfit", "lasso", "lasso_crossfit", "kde_nll",
      "kde_nll_crossfit", "kde_l2", "nmcd",
      "ranger_crossfit", "mlp_crossfit"
    ),
    "reliever"
  )
  family_args <- list(...)
  .reliever_reject_renamed_cpn_max(family_args, "reliever")

  input_spec <- if (cpd_family %in% c("lasso", "lasso_crossfit")) {
    .reliever_response_input_spec(X, y)
  } else {
    NULL
  }
  if (cpd_family %in% c("lm", "glm")) {
    prepared <- .reliever_prepare_parametric_regression(
      X, y, cpd_family, family_args
    )
    data <- prepared$data
    input_spec <- prepared$input_spec
    family_args <- prepared$family_args
  } else {
    data <- .reliever_prepare_input(
      X, y, cpd_family, response_families = c("lasso", "lasso_crossfit")
    )
  }
  if (cpd_family %in% c("glm", "em") &&
      is.null(family_args$family)) {
    stop(
      "family is required for cpd_family = \"", cpd_family, "\".",
      call. = FALSE
    )
  }

  fun <- switch(
    cpd_family,
    var = reliever_var,
    meanvar = reliever_meanvar,
    lm = reliever_lm,
    glm = reliever_glm,
    em = reliever_em,
    lasso = reliever_lasso,
    lasso_crossfit = reliever_lasso_crossfit,
    mean = reliever_mean,
    mean_crossfit = reliever_mean_crossfit,
    kde_nll = reliever_kde_nll,
    kde_nll_crossfit = reliever_kde_nll_crossfit,
    kde_l2 = reliever_kde_l2,
    nmcd = reliever_nmcd,
    ranger_crossfit = reliever_ranger_crossfit,
    mlp_crossfit = reliever_mlp_crossfit
  )

  out <- do.call(fun, c(
    list(
      data = data, cpn_max = cpn_max, dm = dm, cov_rate = cov_rate,
      method = method, cpn_crit = cpn_crit,
      pen_val = pen_val, prune_value = prune_value,
      M = M, wbs_seed = wbs_seed, wbs_stop_crit = wbs_stop_crit,
      detail = detail, cache_backend = cache_backend,
      owner_key = owner_key, echo = echo,
      dc_grid_size = dc_grid_size
    ),
    family_args,
    list(dc_grid = dc_grid)
  ))
  out$settings <- c(list(cpd_family = cpd_family), out$settings)
  if (is.null(input_spec)) {
    out
  } else {
    .reliever_set_input_spec(out, input_spec)
  }
}
