#' relieverChangepoint: Efficient Changepoint Detection and Model Selection
#'
#' `relieverChangepoint` provides Reliever-style changepoint detection for
#' settings where fitting the interval model is the expensive step. Built-in
#' functions cover mean, regression, kernel, nonparametric, and classifier-based
#' changes. Information criteria can select the number of changepoints.
#' Hold-out data, recycled cross-validation, and outer cross-validation can
#' also select model hyperparameters when comparable candidates are available.
#'
#' @section Primary fitting interfaces:
#' The two primary fitting functions are:
#' \itemize{
#'   \item `reliever()` for built-in losses. It fits mean changes by default;
#'     choose another method with `cpd_family`, for example
#'     `reliever(X, y, cpd_family = "lasso_crossfit")` for regression changes
#'     or `reliever(X, cpd_family = "kde_nll_crossfit")` for
#'     distribution changes.
#'   \item `reliever_generic(data, reg_fun = ...)` for a user-supplied
#'     interval-loss function.
#' }
#' Focused `reliever_*()` functions expose the complete arguments for the same
#' built-in methods.
#'
#' @section Fitting and selection:
#' A `reliever()` or `reliever_generic()` call fits candidate changepoint paths
#' with `cpn_crit = "none"` by default. Select a path with one of the following
#' workflows.
#'
#' For ordinary mean and lasso losses, `cv.reliever()` implements the
#' sample-efficient cross-fitted version of CPSS. CPSS is the sample-splitting
#' selector introduced by Zou, Wang, and Li (2020):
#' \preformatted{
#' mean_cv <- cv.reliever(
#'   X = x, cpd_family = "mean",
#'   cpn_max = 5, dm = 30, cov_rate = 0.8
#' )
#' plot(mean_cv)
#' summary(mean_cv)
#'
#' lasso_cv <- cv.reliever(
#'   X = X, y = y, cpd_family = "lasso",
#'   cpn_max = 5, dm = 30, cov_rate = 0.8
#' )
#' plot(lasso_cv)
#' summary(lasso_cv)
#' }
#' Each `cv.reliever()` call automatically performs the final full-data
#' `reliever()` fit after outer-CV selection and returns it in
#' `full_data_fit`. Do not call `reliever()` again merely to refit the selected
#' result.
#'
#' Use `cv.reliever_generic()` for a compatible custom loss. Fixed-kernel
#' KDE-L2 and NMCD paths support explicit information criteria or independent
#' hold-out selection. KDE bandwidth selection is available through
#' `cpd_family = "kde_nll_crossfit"`.
#'
#' A cross-fitted fit returns only the interval-adaptive \code{recv} loss by
#' default. Request homogeneous-hyperparameter paths when they are needed:
#' \preformatted{
#' fit <- reliever(
#'   X = X, y = y, cpd_family = "lasso_crossfit",
#'   cpn_crit = "none",
#'   loss_output_types = c("recv", "crossfit_homo_hyper")
#' )
#' plot(fit, run_type = "recv")
#'
#' recv <- select_by_run(
#'   result = fit, run_type = "recv", cpn_crit = "loss"
#' )
#'
#' recv_homo_hyper <- select_across_runs(
#'   result = fit, run_type = "crossfit_homo_hyper", cpn_crit = "loss"
#' )
#' }
#' `run_type = "recv"` lets the hyperparameter adapt by interval and selects K
#' and changepoints. `run_type = "crossfit_homo_hyper"` compares cross-fitted paths
#' with one global hyperparameter and jointly selects that setting, K, and
#' changepoints. For lasso this is the normalized `lam_set`. Both selectors
#' reuse `fit`. Add `"incv"` to `loss_output_types` only when the matching
#' interval-adaptive in-sample loss is also required.
#'
#' Information criteria are explicit alternatives applied with
#' `select_by_run()`; see `reliever()` for the available criteria.
#'
#' With `cpn_crit = "none"`, `summary(fit)` is usually empty and
#' `plot(fit, run_type = "recv")` shows the interval-adaptive fitted path.
#' A fixed-K hyperparameter curve uses:
#' \preformatted{
#' plot(
#'   x = fit,
#'   x_axis = "hyperparameter",
#'   K = recv$K_hat[[1L]],
#'   run_type = "crossfit_homo_hyper"
#' )
#' }
#' Display the selected changepoints with
#' `plot_reliever_data(result = fit, data = y, selection = recv)`.
#'
#' @section Post-fit selectors:
#' \describe{
#'   \item{\code{\link{select_by_run}()}}{Select K separately in each requested
#'   loss path.}
#'   \item{\code{\link{select_across_runs}()}}{Select one comparable loss path
#'   and K jointly.}
#'   \item{\code{\link{select_holdout}()}}{Select a fitted path using
#'   independent evaluation data.}
#' }
#'
#' @section Extending Reliever:
#' `reliever_generic()` accepts a custom interval-loss function (`reg_fun`).
#' Its help page gives the loss-output contract and a minimal implementation.
#' Compatible functions can also be used with `cv.reliever_generic()`.
#'
#' Search algorithms include segment neighbourhood, optimal partitioning,
#' PELT, globally greedy and recursive wild binary segmentation, seeded binary
#' segmentation, and binary segmentation. `twostep()` and `local_refine()`
#' provide alternative low-model-fit and refinement workflows.
#'
#' @references
#' Qian, C., Wang, G., and Zou, C. (2025). Reliever: Relieving the burden of
#' costly model fits for changepoint detection. \emph{Journal of Machine
#' Learning Research}, 26(203), 1--57.
#'
#' Qian, C., Wang, G., Wang, Z., and Zou, C. (2024). Changepoint detection in
#' complex models: Cross-fitting is needed. arXiv:2411.07874.
#'
#' Zou, C., Wang, G., and Li, R. (2020). Consistent selection of the number of
#' change-points via sample-splitting. \emph{The Annals of Statistics}, 48(1),
#' 413--439.
#'
#' @seealso \code{\link{reliever}()}, \code{\link{cv.reliever}()}, and
#'   \code{\link{reliever_generic}()}.
#' @keywords package
#' @importFrom Rcpp evalCpp
#' @importFrom stats predict
#' @importFrom stats toeplitz
#' @importFrom utils tail
#' @useDynLib relieverChangepoint, .registration = TRUE
#'
#' @examples
#' set.seed(2026)
#' n_seg <- 300
#' x <- c(
#'   rnorm(n_seg, mean = 0, sd = 0.5),
#'   rnorm(n_seg, mean = 4, sd = 0.5),
#'   rnorm(n_seg, mean = -4, sd = 0.5)
#' )
#' fit_cv <- cv.reliever(
#'   X = x, cpn_max = 5, dm = 30, cov_rate = 0.8,
#'   nfolds = 3, method = "SN"
#' )
#' plot(fit_cv)
#' selected <- summary(fit_cv)
#' selected
#' stopifnot(identical(selected$K_hat, 2L))
"_PACKAGE"
