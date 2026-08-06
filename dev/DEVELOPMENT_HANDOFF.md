# relieverChangepoint Development Handoff

Last updated: 2026-07-20

This document records the current package architecture, public workflows,
design decisions that should not be accidentally reversed, and the verification
steps needed before release. It intentionally omits cumulative test counts,
artifact hashes, and superseded private-development decisions.

## 1. Current Snapshot

| Item | Current state |
| --- | --- |
| R package | `relieverChangepoint` |
| GitHub repository | `https://github.com/qianchd/reliever_changepoint` |
| Package version | `0.9.0` |
| Branch | `main` |
| Current Windows R used locally | R 4.1.3, 64 bit |

The package was renamed from `Xeed` to `relieverChangepoint`. The GitHub
repository uses the snake-case name `reliever_changepoint`. The old `Xeed`
name should not be reintroduced into public or internal interfaces.

This is still a pre-release package. We deliberately do not retain deprecated
aliases or compatibility wrappers for APIs that were only used during private
development. New code should use the current names directly.

## 2. Migration Checklist

### 2.1 Preserve Git history

Before moving development to another computer, verify that the intended branch
is committed and synchronized:

```powershell
git status --short --branch
git log -3 --oneline
git rev-parse HEAD
git push origin main
```

If pushing is not possible, create a Git bundle that includes the current
branch:

```powershell
git bundle create ..\relieverChangepoint-main.bundle main
git bundle verify ..\relieverChangepoint-main.bundle
```

On the new computer:

```powershell
git clone relieverChangepoint-main.bundle relieverChangepoint
cd relieverChangepoint
git remote set-url origin https://github.com/qianchd/reliever_changepoint.git
```

Copying the whole project directory is valid only if the hidden `.git`
directory is included.

### 2.2 Preserve ignored diagnostics separately

The following paths are intentionally ignored by Git:

- `dev/diagnostics/tmp*`
- `dev/diagnostics/results/`
- `*.Rcheck`
- `*.tar.gz`
- `*.pdf`

They contain benchmark scripts, statistical experiments, plots, and RDS
results. They are not package source and are not required to build the package,
but they are useful research records. Copy or archive `dev/diagnostics/`
separately if those records are needed on the new computer.

Some old diagnostic scripts still use the former words `sigma` or `precision`
and call removed internal helpers such as `.kde_resolve_sigma_vec()`. They are
historical experiments, not authoritative package code. Update them to
bandwidth semantics before rerunning them.

### 2.3 Recreate the R toolchain

On Windows, install a current R release and the matching Rtools. R 4.1.3 was
used for the latest local check, but a current R release should be used before
CRAN submission.

Required package dependencies:

```r
install.packages(c(
  "Rcpp", "RcppArmadillo", "glmnet", "R6",
  "testthat", "roxygen2", "pkgload", "devtools", "rcmdcheck"
))
```

Optional model dependencies:

```r
install.packages(c("nnet", "ranger"))
```

After opening the repository:

```r
Rcpp::compileAttributes()
roxygen2::roxygenise(".", roclets = c("rd", "namespace"))
devtools::test()
```

Do not edit `R/RcppExports.R` or `src/RcppExports.cpp` by hand.
The generated `man/*.Rd`, `NAMESPACE`, and Rcpp export files are currently
tracked. Commit their regenerated changes together with the source that caused
them. PDF manuals, source tarballs, and check directories remain untracked.

### 2.4 Line endings

The root `.gitattributes` makes every detected text file use LF in both the
Git index and working tree, while common image, PDF, and serialized R formats
are explicitly binary. This is the repository-level rule on macOS, Linux, and
Windows; do not rely on a developer's global `core.autocrlf` setting. The
current repository-local configuration is `core.autocrlf = false` and
`core.eol = lf`. After introducing the rule, all 133 tracked text files report
`i/lf w/lf`, and the only index content that required normalization was the
metrics CSV previously stored with CRLF.

### 2.5 GitHub Actions and secrets

The repository contains:

- `.github/workflows/R-CMD-check.yaml`;
- `.github/workflows/update-usage-metrics.yml`.

The usage workflow runs on days 1, 9, 17, and 25 and accumulates delayed GitHub
traffic counts rather than displaying only the latest 14-day window. It needs
the repository secret `GH_TRAFFIC_TOKEN`. Secrets stay with the existing
GitHub repository but are not carried by a Git clone, bundle, or copied working
tree. Recreate the secret if development moves to a new GitHub repository.

The metrics history is tracked under `dev/metrics/data/`: daily traffic in
`github_traffic_daily.csv` and collection-date summaries in
`github_usage_snapshots.csv`. These records are not written into the package
README and should not be presented as user downloads.

## 3. Package Goals and User Groups

The package is being designed for two user groups.

### 3.1 Applied users

Examples include biomedical or interdisciplinary users who mainly want one
estimated changepoint number and one changepoint set. They should normally use:

- `reliever()` for built-in loss families;
- `cv.reliever()` when outer cross-validation should select K and possibly a
  hyperparameter;
- `summary(fit)` for the compact statistical answer;
- `plot(fit)` to inspect the selection path.

Their workflow should not require understanding cache rows, C++ result
structures, or internal loss-output identifiers.

### 3.2 Statistical-method developers

Method developers may need:

- `reliever_generic()` with a custom R `reg_fun`;
- `cv.reliever_generic()` for outer-CV selection with a custom loss;
- `reg_fun_crossfit_template()` or `reg_fun_clf_crossfit_template()` for
  interval-level cross-fitting, including adaptive ReCV and fixed-hyperparameter
  outputs;
- `select_by_run()` for K selection separately within arbitrary stored loss
  paths, including outputs from a custom `reg_fun`;
- `select_across_runs()` for joint selection across arbitrary comparable
  stored loss paths;
- `run_meta`, `cpd_path`, cache profiles, and the R/C++ extension points.

The applied layer should stay simple without hiding the lower-level contract
needed by this second group.

## 4. Recommended Public Workflows

| Statistical task | Recommended first call | Recommended selection |
| --- | --- | --- |
| Mean changes, unknown K | `cv.reliever(X = X)` | Global outer CV selects K |
| Mean path without outer CV | `reliever(X = X, cpd_family = "mean")` | Apply `select_by_run(..., cpn_crit = "rss_sic")` |
| Lasso changes with interval-adaptive lambda | `reliever(X = X, y = y, cpd_family = "lasso_crossfit")` | Apply `select_by_run(..., run_type = "recv", cpn_crit = "loss")`; lambda may differ between intervals |
| Lasso with one globally selected lambda | Apply `select_across_runs()` to the lasso crossfit fit's `crossfit_homo_hyper` runs | Homogeneous-lambda cross-fitted loss jointly selects lambda and K without outer CV or refitting |
| Lasso path without CV | `reliever(X = X, y = y, cpd_family = "lasso")` | `select_by_run(..., cpn_crit = "rss_sic")` selects K separately within each lambda; it does not select lambda across runs |
| KDE NLL distribution changes | `reliever(X = X, cpd_family = "kde_nll_crossfit")` | Apply `select_by_run(..., run_type = "recv", cpn_crit = "loss")`; `sic` is a more conservative K criterion |
| Fixed-kernel KDE L2 distribution changes | `reliever(X = X, cpd_family = "kde_l2", kernel = ..., bandwidth = ...)` | Fix the density kernel and bandwidth first, then apply `select_by_run(..., cpn_crit = "rss_sic")` |
| Univariate nonparametric changes | `reliever(X = x, cpd_family = "nmcd")` | Apply `select_by_run(..., cpn_crit = "sic")` |
| Classifier-based distribution changes | `reliever(X = X, cpd_family = "ranger_crossfit")` or `"mlp_crossfit"` | Apply `select_by_run(..., run_type = "recv", cpn_crit = "loss")` |
| Custom interval loss | `reliever_generic(data = data, reg_fun = reg_fun)` | No automatic rule unless the caller supplies an appropriate `cpn_crit` |

Important interpretation rules:

- A `reliever*()` function fits candidate paths. The public fitting defaults
  use `cpn_crit = "none"` so K and hyperparameter selection remain explicit
  through outer CV or the post-fit selectors.
- `cv.reliever()` or `cv.reliever_generic()` reruns the fit in outer folds
  and explicitly selects among comparable model-setting/path combinations by
  held-out loss.
- Full-data lasso training loss must not be minimized across lambda. The
  interval-adaptive choice is
  `select_by_run(..., run_type = "recv", cpn_crit = "loss")`. When one lambda
  must be shared across intervals, pass the already produced
  `crossfit_homo_hyper` runs to `select_across_runs()` for joint lambda/K selection.
  Outer CV and an independent hold-out set are optional alternatives.
- An RSS information criterion can select K separately inside each
  fixed-lambda run, but it must not compare training losses across lambda.
- KDE NLL can compare bandwidths through cross-fitted NLL.
- KDE-L2 raw losses from different RBF bandwidths are not comparable. As
  bandwidth tends to infinity, the feature map becomes constant and the L2
  loss tends to zero even without a meaningful fit. `reliever_kde_l2()` thus
  intentionally accepts one fixed feature/kernel matrix and has no automatic
  bandwidth path.

## 5. Main R API

### 5.1 Dispatchers

`reliever()` is the simple built-in dispatcher:

```r
reliever(
  X = X,
  y = NULL,
  cpd_family = "mean",
  cpn_max = 3,
  dm = 50,
  cov_rate = 0.8,
  method = "SN",
  cpn_crit = "none",
  ...
)
```

Supported `cpd_family` values:

- `mean`
- `mean_crossfit`
- `lasso`
- `lasso_crossfit`
- `kde_nll`
- `kde_nll_crossfit`
- `kde_l2`
- `nmcd`
- `ranger_crossfit`
- `mlp_crossfit`

`reliever_generic()` is the lower-level custom-loss interface. Its default
`cpn_crit = "none"` is intentional because it cannot infer whether an arbitrary
loss is RSS, negative log-likelihood, or cross-fitted evaluation loss.

`cv.reliever()` currently supports built-in `mean` and `lasso` families. Its
mean branch uses the faster native C++ mean-loss path.
`cv.reliever_generic()` supports a custom `reg_fun`.

### 5.2 Focused wrappers

Focused wrappers remain public because they document model-specific arguments
more clearly than `...` in the dispatcher:

- `reliever_mean()`
- `reliever_mean_crossfit()`
- `reliever_lasso()`
- `reliever_lasso_crossfit()`
- `reliever_kde_nll()`
- `reliever_kde_nll_crossfit()`
- `reliever_kde_l2()`
- `reliever_nmcd()`
- `reliever_ranger_crossfit()`
- `reliever_mlp_crossfit()`

Every focused wrapper defaults to `cpn_crit = "none"`, matching `reliever()`
and keeping fitting separate from selection. Model-specific criteria remain
explicit: RSS-SIC for mean/lasso/KDE-L2 paths, additive SIC for
likelihood/CDF paths, and ReCV through `run_type = "recv"` with
`cpn_crit = "loss"`. Do not add a hidden wrapper-specific default resolver.

### 5.3 Search methods

Supported methods:

- `SN`
- `WBS`
- `WBS_recursive`
- `SeedBS`
- `BS`
- `PELT`
- `OP`

`WBS` uses the package's global greedy path: after each split, compare the best
gain over all current child segments and split the largest one.

`WBS_recursive` follows the original WBS-style recursive tree traversal.

For WBS-family methods, `wbs_stop_crit` maps each threshold to the prefix before
the first proposed split whose gain is no larger than the threshold. The
current implementation remains capped at `cpn_max` splits. Without
`wbs_stop_crit`, the path is indexed directly by K from 0 to the effective
`cpn_max`.

For `PELT` and `OP`, `pen_val` is the search selector. Each penalty returns one
segmentation. The result keeps the penalty-to-candidate mapping even when
several penalties produce the same changepoint set.

`twostep()` is parallel to Reliever, not a post-processing function. It uses a
coarse initial split plus local split refinement within WBS-family searches.
`local_refine()` is the separate post-estimation refinement function.

## 6. Result Objects

### 6.1 `reliever_result`

A normal fit has these user-facing components:

```text
summary
cpd_path
run_meta
settings
timing
diagnostics
cache_profile
```

Their roles are:

- `summary`: compact statistical result returned by `summary(fit)` and shown by
  `print(fit)`;
- `cpd_path$candidates`: all distinct candidate segmentations, with `run_id`,
  K, changepoints, loss, and an internal `candidate_id`;
- `cpd_path$selector_map`: optional mapping from `wbs_stop_crit` or `pen_val`
  values to candidate rows;
- `run_meta`: one row per searched loss output, including model metadata such
  as lambda, bandwidth, row type, and hyperparameter identifier when supplied
  by `reg_fun`;
- `settings`: effective arguments after validation and automatic resolution,
  including the resolved lambda/bandwidth path in `loss_args`;
- `timing`: model-fit counts and elapsed times by run;
- `diagnostics`: debug/search details only when requested;
- `cache_profile`: reusable cache data only when requested.

The compact result intentionally does not duplicate `K_hat` and `cpd_hat` as
top-level scalar fields. Use `summary(fit)` rather than manually inspecting the
internal path for ordinary use.

### 6.2 `cv_reliever_result`

An outer-CV fit contains:

- `summary`: one selected full-data result;
- `cv_loss`: every candidate compared by held-out loss, including `cv_mean` and
  `cv_se`;
- `full_data_fit`: the complete final Reliever path fitted to all observations;
- `settings`;
- `diagnostics` when requested.

`summary(fit)` is the final answer. `fit$cv_loss` is retained for plotting and
alternative selection rules.

### 6.3 Printing, summaries, and plotting

Implemented S3 methods:

- `print.reliever_result()`
- `summary.reliever_result()`
- `plot.reliever_result()`
- `print.cv_reliever_result()`
- `summary.cv_reliever_result()`
- `plot.cv_reliever_result()`
- `print.reliever_model_selection()`

`plot(fit)` shows the K-selection path. Other supported views:

```r
plot(fit, x_axis = "hyperparameter", K = 2)
plot(fit, x_axis = "search_value")
```

For comparable fixed-hyperparameter cross-fitted or outer-CV losses, an
endpoint minimum triggers a warning that the lambda/bandwidth range may be too
narrow. Ordinary in-sample solution-path loss is not treated as a valid
hyperparameter selector.

### 6.4 Post-fit selection

- A run is one fitted loss path, including one output of a custom `reg_fun`.
- `select_by_run()` selects K independently inside each requested run and
  does not compare runs or hyperparameters.
- `select_across_runs()` pools candidates across runs and selects one
  run/model setting together with K; callers choose the intended paths with
  either exact `run_ids` or a statistical `run_type`.
- `evaluate_reliever_segments()` computes hold-out segment losses for an
  existing candidate path.
- `select_holdout()` evaluates and selects inside one call.

These public functions are the general selection interface for built-in and
custom loss outputs. They operate on an already fitted candidate path and
never rerun the changepoint search. For lasso crossfit,
`select_by_run(..., run_type = "recv", cpn_crit = "loss")` performs
interval-adaptive tuning; `select_across_runs()` on `crossfit_homo_hyper` performs
optional homogeneous tuning.

`evaluate_reliever_segments()` and `select_holdout()` require independent
`eval_data`. `eval_index` may be omitted only when evaluation rows map
one-to-one to the original observations; sparse, repeated, or reordered
evaluation samples require it explicitly. They must never silently fall back
to in-sample loss.

## 7. Loss and Model-Selection Semantics

Supported K-selection criteria:

| Criterion | Score |
| --- | --- |
| `loss` | stored loss |
| `aic` | `loss + 2 K` |
| `hqc` | `loss + 2 log(log(n)) K` |
| `sic` | `loss + log(n) K` |
| `rss_aic` | `n / 2 * log(RSS / n) + 2 K` |
| `rss_hqc` | `n / 2 * log(RSS / n) + 2 log(log(n)) K` |
| `rss_sic` | `n / 2 * log(RSS / n) + log(n) K` |
| non-negative numeric value | `loss + value * K` |
| `none` | no loss-based K selection |

The `rss_*` rules apply a log-RSS transformation and therefore require a
non-negative loss. The additive rules use the stored loss directly. The
software documents these interpretations but does not infer or restrict the
criterion from `loss_kind`.

## 8. Custom `reg_fun` Contract

An R interval-loss function has the conceptual signature:

```r
reg_fun <- function(
  data,
  l,
  r,
  l_end = l,
  r_end = r,
  save_model = FALSE,
  is_virtual_run = FALSE,
  ...
) {
  # ...
}
```

The fitted core interval is `l:r`; losses are returned for `l_end:r_end`.

For `is_virtual_run = TRUE`, return either:

```r
n_loss_outputs
```

or:

```r
list(
  n_loss_outputs = ...,
  loss_output_meta = data.frame(...)
)
```

`loss_output_meta` should have one row per loss column. Useful optional fields
include:

- `loss_output_id`
- `row_type`
- `hyper_id`
- `hyper_value`
- `hyper_name`
- `default_selection`

For a real fit, return:

```r
list(
  loss = loss_matrix,
  model = optional_model
)
```

The loss matrix must have one row per evaluated observation and one column per
declared loss output. Hyperparameters belong in metadata, not inside the loss
matrix.

`twostep()` requires observation-level `individual_loss`; Reliever can use the
same interface and aggregate it into block losses. Native C++ losses may
override block-loss computation for a faster path.

## 9. Crossfit Outputs and Outer CV

The `*_crossfit` functions are containers for several interval-level loss
types. ReCV is one output within that container. Global outer CV is a distinct
procedure and should not be renamed into either one.

### 9.1 Interval-level cross-fitting

Implemented by `reg_fun_crossfit_template()` and the `reliever_*_crossfit()` wrappers.
Folds are reconstructed inside each fitted interval. The output can include:

- interval-adaptive ReCV (`recv`);
- corresponding in-sample CV (`incv`);
- one fixed-hyperparameter cross-fitted row per hyperparameter;
- one matching fixed-hyperparameter in-sample row per hyperparameter.

The adaptive `recv` row chooses a hyperparameter separately in every interval.
It therefore does not represent one global lambda or bandwidth.

For lasso, `select_by_run(..., run_type = "recv", cpn_crit = "loss")` uses the
interval-adaptive row and therefore does not estimate one global lambda. When
homogeneous tuning is required, use `run_type = "crossfit_homo_hyper"` and
`cpn_crit = "loss"` in `select_across_runs()`. This jointly selects one lambda
and K from the stored cross-fitted losses; no second fit or outer CV is
required.

### 9.2 Global outer CV

Implemented by `cv.reliever()` for built-in families and
`cv.reliever_generic()` for custom `reg_fun` losses.

Folds are constructed once on the complete time axis. Training rows retain
their original order. A boundary after training row `t` maps back to original
time `train_id[t]`; held-out rows are assigned to the resulting original-time
segments.

Outer CV can select:

- K for SN/K-indexed paths;
- `wbs_stop_crit` for WBS-family threshold paths;
- `pen_val` for PELT/OP;
- a model hyperparameter such as lasso lambda when `reg_fun` has several loss
  outputs.

## 10. Relief Intervals and Cache Backends

### 10.1 Interval mapping

The current names are:

- `create_relief_itv()`
- `exact2relief_itv_routine_c()`
- `relief2exact_itv_routine_c()`
- `create_wbs_itv()`
- `create_seed_itv()`

The two mapping routines are intended to be algebraic inverses for the owned
sawtooth region. Their exhaustive tests are expensive and are isolated in the
interval test mode.

### 10.2 Cache backends

`cache_backend = "by_loss_block"` is the default for non-full Reliever search.
It stores fitted raw relief blocks and reconstructs exact interval costs.

`cache_backend = "by_cost_mat"` stores expanded exact interval costs in a
dense matrix. It remains useful for full search and for explicitly supplied
cost-matrix cache reuse.

When `cov_rate` is effectively 1, a requested loss-block backend warns and
automatically switches to `by_cost_mat`.

Both backends support `dc_grid`. `dm` is specified on the original time scale
and is converted to a grid-scale minimum distance internally.

`owner_key` applies only to `by_loss_block`. It trades an integer owner lookup
array for faster repeated exact-to-relief lookup. It must not affect statistical
results.

`detail = FALSE` does not serialize a cache profile. `detail = TRUE` or
`"cache"` returns reusable cache objects. Expensive diagnostic counters belong
to that explicit cache profile, not the default result. The current accepted
detail modes are `FALSE`/`"none"` and `TRUE`/`"cache"`; there is no separate
public debug mode.

Cache reuse is validated against data/loss context:

- `by_loss_block` reuse requires the same relief interval set;
- `by_cost_mat` reuses the double cost matrix;
- full search ignores an incompatible loss-block cache and creates a
  cost-matrix cache.

The R-visible loss-block cache is an inspectable R list rather than an external
pointer. This was a deliberate choice. Do not replace it with an opaque C++
pointer merely to avoid serialization.

## 11. C++ Architecture

The intended layers are:

```text
RegLossFunction
  - RRegLossFunction (R callback, bridge layer)
  - MeanSquareRegLossFunction (native fast mean loss)
  - future native model losses

CostEngine
  - CostEngineByCostMat
  - CostEngineByLossBlock

Search algorithms
  - SN
  - global-greedy WBS / SeedBS / BS
  - recursive WBS
  - OP / PELT
  - TwoStepSearch

C++ result model
  - CpdCandidates
  - SingleCpdResult
  - AllCpdResults

R bridge
  - flatten C++ results
  - convert cache state
  - call R callbacks
  - format final R objects in R
```

Important files:

| File | Responsibility |
| --- | --- |
| `src/cost_engine.h/.cpp` | Loss interface, native mean loss, both cache engines |
| `src/cpd_algorithms.h/.cpp` | Changepoint search algorithms and split evaluators |
| `src/cpd_result.h/.cpp` | R-independent result structures |
| `src/twostep.h/.cpp` | Pure C++ two-step search core |
| `src/relief_interval.h/.cpp` | Exact/relief interval mapping |
| `src/reg_fun.h/.cpp` | Native NMCD loss primitive |
| `src/*_r_bridge.*` | Rcpp-only conversion and callback boundaries |

The core classes expose `arma::*`, `std::*`, and primitive C++ types rather
than `Rcpp::List` or `Rcpp::DataFrame`. Rcpp data structures should remain in
the bridge files. Some core headers still include `RcppArmadillo.h` because
this is an R package build; a future standalone library could replace those
includes with Armadillo-only headers without redesigning the class APIs.

The current public `twostep()` path still reaches the C++ core through
`RRegLossFunction`. The C++ core itself accepts any `RegLossFunction`, so a
future native mean two-step wrapper can reuse it without rewriting the search.

## 12. Major Work Already Completed

### 12.1 API and naming

- Renamed the package and internal `xeed_*` names to `reliever*`.
- Standardized interval names around `itv`.
- Standardized cache backend names to `by_cost_mat` and `by_loss_block`.
- Standardized KDE public APIs and metadata on `bandwidth`; removed public
  precision semantics.
- Renamed interval CV methods to ReCV to distinguish them from outer CV.
- Removed obsolete API compatibility aliases because no public release exists.
- Removed `data_te` and post-fit model evaluation from the core Reliever call.
- Removed the old `save_hyper` workflow.
- Moved hold-out evaluation and selection to explicit post-fit functions.

### 12.2 Result and post-fit selection API

- Replaced dense/duplicated changepoint output tables with a compact candidate
  path plus selector map.
- Added `run_meta` as the canonical loss-output lookup table.
- Added compact `summary`, `print`, and `plot` methods.
- Added separate within-run and across-run post-fit selection helpers.
- Added outer-CV result objects with one final selected full-data result.
- Preserved complete paths for statistical-method developers without making
  applied users read them.

### 12.3 Search methods

- Unified SN, WBS, SeedBS, BS, OP, and PELT around common C++ result structures.
- Added `WBS_recursive` without changing the existing global-greedy WBS.
- Added `wbs_stop_crit` path mapping.
- Added the same WBS-family semantics to `twostep()`.
- Added `dc_grid` tests and corrected grid-scale `dm`.
- Added protection for too-short or non-finite split intervals.

### 12.4 Cache and performance work

- Implemented symmetric `CostEngineByCostMat` and
  `CostEngineByLossBlock` interfaces.
- Implemented R-visible cache reuse for both backends.
- Added optional lazy `owner_key` acceleration for the loss-block backend.
- Kept full search on the cost-matrix backend.
- Reduced default result/cache diagnostics and separated them from normal
  cache operation.
- Investigated explicit periodic R garbage collection, but did not retain a
  public `gc_every` option because frequent collection substantially slowed
  model fitting and did not explain the stable final memory results.
- Benchmarked cost-matrix/loss-block timing, peak RSS, cache size, and result
  parity across several settings.

### 12.5 R/C++ boundary

- Removed Rcpp containers from C++ result and cache core APIs.
- Moved R list/data-frame conversion to explicit bridge files.
- Replaced nested lambda-heavy runners with named search functions and simpler
  result ownership.
- Unified native mean and R callback paths through `RegLossFunction`.
- Added a C++ `TwoStepSearch` core that is not structurally tied to R lists.
- Removed unused one-sided/null/legacy C++ paths and stale runtime state.

### 12.6 Built-in models

- Native fast multidimensional mean loss with live interval updates.
- Lasso solution path with automatic glmnet-based lambda grid.
- Cross-fitted lasso outputs, including fixed-lambda CF and adaptive ReCV.
- KDE NLL solution path and crossfit outputs with a scale-adaptive bandwidth
  grid.
- Fixed-kernel KDE-L2 loss evaluated from one kernel-density contribution
  matrix; generic fixed-feature inputs remain a computational extension.
- Native NMCD loss primitive and public wrapper.
- Ranger and MLP classifier crossfit wrappers.
- Generic crossfit and classifier-crossfit templates for extension.

### 12.7 Documentation and packaging

- Rewrote README for applied users first and developers second.
- Added recommended K-selection and model-selection workflows.
- Added statistical examples with two changepoints, generally `n = 900`,
  `set.seed(2026)`, and `cov_rate = 0.8`.
- Added `See Also` links and clarified ambiguous internal terminology.
- Added package-level references to the Reliever and cross-fitting papers.
- Added GitHub usage-metric automation.
- Added CRAN-oriented DESCRIPTION, build, manual, and check workflows.

## 13. Testing Structure and Current Evidence

Tests are organized in `tests/testthat/README.md`.

### 13.1 Test modes

Quick tests:

```r
devtools::test()
```

Full model/API/backend tests:

```r
Sys.setenv(RELIEVER_TEST_MODE = "full")
devtools::test()
```

Exhaustive relief ownership tests:

```r
Sys.setenv(RELIEVER_TEST_MODE = "interval")
devtools::test(filter = "interval-routines")
```

Everything:

```r
Sys.setenv(RELIEVER_TEST_MODE = "all")
devtools::test()
```

Individual tests use `with_test_timeout()` and should complete within two
minutes. If a test exceeds that limit, reduce its computational settings
without weakening the statistical assertion.

### 13.2 What is tested

- API and wrapper contracts;
- all search methods;
- both cache backends and cache reuse;
- native mean versus generic R loss;
- `dc_grid`;
- WBS and recursive WBS;
- TwoStep and local refinement;
- outer CV for K, WBS thresholds, and PELT/OP penalties;
- ReCV metadata and selection;
- two-changepoint statistical accuracy;
- lasso and KDE bandwidth paths;
- optional ranger/MLP wrappers;
- interval ownership algebra.

Tests should assert statistical content such as K, changepoint locations, loss
ordering, or backend equality. Merely checking that code runs without an error
is not sufficient for core methods. For a stochastic DGP, assert K separately
and bound changepoint-location error with `cp_error()`; do not require the
estimated locations to equal the truth exactly. A tolerance of 10 observations
is the standard acceptance threshold for the documented `n = 900` workflows.
Exact equality remains appropriate for deterministic contracts such as backend
parity, cache reuse, index conversion, and artificial selector tables.

### 13.3 Required verification evidence

Do not preserve cumulative expectation counts, tarball sizes, or artifact
hashes here; they become stale as soon as the next change lands. For each
release candidate, record the commit hash in the release notes or CI run and
verify:

1. targeted tests for every changed subsystem;
2. `RELIEVER_TEST_MODE = "all"` with no unexpected failure, warning, or skip;
3. a clean source build;
4. `R CMD check --as-cran` on the source tarball;
5. generated Rd files and `NAMESPACE` match their roxygen sources.

Generated tarballs, check directories, manuals, and diagnostic results remain
untracked. While the GitHub repository is private, an incoming URL check may
report that the project URLs are unavailable; make the repository public and
rerun that check before CRAN submission.

## 14. Current Decisions and Remaining Work

### 14.1 Public API state

The package is still unpublished, so clean API changes are preferred over
private-development compatibility layers. Current fitting signatures put the
primary data, model family, common search controls, and selection controls
before specialized or diagnostic arguments. Keep this ordering unless a clear
user workflow justifies changing it.

Tests and examples should name every argument after the primary data or result
object. In particular, the current `y` and `eval_y` positions in
`evaluate_reliever_segments()` and `select_holdout()` are intentional; call
sites must not rely on the positions of later arguments.

`select_by_run()` returns one row for every requested run. If a run has no
finite candidate score, its row is retained with `NA` selection fields.
`select_across_runs()` may warn when metadata suggest that losses could be
incomparable, but it does not forbid deliberate advanced comparisons.
`loss_kind` is diagnostic metadata, not a gate on the user's `cpn_crit` choice.

### 14.2 Statistical implementation state

- Native mean costs use the optimized C++ path and stabilize numerically risky
  segment calculations without changing the ordinary fast path.
- `prune_value = 0` is the ordinary PELT pruning constant for the supported
  squared-error costs; choosing PELT does not silently switch to OP.
- CV and outer-CV folds may differ in size when `n` is not divisible by
  `nfolds`; every observation must still appear in exactly one evaluation fold.
- KDE-NLL reuses a bandwidth-independent pairwise distance representation.
  KDE-L2 uses one fixed kernel-evaluation matrix at its supplied bandwidth.
- Categorical hyperparameters are supported in profile plots and follow their
  declared run order.

### 14.3 Remaining release work

1. Run targeted tests, all-mode tests, and a source-tarball CRAN-style check for
   the final commit.
2. Push the verified commit and inspect the GitHub Actions matrix.
3. Make the repository public before CRAN submission, then rerun incoming URL
   checks.

## 15. Development Rules to Preserve

- Use current disk files as authoritative; do not review only an old commit.
- Do not revert unrelated user changes in a dirty worktree.
- Keep public interfaces explicit and put common user controls before advanced
  or diagnostic controls.
- Use named arguments after the primary data or result object in tests,
  examples, and internal calls.
- Use `select_by_run()`, `select_across_runs()`, and `select_holdout()`; do not
  reintroduce the former `reselect_*()` aliases.
- Preserve statistical content in tests. A test that merely runs is not enough.
- Keep every individual test under two minutes.
- Keep generated Rd and `NAMESPACE` changes with their roxygen sources; do not
  edit generated Rcpp export files manually.
- Use `bandwidth`, not kernel precision, in public kernel/KDE interfaces.
- Do not store hyperparameters inside cost matrices.
- Do not make normal `detail = FALSE` runs serialize debug cache diagnostics.
- Keep the inspectable R cache rather than replacing it with an opaque external
  pointer.

## 16. Verification Workflow

```r
Rcpp::compileAttributes()
roxygen2::roxygenise(".", roclets = c("rd", "namespace"))
devtools::test()

Sys.setenv(RELIEVER_TEST_MODE = "all")
devtools::test()
```

Then build and check the source tarball:

```powershell
R CMD build .
R CMD check --as-cran relieverChangepoint_0.9.0.tar.gz
```

Before transferring or releasing the repository, confirm the branch and remote
state with `git status --short --branch`, `git log -3 --oneline`, and
`git rev-parse HEAD`.
