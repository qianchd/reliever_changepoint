# relieverChangepoint

`relieverChangepoint` implements efficient changepoint detection and
data-driven selection of tuning parameters and the number of changepoints. Its
methodology is mainly based on:

- Qian, Wang, and Zou (2025). Reliever: Relieving the burden of costly model
  fits for changepoint detection. *Journal of Machine Learning Research*,
  26(203), 1--57. <https://www.jmlr.org/papers/v26/24-1108.html>
- Qian, Wang, Wang, and Zou (2024). Changepoint detection in complex models:
  Cross-fitting is needed. <https://arxiv.org/abs/2411.07874>
- Zou, Wang, and Li (2020). Consistent selection of the number of change-points
  via sample-splitting. *The Annals of Statistics*, 48(1), 413--439.

The efficient algorithm, *reliever*, is designed for settings where repeatedly
fitting an interval model is expensive. Built-in changepoint families in the
package cover mean, covariance, regression, kernel, nonparametric,
exponential-family, and classifier-based changes. Search algorithms include
SN, OP, PELT, WBS, SeedBS, BS, and two-step/local-refinement workflows.

The package separates two tasks deliberately:

1. `reliever()` and `reliever_generic()` fit candidate changepoint models. By
   default they retain the complete path and do not select a model or K.
2. Tuning parameter and K selections are then performed by outer CV,
   cross-fitting ReCV, an information criterion, or independent hold-out data.

Documentation is available as an
[online reference and set of guides](https://qianchd.github.io/reliever_changepoint/)
or as a complete [PDF reference manual](docs/relieverChangepoint-manual.pdf).

## Installation

Install the development version after the GitHub repository becomes public:

```r
remotes::install_github("qianchd/reliever_changepoint")
```

After a CRAN release, install it with:

```r
install.packages("relieverChangepoint")
```

## Four main entry points

If you are new to the package, start with `cv.reliever()`: for ordinary
built-in losses it selects K and returns the final full-data fit in one call.

| Entry point | Purpose | Returned result |
| --- | --- | --- |
| `reliever()` | Fit a built-in changepoint loss family | Candidate paths; no model or K selection by default |
| `cv.reliever()` | Select K and, for lasso, its normalized lambda setting by CPSS outer CV | Selected result plus an automatic final full-data `reliever()` fit |
| `reliever_generic()` | Fit a user-defined interval loss | Candidate paths for all declared loss outputs |
| `cv.reliever_generic()` | Apply outer CV to a compatible custom loss | Selected result plus an automatic final full-data `reliever_generic()` fit |

`cv.reliever()` and `cv.reliever_generic()` already perform the final fit on
all observations. Their `full_data_fit` component contains that fit, so an
additional external `reliever()` call is not needed.

As a quick workflow guide:

| Goal | Recommended workflow |
| --- | --- |
| Select K for an ordinary built-in loss | `cv.reliever()` |
| Fit and retain a candidate path without selecting K | `reliever()` |
| Select K from an interval-level cross-fitted loss | `reliever()` followed by `select_by_run(..., run_type = "recv")` |
| Supply a custom interval loss | `reliever_generic()` or `cv.reliever_generic()` |

## Quick start: mean changes

For a built-in single-path parametric loss, `cv.reliever()` is the recommended
way to select K. It implements the sample-efficient cross-fitting version of
the CPSS sample-splitting selector.

```r
library(relieverChangepoint)

set.seed(2026)
n_seg <- 300
x <- c(
  rnorm(n_seg, mean = 0, sd = 0.5),
  rnorm(n_seg, mean = 4, sd = 0.5),
  rnorm(n_seg, mean = -4, sd = 0.5)
)

fit_mean_cv <- cv.reliever(
  X = x,
  cpn_max = 5,
  dm = 30,
  cov_rate = 0.8,
  method = "SN",
  nfolds = 5
)

summary(fit_mean_cv)
plot(fit_mean_cv)
fit_mean_cv$full_data_fit
plot_reliever_data(result = fit_mean_cv, data = x)
```

To fit the candidate path without selecting K:

```r
fit_mean_path <- reliever(
  X = x,
  cpn_max = 5,
  dm = 30,
  cov_rate = 0.8,
  method = "SN"
)

summary(fit_mean_path)  # explains the available selection workflows
plot(fit_mean_path)

# Explicit information-criterion alternative for an RSS path.
mean_sic <- select_by_run(
  result = fit_mean_path,
  cpn_crit = "rss_sic"
)
mean_sic
```

## Cross-fitting model selection

For a cross-fitting changepoint loss, the interval-adaptive `recv` output is the
recommended primary path:

```r
set.seed(2026)
n <- 450
p <- 20
tau <- c(150, 300)
b0 <- c(3, -2.5, 2, -1.5, 1.5, rep(0, p - 5))
delta <- cbind(-2 * b0, 1.8 * b0)
reg_data <- dgp_linear_regression(n, p, tau, b0, delta, sig = 1)$data
reg_y <- reg_data[, 1]
reg_x <- reg_data[, -1, drop = FALSE]

fit_lasso_cf <- reliever(
  X = reg_x,
  y = reg_y,
  cpd_family = "lasso_crossfit",
  cpn_max = 5,
  dm = 15,
  cov_rate = 0.7,
  method = "SN",
  nfolds = 2,
  nlambda = 25,
  loss_output_types = c("recv", "crossfit_homo_hyper")
)

# Omission also defaults to recv, but documentation keeps it explicit.
plot(fit_lasso_cf, run_type = "recv")

lasso_recv <- select_by_run(
  result = fit_lasso_cf,
  run_type = "recv",
  cpn_crit = "loss"
)
lasso_recv
```

The `recv` loss lets the normalized lambda adapt to each fitted interval. To
select one normalized lambda setting jointly with K, request and compare the
homogeneous-hyperparameter paths:

```r
lasso_homogeneous <- select_across_runs(
  result = fit_lasso_cf,
  run_type = "crossfit_homo_hyper",
  cpn_crit = "loss"
)
lasso_homogeneous

plot(
  fit_lasso_cf,
  x_axis = "hyperparameter",
  K = lasso_homogeneous$K_hat[[1L]],
  run_type = "crossfit_homo_hyper",
  cpn_crit = "loss",
  log = "x"
)
```

If either unpenalized ReCV rule selects too many changepoints, replace
`cpn_crit = "loss"` with `cpn_crit = "sic"` to add a `log(n) * K` penalty.
For ordinary lasso without interval-level cross-fitting, use
`cv.reliever(X = reg_x, y = reg_y, cpd_family = "lasso")` to select the
normalized lambda setting and K by outer CV.

## Built-in cpd families

The following methods are available directly through `cpd_family`:

| Analysis | `cpd_family` | Recommended entry point | Selection |
| --- | --- | --- | --- |
| Mean changes | `"mean"` | `cv.reliever()` | CPSS outer CV for K |
| Covariance or joint mean/covariance changes | `"var"`, `"meanvar"` | `cv.reliever()` | CPSS outer CV for K |
| Linear or generalized linear model changes | `"lm"`, `"glm"` | `cv.reliever()` | CPSS outer CV for K |
| Exponential-family distribution changes | `"em"` | `cv.reliever()` | CPSS outer CV for K |
| Ordinary lasso path | `"lasso"` | `cv.reliever()` | CPSS jointly selects lambda and K |
| Cross-fitting lasso | `"lasso_crossfit"` | `reliever()` | `recv`, or `crossfit_homo_hyper` when one global lambda is required |
| KDE negative log-likelihood | `"kde_nll"`, `"kde_nll_crossfit"` | `reliever()` | ReCV for data-driven bandwidth selection |
| Fixed-kernel KDE-L2 | `"kde_l2"` | `cv.reliever()` | CPSS outer CV for K with fixed kernel and bandwidth |
| Univariate nonparametric CDF loss | `"nmcd"` | `cv.reliever()` | CPSS outer CV for K |
| Optional classifier losses | `"ranger_crossfit"`, `"mlp_crossfit"` | `reliever()` | `recv`, or `crossfit_homo_hyper` for one global model setting |

Cross-fitting fits return only `recv` by default. Include
`loss_output_types = c("recv", "crossfit_homo_hyper")` when homogeneous
hyperparameter selection will be needed. The order supplied to
`loss_output_types` does not affect the output ordering, and `recv` is required.

## Selection reference

- `cv.reliever()` and `cv.reliever_generic()` use outer held-out loss and return
  an automatically fitted full-data result.
- `select_by_run()` selects K separately within each requested stored path.
- `select_across_runs()` jointly selects one path and K across comparable
  losses.
- `select_holdout()` refits segments and selects with independent evaluation
  data.

The most common criteria are:

| `cpn_crit` | Score |
| --- | --- |
| `"loss"` | Stored loss, without a K penalty |
| `"sic"` | `loss + log(n) * K` |
| `"rss_sic"` | `n / 2 * log(RSS / n) + log(n) * K` |
| non-negative number | `loss + penalty * K` |
| `"none"` | Keep the complete path without selecting K |

`"sic"` applies directly to any stored loss, including RSS. The alternative
`"rss_sic"` first applies the log-RSS transformation and may be preferable
when residual variance is unknown or potentially heterogeneous.

## Plotting and result structure

For a cross-fitting result, omitting `run_type` defaults to `recv`. Examples
still specify the run explicitly:

```r
plot(fit_lasso_cf, run_type = "recv")
```

Hyperparameter plots must identify the homogeneous paths explicitly:

```r
plot(
  fit_lasso_cf,
  x_axis = "hyperparameter",
  K = 2,
  run_type = "crossfit_homo_hyper",
  cpn_crit = "loss",
  log = "x"
)
```

The main result components are:

| Component | Meaning |
| --- | --- |
| `summary` | Explicitly selected result; empty by default for a K-indexed fit |
| `cpd_path$candidates` | Candidate segmentations with K, loss, and changepoints |
| `run_meta` | Loss-output type and hyperparameter metadata |
| `settings` | Resolved search, cache, grid, and loss settings |
| `full_data_fit` | Final complete-data fit returned by a `cv.*` entry point |

## Main search controls

- `cpn_max`: largest changepoint count considered.
- `dm`: minimum segment length on the original time scale.
- `cov_rate`: interval coverage; `1` gives a full search and smaller values
  reduce model fits.
- `method`: `"SN"`, `"OP"`, `"PELT"`, `"WBS"`, `"WBS_recursive"`,
  `"SeedBS"`, or `"BS"`.

See `?reliever` for the full workflow and search definitions, and
`?relieverChangepoint` for a package-level overview. The
[online guides](https://qianchd.github.io/reliever_changepoint/articles/)
provide fuller workflow explanations and examples.

## Advanced features: custom losses

`reliever_generic()` accepts an interval-loss function. A virtual call declares
the number of outputs; an ordinary call returns one observation-level loss per
evaluation row. Advanced custom losses can declare several outputs and
metadata; use `run_cpd_ids` to identify deliberately comparable paths.

See the online guide
[Writing a custom `reg_fun`](https://qianchd.github.io/reliever_changepoint/articles/extending.html)
for a complete worked example, and `?reliever_generic` plus
`?cv.reliever_generic` for the precise interface contracts.
