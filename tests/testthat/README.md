# Test layout

The package tests are grouped by responsibility rather than by development
history:

| Group | Files | Purpose |
| --- | --- | --- |
| Public API | `test-api-*`, `test-selection.R`, `test-reg-fun*.R` | Wrappers, result selection, and the custom `reg_fun` contract |
| Search | `test-reliever-core.R`, `test-search-*`, `test-cv-outer.R` | Search paths, grids, TwoStep, outer CV, validation, and hold-out evaluation |
| Backends | `test-backend-*` | Native/generic and cost-matrix/loss-block result equivalence |
| Statistical results | `test-accuracy.R` | Two-changepoint recovery and model-selection behavior |
| Intervals | `test-interval-constructors.R`, `test-interval-routines.R` | Interval construction and relief ownership algebra |
| Utilities | `test-dgp-utils.R`, `test-documentation-coverage.R` | DGP helpers, exported documentation, README examples, and citation metadata |

## Test modes

The default `quick` mode runs all routine unit, integration, and statistical
tests, but skips optional-package checks and exhaustive relief ownership tests.

```r
devtools::test()
```

Use `full` after broad API, backend, CV, or model-wrapper changes:

```r
Sys.setenv(RELIEVER_TEST_MODE = "full")
devtools::test()
```

The two relief mapping directions are tested separately because their
partition checks enumerate exact intervals. Run this block after changing
`create_relief_itv()`, `exact2relief_itv_routine_c()`,
`relief2exact_itv_routine_c()`, or the corresponding ownership logic in
`src/cost_engine.cpp`:

```r
Sys.setenv(RELIEVER_TEST_MODE = "interval")
devtools::test(filter = "interval-routines")
```

`RELIEVER_TEST_MODE = "all"` enables both full and interval-only tests.

## Coverage

For an optional local coverage report, run:

```r
covr::package_coverage(type = "tests")
```

`covr` executes tests from a temporary installed package. Source-only checks
of `R/` and `man/` skip in that environment, while the installed `CITATION`
metadata check still runs. The ordinary development-tree test continues to
check the generated source documentation directly.

## Statistical assertions

For a stochastic DGP, first assert the selected changepoint count and then
bound the location error with `expect_cpd_error_lte()`. The default tolerance is
10 observations:

```r
expect_equal(selected$K_hat, length(truth))
expect_cpd_error_lte(selected$cpd_hat[[1L]], truth)
```

Do not require estimated changepoints to equal the true locations exactly.
Exact equality remains appropriate when testing deterministic relationships,
such as backend parity, cache reuse, index conversion, or whether a selector
returns a known row from an artificial candidate table.
