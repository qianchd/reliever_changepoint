test_that("DGP helpers return stable dimensions and identities", {
  with_test_timeout({
    set.seed(42)
    n <- 8L
    p <- 2L
    tau <- 4L
    b0 <- c(1, 0)
    delta <- matrix(c(1, -1), nrow = p)
    res <- relieverChangepoint::dgp_linear_regression(n, p, tau, b0, delta, sig = 0)

    expect_equal(dim(res$data), c(n, p + 1L))
    expect_equal(length(res$mu), n)
    expect_equal(res$ep, numeric(n))
    expect_equal(res$data[, 1], res$mu)

    set.seed(42)
    res_ar0 <- relieverChangepoint::dgp_linear_regression(
      n, p, tau, b0, delta, sig = 0, rho_ar1 = 0
    )
    expect_equal(res_ar0, res)

    set.seed(43)
    temporal <- relieverChangepoint::dgp_linear_regression(
      n, p, tau, b0, delta, sig = 0, rho_ar1 = 0.3
    )
    expect_equal(dim(temporal$data), c(n, p + 1L))
    expect_equal(temporal$data[, 1], temporal$mu)

    set.seed(44)
    equi_y <- relieverChangepoint::dgp_linear_regression_equal_response(
      n, p, tau, b0, delta, snr_y = 1, rho_ar1 = 0.2
    )
    expect_equal(dim(equi_y$data), c(n, p + 1L))
    expect_false(anyNA(equi_y$data))

    expect_equal(
      relieverChangepoint::cp_error(
        estimate = c(10, 20), truth = c(12, 18)
      ),
      2
    )
    expect_equal(relieverChangepoint::cp_error(integer(), integer()), 0)
    expect_equal(relieverChangepoint::cp_error(integer(), 10L), Inf)
    expect_equal(relieverChangepoint::cp_error(10L, integer()), Inf)
  })
})

test_that("linear-regression DGPs reject malformed statistical inputs", {
  with_test_timeout({
    n <- 12L
    p <- 2L
    tau <- c(4L, 8L)
    b0 <- c(1, 0)
    delta <- matrix(c(1, -1, -1, 1), nrow = p)

    expect_error(
      relieverChangepoint::dgp_linear_regression(
        n, p, rev(tau), b0, delta
      ),
      "strictly increasing"
    )
    expect_error(
      relieverChangepoint::dgp_linear_regression(
        n, p, tau, b0, delta[, 1L, drop = FALSE]
      ),
      "p by length"
    )
    expect_error(
      relieverChangepoint::dgp_linear_regression(
        n, p, tau, b0, delta, rho = 1
      ),
      "rho must"
    )
    expect_error(
      relieverChangepoint::dgp_linear_regression(
        n, p, tau, b0, delta, sig = -1
      ),
      "sig must"
    )
    expect_error(
      relieverChangepoint::dgp_linear_regression_equal_response(
        n, p, tau, b0, delta, snr_y = 0
      ),
      "snr_y must"
    )
  })
})

test_that("AR(1) generation handles a single observation", {
  with_test_timeout({
    generated <- relieverChangepoint::dgp_linear_regression(
      n = 1L, p = 1L, tau = integer(), b0 = 1,
      delta = matrix(numeric(), nrow = 1L), rho_ar1 = 0.3
    )
    expect_equal(dim(generated$data), c(1L, 2L))
  })
})
