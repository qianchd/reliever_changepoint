.kernel_validate_args <- function(kernel_args) {
  if (is.null(kernel_args)) {
    return(list())
  }
  if (!is.list(kernel_args)) {
    stop("kernel_args must be a named list.", call. = FALSE)
  }
  if (length(kernel_args) > 0L) {
    arg_names <- names(kernel_args)
    if (is.null(arg_names) || anyNA(arg_names) || any(!nzchar(arg_names)) ||
        anyDuplicated(arg_names)) {
      stop(
        "Every element of kernel_args must have a unique non-empty name.",
        call. = FALSE
      )
    }
  }
  kernel_args
}

.kernel_reject_unused_args <- function(kernel_args, allowed, kernel) {
  unused <- setdiff(names(kernel_args), allowed)
  if (length(unused) > 0L) {
    stop(
      "Unused kernel_args for kernel = ", dQuote(kernel), ": ",
      paste(unused, collapse = ", "), ".",
      call. = FALSE
    )
  }
  invisible(kernel_args)
}

.kernel_positive_scalar <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) ||
      !is.finite(x) || x <= 0) {
    stop(name, " must be a positive finite number.", call. = FALSE)
  }
  as.numeric(x)
}

.kernel_nonnegative_scalar <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) ||
      !is.finite(x) || x < 0) {
    stop(name, " must be a finite non-negative number.", call. = FALSE)
  }
  as.numeric(x)
}

.kernel_data_matrix <- function(data, name = "data") {
  if (is.vector(data) && is.null(dim(data))) {
    data <- matrix(data, ncol = 1L)
  } else if (is.data.frame(data)) {
    if (!all(vapply(data, is.numeric, logical(1L)))) {
      stop(name, " must contain only numeric columns.", call. = FALSE)
    }
    data <- as.matrix(data)
  } else if (is.matrix(data)) {
    data <- as.matrix(data)
  } else {
    stop(
      name, " must be a numeric vector, matrix, or data frame.",
      call. = FALSE
    )
  }
  if (!is.numeric(data) || nrow(data) < 1L || ncol(data) < 1L ||
      anyNA(data) || any(!is.finite(data))) {
    stop(
      name, " must contain finite numeric observations in rows.",
      call. = FALSE
    )
  }
  data
}

.kernel_dist_sq <- function(x, y = NULL) {
  x <- .kernel_data_matrix(x, "x")
  if (is.null(y)) {
    y <- x
  } else {
    y <- .kernel_data_matrix(y, "y")
    if (ncol(x) != ncol(y)) {
      stop("x and y must have the same number of columns.", call. = FALSE)
    }
  }
  same_points <- identical(x, y)
  if (same_points) {
    return(unname(as.matrix(stats::dist(x))^2))
  }

  # The usual norm expansion can lose every meaningful digit when a data set
  # contains both a large outlier and a tightly clustered group. Direct
  # coordinate differences retain translation invariance and use only one
  # additional cross-distance-sized work matrix.
  out <- matrix(0, nrow(x), nrow(y))
  for (j in seq_len(ncol(x))) {
    difference <- outer(
      unname(x[, j]), unname(y[, j]), "-"
    )
    out <- out + difference^2
  }
  dimnames(out) <- NULL
  out
}

.kde_validate_var_dim <- function(var_dim) {
  var_dim_int <- suppressWarnings(as.integer(var_dim))
  if (!is.numeric(var_dim) || length(var_dim) != 1L || is.na(var_dim) ||
      is.na(var_dim_int) || !is.finite(var_dim) ||
      var_dim != var_dim_int || var_dim_int < 1L) {
    stop("var_dim must be a positive integer.", call. = FALSE)
  }
  var_dim_int
}

.kde_validate_distance_power <- function(distance_power) {
  if (!is.numeric(distance_power) || length(distance_power) != 1L ||
      is.na(distance_power) || !distance_power %in% c(1, 2)) {
    stop("distance_power must be 1 or 2.", call. = FALSE)
  }
  as.integer(distance_power)
}

.kde_kernel_settings <- function(kernel, kernel_args = list()) {
  if (!is.character(kernel) || length(kernel) != 1L ||
      is.na(kernel) || !nzchar(kernel)) {
    stop("kernel must be one non-empty character string.", call. = FALSE)
  }
  kernel <- tolower(kernel)
  kernel <- match.arg(kernel, c("gaussian", "laplace", "student"))
  kernel_args <- .kernel_validate_args(kernel_args)
  if (kernel == "student") {
    .kernel_reject_unused_args(kernel_args, "df", kernel)
    df <- if (is.null(kernel_args$df)) 5 else kernel_args$df
    kernel_args <- list(df = .kernel_positive_scalar(df, "kernel_args$df"))
  } else {
    .kernel_reject_unused_args(kernel_args, character(), kernel)
    kernel_args <- list()
  }
  list(kernel = kernel, kernel_args = kernel_args)
}

.kde_is_distance_like <- function(data) {
  if (!is.matrix(data) || !is.numeric(data) ||
      nrow(data) != ncol(data) || anyNA(data) ||
      any(!is.finite(data)) || any(data < 0)) {
    return(FALSE)
  }
  scale <- max(1, max(abs(data)))
  tolerance <- 1e-8 * scale
  max(abs(diag(data))) <= tolerance &&
    max(abs(data - t(data))) <= tolerance
}

.kde_prepare_input <- function(data, var_dim = NULL,
                               input_type = c("auto", "data", "distance")) {
  input_type <- match.arg(input_type)
  distance_like <- .kde_is_distance_like(data)
  if (input_type == "auto") {
    if (distance_like && !is.null(var_dim)) {
      input_type <- "distance"
    } else if (distance_like) {
      stop(
        "data looks like a squared-distance matrix. Supply var_dim, or set ",
        "input_type = \"data\" if these are raw observations.",
        call. = FALSE
      )
    } else {
      input_type <- "data"
    }
  }

  if (input_type == "distance") {
    if (!.kde_is_distance_like(data)) {
      stop(
        "A precomputed squared-distance matrix must be symmetric with a ",
        "zero diagonal.",
        call. = FALSE
      )
    }
    var_dim <- .kde_validate_var_dim(var_dim)
    return(list(
      data = data,
      var_dim = var_dim,
      input_spec = list(
        type = "kde_nll",
        original_form = "squared_distance",
        var_dim = var_dim
      )
    ))
  }

  x <- .kernel_data_matrix(data)
  inferred_dim <- ncol(x)
  if (!is.null(var_dim) &&
      .kde_validate_var_dim(var_dim) != inferred_dim) {
    stop(
      "var_dim must equal ncol(data) when input_type = \"data\".",
      call. = FALSE
    )
  }
  list(
    data = .kernel_dist_sq(x),
    var_dim = inferred_dim,
    input_spec = list(
      type = "kde_nll",
      original_form = "raw_data",
      var_dim = inferred_dim,
      n_features = inferred_dim
    )
  )
}

.kernel_l2_settings <- function(kernel, bandwidth = NULL,
                                kernel_args = list()) {
  kernel_args <- .kernel_validate_args(kernel_args)
  if (is.function(kernel)) {
    reserved <- intersect(names(kernel_args), c("x", "y", "bandwidth"))
    if (length(reserved) > 0L) {
      stop(
        "Custom kernel_args must not use reserved names: ",
        paste(reserved, collapse = ", "), ".",
        call. = FALSE
      )
    }
    if (!is.null(bandwidth)) {
      bandwidth <- .kernel_positive_scalar(bandwidth, "bandwidth")
    }
    return(list(
      kernel = kernel,
      kernel_name = "custom",
      bandwidth = bandwidth,
      kernel_args = kernel_args
    ))
  }
  if (!is.character(kernel) || length(kernel) != 1L ||
      is.na(kernel) || !nzchar(kernel)) {
    stop(
      "kernel must be NULL, a supported kernel name, or a function.",
      call. = FALSE
    )
  }
  kernel <- tolower(kernel)
  if (kernel == "rbf") {
    kernel <- "gaussian"
  }
  kernel <- match.arg(kernel, c(
    "gaussian", "laplace", "linear", "polynomial",
    "matern32", "matern52", "rational_quadratic"
  ))

  radial <- kernel %in% c(
    "gaussian", "laplace", "matern32", "matern52",
    "rational_quadratic"
  )
  if (radial) {
    bandwidth <- .kernel_positive_scalar(bandwidth, "bandwidth")
  } else if (!is.null(bandwidth)) {
    stop(
      "bandwidth is used only by radial kernels; leave it NULL for kernel = ",
      dQuote(kernel), ".",
      call. = FALSE
    )
  }

  if (kernel == "polynomial") {
    .kernel_reject_unused_args(
      kernel_args, c("degree", "scale", "offset"), kernel
    )
    degree <- if (is.null(kernel_args$degree)) 2L else kernel_args$degree
    degree_int <- suppressWarnings(as.integer(degree))
    if (!is.numeric(degree) || length(degree) != 1L || is.na(degree) ||
        is.na(degree_int) || degree != degree_int || degree_int < 1L) {
      stop("kernel_args$degree must be a positive integer.", call. = FALSE)
    }
    scale <- if (is.null(kernel_args$scale)) 1 else kernel_args$scale
    offset <- if (is.null(kernel_args$offset)) 1 else kernel_args$offset
    kernel_args <- list(
      degree = degree_int,
      scale = .kernel_nonnegative_scalar(scale, "kernel_args$scale"),
      offset = .kernel_nonnegative_scalar(offset, "kernel_args$offset")
    )
  } else if (kernel == "rational_quadratic") {
    .kernel_reject_unused_args(kernel_args, "alpha", kernel)
    alpha <- if (is.null(kernel_args$alpha)) 1 else kernel_args$alpha
    kernel_args <- list(
      alpha = .kernel_positive_scalar(alpha, "kernel_args$alpha")
    )
  } else {
    .kernel_reject_unused_args(kernel_args, character(), kernel)
    kernel_args <- list()
  }

  list(
    kernel = kernel,
    kernel_name = kernel,
    bandwidth = bandwidth,
    kernel_args = kernel_args
  )
}

.kernel_l2_matrix <- function(x, y = NULL, kernel,
                              bandwidth = NULL, kernel_args = list()) {
  x <- .kernel_data_matrix(x, "x")
  y_missing <- is.null(y)
  y <- if (y_missing) x else .kernel_data_matrix(y, "y")
  if (ncol(x) != ncol(y)) {
    stop("x and y must have the same number of columns.", call. = FALSE)
  }
  settings <- .kernel_l2_settings(kernel, bandwidth, kernel_args)
  if (is.function(settings$kernel)) {
    call_args <- c(list(x = x, y = y), settings$kernel_args)
    if (!is.null(settings$bandwidth)) {
      call_args$bandwidth <- settings$bandwidth
    }
    out <- do.call(settings$kernel, call_args)
  } else if (settings$kernel %in% c("linear", "polynomial")) {
    inner_product <- tcrossprod(x, y)
    out <- if (settings$kernel == "linear") {
      inner_product
    } else {
      (
        settings$kernel_args$scale * inner_product +
          settings$kernel_args$offset
      )^settings$kernel_args$degree
    }
  } else {
    dist_sq <- .kernel_dist_sq(x, y)
    bandwidth <- settings$bandwidth
    out <- switch(
      settings$kernel,
      gaussian = {
        value <- exp(-dist_sq / (2 * bandwidth^2))
        value[dist_sq == 0] <- 1
        value
      },
      laplace = exp(-sqrt(dist_sq) / bandwidth),
      matern32 = {
        scaled <- sqrt(3 * dist_sq) / bandwidth
        value <- exp(log1p(scaled) - scaled)
        value[is.infinite(scaled)] <- 0
        value[dist_sq == 0] <- 1
        value
      },
      matern52 = {
        scaled <- sqrt(5 * dist_sq) / bandwidth
        log_polynomial <- numeric(length(scaled))
        dim(log_polynomial) <- dim(scaled)
        small <- is.finite(scaled) & scaled <= 1
        log_polynomial[small] <- log1p(
          scaled[small] + scaled[small]^2 / 3
        )
        large <- is.finite(scaled) & !small
        log_polynomial[large] <-
          2 * log(scaled[large]) - log(3) +
          log1p(3 / scaled[large] + 3 / scaled[large]^2)
        value <- exp(log_polynomial - scaled)
        value[is.infinite(scaled)] <- 0
        value[dist_sq == 0] <- 1
        value
      },
      rational_quadratic = {
        alpha <- settings$kernel_args$alpha
        value <-
          (1 + dist_sq / (2 * alpha * bandwidth^2))^(-alpha)
        value[dist_sq == 0] <- 1
        value
      }
    )
  }
  if (!is.matrix(out) || !is.numeric(out) ||
      nrow(out) != nrow(x) || ncol(out) != nrow(y) ||
      anyNA(out) || any(!is.finite(out))) {
    stop(
      "kernel must return a finite numeric matrix with nrow(x) rows and ",
      "nrow(y) columns.",
      call. = FALSE
    )
  }
  if (y_missing) {
    scale <- max(1, max(abs(out)))
    if (max(abs(out - t(out))) > 1e-8 * scale) {
      stop("A Gram matrix returned by kernel must be symmetric.", call. = FALSE)
    }
    out <- (out + t(out)) / 2
  }
  out
}

.kernel_l2_prepare_input <- function(data, kernel = NULL, bandwidth = NULL,
                                     kernel_args = list()) {
  if (is.null(kernel)) {
    if (!is.null(bandwidth)) {
      stop(
        "bandwidth requires kernel; leave both NULL for precomputed features.",
        call. = FALSE
      )
    }
    kernel_args <- .kernel_validate_args(kernel_args)
    if (length(kernel_args) > 0L) {
      stop(
        "kernel_args requires kernel; leave both empty for precomputed features.",
        call. = FALSE
      )
    }
    if (is.vector(data) && is.null(dim(data))) {
      prepared <- data
      n_features <- 1L
    } else {
      prepared <- .kernel_data_matrix(data)
      n_features <- ncol(prepared)
    }
    return(list(
      data = prepared,
      input_spec = list(
        type = "kernel_features",
        original_form = "precomputed_features",
        n_features = n_features
      )
    ))
  }

  x <- .kernel_data_matrix(data)
  settings <- .kernel_l2_settings(kernel, bandwidth, kernel_args)
  list(
    data = .kernel_l2_matrix(
      x,
      kernel = settings$kernel,
      bandwidth = settings$bandwidth,
      kernel_args = settings$kernel_args
    ),
    input_spec = list(
      type = "kernel_features",
      original_form = "raw_data",
      n_features = ncol(x),
      kernel = settings$kernel,
      kernel_name = settings$kernel_name,
      bandwidth = settings$bandwidth,
      kernel_args = settings$kernel_args
    )
  )
}

.kde_resolve_bandwidth_vec <- function(bandwidth_vec, var_dim,
                                       dist_sq = NULL,
                                       kernel = "gaussian",
                                       kernel_args = list()) {
  var_dim <- .kde_validate_var_dim(var_dim)
  kernel_settings <- .kde_kernel_settings(kernel, kernel_args)
  if (is.null(bandwidth_vec)) {
    .kde_validate_dist_sq(dist_sq)
    n <- nrow(dist_sq)
    if (n < 2L) {
      stop("At least two observations are needed to create bandwidth_vec.",
           call. = FALSE)
    }
    max_lag <- min(n - 1L, ceiling(sqrt(n)))
    pair_dist_sq <- unlist(lapply(seq_len(max_lag), function(lag) {
      dist_sq[cbind(seq_len(n - lag), seq.int(lag + 1L, n))]
    }), use.names = FALSE)
    if (any(pair_dist_sq < 0, na.rm = TRUE)) {
      stop("dist_sq cannot contain negative squared distances.", call. = FALSE)
    }
    pair_dist_sq <- pair_dist_sq[is.finite(pair_dist_sq) & pair_dist_sq > 0]
    if (!length(pair_dist_sq)) {
      stop(
        "Cannot create a default bandwidth_vec without positive finite distances.",
        call. = FALSE
      )
    }
    distance_scale <- stats::median(pair_dist_sq)
    bandwidth_vec <- exp(seq(
      log(sqrt(distance_scale / 120)),
      log(sqrt(distance_scale)),
      length.out = 20L
    ))
    if (kernel_settings$kernel == "laplace") {
      bandwidth_vec <- bandwidth_vec / sqrt(var_dim + 1)
    } else if (kernel_settings$kernel == "student") {
      df <- kernel_settings$kernel_args$df
      if (df <= 2) {
        stop(
          "bandwidth_vec must be supplied for Student kernels with df <= 2 ",
          "because their marginal variance is not finite.",
          call. = FALSE
        )
      }
      bandwidth_vec <- bandwidth_vec * sqrt((df - 2) / df)
    }
  }
  bandwidth_vec <- as.numeric(unlist(bandwidth_vec, use.names = FALSE))
  if (length(bandwidth_vec) < 1L || anyNA(bandwidth_vec) ||
      any(!is.finite(bandwidth_vec)) || any(bandwidth_vec <= 0)) {
    stop("bandwidth_vec must contain positive finite values.", call. = FALSE)
  }
  bandwidth_vec
}

.kde_prepare_fit_input <- function(data, var_dim, input_type,
                                   kernel, kernel_args, bandwidth_vec) {
  prepared <- .kde_prepare_input(data, var_dim, input_type)
  kernel_settings <- .kde_kernel_settings(kernel, kernel_args)
  bandwidth_vec <- .kde_resolve_bandwidth_vec(
    bandwidth_vec = bandwidth_vec,
    var_dim = prepared$var_dim,
    dist_sq = prepared$data,
    kernel = kernel_settings$kernel,
    kernel_args = kernel_settings$kernel_args
  )
  distance_power <- 2L
  if (kernel_settings$kernel == "laplace") {
    prepared$data <- sqrt(prepared$data)
    distance_power <- 1L
  }
  prepared$input_spec$kernel <- kernel_settings$kernel
  prepared$input_spec$kernel_args <- kernel_settings$kernel_args
  prepared$input_spec$distance_power <- distance_power
  c(
    prepared,
    list(
      kernel = kernel_settings$kernel,
      kernel_args = kernel_settings$kernel_args,
      bandwidth_vec = bandwidth_vec,
      distance_power = distance_power
    )
  )
}

.kde_bandwidth_reference <- function(data, bandwidth_vec, distance_power) {
  if (!is.null(bandwidth_vec)) {
    return(data)
  }
  n_train <- attr(data, ".reliever_external_train_n", exact = TRUE)
  if (!is.null(n_train)) {
    n_train_int <- suppressWarnings(as.integer(n_train))
    if (!is.numeric(n_train) || length(n_train) != 1L || is.na(n_train) ||
        is.na(n_train_int) || n_train != n_train_int || n_train_int < 1L ||
        !is.matrix(data) || nrow(data) < n_train_int ||
        ncol(data) < n_train_int) {
      stop("External KDE evaluation has an invalid training-row marker.",
           call. = FALSE)
    }
    data <- data[
      seq_len(n_train_int), seq_len(n_train_int), drop = FALSE
    ]
  }
  if (distance_power == 1L) data^2 else data
}

.kde_nll_loss <- function(dist_sub, bandwidth_vec, var_dim,
                          kernel = "gaussian", kernel_args = list(),
                          distance_power = 2L) {
  if (any(!is.finite(dist_sub)) || any(dist_sub < 0)) {
    stop("dist_sq must contain finite non-negative values.", call. = FALSE)
  }
  distance_power <- .kde_validate_distance_power(distance_power)
  var_dim <- .kde_validate_var_dim(var_dim)
  kernel_settings <- .kde_kernel_settings(kernel, kernel_args)
  kernel <- kernel_settings$kernel
  kernel_args <- kernel_settings$kernel_args
  bandwidth_vec <- .kde_resolve_bandwidth_vec(
    bandwidth_vec, var_dim, kernel = kernel, kernel_args = kernel_args
  )
  loss <- matrix(0, nrow(dist_sub), length(bandwidth_vec))
  distance <- if (kernel == "laplace") {
    if (distance_power == 1) dist_sub else sqrt(dist_sub)
  } else {
    NULL
  }
  dist_sq <- if (kernel != "laplace") {
    if (distance_power == 2) dist_sub else dist_sub^2
  } else {
    NULL
  }
  radial_value <- if (kernel == "laplace") distance else dist_sq
  min_col <- max.col(-radial_value, ties.method = "first")
  row_min <- radial_value[cbind(seq_len(nrow(radial_value)), min_col)]
  centered_value <- radial_value - row_min
  max_radial_value <- max(radial_value)

  for (j in seq_along(bandwidth_vec)) {
    bandwidth <- bandwidth_vec[j]
    if (kernel == "gaussian") {
      log_relative <-
        -0.5 * ((centered_value / bandwidth) / bandwidth)
      log_shape_max <- -0.5 * ((row_min / bandwidth) / bandwidth)
      log_normalizer <-
        -(var_dim / 2) * log(2 * pi) - var_dim * log(bandwidth)
    } else if (kernel == "laplace") {
      log_relative <- -centered_value / bandwidth
      log_shape_max <- -row_min / bandwidth
      log_normalizer <-
        lgamma(var_dim / 2) - log(2) -
        (var_dim / 2) * log(pi) - lgamma(var_dim) -
        var_dim * log(bandwidth)
    } else {
      df <- kernel_args$df
      tail_power <- (df + var_dim) / 2
      scale <- df * bandwidth^2
      use_log_domain <-
        !is.finite(tail_power) || !is.finite(scale) || scale <= 0 ||
        !is.finite(max_radial_value / scale) ||
        any(!is.finite(scale + row_min))
      if (use_log_domain) {
        log_scale <- log(df) + 2 * log(bandwidth)
        log_ratio <- log(radial_value) - log_scale
        log_ratio_min <- log(row_min) - log_scale
        softplus_ratio <- (
          pmax(log_ratio, 0) + log1p(exp(-abs(log_ratio)))
        )
        softplus_ratio_min <- (
          pmax(log_ratio_min, 0) +
            log1p(exp(-abs(log_ratio_min)))
        )
        log_relative <-
          -tail_power * (softplus_ratio - softplus_ratio_min)
        log_relative[centered_value == 0] <- 0
        log_shape_max <- -tail_power * softplus_ratio_min
        log_shape_max[row_min == 0] <- 0
      } else {
        relative_arg <- centered_value / (scale + row_min)
        log_relative <- -tail_power * log1p(relative_arg)
        log_shape_max <- -tail_power * log1p(row_min / scale)
      }
      log_normalizer <-
        lgamma(var_dim / 2) - lbeta(df / 2, var_dim / 2) -
        (var_dim / 2) * (log(df) + log(pi)) -
        var_dim * log(bandwidth)
    }
    log_kernel_mean <-
      log_shape_max + log(rowMeans(exp(log_relative)))
    loss[, j] <- -(log_normalizer + log_kernel_mean)
  }
  loss
}

#' Radial-kernel KDE negative-log-likelihood interval loss
#'
#' Fit an isotropic kernel density estimate on rows \code{l:r} and return the
#' negative log density of rows \code{l_end:r_end}, once for each bandwidth in
#' \code{bandwidth_vec}. The low-level input is an \eqn{n} by \eqn{n}
#' squared-Euclidean-distance matrix, so pairwise distances are computed only
#' once. Most users can instead pass the original observations to
#' \code{\link{reliever}()}.
#'
#' In dimension \eqn{p}, Gaussian uses the normalized density
#' \deqn{(2\pi)^{-p/2}h^{-p}\exp\{-r^2/(2h^2)\},}
#' radial Laplace uses
#' \deqn{\frac{\Gamma(p/2)}{2\pi^{p/2}\Gamma(p)}
#' h^{-p}\exp\{-r/h\},}
#' and Student with degrees of freedom \eqn{\nu} uses
#' \deqn{\frac{\Gamma((\nu+p)/2)}
#' {\Gamma(\nu/2)(\nu\pi)^{p/2}h^p}
#' \{1+r^2/(\nu h^2)\}^{-(\nu+p)/2},}
#' where \eqn{r=\|x-y\|_2}. All calculations use row-wise log-sum-exp rather
#' than materializing one full kernel matrix per bandwidth.
#' One Euclidean-distance matrix can therefore be reused for any isotropic
#' radial kernel. Anisotropic, coordinate-product, or observation-adaptive
#' kernels require a different representation or a custom \code{reg_fun}.
#'
#' @param data Square squared-distance matrix.
#' @param l,r Rows used to fit the interval KDE.
#' @param l_end,r_end Rows whose negative log densities are returned.
#' @param save_model Accepted so this function can be used as a \code{reg_fun};
#'   KDE models are not retained.
#' @param is_virtual_run Query flag used by \code{reliever_generic()}. When
#'   \code{TRUE}, return only the number and metadata of loss outputs.
#' @param bandwidth_vec Positive KDE bandwidths. \code{NULL} uses the
#'   scale-adaptive grid documented in \code{\link{reliever_kde_nll}()}.
#' @param var_dim Dimension of the original observations used to construct the
#'   squared distances. It appears in the density normalization constant and
#'   must be a positive integer because it cannot be recovered from a distance
#'   matrix.
#' @param kernel Isotropic density kernel: \code{"gaussian"},
#'   \code{"laplace"}, or \code{"student"}.
#' @param kernel_args Named list of kernel-specific settings; Student accepts
#'   positive \code{df} (default 5).
#' @param distance_power Advanced representation flag. The ordinary value 2
#'   means that \code{data} contains squared distances. The high-level Laplace
#'   wrapper internally uses 1 after taking the square root once for the whole
#'   data set; users normally should not change this argument.
#'
#' @return With \code{is_virtual_run = TRUE}, metadata linking each loss output
#'   to its bandwidth. Otherwise, a list containing an observation-by-bandwidth
#'   \code{loss} matrix and \code{model = NULL}. See the
#'   \emph{Writing a custom reg_fun} section of
#'   \code{\link{reliever_generic}()} for the common return convention.
#' @seealso \code{\link{reliever_kde_nll}()},
#'   \code{\link{reg_fun_kde_nll_crossfit}()}, and
#'   \code{\link{reliever_generic}()}.
#' @references Qian, C., Wang, G., Wang, Z., and Zou, C. (2024).
#'   Changepoint Detection in Complex Models: Cross-Fitting Is Needed.
#'   arXiv:2411.07874.
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
#' dist_sq <- as.matrix(dist(x))^2
#' out <- reg_fun_kde_nll_solpath(dist_sq, 241, 660, 211, 690,
#'                        var_dim = 5)
#' dim(out$loss)
reg_fun_kde_nll_solpath <- function(data, l, r, l_end = l, r_end = r,
                            save_model = FALSE, is_virtual_run = FALSE,
                            bandwidth_vec = NULL, var_dim = NULL,
                            kernel = "gaussian", kernel_args = list(),
                            distance_power = 2L) {
  distance_power <- .kde_validate_distance_power(distance_power)
  bandwidth_data <- .kde_bandwidth_reference(
    data, bandwidth_vec, distance_power
  )
  bandwidth_vec <- .kde_resolve_bandwidth_vec(
    bandwidth_vec, var_dim, bandwidth_data, kernel = kernel,
    kernel_args = kernel_args
  )
  kernel_settings <- .kde_kernel_settings(kernel, kernel_args)
  if (is_virtual_run) {
    return(list(
      n_loss_outputs = length(bandwidth_vec),
      loss_output_meta = data.frame(
        loss_output_id = seq_along(bandwidth_vec),
        row_type = "kde_nll",
        hyper_id = seq_along(bandwidth_vec),
        hyper_value = I(as.list(bandwidth_vec)),
        hyper_name = "bandwidth",
        loss_kind = "negative_log_likelihood"
      )
    ))
  }
  full_indices <- l_end:r_end
  core_indices <- l:r
  dist_sub <- data[full_indices, core_indices, drop = FALSE]
  list(
    loss = .kde_nll_loss(
      dist_sub, bandwidth_vec, var_dim,
      kernel = kernel_settings$kernel,
      kernel_args = kernel_settings$kernel_args,
      distance_power = distance_power
    ),
    model = NULL
  )
}

#' Empirical-grid KDE L2 interval loss
#'
#' This function evaluates a kernel-density \eqn{L_2} loss on a fixed common
#' grid. Let \eqn{z_1,\ldots,z_m} be that grid and let
#' \eqn{G_{ij}=k_h(X_i,z_j)} be row \eqn{i} of \code{data}. For a segment
#' \eqn{I}, its KDE on the grid is
#' \deqn{\widehat f_I(z_j)=|I|^{-1}\sum_{i\in I}G_{ij}.}
#' The observation-level interval loss returned for row \eqn{i} is
#' \deqn{\ell_i(I)=m^{-1}\sum_{j=1}^m
#' \{G_{ij}-\widehat f_I(z_j)\}^2.}
#' Summing these losses over the evaluation rows gives the empirical-grid
#' discretization of the KDE \eqn{L_2} criterion used for changepoint
#' detection.
#'
#' With a Gram matrix \eqn{K}, the common evaluation locations are the observed
#' sample points and \eqn{K_{ij}=k_h(X_i,X_j)}. The matrix is computed once and
#' reused to construct each segment KDE and evaluate its \eqn{L_2} loss; it is
#' a computational representation rather than the statistical target. A
#' generic precomputed feature matrix is also accepted; it represents KDE L2
#' when its columns are evaluations of one fixed density kernel on a common
#' grid.
#' Omitting a normalizing constant shared by that fixed kernel rescales all
#' candidate losses by the same positive factor.
#'
#' The kernel and bandwidth are fixed for this loss. For data-driven bandwidth
#' selection, use the \code{"kde_nll_crossfit"} family. For fixed KDE-L2, use
#' the \code{"kde_l2"} family or \code{\link{reliever_kde_l2}()}.
#' KDE-NLL constructs and reuses a bandwidth-independent pairwise distance
#' matrix. KDE-L2 instead constructs one fixed kernel-evaluation matrix at the
#' supplied bandwidth.
#'
#' @param data Numeric matrix whose rows are kernel-density contributions
#'   evaluated at common locations. A numeric vector is treated as a
#'   one-column matrix. A generic fixed feature matrix is also accepted as a
#'   computational extension.
#' @param l,r Rows used to estimate the segment KDE.
#' @param l_end,r_end Rows whose losses are returned.
#' @param save_model Return the fitted segment KDE as
#'   \code{model = list(center = ...)}.
#' @param is_virtual_run Query flag used by \code{reliever_generic()}. When
#'   \code{TRUE}, report one loss output without fitting.
#'
#' @return With \code{is_virtual_run = TRUE}, metadata for one KDE-L2 loss
#'   output. Otherwise, a list containing a one-column \code{loss} matrix and
#'   \code{model}. When
#'   \code{save_model = TRUE}, \code{model} is
#'   \code{list(center = ...)}, where \code{center} contains the fitted segment
#'   KDE values on the common grid, or the segment mean for a generic
#'   fixed-feature extension; otherwise it is \code{NULL}. See the
#'   \emph{Writing a custom reg_fun} section of
#'   \code{\link{reliever_generic}()} for the common return convention.
#' @references Padilla, O. H. M., Yu, Y., Wang, D., and Rinaldo, A. (2021).
#'   Optimal nonparametric multivariate change point detection and localization.
#'   \emph{IEEE Transactions on Information Theory}, 68(3), 1922--1944.
#' @seealso \code{\link{reliever_kde_l2}()},
#'   \code{\link{reliever_kde_nll_crossfit}()},
#'   \code{\link{reliever_generic}()}
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
#' dist_sq <- as.matrix(dist(x))^2
#' bandwidth <- 1
#' kernel_mat <- exp(-dist_sq / (2 * bandwidth^2))
#' out <- reg_fun_kde_l2(kernel_mat, 241, 660, 211, 690)
#' dim(out$loss)
reg_fun_kde_l2 <- function(data, l, r, l_end = l, r_end = r,
                      save_model = FALSE, is_virtual_run = FALSE) {
  if (is_virtual_run) {
    return(list(
      n_loss_outputs = 1L,
      loss_output_meta = data.frame(
        loss_output_id = 1L,
        row_type = "kde_l2",
        loss_kind = "rss"
      )
    ))
  }
  .reg_fun_mean_loss(
    data = data, l = l, r = r, l_end = l_end, r_end = r_end,
    save_model = save_model
  )
}
