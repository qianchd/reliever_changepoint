.reliever_nobs <- function(data) {
  if (is.vector(data)) {
    return(length(data))
  }
  if (requireNamespace("R6", quietly = TRUE) && R6::is.R6(data)) {
    if (!is.function(data$.total_size)) {
      stop("R6 data objects must provide a .total_size() method.", call. = FALSE)
    }
    return(data$.total_size())
  }
  n <- nrow(data)
  if (is.null(n)) {
    stop("data must be a vector, matrix/data.frame-like object, or supported R6 data object.", call. = FALSE)
  }
  n
}

.reliever_validate_positive_integer <- function(x, name, allow_zero = FALSE) {
  lower <- if (allow_zero) 0 else 1
  x_int <- suppressWarnings(as.integer(x))
  if (!is.numeric(x) || length(x) != 1 || is.na(x) || is.na(x_int) ||
      x != x_int || x_int < lower) {
    stop(name, " must be a single integer >= ", lower, ".", call. = FALSE)
  }
  x_int
}

.reliever_validate_positive_integer_vector <- function(x, name) {
  x_int <- suppressWarnings(as.integer(x))
  if (!is.numeric(x) || length(x) == 0L ||
      anyNA(x) || anyNA(x_int) || any(x != x_int) || any(x_int < 1L)) {
    stop(name, " must contain positive integer values.", call. = FALSE)
  }
  x_int
}

.reliever_reject_renamed_cpn_max <- function(args, caller) {
  if ("L" %in% names(args)) {
    stop(
      caller, "(): argument L was renamed to cpn_max.",
      call. = FALSE
    )
  }
  invisible(args)
}

.reliever_grid_from_size <- function(n, dc_grid_size) {
  n <- .reliever_validate_positive_integer(n, "n")
  dc_grid_size <- .reliever_validate_positive_integer(
    dc_grid_size, "dc_grid_size"
  )
  interior <- if (dc_grid_size < n) {
    seq.int(dc_grid_size, n - 1L, by = dc_grid_size)
  } else {
    integer()
  }
  as.integer(c(interior, n))
}

.reliever_resolve_dc_grid <- function(data, dc_grid_size = NULL,
                                      dc_grid = NULL) {
  if (!is.null(dc_grid_size) && !is.null(dc_grid)) {
    stop("Supply only one of dc_grid_size and dc_grid.", call. = FALSE)
  }
  if (is.null(dc_grid_size)) {
    return(dc_grid)
  }
  .reliever_grid_from_size(
    n = .reliever_nobs(data),
    dc_grid_size = dc_grid_size
  )
}

.reliever_with_seed <- function(seed, expr) {
  if (is.null(seed)) {
    return(force(expr))
  }
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) {
    old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  }
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed)
  force(expr)
}
