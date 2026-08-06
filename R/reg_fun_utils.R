.reliever_loss_matrix <- function(loss) {
  if (is.null(loss)) {
    stop("reg_fun must return a list containing `loss`.", call. = FALSE)
  }
  if (is.vector(loss)) {
    matrix(loss, ncol = 1)
  } else {
    loss
  }
}

.reliever_uses_builtin_data_preparation <- function(reg_fun) {
  any(vapply(
    list(
      reg_fun_var,
      reg_fun_meanvar,
      reg_fun_lm,
      reg_fun_glm,
      reg_fun_em,
      reg_fun_lasso_solpath,
      reg_fun_lasso_crossfit,
      reg_fun_kde_nll_solpath,
      reg_fun_kde_nll_crossfit
    ),
    identical,
    logical(1L),
    y = reg_fun
  ))
}

.reliever_prepare_builtin_reg_fun_data <- function(data, reg_fun, args) {
  if (identical(reg_fun, reg_fun_var) ||
      identical(reg_fun, reg_fun_meanvar) ||
      identical(reg_fun, reg_fun_em)) {
    return(.reliever_prepare_parametric_data(data))
  }
  if (identical(reg_fun, reg_fun_lm)) {
    intercept <- if (is.null(args$intercept)) TRUE else args$intercept
    return(.reliever_prepare_parametric_data(
      data, response_ncol = 1L, intercept = intercept
    ))
  }
  if (identical(reg_fun, reg_fun_glm)) {
    response_ncol <- if (is.null(args$response_ncol)) {
      1L
    } else {
      args$response_ncol
    }
    intercept <- if (is.null(args$intercept)) TRUE else args$intercept
    family <- if (is.null(args$family)) stats::gaussian() else args$family
    return(.reliever_prepare_parametric_data(
      data,
      response_ncol = response_ncol,
      intercept = intercept,
      family = family
    ))
  }
  if (identical(reg_fun, reg_fun_lasso_solpath) ||
      identical(reg_fun, reg_fun_lasso_crossfit)) {
    family <- if (is.null(args$family)) "gaussian" else args$family
    data <- .reliever_lasso_data_matrix(data)
    family <- .reliever_validate_glmnet_family(family)
    .reliever_validate_lasso_response(data[, 1L], family)
    return(data)
  }
  if (identical(reg_fun, reg_fun_kde_nll_solpath) ||
      identical(reg_fun, reg_fun_kde_nll_crossfit)) {
    square <- is.null(attr(
      data, ".reliever_external_train_n", exact = TRUE
    ))
    .kde_validate_dist_sq(data, square = square)
    return(data)
  }
  data
}

.reliever_make_individual_loss_fun <- function(data, reg_fun, dc_grid = NULL,
                                               para_list,
                                               n_loss_outputs) {
  data <- .reliever_prepare_builtin_reg_fun_data(data, reg_fun, para_list)
  grid_left <- if (is.null(dc_grid)) {
    NULL
  } else {
    c(1L, dc_grid[-length(dc_grid)] + 1L)
  }

  function(l, r, l_end, r_end) {
    if (r < l) {
      stop("r must be no smaller than l.", call. = FALSE)
    }

    fit_l <- l
    fit_r <- r
    eval_l <- l_end
    eval_r <- r_end
    if (!is.null(dc_grid)) {
      fit_l <- grid_left[l]
      fit_r <- dc_grid[r]
      eval_l <- grid_left[l_end]
      eval_r <- dc_grid[r_end]
    }

    fit <- do.call(
      reg_fun,
      c(
        list(
          data = data,
          l = fit_l,
          r = fit_r,
          l_end = eval_l,
          r_end = eval_r,
          save_model = FALSE
        ),
        para_list
      )
    )
    observation_loss <- .reliever_loss_matrix(fit[["loss"]])

    if (is.null(dc_grid)) {
      loss <- observation_loss
    } else {
      loss <- matrix(0, r_end - l_end + 1L, ncol(observation_loss))
      for (grid_id in l_end:r_end) {
        obs_l <- grid_left[grid_id] - eval_l + 1L
        obs_r <- dc_grid[grid_id] - eval_l + 1L
        loss[grid_id - l_end + 1L, ] <- colSums(
          observation_loss[obs_l:obs_r, , drop = FALSE]
        )
      }
    }

    expected_n <- r_end - l_end + 1L
    if (nrow(loss) != expected_n) {
      stop(
        sprintf(
          "reg_fun returned %d observation losses, but reliever expected %d.",
          nrow(loss), expected_n
        ),
        call. = FALSE
      )
    }
    if (ncol(loss) != n_loss_outputs) {
      stop(
        sprintf(
          paste0(
            "reg_fun returned %d loss outputs, but its virtual call reported ",
            "%d. Return metadata such as lambda separately from the loss matrix."
          ),
          ncol(loss), n_loss_outputs
        ),
        call. = FALSE
      )
    }
    loss
  }
}

.reliever_stable_col_means <- function(x) {
  center <- as.vector(stable_col_means_cpp(x))
  names(center) <- colnames(x)
  center
}

.reg_fun_mean_loss <- function(data, l, r, l_end, r_end, save_model) {
  x <- if (is.vector(data)) matrix(data, ncol = 1L) else as.matrix(data)
  if (!is.numeric(x)) {
    stop("data must be numeric.", call. = FALSE)
  }
  center <- .reliever_stable_col_means(x[l:r, , drop = FALSE])
  loss <- matrix(
    rowMeans(sweep(x[l_end:r_end, , drop = FALSE], 2L, center, "-")^2),
    ncol = 1L
  )
  list(
    loss = loss,
    model = if (isTRUE(save_model)) list(center = center) else NULL
  )
}

.kde_validate_dist_sq <- function(dist_sq, square = TRUE) {
  if (!is.matrix(dist_sq) || !is.numeric(dist_sq)) {
    stop("KDE reg_fun functions require a numeric distance matrix.",
         call. = FALSE)
  }
  if (isTRUE(square) && nrow(dist_sq) != ncol(dist_sq)) {
    stop("dist_sq must be a square squared-distance matrix.", call. = FALSE)
  }
  if (anyNA(dist_sq) || any(!is.finite(dist_sq)) || any(dist_sq < 0)) {
    stop(
      "The distance matrix must contain finite non-negative values.",
      call. = FALSE
    )
  }
  invisible(dist_sq)
}

.reliever_glmnet_loss <- function(model, x, y, family, s = NULL) {
  x <- as.matrix(x)
  predict_args <- list(object = model, newx = x)
  if (!is.null(s)) {
    predict_args$s <- s
  }
  if (family %in% c("binomial", "poisson")) {
    predict_args$type <- "response"
  }
  preds <- do.call(predict, predict_args)

  if (family == "gaussian") {
    loss <- (y - preds)^2
  } else if (family == "binomial") {
    preds <- pmax(pmin(preds, 1 - 1e-15), 1e-15)
    loss <- -(y * log(preds) + (1 - y) * log(1 - preds))
  } else if (family == "poisson") {
    preds <- pmax(preds, 1e-15)
    loss <- -(y * log(preds) - preds)
  } else {
    stop("Unsupported glmnet family: ", family)
  }
  .reliever_loss_matrix(loss)
}
