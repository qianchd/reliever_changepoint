#' Univariate nonparametric CDF interval loss
#'
#' Fit the empirical distribution on rows \code{l:r} and score each row in
#' \code{l_end:r_end} by aggregating binary log loss across empirical-CDF
#' cutpoints. Because it models the full univariate distribution rather than
#' only its mean, the loss can detect changes in location, scale, or shape. The
#' calculation is implemented in C++. For changepoint fitting, most users
#' should call \code{reliever(x, cpd_family = "nmcd")}; this function is the
#' lower-level interval-loss implementation.
#'
#' Let \eqn{z_1,\ldots,z_{m-1}} be the ordered reference cutpoints,
#' \eqn{\widehat F_I} the empirical CDF fitted on interval \eqn{I}, and
#' \eqn{B_{ij}=1\{X_i\le z_j\}}. Apart from optional tail truncation, the
#' observation loss is
#' \deqn{-\sum_{j=1}^{m-1}\frac{m}{j(m-j)}
#' [B_{ij}\log\widehat F_I(z_j)
#' +(1-B_{ij})\log\{1-\widehat F_I(z_j)\}].}
#' Fitted probabilities equal to 0 or 1 are replaced by
#' \eqn{1/(2|I|)} or \eqn{1-1/(2|I|)}. This is a Reliever-compatible weighted
#' empirical-CDF likelihood based on the NMCD construction.
#'
#' @param data Numeric vector, or one-column matrix, of observations.
#' @param l,r Rows used to estimate the interval distribution.
#' @param l_end,r_end Rows whose distribution losses are returned.
#' @param save_model Accepted so this function can be used as a \code{reg_fun};
#'   no fitted object is retained.
#' @param is_virtual_run Query flag used by \code{reliever_generic()}. When
#'   \code{TRUE}, report one loss output without fitting.
#' @param w_trunc Fraction removed from each tail of the ordered full-sample CDF
#'   cutpoints. If the reference vector has length \eqn{m},
#'   \code{floor(w_trunc * m)} cutpoints are omitted from both the lower and
#'   upper tails. The default 0 keeps every cutpoint used by the criterion.
#' @param sort_X Optional sorted reference vector that defines the empirical-CDF
#'   cutpoints. The default sorts the training observations in \code{data};
#'   externally appended evaluation rows are excluded automatically. Model
#'   wrappers may resolve this reference once and reuse it.
#'
#' @return With \code{is_virtual_run = TRUE}, metadata for one additive
#'   empirical-CDF negative-log-likelihood output. Otherwise, a list containing
#'   a one-column \code{loss} matrix and \code{model = NULL}. See the
#'   \emph{Writing a custom reg_fun} section of
#'   \code{\link{reliever_generic}()} for the common return convention.
#' @references Zou, C., Yin, G., Feng, L., and Wang, Z. (2014).
#'   Nonparametric maximum likelihood approach to multiple change-point
#'   problems. \emph{The Annals of Statistics}, 42(3), 970--1002.
#'   DOI: 10.1214/14-AOS1210.
#' @seealso \code{\link{reliever_nmcd}()},
#'   \code{\link{reliever_generic}()},
#'   \code{\link{reg_fun_kde_l2}()},
#'   \code{\link{reg_fun_kde_nll_solpath}()}
#' @export
#'
#' @examples
#' set.seed(2026)
#' n_seg <- 300
#' x <- c(
#'   rnorm(n_seg, mean = 0, sd = 0.2),
#'   rnorm(n_seg, mean = 4, sd = 0.2),
#'   rnorm(n_seg, mean = -4, sd = 0.2)
#' )
#' out <- reg_fun_nmcd(x, 241, 660, 211, 690)
#' dim(out$loss)
reg_fun_nmcd <- function(data, l, r, l_end = l, r_end = r, save_model = FALSE,
                 is_virtual_run = FALSE, w_trunc = 0, sort_X = NULL) {
  if (is_virtual_run) {
    return(list(
      n_loss_outputs = 1L,
      loss_output_meta = data.frame(
        loss_output_id = 1L,
        row_type = "nmcd",
        loss_kind = "negative_log_likelihood"
      )
    ))
  }
  if (is.matrix(data) && ncol(data) > 1L) {
    stop("reg_fun_nmcd() supports only univariate data.", call. = FALSE)
  }
  external_train_n <- attr(
    data, ".reliever_external_train_n", exact = TRUE
  )
  x <- as.numeric(data)
  if (any(!is.finite(x))) {
    stop("data must contain only finite values.", call. = FALSE)
  }
  if (!is.numeric(w_trunc) || length(w_trunc) != 1L || is.na(w_trunc) ||
      !is.finite(w_trunc) || w_trunc < 0 || w_trunc >= 0.5) {
    stop("w_trunc must be a single number in [0, 1/2).", call. = FALSE)
  }
  if (is.null(sort_X)) {
    reference_x <- x
    if (!is.null(external_train_n)) {
      external_train_n_int <- suppressWarnings(as.integer(external_train_n))
      if (!is.numeric(external_train_n) || length(external_train_n) != 1L ||
          is.na(external_train_n) || is.na(external_train_n_int) ||
          external_train_n != external_train_n_int ||
          external_train_n_int < 1L || external_train_n_int >= length(x)) {
        stop(
          "data has an invalid external training-row marker.",
          call. = FALSE
        )
      }
      reference_x <- x[seq_len(external_train_n_int)]
    }
    sort_X <- sort(reference_x)
  } else {
    sort_X <- as.numeric(sort_X)
    if (length(sort_X) < 2L || any(!is.finite(sort_X)) ||
        is.unsorted(sort_X)) {
      stop(
        "sort_X must contain at least two finite values in nondecreasing order.",
        call. = FALSE
      )
    }
  }
  reference_n <- length(sort_X)
  m <- floor(w_trunc * reference_n)
  if (m >= (reference_n - 1L) / 2L) {
    stop("w_trunc leaves too few empirical-CDF cutpoints.", call. = FALSE)
  }

  loss <- .reliever_loss_matrix(
    nmcd_individual_loss(
      x, l, r, l_end, r_end,
      sorted_reference = sort_X,
      tail_truncation = m
    )
  )
  list(loss = loss, model = NULL)
}
