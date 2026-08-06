.documentation_source_root <- function(start = getwd()) {
  root <- normalizePath(start, winslash = "/", mustWork = TRUE)
  required <- c(
    "DESCRIPTION", "NAMESPACE", "R", "man", file.path("inst", "CITATION")
  )
  repeat {
    if (all(file.exists(file.path(root, required)))) {
      return(root)
    }
    parent <- dirname(root)
    if (identical(parent, root)) {
      return(NULL)
    }
    root <- parent
  }
}

.readme_r_blocks <- function(path) {
  lines <- readLines(path, warn = FALSE)
  start_ids <- grep("^```[rR][[:space:]]*$", lines)
  closing_ids <- grep("^```[[:space:]]*$", lines)
  end_ids <- vapply(start_ids, function(start_id) {
    following <- closing_ids[closing_ids > start_id]
    if (length(following) == 0L) {
      stop("README contains an unclosed R code block.", call. = FALSE)
    }
    following[[1L]]
  }, integer(1L))
  Map(function(start_id, end_id) {
    lines[seq.int(start_id + 1L, end_id - 1L)]
  }, start_ids, end_ids)
}

.rd_examples_text <- function(path) {
  parsed <- tools::parse_Rd(path)
  tags <- vapply(parsed, function(node) {
    tag <- attr(node, "Rd_tag")
    if (is.null(tag)) "" else tag
  }, character(1L))
  paste(
    unlist(parsed[tags == "\\examples"], recursive = TRUE, use.names = FALSE),
    collapse = ""
  )
}

.focused_wrapper_topics <- function() {
  paste0(
    "reliever_",
    c(
      "mean", "var", "meanvar", "lm", "glm", "em",
      "lasso", "lasso_crossfit", "mean_crossfit",
      "kde_nll", "kde_nll_crossfit", "kde_l2", "nmcd",
      "ranger_crossfit", "mlp_crossfit"
    ),
    ".Rd"
  )
}

.focused_wrapper_call_pattern <- function() {
  wrappers <- sub("\\.Rd$", "", .focused_wrapper_topics())
  paste0("\\b(", paste(wrappers, collapse = "|"), ")[[:space:]]*\\(")
}

test_that("public documentation has complete topics and usable examples", {
  with_test_timeout({
    root <- .documentation_source_root()
    skip_if(is.null(root), "source documentation is unavailable")

    namespace <- readLines(file.path(root, "NAMESPACE"), warn = FALSE)
    exports <- sub(
      "^export\\((.*)\\)$", "\\1",
      grep("^export\\(", namespace, value = TRUE)
    )
    rd_files <- list.files(
      file.path(root, "man"), pattern = "\\.Rd$", full.names = TRUE
    )
    rd_source <- setNames(vapply(rd_files, function(path) {
      paste(readLines(path, warn = FALSE), collapse = "\n")
    }, character(1L)), basename(rd_files))

    topics <- lapply(exports, function(export) {
      alias <- paste0("\\alias{", export, "}")
      names(rd_source)[vapply(
        rd_source,
        function(source) grepl(alias, source, fixed = TRUE),
        logical(1L)
      )]
    })
    names(topics) <- exports
    undocumented <- names(topics)[lengths(topics) == 0L]
    expect_true(
      length(undocumented) == 0L,
      info = paste(undocumented, collapse = ", ")
    )

    has_examples <- vapply(topics, function(topic) {
      any(grepl("\\examples{", rd_source[topic], fixed = TRUE))
    }, logical(1L))
    expect_true(
      all(has_examples),
      info = paste(names(has_examples)[!has_examples], collapse = ", ")
    )

    rd_examples <- vapply(rd_files, .rd_examples_text, character(1L))
    k_selection_pattern <- paste0(
      "\\b(cv\\.reliever(?:_generic)?|select_by_run|select_across_runs|",
      "select_holdout)[[:space:]]*\\("
    )
    k_selection_topics <- basename(rd_files)[
      grepl(k_selection_pattern, rd_examples, perl = TRUE)
    ]
    stopifnot_topics <- basename(rd_files)[
      grepl("\\bstopifnot[[:space:]]*\\(", rd_examples)
    ]
    expect_true(
      all(k_selection_topics %in% stopifnot_topics),
      info = paste(setdiff(k_selection_topics, stopifnot_topics))
    )
    expect_false(any(grepl("\\dontrun{", rd_source, fixed = TRUE)))
    expect_setequal(
      names(rd_source)[grepl("\\donttest{", rd_source, fixed = TRUE)],
      c(
        "cv.reliever.Rd", "evaluate_reliever_segments.Rd", "reliever.Rd",
        "reliever_kde_l2.Rd", "reliever_kde_nll_crossfit.Rd",
        "reliever_lasso.Rd", "reliever_lasso_crossfit.Rd",
        "reliever_mlp_crossfit.Rd", "reliever_ranger_crossfit.Rd",
        "select_across_runs.Rd", "select_holdout.Rd"
      )
    )

    package_doc <- rd_source[["relieverChangepoint-package.Rd"]]
    package_title <- read.dcf(
      file.path(root, "DESCRIPTION"), fields = "Title"
    )[1L, 1L]
    expect_match(package_doc, package_title, fixed = TRUE)
    expect_match(package_doc, "\\keyword{package}", fixed = TRUE)
    expect_false(grepl("\\keyword{internal}", package_doc, fixed = TRUE))
    primary_links <- c(
      "reliever", "reliever_generic", "cv.reliever",
      "select_by_run", "select_across_runs", "select_holdout"
    )
    missing_links <- primary_links[!vapply(primary_links, function(topic) {
      grepl(paste0("\\link{", topic, "}"), package_doc, fixed = TRUE)
    }, logical(1L))]
    expect_true(length(missing_links) == 0L, info = paste(missing_links))

    focused_topics <- .focused_wrapper_topics()
    expect_true(all(focused_topics %in% names(rd_source)))
    focused_help <- rd_source[focused_topics]
    expect_true(all(grepl("\\code{reliever(", focused_help, fixed = TRUE)))
    expect_true(all(grepl('cpn_crit = "none"', focused_help, fixed = TRUE)))

    ordinary_examples <- rd_examples[
      !basename(rd_files) %in% focused_topics
    ]
    unexpected_calls <- basename(rd_files)[
      !basename(rd_files) %in% focused_topics
    ][grepl(.focused_wrapper_call_pattern(), ordinary_examples)]
    expect_true(
      length(unexpected_calls) == 0L,
      info = paste(unexpected_calls)
    )

    stale_topics <- c(
      "cv.reliever_mean.Rd", "reliever_result.Rd", "cv_reliever_result.Rd",
      "reg_fun_contract.Rd", "reg_fun_crossfit_contract.Rd",
      "reselect_by_run.Rd", "reselect_across_runs.Rd", "reselect_holdout.Rd"
    )
    expect_false(any(file.exists(file.path(root, "man", stale_topics))))
    expect_false(any(grepl("reselect_", rd_source, fixed = TRUE)))
  })
})

test_that("README R examples are syntactically valid", {
  with_test_timeout({
    root <- .documentation_source_root()
    skip_if(is.null(root), "source documentation is unavailable")

    readme <- file.path(root, "README.md")
    readme_text <- paste(readLines(readme, warn = FALSE), collapse = "\n")
    blocks <- .readme_r_blocks(readme)
    expect_gt(length(blocks), 0L)
    expect_false(any(vapply(blocks, function(block) {
      grepl("\\bstopifnot[[:space:]]*\\(", paste(block, collapse = "\n"))
    }, logical(1L))))
    expect_false(grepl(.focused_wrapper_call_pattern(), readme_text))
    expect_false(grepl("reselect_", readme_text, fixed = TRUE))

    primary_calls <- c(
      "reliever(", "reliever_generic(", "cv.reliever(",
      "cv.reliever_generic(",
      "select_by_run(", "select_across_runs("
    )
    missing_calls <- primary_calls[!vapply(primary_calls, function(call) {
      grepl(call, readme_text, fixed = TRUE)
    }, logical(1L))]
    expect_true(length(missing_calls) == 0L, info = paste(missing_calls))
    expect_true(grepl(
      'plot(fit_lasso_cf, run_type = "recv")', readme_text, fixed = TRUE
    ))
    expect_true(grepl(
      'run_type = "crossfit_homo_hyper"', readme_text, fixed = TRUE
    ))

    parse_errors <- vapply(seq_along(blocks), function(i) {
      tryCatch(
        {
          parse(text = blocks[[i]])
          ""
        },
        error = function(error) {
          sprintf("README R block %d: %s", i, conditionMessage(error))
        }
      )
    }, character(1L))
    expect_false(
      any(nzchar(parse_errors)),
      info = paste(parse_errors[nzchar(parse_errors)], collapse = "\n")
    )
  })
})

test_that("publication citations match package metadata", {
  with_test_timeout({
    root <- .documentation_source_root()
    checking_installed <- is.null(root)
    if (checking_installed) {
      citations <- citation(package = "relieverChangepoint")
    } else {
      meta <- as.list(read.dcf(file.path(root, "DESCRIPTION"))[1L, ])
      citations <- utils::readCitationFile(
        file.path(root, "inst", "CITATION"), meta = meta
      )
    }

    expect_length(citations, 2L)
    expect_equal(citations[[1L]]$number, "203")
    expect_equal(citations[[1L]]$pages, "1--57")
    expect_equal(
      as.character(citations[[2L]]$author),
      c("Chengde Qian", "Guanghui Wang", "Zhaojun Wang", "Changliang Zou")
    )

    if (!checking_installed) {
      references <- paste(
        readLines(
          file.path(root, "man", "relieverChangepoint-package.Rd"),
          warn = FALSE
        ),
        collapse = "\n"
      )
      expect_match(references, "26(203), 1--57", fixed = TRUE)
      expect_match(
        references,
        "Qian, C., Wang, G., Wang, Z., and Zou, C.",
        fixed = TRUE
      )
    }
  })
})
