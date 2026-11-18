.default_diagnostics <- c(
  "gelman.diag", "effectiveSize", "geweke.diag"
)

.available_diagnostics <- c(
  .default_diagnostics, "heidel.diag", "raftery.diag"
)

.check_diagnostic <- function(
  result, test, gelman.threshold = 1.01, ess.threshold = 400,
  geweke.threshold = stats::qnorm(0.975), raftery.threshold = 5
) {
  switch(test,
    "gelman.diag" = {
      psrf <- result$psrf[, "Point est."]

      return(names(psrf)[!is.finite(psrf) | psrf > gelman.threshold])
    },
    "effectiveSize" = {
      ess <- unlist(result)

      return(names(ess)[!is.finite(ess) | ess < ess.threshold])
    },
    "geweke.diag" = {
      chain_results <- if (!is.null(result$z)) {
        list(result)
      } else {
        result
      }

      problems <- lapply(
        seq_along(chain_results),
        function(chain) {
          z <- chain_results[[chain]]$z

          bad <- !is.finite(z) | abs(z) > geweke.threshold

          if (!any(bad)) {
            return(NULL)
          }

          return(
            data.frame(
              parameter = names(z)[bad],
              chain = chain,
              stringsAsFactors = FALSE
            )
          )
        }
      )

      problems <- Filter(Negate(is.null), problems)

      if (length(problems) == 0L) {
        return(NULL)
      } else {
        return(unique(do.call(rbind, problems)))
      }
    },
    "heidel.diag" = {
      chain_results <- if (is.matrix(result)) {
        list(result)
      } else {
        result
      }

      problems <- lapply(
        seq_along(chain_results),
        function(chain) {
          results <- chain_results[[chain]]

          failed_stationarity <- rep(FALSE, nrow(results))
          failed_halfwidth <- rep(FALSE, nrow(results))

          if ("stest" %in% colnames(results)) {
            stationarity <- suppressWarnings(as.numeric(results[, "stest"]))
            failed_stationarity <- !is.finite(stationarity) | stationarity == 0
          }

          if ("htest" %in% colnames(results)) {
            halfwidth <- suppressWarnings(as.numeric(results[, "htest"]))
            failed_halfwidth <- !is.finite(halfwidth) | halfwidth == 0
          }

          bad <- failed_stationarity | failed_halfwidth

          if (!any(bad)) {
            return(NULL)
          }

          return(
            data.frame(
              parameter = rownames(results)[bad],
              chain = chain,
              stringsAsFactors = FALSE
            )
          )
        }
      )

      problems <- Filter(Negate(is.null), problems)

      if (length(problems) == 0L) {
        return(NULL)
      } else {
        return(unique(do.call(rbind, problems)))
      }
    },
    "raftery.diag" = {
      chain_results <- if (!is.null(result$resmatrix)) {
        list(result)
      } else {
        result
      }

      problems <- lapply(
        seq_along(chain_results),
        function(chain) {
          results <- chain_results[[chain]]$resmatrix

          dependence_column <- grep(
            "^I$|dependence", colnames(results),
            value = TRUE, ignore.case = TRUE
          )

          if (length(dependence_column) == 0L) {
            return(NULL)
          }

          dependence <- suppressWarnings(
            as.numeric(results[, dependence_column[1L]])
          )

          bad <- !is.finite(dependence) | dependence > raftery.threshold

          if (!any(bad)) {
            return(NULL)
          }

          return(
            data.frame(
              parameter = rownames(results)[bad],
              chain = chain,
              stringsAsFactors = FALSE
            )
          )
        }
      )

      problems <- Filter(Negate(is.null), problems)

      if (length(problems) == 0L) {
        return(NULL)
      } else {
        return(unique(do.call(rbind, problems)))
      }
    },
    return(character())
  )
}

.format_diagnostic_warning <- function(problematic, test, n.chains) {
  if (is.null(problematic) || length(problematic) == 0L) {
    return(NULL)
  }

  if (
    is.data.frame(problematic) &&
      all(c("parameter", "chain") %in% names(problematic))
  ) {
    problematic <- unique(
      problematic[, c("parameter", "chain"), drop = FALSE]
    )

    parameter_order <- unique(problematic$parameter)

    by_parameter <- split(
      problematic$chain,
      factor(problematic$parameter, levels = parameter_order)
    )

    by_parameter <- lapply(by_parameter, function(x) sort(unique(x)))

    details <- paste(
      sprintf(
        "%s (%s)", names(by_parameter),
        vapply(by_parameter, paste, collapse = ", ", FUN.VALUE = character(1))
      ),
      collapse = "; "
    )

    warning(
      "Diagnostic `", test, "` identified a possible convergence issue for ",
      length(by_parameter), " parameter",
      if (length(by_parameter) != 1L) "s" else "",
      " in ", length(unique(problematic$chain)), " chain",
      if (length(unique(problematic$chain)) != 1L) "s" else "",
      ": ", details, ".",
      call. = FALSE
    )

    return(invisible(NULL))
  }

  problematic <- unique(problematic)

  warning(
    "Diagnostic `", test, "` identified a possible convergence issue for ",
    length(problematic), " parameter",
    if (length(problematic) != 1L) "s" else "",
    ": ", paste(problematic, collapse = ", "), ".",
    call. = FALSE
  )

  return(invisible(NULL))
}

.run_diagnostics <- function(
  samples, do.diagnostics = .default_diagnostics, gelman.threshold = 1.01,
  ess.threshold = 400, geweke.threshold = stats::qnorm(0.975), raftery.threshold = 5
) {
  if (is.null(do.diagnostics) || identical(do.diagnostics, FALSE)) {
    warning("No diagnostic tests run. Consider checking convergence manually.")
    return(NULL)
  }

  if (identical(do.diagnostics, TRUE)) {
    do.diagnostics <- .default_diagnostics
  }

  if (!is.character(do.diagnostics)) {
    stop(
      "`do.diagnostics` must be TRUE, FALSE, NULL, or a ",
      "character vector of coda diagnostic function names.",
      call. = FALSE
    )
  }

  do.diagnostics <- unique(do.diagnostics)

  unknown <- setdiff(do.diagnostics, .available_diagnostics)

  if (length(unknown) > 0L) {
    stop(
      "Unsupported diagnostic", if (length(unknown) > 1L) "s" else "",
      ": ", paste(unknown, collapse = ", "), ". Supported diagnostics are: ",
      paste(.available_diagnostics, collapse = ", "), ".",
      call. = FALSE
    )
  }

  diagnostics <- stats::setNames(
    vector("list", length(do.diagnostics)), do.diagnostics
  )

  for (test in do.diagnostics) {
    fun <- getExportedValue("coda", test)

    result <- tryCatch(
      {
        if (identical(test, "gelman.diag")) {
          fun(samples, multivariate = FALSE)
        } else {
          fun(samples)
        }
      },
      error = function(e) {
        e
      }
    )

    if (inherits(result, "error")) {
      warning(
        "Diagnostic `", test, "` could not be calculated: ",
        conditionMessage(result),
        call. = FALSE
      )

      diagnostics[[test]] <- structure(
        list(message = conditionMessage(result)),
        class = "dodo_diagnostic_error"
      )

      next
    }

    diagnostics[[test]] <- result

    problematic <- tryCatch(
      .check_diagnostic(
        result = result,
        test = test,
        gelman.threshold = gelman.threshold,
        ess.threshold = ess.threshold,
        geweke.threshold = geweke.threshold,
        raftery.threshold = raftery.threshold
      ),
      error = function(e) {
        warning(
          "Diagnostic `",
          test,
          "` was calculated, but its result could not be checked: ",
          conditionMessage(e),
          call. = FALSE
        )

        character()
      }
    )

    if (!is.null(problematic) && length(problematic) > 0L) {
      .format_diagnostic_warning(
        problematic = problematic,
        test = test,
        n.chains = length(samples)
      )
    }
  }

  return(diagnostics)
}
