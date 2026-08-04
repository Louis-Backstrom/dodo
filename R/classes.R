dodo_fit <- function(x, model = NULL) {
  if (!is.list(x)) {
    stop("x must be a list")
  }

  if (!is.null(model)) {
    x$model <- model
    class(x) <- c(paste0(model, "_fit"), "dodo_fit", "list")
  } else {
    class(x) <- c("dodo_fit", "list")
  }

  return(x)
}

#' @export
print.dodo_fit <- function(x, ...) {
  cat("Extinction model fit\n")
  cat("--------------------\n")

  # Print model name
  if (!is.null(x$model)) {
    cat("Model:        ", x$model, "\n", sep = "")
  }

  # Print key model outputs
  model_outputs <- c(
    "p.extant",
    "p.value",
    "estimate",
    "cred.int",
    "conf.int",
    "Bayes.factor"
  )

  model_labels <- c(
    p.extant = "p(Extant)",
    p.value = "p-value",
    estimate = "Estimate",
    cred.int = paste0(100 * (1 - x$alpha), "% CI"),
    conf.int = paste0(100 * (1 - x$alpha), "% CI"),
    Bayes.factor = "Bayes Factor"
  )

  for (output in model_outputs) {
    if (!is.null(x[[output]])) {
      val <- x[[output]]

      if (is.numeric(val) && length(val) == 1) {
        # Format non-interval values
        val <- paste0(format(val, digits = 5))
      } else if (is.numeric(val) && length(val) == 2) {
        # Format standard interval values
        val <- paste0("[", paste(format(val, digits = 5), collapse = ", "), "]")
      } else if (is.numeric(val) && length(val) == 3) {
        # Format non-standard interval values
        val <- paste0(
          "[", format(val[1], digits = 5), ", [",
          paste(format(val[2:3], digits = 5), collapse = ", "), "]]"
        )
      }

      label <- model_labels[[output]]

      cat(sprintf("%-13s %s\n", paste0(label, ":"), val))
    }
  }

  invisible(x)
}
