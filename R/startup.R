.onAttach <- function(libname, pkgname) {
  # Check Suggests
  suggests <- c(
    "readxl",
    "rjags"
  )

  missing_suggests <- suggests[
    !vapply(suggests, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
  ]

  if (length(missing_suggests) > 0) {
    packageStartupMessage(
      paste0(
        "Optional package(s) not installed: ",
        paste(missing_suggests, collapse = ", "),
        ". Some dodo functionality may be unavailable."
      )
    )
  }

  # Check JAGS
  jags_available <- FALSE

  if (requireNamespace("rjags", quietly = TRUE)) {
    jags_available <- tryCatch(
      {
        rjags::jags.version()
        TRUE
      },
      error = function(e) FALSE
    )
  }

  if (!jags_available) {
    packageStartupMessage(
      paste(
        "JAGS was not detected.",
        "Model fitting functions requiring JAGS will not work."
      )
    )
  }
}
