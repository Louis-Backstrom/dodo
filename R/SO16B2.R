#' @title Solow's (2016) "Sequential" model
#'
#' @description
#' Equation 9 from Solow 2016. Estimates a posterior probability that the
#' species is extant.
#'
#' @param records sighting records in `cbin` format (see
#' \code{\link{convert_dodo}} for details).
#' @param init.time start of the observation period. Defaults to `NULL`; this
#' parameter is only necessary if specifying a different `test.time` to the end
#' of the observation period.
#' @param test.time time point to retrospectively calculate extinction
#' probability at. Defaults to `NULL`, in which case the probability is
#' estimated for the end of the observation period.
#' @param curr.time end of the observation period. Defaults to `NULL`; this
#' parameter is only necessary if specifying a different `test.time` to the end
#' of the observation period. If specified, `curr.time` - `init.time` must equal
#' `length(records)`.
#'
#' @returns a `list` object with the original parameters and p(extant) included
#' as elements.
#'
#' @note
#' All sighting records are assumed to be certain and sampling effort is assumed
#' to be constant. Uses the uniform prior (1 / (m + T)) specified on p. 797 of
#' Solow (2016).
#'
#' @references
#' **Key Reference**
#'
#' Solow, A. R. (2016). A simple Bayesian method of inferring extinction:
#' comment. *Ecology*, 97(3), 796-798. \doi{10.1890/15-0336.1}
#'
#' **Other References**
#'
#' Alroy, J. (2014). A simple Bayesian method of inferring extinction.
#' *Paleobiology*, 40(4), 584-607. \doi{10.1666/13074}
#'
#' Alroy, J. (2016). A simple Bayesian method of inferring extinction: reply.
#' *Ecology*, 97(3), 798-800. \doi{10.1002/ecy.1321}
#'
#' @seealso [AL14B1()]; [SO16B1()]
#'
#' @examples
#' # Run the Dodo analysis from Solow 2016
#' SO16B2(
#'   records = as.integer((min(dodos) + 1):2015 %in% dodos),
#'   init.time = min(dodos), test.time = 1672, curr.time = 2015
#' )
#' # NB: Solow 2016 presents p(extinct), which is 1 - p(extant) used here.
#' \dontrun{
#' # Run an example analysis using the Slender-billed Curlew data
#' SO16B2(curlew$cbin)
#' }
#'
#' @export

SO16B2 <- function(records, init.time = NULL, test.time = NULL,
                   curr.time = NULL) {
  # Determine number of records
  n <- sum(records)

  # Determine length of sighting record
  m <- max(which(records != 0))

  # If init.time, test.time, and curr.time are specified, check they are valid
  if (!is.null(init.time) | !is.null(test.time) | !is.null(curr.time)) {
    if (curr.time - init.time != length(records)) {
      stop("curr.time - init.time != length(records)")
    }
  }

  # Determine time since last record
  r <- rle(rev(records))
  if (r$values[1] == 0) {
    bigT <- r$lengths[1]
  } else {
    bigT <- 0
  }
  rm(r)

  # Determine test time
  if (!is.null(test.time)) {
    j <- test.time - (init.time + m)
  } else {
    j <- bigT
  }

  if (j > bigT | j < 0) {
    stop("Invalid test.time!")
  }

  # Calculate p(extinct)
  p.extinct <- Eq9(m, j, n, bigT)

  # Output
  output <- list(
    records = records,
    init.time = init.time,
    test.time = test.time,
    curr.time = curr.time,
    p.extant = 1 - p.extinct
  )

  return(output)
}

#' @title Equation 9 from Solow (2016).
#'
#' @description
#' Helper function. From Equation 7 in Solow 2016.
#'
#' @param m a number.
#' @param j a number.
#' @param n a number.
#' @param bigT a number
#'
#' @returns a number.
#'
#' @references
#' **Key Reference**
#'
#' Solow, A. R. (2016). A simple Bayesian method of inferring extinction:
#' comment. *Ecology*, 97(3), 796-798. \doi{10.1890/15-0336.1}
#'
#' @noRd

Eq9 <- function(m, j, n, bigT) {
  numerator <- sum(choose(m - 1 + 1:j, n)^-1 * 1 / (bigT + m))
  denominator_a <- numerator
  denominator_b <- choose(m + j, n)^-1
  denominator_c <- 1 - (j / (bigT + m))

  return(numerator / (denominator_a + denominator_b * denominator_c))
}
