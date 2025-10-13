#' @title Solow's (2016) "Retrospective" model
#'
#' @description
#' Equation 7 from Solow 2016. Estimates a posterior probability that the
#' species is/was extant at the test time.
#'
#' @param records numeric vector object containing all sighting records of the
#' taxon of interest.
#' @param init.time start of the observation period. Defaults to the time of
#' the first sighting, in which case this sighting is removed from the record.
#' @param test.time time point to retrospectively calculate extinction
#' probability at. Must not be earlier than the time of the most recent
#' sighting nor later than `curr.time`.
#' @param curr.time end of the observation period, typically the present day
#' (defaults to the current year).
#'
#' @returns a `list` object with the original parameters and p(extant) included
#' as elements.
#'
#' @note
#' All sighting records are assumed to be certain and sampling effort is assumed
#' to be constant. Uses the uniform prior (1 / (m + T)) specified on p. 797.
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
#' @seealso [AL14B1()]; [SO16B2()]
#'
#' @examples
#' # Run the Dodo analysis from Solow 2016
#' SO16B1(dodos, test.time = 1672, curr.time = 2015)
#' # Note that Solow 2016 presents p(extinct), which is 1 - p(extant) used here.
#' # Also note that Solow 2016 rounds p to 2 d.p., which explains the slight
#' # difference between the two final p(extinct) results (0.688 vs 0.692).
#'
#' @export

SO16B1 <- function(records, init.time = min(records), test.time,
                   curr.time = as.numeric(format(Sys.Date(), "%Y"))) {
  # Sort records
  records <- sort(records)

  # Determine number of records
  n <- length(records)

  # If using first record as init.time, remove this from the record sequence
  if (init.time == min(records)) {
    records <- tail(records, -1)
    n <- n - 1
  }

  # Determine length of sighting record
  m <- max(records) - init.time

  # Determine test time
  j <- test.time - max(records)

  # Determine sighting rate
  p <- (n - 2) / (m - 1) # this is rounded to 2 d.p. in the original paper

  # Determine time since last record
  bigT <- curr.time - max(records)

  if (j > bigT | j < 0) {
    stop("Invalid test.time!")
  }

  # Calculate p(extinct)
  p.extinct <- Eq7(m, j, p, bigT)

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

#' @title Equation 7 from Solow (2016).
#'
#' @description
#' Helper function. From Equation 7 in Solow 2016.
#'
#' @param m a number.
#' @param j a number.
#' @param p a number.
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

Eq7 <- function(m, j, p, bigT) {
  numerator <- sum((1 - p)^((1:j) - 1) * (1 / (bigT + m)))
  denominator_a <- sum((1 - p)^((1:bigT) - 1) * (1 / (bigT + m)))
  denominator_b <- (1 - p)^bigT
  denominator_c <- 1 - (j / (bigT + m))

  return(numerator / (denominator_a + denominator_b * denominator_c))
}
