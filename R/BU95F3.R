#' @title Burgman et al.'s (1995) "Empty Cells" Model
#'
#' @description
#' Equation 4 from Burgman et al. 1995. Estimates a p-value for
#' testing competing hypotheses of extinction/non-extinction.
#'
#' @param records `data.frame` with two columns: `time` and `records`. The
#' `time` column must extend from the start of the observation period (which
#' may be prior to the first sighting) to the end (typically the present day),
#' with evenly-spaced temporal intervals (typically years).
#'
#' @returns a `list` object with the original parameters and the p-value
#' included as elements.
#'
#' @note
#' All sighting records are assumed to be certain and sampling effort is assumed
#' to be constant.
#'
#' Uses code from the `sExtinct` package
#' (\href{https://github.com/ExperimentalConservation/sExtinct}{GitHub})
#'
#' @references
#' **Key Reference**
#'
#' Burgman, M. A., Grimson, R. C., & Ferson, S. (1995). Inferring Threat from
#' Scientific Collections. *Conservation Biology*, 9(4), 923-928.
#' \doi{0.1046/j.1523-1739.1995.09040923.x}
#'
#' **Other References**
#'
#' Clements, C. F., Collen, B., Blackburn, T. M., & Petchey, O. L. (2014).
#' Effects of recent environmental change on accuracy of inferences of
#' extinction status. *Conservation Biology*, 28(4), 971-981.
#' \doi{10.1111/cobi.12329}
#'
#' Grimson, R. C., Aldrich, T. E., & Drane, J. W. (1992). Clustering in sparse
#' data and an analysis of rhabdomyosarcoma incidence.
#' *Statistics in Medicine*, 11(6), 761-768. \doi{10.1002/sim.4780110607}
#'
#' @export

BU95F3 <- function(records) {

  # Adapted from sExtinct (by Christopher Clements) package code!

  # Determine number of records
  N <- sum(records$records)

  # Determine the total number of sighting intervals
  CT <- nrow(records)

  # Determine the length of the longest run of empty cells
  r <- max(rle(records$records == 0)$lengths[which(
    rle(records$records == 0)$values == TRUE)])

  # Calculate sum in for loop (from sExtinct code)
  full_sum <- 0

  for (j in 1:N) {

    TT <- seq(1, j)
    KK <- seq(1, j + 1)

    for (k in 1:(j + 1)) {

      if (k <= (CT / r)) {

        ff <- 1

        for(indexn in 1:length(TT)) {
          n <- indexn - 1
          ff <-  ff * (CT - (r * k) - n)
        }

        dummy <- 0

        for(indexi in 1:length(KK)) {
          i <- indexi - 1
          dummy <- dummy + (((-1) ^ i) * choose(j, i) * ((j - i) ^ N))
        }

        Sn <-  (1 / factorial(j)) * dummy
        calc <- ((-1) ^ (k+1)) * choose(j + 1, k) * ff * Sn

      } else {calc <- 0}

      full_sum <- full_sum + calc
    }
  }

  # Calculate p-value
  p.value <- CT ^ -N * full_sum

  # Output
  output <- list(
    records = records,
    p.value = p.value
  )

  return(output)

}
