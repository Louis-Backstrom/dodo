#' Thompson et al.'s (2017) "Iterative" model
#'
#' @description
#' The records and surveys model from Thompson et al. 2017. Estimates a
#' probability that the species is extant at the end of the observation period.
#'
#' @param records sighting records in `iucn` format (see
#' \code{\link{convert_dodo}} for details).
#' @param surveys survey effort data (incorporating both dedicated and passive
#' survey effort). A `data.frame` with 11 columns: `time`, `survey`, `epsilon`,
#' `epsilon_lower`, `epsilon_upper`, `p_i`, `p_i_lower`, `p_i_upper`, `p_r`,
#' `p_r_lower`, and `p_r_upper`. The number of rows must match `records`.
#'
#' @returns a `list` object with the original parameters and the iterative and
#  final estimates of p(extant) (\eqn{P(X_t)}) included as elements.
#'
#' @note
#' This model incorporates both sighting uncertainty and variable survey effort.
#'
#' @references
#' **Key Reference**
#'
#' Thompson, C. J., Koshkina, V., Burgman, M. A., Butchart, S. H. M., & Stone,
#' L. (2017). Inferring extinctions II: A practical, iterative model based on
#' records and surveys. *Biological Conservation*, 214, 328-335.
#' \doi{10.1016/j.biocon.2017.07.029}
#'
#' @examples
#' # Run the Alaotra Grebe analysis from Thompson et al. 2017
#' TH17I1(records = grebe, surveys = grebe_surveys)
#'
#' @export

TH17I1 <- function(records, surveys) {
  # Set up merged records and surveys data.frame
  data <- merge(records, surveys, by = "time")
  data$PXt <- NA
  data$PXt_lower <- NA
  data$PXt_upper <- NA

  PXt <- 1
  PXt_lower <- 1
  PXt_upper <- 1

  # Iteratively calculate PXt values
  for (i in 1:nrow(data)) {
    if (data[i, ]$record == TRUE) {
      PXt1 <- PXt + data[i, ]$record *
        data[i, ]$certainty * (1 - PXt)
      PXt1_lower <- PXt_lower + data[i, ]$record *
        data[i, ]$certainty_lower * (1 - PXt_lower)
      PXt1_upper <- PXt_upper + data[i, ]$record *
        data[i, ]$certainty_upper * (1 - PXt_upper)
    } else {
      PXt1 <- (1 - data[i, ]$epsilon * data[i, ]$p_i *
        data[i, ]$p_r) * PXt
      PXt1_lower <- (1 - data[i, ]$epsilon_lower * data[i, ]$p_i_lower *
        data[i, ]$p_r_lower) * PXt_lower
      PXt1_upper <- (1 - data[i, ]$epsilon_upper * data[i, ]$p_i_upper *
        data[i, ]$p_r_upper) * PXt_upper
    }

    data[i, ]$PXt <- PXt1
    data[i, ]$PXt_lower <- PXt1_lower
    data[i, ]$PXt_upper <- PXt1_upper

    PXt_lower <- PXt1_lower
    PXt_upper <- PXt1_upper
    PXt <- PXt1
  }

  # Swap lower and upper names
  n <- ncol(data)
  temp <- data[[n - 1]]
  data[[n - 1]] <- data[[n]]
  data[[n]] <- temp

  # Output
  output <- list(
    records = records,
    surveys = surveys,
    estimates = data[, c("time", "PXt", "PXt_lower", "PXt_upper")],
    p.extant = PXt
  )

  return(output)
}
