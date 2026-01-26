#' @title Convert basic sightings spreadsheet into all `dodo` formats
#'
#' @description
#' Converts a sightings spreadsheet into all the data formats used in `dodo`
#' models.
#'
#' @param x sightings spreadsheet to be converted. Must be a `data.frame` with
#' all of the requisite columns included. See Table S1 in Buchanan et al. (2025)
#' for an example.
#' @param init.time start of the observation period.
#' @param test.time end of the observation period, typically the present day
#' (defaults to the current year).
#' @param threshold cutoff certainty value for which sightings to consider
#' certain (for binary certain/uncertain models). Defaults to `0.9`.
#' @param unique whether to deduplicate the dataset (to mitigate nonidependence
#' issues). Defaults to `TRUE`, in which case the dataset is grouped by every
#' column except `time`, with the most recent row kept for each group.
#' @param aggregate what function to use when aggregating records into `iucn`
#' format. Default is `"pci_prod"`, which returns \eqn{1 - \prod(1 - x)} where
#' \eqn{x} is the vector of certainty values for each time step. Other sensible
#' functions include `"min`, `"mean"`, or `"max"`.
#' @param time name of the column with time values.
#' @param certainty name of the column with certainty values.
#' @param certainty_lower name of the column with certainty lower bound values.
#' @param certainty_upper name of the column with certainty upper bound values.
#'
#' @returns a `list` object with 11 elements: the original parameters
#' `init.time`, `test.time`, `threshold`, and `unique`, and the 7 converted
#' datasets.
#'
#' @section Data Formats:
#' **`ccon`**: *continuous certain sightings*. The simplest format, a vector of
#' all the (certain) sighting times of the taxon of interest. See e.g.
#' \code{\link{SO93F1}}.
#'
#' **`cbin`**: *binary certain sightings*. A vector of 0/1 (certain) sighting
#' states at all discrete (integer) time intervals between `init.time` and
#' `test.time` for the taxon of interest. Similar to `ccon`, but note that
#' multiple records in an interval collapse to 1. See e.g. \code{\link{CB14B1}}.
#'
#' **`cdis`**: *discrete certain sightings*. A vector of (certain) sighting
#' counts at all discrete (integer) time intervals between `init.time` and
#' `test.time` for the taxon of interest. Similar to `cbin`, but note that
#' multiple records in an interval *do not* collapse to 1. See e.g.
#' \code{\link{BU95F1}}.
#'
#' **`ucon`**: *continuous uncertain sightings*. A `data.frame` of all the
#' sighting times of the taxon of interest, with certainty scores between 0 and
#' 1 associated with each record. See e.g. \code{\link{JR14F1}}.
#'
#' **`ubin`**: *binary uncertain sightings*. A `data.frame` with two columns,
#' the first being the vector of 0/1 certain sighting states at all discrete
#' (integer) time intervals between `init.time` and `test.time` for the taxon
#' of interest, and the second being the equivalent vector for uncertain
#' sightings. See e.g. \code{\link{KO20B2}}.
#'
#' **`umcd`**: *multi-class discrete uncertain sightings*. A `data.frame` of
#' sighting counts at all discrete (integer) time intervals between `init.time`
#' and `test.time` for the taxon of interest. Sightings are split up into
#' classes based on their `certainty` score. See e.g. \code{\link{TH13B1}}.
#'
#' **`iucn`**: *IUCN-format sightings*. A `data.frame` of sighting certainty
#' (\eqn{p(ci)}) bounds (lower and upper) at all discrete (integer) time
#' intervals between `init.time` and `test.time` for the taxon of interest. If
#' there are multiple sightings in a single time interval, they are assumed
#' to be independent of one another and the overall \eqn{p(ci)} score for that
#' period is defined as \eqn{1 - \prod 1 - p(ci)}.
#' See e.g. \code{\link{TH17I1}}.
#'
#' @examples
#' # Convert the raw Slender-billed Curlew data
#' convert_dodo(
#'   x = curlew_raw, init.time = 1817, test.time = 2022, threshold = 0.9,
#'   unique = TRUE, aggregate = "pci_prod", time = "year", certainty = "p_ci",
#'   certainty_lower = "p_ci_min", certainty_upper = "p_ci_max"
#' )
#'
#' @export

convert_dodo <- function(x, init.time,
                         test.time = as.numeric(format(Sys.Date(), "%Y")),
                         threshold = 0.9, unique = TRUE, aggregate = "pci_prod",
                         time, certainty, certainty_lower, certainty_upper) {
  # Check if x is a data.frame
  if (!identical(class(x), "data.frame")) {
    stop("Sightings spreadsheet must be a data.frame object!")
  }

  # Remove any sightings before init.time or after test.time
  x <- x[x[[time]] >= init.time, ]
  x <- x[x[[time]] <= test.time, ]

  # Remove duplicates
  if (unique == TRUE) {
    x <- aggregate(x[[time]], by = x[setdiff(names(x), time)], FUN = max)
    names(x)[names(x) == "x"] <- time
    x <- x[c(time, setdiff(names(x), time))]
    x <- x[order(x[, time]), ]
  }

  # Continuous certain sightings
  x_ccon <- sort(x[x[[certainty]] >= threshold, time])

  # Binary certain sightings from init.time to test.time
  x_cbin <- as.integer(init.time:test.time %in% x_ccon)

  # Discrete certain sightings from init.time to test.time
  x_cdis <- as.integer(table(factor(x_ccon, levels = init.time:test.time)))

  # Continuous uncertain sightings from init.time to test.time
  x_ucon <- x[, c(time, certainty)]
  names(x_ucon) <- c("time", "certainty")
  x_ucon <- sort_by(x_ucon, ~time)

  # Binary uncertain sightings from init.time to test.time
  x_ubin <- data.frame(
    certain = x_cbin,
    uncertain = as.integer(init.time:test.time %in% sort(
      x[x[[certainty]] < threshold, time]
    ))
  )

  # Multi-class discrete uncertain sightings from init.time to test.time
  x_umcd <- data.frame(time = init.time:test.time)
  for (certainty_class in sort(unique(x[[certainty]]), decreasing = TRUE)) {
    x_umcd[paste0("records_", certainty_class)] <- as.integer(table(factor(
      x[x[[certainty]] == certainty_class, time],
      levels = init.time:test.time
    )))
  }

  # Helper function for IUCN records
  pci_prod <- function(x) {
    return(1 - prod(1 - x))
  }

  # IUCN format records from init.time to test.time
  x_iucn <- data.frame(time = init.time:test.time)
  x_iucn$record <- x_iucn$time %in% x[[time]]
  x_iucn <- merge(x_iucn, aggregate(
    x[, c(certainty, certainty_lower, certainty_upper)],
    by = x[time], FUN = aggregate
  ),
  by.x = "time", by.y = time, all.x = TRUE
  )
  x_iucn[is.na(x_iucn)] <- 0
  names(x_iucn) <- c(
    "time", "record", "certainty",
    "certainty_lower", "certainty_upper"
  )

  # Output
  output <- list(
    init.time = init.time,
    test.time = test.time,
    threshold = threshold,
    unique = unique,
    ccon = x_ccon,
    cbin = x_cbin,
    cdis = x_cdis,
    ucon = x_ucon,
    ubin = x_ubin,
    umcd = x_umcd,
    iucn = x_iucn
  )

  return(output)
}
