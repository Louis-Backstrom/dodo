#' @title Lee's (2014) "Simple" model
#'
#' @description
#' The simple numerical tool presented in Lee 2014, including both uncertain
#' sightings and surveys. Estimates a posterior probability that the species is
#' extant at the end of the observation period.
#'
#' @param records sighting records in `iucn` format (see
#' \code{\link{convert_dodo}} for details).
#' @param surveys survey effort data (passive survey effort is ignored). A
#' `data.frame` with 11 columns: `time`, `survey`, `epsilon`, `epsilon_lower`,
#' `epsilon_upper`, `p_i`, `p_i_lower`, `p_i_upper`, `p_r`, `p_r_lower`, and
#' `p_r_upper`. The number of rows must match `records`.
#' @param threshold cutoff certainty value for which sightings to consider
#' certain. Defaults to `0.9`.
#' @param prior a two-element `numeric` vector with the lower and upper bounds
#' of the (uniform) prior for extant (P(X_T)). Defaults to `c(0, 1)`.
#' @param n.iter number of iterations to run (defaults to 10000).
#'
#' @returns a `list` object with the original parameters and p(extant) included
#' as elements.
#'
#' @note
#' Passive survey effort (rows where `survey == FALSE`) are ignored, unlike in
#' [TH17I1()].
#'
#' @references
#' **Key Reference**
#'
#' Lee, T. E. (2014). A simple numerical tool to infer whether a species is
#' extinct. *Methods in Ecology and Evolution*, 5(8), 791-796.
#' \doi{10.1111/2041-210x.12227}
#'
#' @examples
#' # Run the Pohnpei Mountain Starling analysis from Lee 2014.
#' LE14B1(
#'   records = starling, surveys = starling_surveys, threshold = 0.9,
#'   prior = c(0.24, 0.48), n.iter = 1e4
#' )
#' # Run an example analysis using the Alaotra Grebe data
#' \dontrun{
#' LE14B1(
#'   records = grebe, surveys = grebe_surveys, threshold = 0.9,
#'   prior = c(0, 1), n.iter = 1e4
#' )
#' }
#'
#' @export

LE14B1 <- function(records, surveys, threshold = 0.9, prior = c(0, 1),
                   n.iter = 1e4) {
  # NB terms: W is certain sightings; U is uncertain sightings; S is surveys

  # Binary sighting history of certain sightings
  W <- as.integer(records$certainty >= threshold)

  # Key model parameters
  bigT <- length(W)
  tn <- max(which(W == 1))
  S <- sum(W)
  N <- bigT - tn

  # True positive detectability bounds for certain sightings
  pWL <- 0
  pWU <- 2 * S / tn

  # False positive detectability bounds for certain sightings
  qWL <- 0
  qWU <- 0

  # Prior bounds for extant (P(X_T))
  pieL <- min(prior)
  pieU <- max(prior)

  # Prior bounds for extinct in a given year (P(E))
  PEEL <- (1 - pieU) / (bigT - tn)
  PEEU <- (1 - pieL) / (bigT - tn)

  # Number of uncertain sighting types/qualities
  UQ <- records
  UQ_unique <- UQ[UQ$record == TRUE & UQ$certainty < threshold, ]

  # Create uncertain sighting quality list
  US <- list()
  for (h in 1:nrow(UQ_unique)) {
    US[[paste0("U", h)]] <- ifelse(UQ$time %in% merge(
      UQ_unique[h, ], UQ
    )$time, 1, 0)
    US[[paste0("SumU", h)]] <- sum(US[[paste0("U", h)]])
    US[[paste0("pU", h, "L")]] <- (S / tn) / UQ_unique[h, ]$certainty_upper
    US[[paste0("pU", h, "U")]] <- (S / tn) / UQ_unique[h, ]$certainty_lower
    US[[paste0("qU", h, "L")]] <- US[[paste0("pU", h, "L")]] - (S / tn)
    US[[paste0("qU", h, "U")]] <- US[[paste0("pU", h, "U")]] - (S / tn)
  }

  # Survey quality bounds
  SQ <- surveys
  SQ$pSL <- SQ$epsilon_lower * SQ$p_i_lower * SQ$p_r_lower
  SQ$pSU <- SQ$epsilon_upper * SQ$p_i_upper * SQ$p_r_upper
  SQ$qSL <- 0
  SQ$qSU <- SQ$pSL

  # Number of survey types/qualities (incidental survey effort is dropped)
  SQ_unique <- unique(SQ[SQ$survey == TRUE, c("pSL", "pSU", "qSL", "qSU")])

  # Create survey quality list
  SS <- list()
  for (i in 1:nrow(SQ_unique)) {
    SS[[paste0("S", i)]] <- ifelse(SQ$time %in% merge(
      SQ_unique[i, ], SQ
    )$time, 0, NA)
    SS[[paste0("pS", i, "L")]] <- SQ_unique$pSL[i]
    SS[[paste0("pS", i, "U")]] <- SQ_unique$pSU[i]
    SS[[paste0("qS", i, "L")]] <- SQ_unique$qSL[i]
    SS[[paste0("qS", i, "U")]] <- SQ_unique$qSU[i]
  }

  Q <- numeric(n.iter)

  for (k in 1:n.iter) {
    # Random samples
    pie <- runif(1, pieL, pieU)
    PEE <- runif(1, min(c(PEEL, PEEU)), max(c(PEEL, PEEU)))
    pWk <- runif(1, pWL, pWU)
    pW <- ifelse(W == 0, 1 - pWk, pWk)
    qW <- ifelse(W == 0, 1, 0)

    for (i in 1:nrow(SQ_unique)) {
      SS[[paste0("pS", i, "k")]] <- runif(
        1, SS[[paste0("pS", i, "L")]],
        SS[[paste0("pS", i, "U")]]
      )
      SS[[paste0("qS", i, "k")]] <- runif(
        1, SS[[paste0("qS", i, "L")]],
        SS[[paste0("qS", i, "U")]]
      )
      SS[[paste0("pS", i)]] <- ifelse(is.na(SS[[paste0("S", i)]]), 1,
        ifelse(SS[[paste0("S", i)]] == 0,
          1 - SS[[paste0("pS", i, "k")]],
          SS[[paste0("pS", i, "k")]]
        )
      )
      SS[[paste0("qS", i)]] <- ifelse(is.na(SS[[paste0("S", i)]]), 1,
        ifelse(SS[[paste0("S", i)]] == 0,
          1 - SS[[paste0("qS", i, "k")]],
          SS[[paste0("qS", i, "k")]]
        )
      )
    }

    for (h in 1:nrow(UQ_unique)) {
      US[[paste0("pU", h, "k")]] <- runif(1, US[[paste0("pU", h, "L")]],
                                          US[[paste0("pU", h, "U")]])
      US[[paste0("qU", h, "k")]] <- runif(1, US[[paste0("qU", h, "L")]],
                                          US[[paste0("qU", h, "U")]])
      if (US[[paste0("SumU", h)]] > 0) {
        US[[paste0("pU", h)]] <- ifelse(US[[paste0("U", h)]] == 0,
                                        1 - US[[paste0("pU", h, "k")]],
                                        US[[paste0("pU", h, "k")]])
        US[[paste0("qU", h)]] <- ifelse(US[[paste0("U", h)]] == 0,
                                        1 - US[[paste0("qU", h, "k")]],
                                        US[[paste0("qU", h, "k")]])
      } else {
        US[[paste0("pU", h)]] <- rep(1, bigT)
        US[[paste0("qU", h)]] <- rep(1, bigT)
      }
    }

    # Multiply probabilities
    pvec <- pW
    qvec <- qW
    for (h in 1:nrow(UQ_unique)) {
      pvec <- pvec * US[[paste0("pU", h)]]
      qvec <- qvec * US[[paste0("qU", h)]]
    }
    for (i in 1:nrow(SQ_unique)) {
      pvec <- pvec * SS[[paste0("pS", i)]]
      qvec <- qvec * SS[[paste0("qS", i)]]
    }

    PX <- prod(pvec) * pie

    # PEj calculation
    PEj <- 0
    for (j in (tn + 1):bigT) {
      partpPEj <- prod(pvec[1:(j - 1)], na.rm = TRUE)
      partqPEj <- prod(qvec[j:bigT], na.rm = TRUE)
      PEj <- PEj + (partpPEj * partqPEj * PEE)
    }

    Q[k] <- PX / (PX + PEj)
  }

  # Output
  output <- list(
    records = records,
    surveys = surveys,
    threshold = threshold,
    prior = prior,
    n.iter = n.iter,
    p.extant = mean(Q)
  )

  return(output)
}
