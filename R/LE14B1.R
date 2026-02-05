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
#' @param n.iter number of iterations to run (defaults to 10,000).
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
#' \dontrun{
#' # Run an example analysis using the Alaotra Grebe data
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
  pWU <- 2 * S / tn # can exceed 1!
  if (pWU > 1) {
    warning("S / tn > 0.5, setting upper detectability bound at 1")
    pWU <- 1
  }

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
  UQ_unique <- unique(UQ[
    UQ$record == TRUE & UQ$certainty < threshold,
    !names(UQ) %in% c("record", "time")
  ])

  if (nrow(UQ_unique) > 0) {
    # Precompute names
    pU_L_names <- sprintf("pU%dL", 1:nrow(UQ_unique))
    pU_U_names <- sprintf("pU%dU", 1:nrow(UQ_unique))
    qU_L_names <- sprintf("qU%dL", 1:nrow(UQ_unique))
    qU_U_names <- sprintf("qU%dU", 1:nrow(UQ_unique))
    U_names <- sprintf("U%d", 1:nrow(UQ_unique))
    SumU_names <- sprintf("SumU%d", 1:nrow(UQ_unique))
    pU_k_names <- sprintf("pU%dk", 1:nrow(UQ_unique))
    qU_k_names <- sprintf("qU%dk", 1:nrow(UQ_unique))

    # Create uncertain sighting quality list; NB: pU and qU can exceed 1!
    US <- list()
    for (h in 1:nrow(UQ_unique)) {
      US[[U_names[h]]] <- ifelse(UQ$time %in% merge(
        UQ_unique[h, ], UQ
      )$time, 1, 0)
      US[[SumU_names[h]]] <- sum(US[[U_names[h]]])
      US[[pU_L_names[h]]] <- (S / tn) / UQ_unique[h, ]$certainty_upper
      US[[pU_U_names[h]]] <- (S / tn) / UQ_unique[h, ]$certainty_lower
      US[[qU_L_names[h]]] <- US[[pU_L_names[h]]] - (S / tn)
      US[[qU_U_names[h]]] <- US[[pU_U_names[h]]] - (S / tn)
    }

    if (any(US[grepl("pU|qU", names(US))] > 1) == TRUE) {
      stop("At least one uncertain sighting detectability bound > 1!")
    }
  }

  # Survey quality bounds
  SQ <- surveys
  SQ$pSL <- SQ$epsilon_lower * SQ$p_i_lower * SQ$p_r_lower
  SQ$pSU <- SQ$epsilon_upper * SQ$p_i_upper * SQ$p_r_upper
  SQ$qSL <- 0
  SQ$qSU <- SQ$pSL

  # Number of survey types/qualities (incidental survey effort is dropped)
  SQ_unique <- unique(SQ[SQ$survey == TRUE, c("pSL", "pSU", "qSL", "qSU")])

  if (nrow(SQ_unique) > 0) {
    # Precompute names
    S_names <- sprintf("S%d", 1:nrow(SQ_unique))
    pS_L_names <- sprintf("pS%dL", 1:nrow(SQ_unique))
    pS_U_names <- sprintf("pS%dU", 1:nrow(SQ_unique))
    qS_L_names <- sprintf("qS%dL", 1:nrow(SQ_unique))
    qS_U_names <- sprintf("qS%dU", 1:nrow(SQ_unique))
    pS_k_names <- sprintf("pS%dk", 1:nrow(SQ_unique))
    qS_k_names <- sprintf("qS%dk", 1:nrow(SQ_unique))

    # Create survey quality list
    SS <- list()
    for (i in 1:nrow(SQ_unique)) {
      SS[[S_names[i]]] <- ifelse(SQ$time %in% merge(
        SQ_unique[i, ], SQ
      )$time, 0, NA)
      SS[[pS_L_names[i]]] <- SQ_unique$pSL[i]
      SS[[pS_U_names[i]]] <- SQ_unique$pSU[i]
      SS[[qS_L_names[i]]] <- SQ_unique$qSL[i]
      SS[[qS_U_names[i]]] <- SQ_unique$qSU[i]
    }
  }

  Q <- numeric(n.iter)

  for (k in 1:n.iter) {
    # Random samples
    pie <- runif(1, pieL, pieU)
    PEE <- runif(1, min(c(PEEL, PEEU)), max(c(PEEL, PEEU)))
    pWk <- runif(1, pWL, pWU)
    pW <- ifelse(W == 0, 1 - pWk, pWk)
    qW <- ifelse(W == 0, 1, 0)

    if (nrow(SQ_unique) > 0) {
      for (i in 1:nrow(SQ_unique)) {
        SS[[pS_k_names[i]]] <- runif(
          1, SS[[pS_L_names[i]]],
          SS[[pS_U_names[i]]]
        )
        SS[[qS_k_names[i]]] <- runif(
          1, SS[[qS_L_names[i]]],
          SS[[qS_U_names[i]]]
        )
        SS[[sprintf("pS%d", i)]] <- ifelse(is.na(SS[[S_names[i]]]), 1,
          ifelse(SS[[S_names[i]]] == 0,
            1 - SS[[pS_k_names[i]]],
            SS[[pS_k_names[i]]]
          )
        )
        SS[[sprintf("qS%d", i)]] <- ifelse(is.na(SS[[S_names[i]]]), 1,
          ifelse(SS[[S_names[i]]] == 0,
            1 - SS[[qS_k_names[i]]],
            SS[[qS_k_names[i]]]
          )
        )
      }
    }

    if (nrow(UQ_unique) > 0) {
      for (h in 1:nrow(UQ_unique)) {
        US[[pU_k_names[h]]] <- runif(
          1, US[[pU_L_names[h]]],
          US[[pU_U_names[h]]]
        )
        US[[qU_k_names[h]]] <- runif(
          1, US[[qU_L_names[h]]],
          US[[qU_U_names[h]]]
        )
        if (US[[SumU_names[h]]] > 0) {
          US[[sprintf("pU%d", h)]] <- ifelse(US[[U_names[h]]] == 0,
            1 - US[[pU_k_names[h]]],
            US[[pU_k_names[h]]]
          )
          US[[sprintf("qU%d", h)]] <- ifelse(US[[U_names[h]]] == 0,
            1 - US[[qU_k_names[h]]],
            US[[qU_k_names[h]]]
          )
        } else {
          US[[sprintf("pU%d", h)]] <- rep(1, bigT)
          US[[sprintf("qU%d", h)]] <- rep(1, bigT)
        }
      }
    }

    # Multiply probabilities
    pvec <- pW
    qvec <- qW

    if (nrow(UQ_unique) > 0) {
      for (h in 1:nrow(UQ_unique)) {
        pvec <- pvec * US[[sprintf("pU%d", h)]]
        qvec <- qvec * US[[sprintf("qU%d", h)]]
      }
    }

    if (nrow(SQ_unique) > 0) {
      for (i in 1:nrow(SQ_unique)) {
        pvec <- pvec * SS[[sprintf("pS%d", i)]]
        qvec <- qvec * SS[[sprintf("qS%d", i)]]
      }
    }

    PX <- prod(pvec) * pie

    # PEj calculation
    cum_pvec <- cumprod(pvec)
    cum_qvec <- rev(cumprod(rev(qvec)))

    PEj <- 0
    for (j in (tn + 1):bigT) {
      partpPEj <- if (j > 1) {
        cum_pvec[j - 1]
      } else {
        1
      }
      partqPEj <- if (j <= bigT) {
        cum_qvec[j]
      } else {
        1
      }
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
