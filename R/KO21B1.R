#' @title Kodikara et al.'s (2021) "Change-point" model
#'
#' @description
#' The model from Kodikara et al. 2021. Estimates a posterior probability that
#' the species is extant at the test time, and a point estimate and one-tailed
#' \eqn{1 - \alpha} credibile interval on the time of extinction.
#'
#' @param records `data.frame` with two columns: `time` and `certain`. The
#' `time` column contains the date of all sightings, which are each either
#' certain (`certain = TRUE`) or uncertain (`certain = FALSE`).
#' @param alpha desired threshold level (defaults to \eqn{\alpha = 0.05}) of
#' the \eqn{1 - \alpha} credible interval.
#' @param init.time start of the observation period. Defaults to the time of
#' the first sighting, in which case this sighting is removed from the record.
#' @param test.time end of the observation period, typically the present day
#' (defaults to the current year).
#'
#' @returns a `list` object with the original parameters and the p(extant),
#' point estimate, and credible interval included as elements. The credible
#' interval is a two-element numeric vector called `cred.int`.
#'
#' @note
#' Sampling effort is assumed to be constant.
#'
#' @references
#' **Key Reference**
#'
#' Kodikara, S., Demirhan, H., Wang, Y., Solow, A., & Stone, L. (2020).
#' Inferring extinction year using a Bayesian approach. *Methods in Ecology*
#' *and Evolution*, 11(8), 964-973. \doi{10.1111/2041-210x.13408}
#'
#' @examples
#' # Run the Ivory-billed Woodpecker analysis from Kodikara et al. 2021
#' KO21B1(woodpecker2, test.time = 2010)
#'
#' @export

KO21B1 <- function(records, alpha = 0.05, init.time = min(records$time),
                   test.time = as.numeric(format(Sys.Date(), "%Y"))) {

  # Sort records
  records <- sort_by(records, ~time)

  # Create date and certainty vectors
  sighting <- records$time - init.time
  cat <- 1 - as.integer(records$certain)
  if (init.time == min(records$time)) {
    sighting <- sighting[-1]
    cat <- cat[-1]
  }

  # Calculate model parameters
  Tt <- test.time - init.time
  t_n <- max(subset(records, certain == TRUE)$time) - init.time
  N <- length(sighting)
  y <- c()

  # Put the information into a list
  dataList <- list(
    d = sighting,
    y = y,
    t_n = t_n,
    N = N,
    cat = cat,
    Tt = Tt)

  # Define the model
  modelpois <- paste0("
                        data {
                        C <- 10000000000000 # JAGS does not warn if too small!


                        for (i in 1:N) {
                          ones[i] <- 1
                        }


                        }

                        model {


                        for(i in 1:N){
                          y[i]<- ifelse(cat[i]<1,
                          ifelse(tau<=t_n,(10^(-100)), ifelse(tau<=Tt,
                          ((a/(sigma^a))*d[i]^(a-1)*exp(-1*((tau/sigma)^a))),
                          ((a/(sigma^a))*d[i]^(a-1)*exp(-1*((Tt/sigma)^a))))),
                          ifelse(tau<=t_n,(10^(-100)), ifelse(d[i]<=tau,
                          ((a_U1/(sigma_U1^a_U1))*d[i]^(a_U1-1)*
                          exp(-1*((tau/sigma_U1)^a_U1))),
                          ((a_U2/(sigma_U2^a_U2))*d[i]^(a_U2-1)*
                          exp(-1*((Tt/sigma_U2)^a_U2-(tau/sigma_U2)^a_U2))))))/C



                          ones[i]~ dbern( y[i] )
                        }

                        tau ~ dexp(0.005)
                        sigma~dunif(0,1000)
                        a~dunif(0,1000)
                        sigma_U1~dunif(0,1000)
                        a_U1~dunif(0,1000)
                        sigma_U2~dunif(0,1000)
                        a_U2~dunif(0,1000)
                        }

                        ")

  writeLines(modelpois, con = "model_pois.txt")

  # Sink (to suppress hyper-verbose console outputs)
  sink(file = tempfile())

  # Run the chains
  jagsmodelpois <- rjags::jags.model(file = "model_pois.txt", data = dataList,
                                     n.chains = 4, n.adapt = 10000)
  update(jagsmodelpois, n.iter = 10000)
  codaSamplespois <- rjags::coda.samples(jagsmodelpois, variable.names = c(
    "tau", "a", "sigma", "a_U1", "sigma_U1", "a_U2", "sigma_U2"),
    n.iter = 130000, thin = 13)

  sink()

  unlink("model_pois.txt")

  # Extract posteriors
  posterior <- as.data.frame(as.matrix(codaSamplespois))

  # Calculate p(extant)
  p.extant <- mean(posterior$tau + init.time > test.time)

  # Calculate point estimate
  estimate <- median(posterior$tau) + init.time

  # Calculate credible interval bounds
  cred.int.lower <- as.numeric(quantile(posterior$tau, 0)) + init.time
  cred.int.upper <- as.numeric(quantile(posterior$tau, 1 - alpha)) + init.time

  # Output
  output <- list(
    records = records,
    alpha = alpha,
    init.time = init.time,
    test.time = test.time,
    p.extant = p.extant,
    estimate = estimate,
    cred.int = c(cred.int.lower, cred.int.upper)
  )

  return(output)

}

# Declare certain as a known global variable (column in records data.frame):
utils::globalVariables("certain")
