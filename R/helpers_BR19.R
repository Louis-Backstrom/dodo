#' @title rw64 function from Brook et al. 2019
#'
#' @description
#' Helper function. From provided code (ede_functions.r).
#'
#' @param ts the time-series of sighting records.
#'
#' @returns an estimate of the extinction time.
#'
#' @references
#' **Key Reference**
#'
#' Brook, B. W., Buettel, J. C., & Jaric, I. (2019). A fast re-sampling method
#' for using reliability ratings of sightings with extinction-date estimators.
#' *Ecology*, 100(9), e02787. \doi{10.1002/ecy.2787}
#'
#' @noRd

rw64 <- function(ts) {
  n <- length(ts)
  return(ts[n] + (ts[n] - ts[n - 1]))
}

#' @title so93 function from Brook et al. 2019
#'
#' @description
#' Helper function. From provided code (ede_functions.r).
#'
#' @param ts the time-series of sighting records.
#'
#' @returns an estimate of the extinction time.
#'
#' @references
#' **Key Reference**
#'
#' Brook, B. W., Buettel, J. C., & Jaric, I. (2019). A fast re-sampling method
#' for using reliability ratings of sightings with extinction-date estimators.
#' *Ecology*, 100(9), e02787. \doi{10.1002/ecy.2787}
#'
#' @noRd

so93 <- function(ts) {
  sr <- ts - ts[1]
  n <- length(sr) - 1
  return(ts[1] + (n + 1) / n * sr[n + 1])
}

#' @title mc06 function from Brook et al. 2019
#'
#' @description
#' Helper function. From provided code (ede_functions.r).
#'
#' @param ts the time-series of sighting records.
#'
#' @returns an estimate of the extinction time.
#'
#' @references
#' **Key Reference**
#'
#' Brook, B. W., Buettel, J. C., & Jaric, I. (2019). A fast re-sampling method
#' for using reliability ratings of sightings with extinction-date estimators.
#' *Ecology*, 100(9), e02787. \doi{10.1002/ecy.2787}
#'
#' @noRd

mc06 <- function(ts) {
  sr <- ts - ts[1]
  n <- length(sr) - 1
  return(ts[1] + ceiling(sr[n + 1] + log(0.5, 1 - (n / sr[n + 1]))))
}

#' @title ole function from Brook et al. 2019
#'
#' @description
#' Helper function. From provided code (ede_functions.r).
#'
#' @param ts the time-series of sighting records.
#'
#' @returns an estimate of the extinction time.
#'
#' @references
#' **Key Reference**
#'
#' Brook, B. W., Buettel, J. C., & Jaric, I. (2019). A fast re-sampling method
#' for using reliability ratings of sightings with extinction-date estimators.
#' *Ecology*, 100(9), e02787. \doi{10.1002/ecy.2787}
#'
#' @noRd

ole <- function(ts) {
  gam.fit <- function(i, j, v) {
    return((gamma(2 * v + i) * gamma(v + j)) / (gamma(v + i) * gamma(j)))
  }
  sights <- rev(sort(ts))
  k <- length(sights)
  v <- (1 / (k - 1)) * sum(log((sights[1] - sights[k]) /
                                 (sights[1] - sights[2:(k - 1)])))
  lambda <- outer(1:k, 1:k, gam.fit, v = v)
  lambda <- ifelse(lower.tri(lambda), lambda, t(lambda))
  e <- matrix(rep(1, k), ncol = 1)
  a <- as.vector(solve(t(e) %*% solve(lambda) %*% e)) * solve(lambda) %*% e

  return(round(sum(t(a) %*% sights)))
}

#' @title bbj.2018 function from Brook et al. 2019
#'
#' @description
#' Helper function. From provided code (ede_functions.r).
#'
#' @param dd the sightings record dataset.
#' @param iter the number of resampling iterations (defaults to 10,000)/
#' @param ey the test.time parameter (defaults to 2018, typically changed).
#' @param m which model to fit (defaults to "LAD").
#' @param plot whether to produce plots (defaults to `FALSE`).
#' @param alpha desired significance level.
#' @param cores number of cores to use (defaults to `NULL`, in which case
#' parallel::detectCores() is run).
#'
#' @returns a `list` object with the relevant estimates included as elements.
#'
#' @references
#' **Key Reference**
#'
#' Brook, B. W., Buettel, J. C., & Jaric, I. (2019). A fast re-sampling method
#' for using reliability ratings of sightings with extinction-date estimators.
#' *Ecology*, 100(9), e02787. \doi{10.1002/ecy.2787}
#'
#' @noRd

bbj.2018 <- function(dd, iter = 10000, ey = 2018, m = "LAD", plot = FALSE,
                     alpha, cores = NULL) {

  if (is.null(cores)) {
    cores <- parallel::detectCores() - 1
  }

  if (exists(m)) {
    ede <- match.fun(m)
  } else {
    ede <- function(ts) {
      return(max(ts) + 1)
    }
  }

  r.rec <- function(d) {
    sr <- d$year[which(d$prob > runif(length(d$year)))]
    ext <- if(length(sr) > 2) {
      ede(sr[order(sr)])
    }
    return(ifelse(length(ext) == 0, NA, ext))
  }

  if(.Platform$OS.type == "windows") {
    cl <- parallel::makeCluster(cores); dd
    ext.yr <- parallel::parSapply(cl = cl, 1:iter,
                        function(x) {r.rec(dd)})
    parallel::stopCluster(cl)
  }
  else {
    ext.yr <- parallel::mclapply(1:iter, mc.cores = cores,
                       function(x) {r.rec(dd)})
    ext.yr <- do.call(rbind, ext.yr)
  }

  ext.yr <- na.omit(ext.yr)
  min.y <- round(min(ext.yr))
  if(sd(ext.yr) == 0) {
    return(list(MTE = round(ext.yr[1]),
                UCI = round(ext.yr[1]),
                PP = ifelse(round(ext.yr[1]) >= ey, 1, 0)))
  }

  freq.ext <- hist(ext.yr, breaks = round(max(ext.yr)) - min.y, plot = FALSE)
  bbj <- list()
  bbj$Method <- m
  x <- min.y:ey
  pp.vec <- rep(NA, length(x))
  index <- 1

  for(i in min.y:ey) {
    pp.vec[index] <- 1 - (sum(freq.ext$counts[which(freq.ext$mids < i)]) / length(ext.yr))
    index <- index + 1
  }

  if(round(mean(ext.yr)) >= ey) {
    bbj$MTE <- which(pp.vec > 0.5)
    print("Median TE reported because mean TE is later than end year")
  } else {
    bbj$MTE <- round(mean(ext.yr))
  }
  bbj$LCI <- ifelse(x[min(which(pp.vec > alpha / 2))] == ey,
                    NA, x[min(which(pp.vec > alpha / 2))])
  bbj$UCI <- ifelse(x[max(which(pp.vec > alpha / 2))] == ey,
                    NA, x[max(which(pp.vec > alpha / 2))])
  bbj$end_year <- pp.vec[length(x)]
  names(bbj$end_year) <- paste("PP", ey)

  if(plot == TRUE) {
    par(mfrow = c(1, 2))
    hist(ext.yr, main = "Freq. of predicted TE", xlab = "Year")
    abline(v = ey, col = "blue", lty = 2)
    plot(pp.vec ~ x, main = "Cumulative PP", xlab = "Year",
         ylab = "Prob. persistence", xlim = c(min(x), ey - 1))
    lines(x,pp.vec)
    abline(v = bbj$MTE, col = "blue", lty = 2)
    abline(v = bbj$UCI, col = "red", lty = 3)
    par(mfrow = c(1, 1))
  }

  return(bbj)
}
