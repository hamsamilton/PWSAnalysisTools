#' CorrectD
#'
#' Take rawD and convert to the corrected D
#' @param RawD Uncorrected D values
#' @param lmin whatever lmin is, it wasnt well documented
#' @param lmax whatever lmax is, it wasnt well documented
#' @export
#'
acfd <- function(dfs, lmin, lmax) {

  delta <- 0.1
  mid <- 60
  point <- (lmax + lmin) / mid

  d_corrected <- numeric(length(dfs))

  for (idf in 1:length(dfs)) {
    d = dfs[[idf]]

      d_corrected[[idf]] <- 3 + (log(ComputeB_SD(d, lmin, lmax[[idf]], point[[idf]] + delta)) -
                                  log(ComputeB_SD(d, lmin, lmax[[idf]], point[[idf]]))) /
        (log(point[[idf]] + delta) - log(point[[idf]]))
  }
  return(d_corrected)
  }



