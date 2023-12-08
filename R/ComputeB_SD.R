#' Compute the ACF used in the "Size Distribution" Model
#'
#' This function evaluates the model of the Autocorrelation Function (ACF) of chromatin
#' for given input parameters. It is based on the methodology described in the paper
#' "Characterizing chromatin packing scaling in whole nuclei using interferometric microscopy".
#'
#' @param d Numeric. The uncorrected fractal dimension (Db). Must be a scalar.
#' @param lmin Numeric. The minimum size in nanometers that the ACF is expected to be valid for. Must be a scalar.
#' @param lmax Numeric. The maximum size in nanometers that the chromatin is expected to be fractal for. Must be a scalar.
#' @param r Numeric. The size or sizes at which to evaluate the function. Can be a scalar or a 1D array.
#'
#' @return Numeric. Returns the computed value of the ACF (bnr) based on the input parameters.
#' The returned value will have the same length as the input `r`.
#'
#' @examples
#' d_example <- 2.5
#' lmin_example <- 10
#' lmax_example <- 200
#' r_example <- seq(1, 100, by = 1)
#' ComputeB_SD(d_example, lmin_example, lmax_example, r_example)
#'
#' @export
ComputeB_SD <- function(d, lmin, lmax, r) {
  bnr <- (3 - d) * r^(d - 3) / (lmin^(d - 3) - lmax^(d - 3)) * (gamma_inc(r / lmax, 3 - d) - gamma_inc(r / lmin, 3 - d))
    # pgamma needs to be replaced with whatever the original implemented function is
  return(bnr)}





}