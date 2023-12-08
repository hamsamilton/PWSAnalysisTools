#' Upper Incomplete Gamma Function
#'
#' This function evaluates the upper incomplete gamma function \eqn{\Gamma(a, x)}
#' for non-negative values of the argument \code{x}. It extends the computation
#' to include negative values of the parameter \code{a}, which is not supported
#' by the standard gamma functions in R.
#'
#' @param x Numeric vector; function argument which must be non-negative.
#' @param a Numeric scalar; parameter of the incomplete gamma function, can be negative.
#'
#' @return Numeric vector of the same length as \code{x}. 
#'         Returns the computed values of the upper incomplete gamma function.
#'         NaN values are returned where elements of \code{x} are negative.
#'
#' @details The function is an extension of the standard incomplete gamma
#'          functions, allowing for negative values of the parameter \code{a}.
#'          For \code{a > 0}, it uses the \code{pgamma} function from base R. 
#'          For \code{a = 0}, it relies on the exponential integral function
#'          from the \code{expint} package. For \code{a < 0}, the function 
#'          uses recursion to compute the values.
#'
#' @examples
#' x <- seq(0.01, 8, by = 0.01)
#' a_pos <- 1
#' a_neg <- -2.3
#' plot(x, gamma_incomplete(x, a_pos), type = 'l', col = 'blue',
#'      main = 'Upper Incomplete Gamma Function',
#'      xlab = 'x', ylab = expression(Gamma(a, x)))
#' lines(x, gamma_incomplete(x, a_neg), col = 'red')
#'
#' @export
gamma_incomplete <- function(x, a) {
if (length(a) > 1) {
  stop("`a` must be a scalar, not an array")
}

fun <- numeric(length(x))
p <- which(x >= 0)
q <- which(x < 0)

fun[p] <- sapply(x[p], function(xi) gamma_inc(a, xi))
fun[q] <- NaN

return(fun)
}

gamma_inc <- function(x, a) {
  if (a < -10) {
    return(NaN)  # Handling very negative 'a' values
  } else if (a == 0) {
    return(expint::expint_E1(x))  # Using the expint package for the exponential integral
  } else if (a > 0) {
    return(gamma(a) * pgamma(x, a, lower.tail = FALSE))
  } else if (a < 0) {
    return((gamma_inc(a + 1, x) - x^a * exp(-x)) / a)
  } else {
    stop("No case matched conditions.")
  }
}

#' 