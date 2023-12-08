
#'numericalInversion
#'
#'for a function y = f(x) use interpolation and a lookup table to approximate
#'x = g(y)
#'@param func a function that accepts an array of values and returns y values
#'@param xVals an array of x values to pass to func
#'@param y the input to the inverted function for which you want x values
#'@export
numericalInversion <- function(df,xVals, Nf) {
  yVals <- CalcNfs(df, lmin = 1,lmax = xVals)
  x <- spline(yVals, xVals,xout = Nf, method = "fmm")$y  # or use approx for linear interpolation

  if (x < min(xVals) || x > max(xVals)) {
    warning(paste("Numerical Inversion resorting to extrapolation. x:", x, "y:", Nf))
  }

  return(x)
}


"""
df = 2.4476
lmin = 1
lmax_r = seq(10,25000,by = 30)
Code for testing numericalInversion

CalcNfs <- function(df,lmin,lmax_r){
  Nfs = 6*(df-3) / df * (1 - (lmin / lmax_r)^df) / ((lmin / lmax_r)^3 - (lmin / lmax_r)^df)
  return(Nfs)

}
CalcNfs(df,lmin,lmax_r)

xVals = lmax_r
yVals = CalcNfs(df,lmin,xVals)
 x <- spline(yVals, xVals,xout = Nf, method = "fmm")$y
"""
