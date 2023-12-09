#'This file contains an R adaptation of the RMS (root mean squared of reflectance across wavelengths
#'captured by PWS) to fractal dimension (D) code written in matlab by Aya Eid (refer to https://pubmed.ncbi.nlm.nih.gov/32870863) 
#'and later to my deduction refactored by Nick Antony. This implementation of that code is designed to
#'be a trimmed down implementation and is missing some features of the original code that I did
#'not require for my own research.
#'

#'ConvertRMStoD
#'
#'The main function of this software module. It converts RMS to D (Fractal Dimension values)
#'
#'@param RMSvals a vector of RMS values
#'@param noiseRMS the amount of background noise (RMS values of the background image)
#'@param phi Not sure, using defaults in the MatLab script
#'@param Nf size of the packing domain
#'@param ncCenterLambda a constant taken from the script
#'@param na_i  Illumination NA
#'@param na_c Collection NA
#'@param thickness thickness of material in uM
#'@param IsImmersedInOil True/False if so the immersion_ri will be set to the RI of the glass found in ri_def
#'@param ImmersionRI The refractive index of the immersion media used for the objective
#'@param IsCellGlassInterface whether there is a cell glass interface, if true then the fresnel coefficients will reflect
#'@return a vector of D values
#'@export
ConvertRMStoD <- function(RMSvals,noiseRMS,phi = .35, Nf = 5e5, ncCenterLambda = 585,
                          na_i = .52,na_c = 1.49, IsCellGlassInterface = T,IsImmersedInOil = T, thickness = 2){
  
  # calculate other fields required
  sMin       = min(sigmaIn)
  sMax       = max(sigmaIn)
  dfs        = seq(2.1,2.9,length.out = 20)
  sigmaLUT   = numeric(length(dfs))
  lmax_r     = seq(10,10000,by = 20)
  nuSys      = SetSystemConfiguration(RIDefinitions= createFromGladstoneDale(),
                                      na_i                 = na_i,
                                      na_c                 = na_c,
                                      ncCenterLambda       = ncCenterLambda,
                                      IsCellGlassInterface = IsCellGlassInterface,
                                      IsImmersedInOil      = IsImmersedInOil)
  thickness = min(nuSys$DoF, thickness)
  lmaxs = numeric(length(dfs))
  
  # Estimate the function required for estimating D from RMS
  for(idf in 1:length(dfs)){
    df = dfs[[idf]]
    
    lmax = numericalInversion(df,lmax_r, Nf)
    
    lmaxs[idf] = lmax
    sigmaLUT[idf] <- sqrt(integrate(function(lc) Pr(lc,lmin = 1,lmax = lmax,D = df) *
                                      s2_to_int(lc,
                                                nuSys,
                                                thickness,
                                                df),
                                    lower = lmin,
                                    upper = lmax)$value)
  }
  
  
  #Estimate D values
  fit <- lm(dfs ~ poly(sigmaLUT, 3, raw = TRUE))
  dOut <- predict(fit, newdata = data.frame(sigmaLUT = sigmaIn))
  
  #Correct D values
  lmax_corrected = LMaxCorrection(dOut,lmin = 1,Nf = Nf)
  dCorrected <- acfd(dOut, lmin, lmax_corrected)
  
  return(dCorrected)
}
  
#'createFromGladstoneDale
#'Create an RI Definition Object Using the GladStoneGale Equation
#'
#'Original in matlab is an object but again the object use is pretty sparse and 
#'functions primarily as a named list so we will use that structure here is awell.
#'@param ri_media refractive index of the media
#'@param phi Phi the CVC of the nucleus.
#'@return a named list with the following values
#'ri_chromatin: refractive index of chromatin
#'ri_glass:     refractive index of glass
#'ri_media:     refractive index of media
#'sigma_n:      SD of RI fluctutations in chromatin
#'@export
createFromGladstoneDale <- function(ri_media = 1.337, phi = .35) {
  phi_mc <- 0.05 # Phi is CVC - mobile crowders - 5% of remaining volume that is unoccupied
  rho_protein <- 1.35 # Density of protein
  rho_chromatin <- 0.56 # Density of chromatin
  riinc <- 0.1799 # Gladstone equation's alpha
  
  # Calculating ri_chromatin
  ri_chromatin <- ri_media + riinc * (rho_chromatin * phi + rho_protein * phi_mc * (1 - phi))
  
  # Calculating sigma_n
  sigma_n <- sqrt(phi * (1 - phi)) * (riinc * rho_chromatin * (1 - phi_mc) - riinc * rho_protein * phi_mc)
  
  # Returning the result as a list
  return(list(ri_glass = 1.517, ri_media = ri_media, ri_chromatin = ri_chromatin, sigma_n = sigma_n))
}


#'Create an RI Definition Object Using the GladStoneGale Equation
#'
#'The Refractive Index (RI object) is an object that contains RI information on 
#'assorted materials required for D calculations. The "object" also calculates values
#'for other constants which are calculated from tehse values.
#' Default values were pulled from
#'scripts in the original matlab code but I may not have fully understood their intentions
#'
#'Original in matlab is an object but again the object use is pretty sparse and 
#'functions primarily as a named list so we will use that structure here is awell.
#'@param ri_media refractive index of the media
#'@param phi Phi the CVC of the nucleus.
#'@return a named list with the following values
#'ri_chromatin: refractive index of chromatin
#'ri_glass:     refractive index of glass
#'ri_media:     refractive index of media
#'sigma_n:      SD of RI fluctutations in chromatin
#'@export
createRIFromGladstoneDale <- function(ri_media = 1.337, phi = .35) {
  phi_mc <- 0.05 # Phi is CVC - mobile crowders - 5% of remaining volume that is unoccupied
  rho_protein <- 1.35 # Density of protein
  rho_chromatin <- 0.56 # Density of chromatin
  ri_glass <-  1.517
  riinc <- 0.1799 # Gladstone equation's alpha
  
  # Calculating ri_chromatin
  ri_chromatin <- ri_media + riinc * (rho_chromatin * phi + rho_protein * phi_mc * (1 - phi))
  
  # Calculating sigma_n
  sigma_n <- sqrt(phi * (1 - phi)) * (riinc * rho_chromatin * (1 - phi_mc) - riinc * rho_protein * phi_mc)
  
  # Returning the result as a list
  return(list(ri_glass = ri_glass, ri_media = ri_media, ri_chromatin = ri_chromatin, sigma_n = sigma_n))
}

#'SetSystemConfiguration
#'
#'This system configuration object contains fields that are required to calculate
#'D. This instatiation method automatically calculates fields that are calculated from the inptus
#'
#'A function that creates an class-lite kind of object as the original object in matlab
#'is basically just being used as a named list so there's no need to really overcomplicate things
#'here as I'm not really interested in refactoring code just creating something close to a 1to1 copy
#'
#'@param IsCellGlassInterface whether there is a cell glass interface, if true then the fresnel coefficients will reflect
#'that the cell_glass interface is the reference plane. Otherwise the cell/media interface will be used
#'@param RIDefinitions the collection of Refraction indices used refers to RIDefinition object in the matlab code
#'@param IsGlassCellInterface T/F Used for calculating the Fresnel Coefficient.
#'@param ImmersionRI The refractive index of the immersion media used for the objective
#'@param CenterLambda the center wavelenth based on k
#'@param na_c Collection NA
#'@param na_i Illumination NA
#'@param IsImmersedInOil True/False if so the immersion_ri will be set to the RI of the glass found in ri_def
#'@return a named list with the following items
#' RIDefinitions used for input
#' na_i: See params
#' na_c: see params
#' CenterLambda: see params
#' IsImmersedInOil: see params
#' IsCellGlassInterface: see params
#' DoF: Depth of focus
#' CenterWavelengthInVacuum: Used for D calculations
#' FresnelCoefficient: Used for calculating D (Need to learn more)
#' ImmersionRI: RI of immersion media
#'@export
SetSystemConfiguration <- function(RIDefinitions,na_i,na_c,CenterLambda,IsImmersedInOil,IsCellGlassInterface){
  
  #Depth of focus calculation used for calculating effective depth
  CalculateDoF <- function(CenterLambda,ImmersionRI,na_i){
    DoF = CenterLambda * ImmersionRI / (2 * na_i^2) / 1e3
    return(DoF)
  }
  
  # calculates the Fresnel Coefficient used for D calculation
  CalculateFresnelCoefficient <- function(RIDefinitions,IsCellGlassInterface){
    
    # not sure what this means exactly
    R_reference = ( (RIDefinitions$ri_glass - RIDefinitions$ri_media) /
                      (RIDefinitions$ri_glass + RIDefinitions$ri_media))^2
    
    if(IsCellGlassInterface){n0 =RIDefinitions$ri_glass}
    else{n0 = RIDefinitions$ri_media}
    
    n1 = RIDefinitions$ri_chromatin
    
    FresnelCoefficient = 4 * n0 * n1 * (n0 - n1) / (n0 + n1)^3 / R_reference
    
    return(FresnelCoefficient)
  }
  
  # Used in some calculations, not familiar enough with the math to really comment
  CalcCenterWavelengthInVacuum <- function(CenterLambda){
    CenterWavelength <-  2*pi / CenterLambda
    return(CenterWavelength)
  }
  
  # Which RI to use depending on exp design.
  ImmersionRI <- if(IsImmersedInOil){
    RIDefinitions$ri_glass}
  else{RIDefinitions$ri_media}
  
  
  
  SystemConfiguration <- list(RIDefinitions            = RIDefinitions,
                              na_i                     = na_i,
                              na_c                     = na_c,
                              CenterLambda             = CenterLambda,
                              ImmersionRI              = ImmersionRI,
                              IsCellGlassInterface     = IsCellGlassInterface,
                              FresnelCoefficient       = CalculateFresnelCoefficient(RIDefinitions,
                                                                                     IsCellGlassInterface),
                              CenterWavelengthInVacuum = CalcCenterWavelengthInVacuum(CenterLambda),
                              DoF                      = CalculateDoF(CenterLambda,ImmersionRI,na_i))
  
  return(SystemConfiguration)
}


#'Pr
#'
#'The probability distribution used for the integration portion of D calculation
#'@param lc the variable this will be integrating over
#'@param lmin the lower bound of the distribution
#'@param lmax the upper bound of the distribution
#'@param D value of df
#'
#'@export
Pr <- function(lc,lmin = 1,lmax,D){
  ProbDist <- lc^(D-4) * (D-3) / (lmax^(D-3) - lmin^(D-3))
  return(ProbDist)
}



#'s2_to_int
#'
#'The s2 to int equation used within the integrate function, I don't understand the intention of 
#'this math so it's difficult for me to precisely describe
#'
#'# lc is lc L is thickness lambda is center lambda, ric is ri chromatin, na is system_config.na_c,
#'lmin and lmax are lmin and lmax, but arent used in the actual function df is df, k is the center wavelength
#'in vacuum
#'@param lc the variable to integrate over from lmin to lmax
#'@param a systemConfig object
#'@param Thickness the thickness of the material
#'@param d a df value, whatever that means
#'@export
s2_to_int <- function(lc,systemConfig, Thickness, d) {
  L = Thickness
  k = systemConfig$CenterWavelengthInVacuum
  na = systemConfig$na_c
  s2 = systemConfig$FresnelCoefficient^2 * systemConfig$RIDefinitions$sigma_n^2 *
    ((2/pi * (L * (k^4 * lc^3) * na^2) /
        (1 + (k * lc)^2 * (4 + na^2)) /
        (1 + (k * lc)^2 * 4)) +
       1/4 * (1 - 1/sqrt(1 + (k * lc * na)^2)))
  return(s2)
}

#' The function for correcting LMAX values
#'
#' @param dfs uncorrected D values
#' @param lmin not sure for now
#' @param Nf genomic length of the packing domain
#' @export
LMaxCorrection <- function(dfs,lmin,Nf){
  
  lmax_r <- seq(10, 25000, by = 30)
  lmax_corrected <- numeric(length(dfs)) 
  
  for (idf in 1:length(dfs)) {
    df = dfs[[idf]]
    lmax_corrected[[idf]] <- numericalInversion(df, lmax_r, Nf)
  }
  
  return(lmax_corrected)
}

#' acfd
#'
#'With the name of this function what this does your guess is as good as mine
#'Used for correcting d values
#' Take rawD and convert to the corrected D
#' @param dfs Uncorrected D values
#' @param lmin whatever lmin is, it wasnt well documented
#' @param lmax whatever lmax is, it wasnt well documented
#' @export
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


#'numericalInversion
#'
#'for a function y = f(x) use interpolation and a lookup table to approximate
#'x = g(y). The function this function is supposed to use CalcNfs is hardcoded into this function
#'for simplicity. The interpolation method is slightly different between R and MatLab
#'and may cause some small variations
#'@param func a function that accepts an array of values and returns y values
#'@param xVals an array of x values to pass to func
#'@param y the input to the inverted function for which you want x values
#'@export
numericalInversion <- function(df,xVals, Nf) {
  
  
  #Referred to as equation 2 in the matlab code, referring to eqn2 in the paper
  # returns a vector of Nfs for corresponding values of lmax_r % This is equation 2 in the paper.\
  #Used only in the context of numericalInversion from all I can tell.
  
  CalcNfs <- function(df,lmin,lmax_r){
    Nfs = 6*(df-3) / df * (1 - (lmin / lmax_r)^df) / ((lmin / lmax_r)^3 - (lmin / lmax_r)^df)
    return(Nfs)
  }
  
  yVals <- CalcNfs(df, lmin = 1,lmax = xVals)
  x <- spline(yVals, xVals,xout = Nf, method = "fmm")$y  # or use approx for linear interpolation
  
  if (x < min(xVals) || x > max(xVals)) {
    warning(paste("Numerical Inversion resorting to extrapolation. x:", x, "y:", Nf))
  }
  
  return(x)
}



#' ComputeB_SD
#' helper function used to Compute the ACF used in the Size Distribution Model
#' This function evaluates the model of the Autocorrelation Function (ACF) of chromatin
#' for given input parameters. It is based on the methodology described in the pape
#' Characterizing chromatin packing scaling in whole nuclei using interferometric microscopy.
#'
#' @param d Numeric. The uncorrected fractal dimension (Db). Must be a scalar.
#' @param lmin Numeric. The minimum size in nanometers that the ACF is expected to be valid for. Must be a scalar.
#' @param lmax Numeric. The maximum size in nanometers that the chromatin is expected to be fractal for. Must be a scalar.
#' @param r Numeric. The size or sizes at which to evaluate the function. Can be a scalar or a 1D array.
#'
#' @return Numeric. Returns the computed value of the ACF (bnr) based on the input parameters.
#' The returned value will have the same length as the input r.
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
  return(bnr)}


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