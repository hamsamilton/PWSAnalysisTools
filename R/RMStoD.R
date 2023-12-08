#'This file contains an R adaptation of the RMS (root mean squared of reflectance across wavelengths
#'captured by PWS) to fractal dimension (D) code written in matlab by Aya Eid (refer to https://pubmed.ncbi.nlm.nih.gov/32870863) 
#'and later to my deduction refactored by Nick Antony. This implementation of that code is designed to
#'be a trimmed down implementation and is missing some features of the original code that I did
#'not require for my own research.
#'
#'
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
