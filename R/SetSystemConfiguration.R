#'SetSystemConfiguration
#'
#'A function that creates an class-lite kind of object as the original object in matlab
#'is basically just being used as a named list so there's no need to really overcomplicate things
#'here as I'm not really interested in refactoring code just creating something close to a 1to1 copy
#'
#'@param IsCellGlassInterface whether there is a cell glass interface, if true then the fresnel coefficients will reflect
#'that the cell_glass interface is the reference plane. Otherwise the cell/media interface will be used
#'@param RIDefinitions the collection of Refraction indices used refers to RIDefinition object in the matlab code
#'@param ImmersionRI The refractive index of the immersion media used for the objective
#'@param CenterLambda the center wavelenth based on k
#'@param na_c Collection NA
#'@param na_i Illumination NA
#'@param IsImmersedInOil True/False if so the immersion_ri will be set to the RI of the glass found in ri_def
#'@export
SetSystemConfiguration <- function(RIDefinitions,na_i,na_c,CenterLambda,IsImmersedInOil,IsCellGlassInterface){

  CalculateDoF <- function(CenterLambda,ImmersionRI,na_i){
    DoF = CenterLambda * ImmersionRI / (2 * na_i^2) / 1e3
    return(DoF)
  }

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

  CalcCenterWavelengthInVacuum <- function(CenterLambda){
    CenterWavelength <-  2*pi / CenterLambda
    return(CenterWavelength)
  }

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
