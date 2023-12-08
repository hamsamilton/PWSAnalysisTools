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
createFromGladstoneDale()
