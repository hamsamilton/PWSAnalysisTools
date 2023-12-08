
#' @param dfs uncorrected D values
#' @param lmin not sure for now
#' @param Nf genomic length of the oacking domain

LMaxCorrection <- function(dfs,lmin,Nf){
  
  lmax_r <- seq(10, 25000, by = 30)
  
  lmax_corrected <- numeric(length(dfs)) 
  
  for (idf in 1:length(dfs)) {
      df = dfs[[idf]]
      lmax_corrected[[idf]] <- numericalInversion(df, lmax_r, Nf)
  }
  
  return(lmax_corrected)
  }
  

}
