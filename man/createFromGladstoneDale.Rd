% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RMStoD.R
\name{createFromGladstoneDale}
\alias{createFromGladstoneDale}
\title{createFromGladstoneDale
Create an RI Definition Object Using the GladStoneGale Equation}
\usage{
createFromGladstoneDale(ri_media = 1.337, phi = 0.35)
}
\arguments{
\item{ri_media}{refractive index of the media}

\item{phi}{Phi the CVC of the nucleus.}
}
\value{
a named list with the following values
ri_chromatin: refractive index of chromatin
ri_glass:     refractive index of glass
ri_media:     refractive index of media
sigma_n:      SD of RI fluctutations in chromatin
}
\description{
Original in matlab is an object but again the object use is pretty sparse and 
functions primarily as a named list so we will use that structure here is awell.
}
