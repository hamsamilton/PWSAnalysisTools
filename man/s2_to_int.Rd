% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RMStoD.R
\name{s2_to_int}
\alias{s2_to_int}
\title{s2_to_int}
\usage{
s2_to_int(lc, systemConfig, Thickness, d)
}
\arguments{
\item{lc}{the variable to integrate over from lmin to lmax}

\item{Thickness}{the thickness of the material}

\item{d}{a df value, whatever that means}

\item{a}{systemConfig object}
}
\description{
The s2 to int equation used within the integrate function, I don't understand the intention of 
this math so it's difficult for me to precisely describe
}
\details{
# lc is lc L is thickness lambda is center lambda, ric is ri chromatin, na is system_config.na_c,
lmin and lmax are lmin and lmax, but arent used in the actual function df is df, k is the center wavelength
in vacuum
}
