library(dplyr)

filename = "/Users/samhamilton/SigmaConversion/examples/2023_11_17 digoxin exp.xlsx"
data = rio::import(filename)
sigmaIn = data$RMS[!is.na(data$RMS)]
noiseRms = .041644926
sigmaIn = sqrt(sigmaIn^2 - noiseRms^2)
sigmaIn = Re(sigmaIn)
sigmaIn[is.na(sigmaIn)] <- 0 # same at this point
phi = 0.35
Nf = 5e5
ncCenterLambda = 585 #2 * pi / ((2 * pi / 450 + 2 * pi / 700) / 2)
sMin = min(sigma)
sMax = max(sigma)
dfs = seq(2.1,2.9,length.out = 20)
sigmaLUT  <- numeric(length(dfs))
lmax_r = seq(10,10000,by = 20)
liveCellRI = createFromGladstoneDale()
nuSys = SetSystemConfiguration(liveCellRI,na_i = .52,na_c = 1.49,ncCenterLambda,T,T)
thickness = 2 #in u
thickness = min(nuSys$DoF,thickness)* 1e3 # Effective thickness
lmaxs <- numeric(length(dfs))
for(idf in 1:length(dfs)){
  df = dfs[[idf]]

#Nf = CalcNfs(df = df,lmin = 1,lmax_r)
  lmax = numericalInversion(df,lmax_r, Nf)
  lmaxs[idf] = lmax
  lmax = lmaxs[idf]

  sigmaLUT[idf] <- sqrt(integrate(function(lc) Pr(lc,lmin = 1,lmax = lmax,D = df) *
                                    s2_to_int(lc,
                                              nuSys,
                                              thickness,
                                              df),
                                  lower = lmin,
                                  upper = lmax)$value)
}
sigmaLUT
testlmaxs = rio::import("/Users/samhamilton/SigmaConversion/lmaxs.csv")
testlmaxs = testlmaxs[1,] %>% as.vector() %>% unlist()
plot(dfs,sigmaLUT)

## testing just for s2_to_int
for(idf in 1:length(dfs)){


fit <- lm(dfs ~ poly(sigmaLUT, 3, raw = TRUE))
dOut <- predict(fit, newdata = data.frame(sigmaLUT = sigmaIn))

lmax_corrected = LMaxCorrection(dOut,lmin = 1,Nf = Nf)  # i have no idea why the 1/10000 correction is required but it works? weird
dCorrected <- acfd(dOut, lmin, lmax_corrected)

dOut <- predict(fit, newdata = data.frame(sigmaLUT = sigmaIn))
Sigma2DApprox <- function(nuSys,Nf,thickness,sMin,sMax){
  SigmaLUT <-seq(sMin,sMax, length.out = 1000)

  Db = SigmaToDAllInputs(SigmaLUT,nuSys,Nf,2)

}
polyVals = Sigma2DApprox(nuSys,Nf,thickness,sMin,sMax,5)

# Referred to as equation 2 in the matlab code, referring to eqn2 in the paper
# returns a vector of Nfs for corresponding values of lmax_r % This is equation 2 in the paper.
CalcNfs <- function(df,lmin,lmax_r){
  Nfs = 6*(df-3) / df * (1 - (lmin / lmax_r)^df) / ((lmin / lmax_r)^3 - (lmin / lmax_r)^df)
  return(Nfs)

}

Pr <- function(lc,lmin = 1,lmax,D){
  ProbDist <- lc^(D-4) * (D-3) / (lmax^(D-3) - lmin^(D-3))
  return(ProbDist)
}
#s2_to_int(lc, thick, system_config.center_lambda, ri_chromatin, system_config.na_c, lmin, lmax, df)
# lc is lc L is thickness lambda is center lambda, ric is ri chromatin, na is system_config.na_c,
#lmin and lmax are lmin and lmax, but arent used in the actual function df is df, k is the center wavelength
#in vacuum
s2_to_int <- function(lc,systemConfig, Thickness, d) {
  L = Thickness
  k = systemConfig$CenterWavelengthInVacuum
  na = systemConfig$na_c
  s2 = systemConfig$FresnelCoefficient^2 * systemConfig$RIDefinitions$sigma_n^2 *
    ((2/pi * (L * (k^4 * lc^3) * na^2) /
       (1 + (k * lc)^2 * (4 + na^2)) /
       (1 + (k * lc)^2 * 4)) +
      1/4 * (1 - 1/sqrt(1 + (k * lc * na)^2)))
  print(s2)
  return(s2)
}
s2_to_int(lc = 1.5,nuSys,Thickness=thickness)
sigmaLUT[idf] <- sqrt(integrate(function(lc) Pr(lc, lmin, lmax, df) * s2_to_int(lc, thick, system_config$center_lambda, ri_chromatin, system_config$na_c, lmin, lmax, df), lower = lmin, upper = lmax)$value)


