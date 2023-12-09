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
sMin = min(sigmaIn)
sMax = max(sigmaIn)
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

fit <- lm(dfs ~ poly(sigmaLUT, 3, raw = TRUE))
dOut <- predict(fit, newdata = data.frame(sigmaLUT = sigmaIn))

lmax_corrected = LMaxCorrection(dOut,lmin = 1,Nf = Nf)
dCorrected <- acfd(dOut, lmin, lmax_corrected)


