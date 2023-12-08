#setwd("/Users/samhamilton/Library/Mobile Documents/com~apple~CloudDocs/Thesis_Aim1/libraries/PWSAnalysisTools")
# devtools::document() ; devtools::build() ; devtools::install()
#'AdaptPWSData
#'
#'@param df a table produced by the PWS system
#'@param ExpName Name for the experiment
#'@param ExpDate Date of the experiment
#'@param IsolateTimeStrat A function to isolate the time variable from the
#'Path column
#'@param IsolateCondStrat A function to isolate the condition variable from the path column
AdaptPWSData <- function(df,
         ExpName = "ExpName",
         ExpDate = "ExpDate",
         IsolateTimeStrat,
         IsolateCondStrat) {

  df = df %>%
    rename("Image" = "Cell#",
           "ROI"   =  "ROI#") %>%
    mutate(Experiment_Name = ExpName,
           ExpDate         = ExpDate,
           ExpTime         = IsolateTimeStrat(Path),
           ExpCond         = IsolateCondStrat(Path),
           Image_ROI       = str_c(Image,`ROI Name`,ROI),
           RMS             = as.numeric(RMS),
           `Dynamics Reflectance` = as.numeric(`Dynamics Reflectance`),
           .before = 1)
  df = df[,names(df) %!in% c("Path","ROI Name","Image")]
  df = df[!str_detect(df$Image_ROI,"bgd"),]
  return(df)}


#'Plot_SampleWise
#'
#' The goal of this function is to generate a box and whisker plots with a connected
#' line graph showing the change in a value for an item over time.
#'
#'@param df an adapted dataframe ready for analysis
#'@param value_col what value should be plotted
#'@return a list of box and whiskers plots, with samples connected by lines across
#'their trajectories
Plot_SampleWise <- function(df,value_col){
  plts = lapply(unique(df$ExpCond),function(Cond){
    df = df %>% filter(ExpCond == Cond)
    plt =   ggplot(df,aes(x = ExpTime,y = get(value_col))) +
      geom_boxplot(outlier.alpha = NULL) +
      geom_point(aes(color = Image_ROI),alpha = .2) +
      geom_line(aes(x = ExpTime,
                    y = get(value_col),
                    group = Image_ROI,
                    color = Image_ROI),
                alpha = .2) +
      labs(x = "Time after Treatment (min)",
           y = value_col,
           title = Cond) +
      theme_classic() +
      theme(legend.position="none")
    return(plt)})
  return(plts)
}

#'t.testr
#'
#' Perform t.tests between timepoints for each condition
#' @param df an adapted dataframe ready for analysis
#' @param value_col the name of column containing the variable to evalutate
#' @param paired    should the t-test be paired or not.
#'
#'
#'@param df an adapted dataframe ready for analysis
#'@param value_col what value should be plotted
#'@return a list of box and whiskers plots, with samples connected by lines across
#'their trajectories
t.testr = function(df,value_col,paired = T){
  # The helper function that processes each pair of data frames
  process_pair <- function(df1, df2, id_col, value_col,paired = T) {
    if(paired){
      # Step 1,2,3: Find common IDs, filter data frames,and arrange them to match
      common_ids <- intersect(df1[["Image_ROI"]], df2[["Image_ROI"]])
      df1 <- df1 %>%
        filter(Image_ROI %in% common_ids) %>%
        arrange(Image_ROI)
      df2 <- df2 %>%
        filter(Image_ROI %in% common_ids) %>%
        arrange(Image_ROI)
    }
    # Step 4: Perform paired t-test
    t_result <- t.test(df1[[value_col]], df2[[value_col]], paired = paired)
    return(t_result = t_result)
  }
  # For each experimental condition....
  expdfs = lapply(unique(df[["ExpCond"]]),function(condition){
    df_exp = df %>% filter(ExpCond == condition)

    # create a list of dataframes for each timepoint.
    timepoint_dfs = lapply(sort(unique(df_exp[["ExpTime"]])),function(timepoint){
      timepoint_df = df_exp %>% filter(ExpTime == timepoint) %>%
        arrange(Image_ROI)
      return(timepoint_df)
    })
    return(timepoint_dfs)
  })

  # name the list of list of dataframes
  names(expdfs) = unique(df[["ExpCond"]])
  results_table = lapply(expdfs,function(things){

    # Initialize list to store results
    results_list <- list()

    # Loop over each unique pair of data frames
    for (i in 1:(length(things) - 1)) {
      for (j in (i + 1):length(things)) {
        # Perform the paired t-test and store results
        result <- process_pair(things[[i]], things[[j]], "ID", value_col)
        results_list[[
          paste0(sort(unique(df[["ExpTime"]]))[i],
                 "_",
                 sort(unique(df[["ExpTime"]]))[j])]] <- result
      }
    }
    # Extract the pvalues from the t.test results
    pvals = sapply(results_list,function(comparison){
      pval = comparison$p.value
      return(pval)
    })
    return(pvals)
  })
  # Structure data as an easy to read table
  results_table = data.frame(results_table) %>%
    rownames_to_column() %>%
    tidyr::separate(rowname,into=c("Time 1","Time 2"),sep = "_")
  return(results_table)
}

#'
#'@param sigmaIn Sigma "RMS" values to convert
#'@param systemConfig the systemConfig object to obtain information from
#'@param Nf the genomic length of a packing domain. Often values of 5e5 and 1e6 are used
#'@param sampleThickness estimated thickness of the sample. The DOF of a LCPWS2 is 1.1uM
#'
#'@return dOut analgous to D_b in ayas papet he model parameter
#'@return dCorrected the True D we care about
#'@return Nf_expected the genomic length we expect based on D
#'@return lmax_corrected calculation of lmax from nf and db based on equation 2
ConvertRMStoD <- function(sigmaIn,systemConfig,Nf,sampleThickness)

  CalcEffectiveCellThickness <- function(systemConfig,thickIn){

    center_lambda = systemConfig$CenterLamba
    systemConfig$ImmersionRI
    dof <- systemConfig$CenterLamba * systemConfig$ImmersionRI /
          (2 * systemConfig$na_i^2) / 1e3

    EffectiveThickness <- min(dof,thickIn) * 1e3

    return(EffectiveThickness)
  }

  CalcEffectiveCellThickness(nuSys,2)

  dfs = seq(2.1,2.9,length.out = 20)
  sigmaLUT <- matrix(0, nrow = length(dfs), ncol = length(dfs))
    #lc is lc L = thick, center lambda = lambda, ric = richromatin,na = nca.c,lmin, minimal fractal length possible (1)
  #lmax is the output of the numerical conversion
  SolveSigma2 <- function(lc, L, lambda, ric, na, lmin = 1, lmax){ ...

                          Fcoeff ^ 2 * systemConfig$RIDefinitions$sigma_n^2 * (
                                       (2/pi * (L * (k^4 * lc^3) * na^2) ...
                                                                                  ./ (1 + (k.*lc).^2*(4 + na^2)) ...
                                                                                  ./ (1 + (k.*lc).^2 * 4)) ...
                                                                                + 1/4 * (1 - 1./sqrt(1 + (k.*lc * na).^2)))
  )}

  CalcProbDistFunc <- function(lc,lmin = 1,lmax,D){
    ProbDist <- lc^(D-4) * (D-3) / (lmax^(D-3) - lmin^(D-3))
    return(ProbDist)
  }



  fs <- seq(2.1, 2.9, length.out = 20)

  sigmaLUT <- matrix(0, nrow = length(dfs), ncol = length(dfs))

  lmax_r <- seq(10, 10000, by = 20)

  }
# k is AverageWavelength
#
s2_to_int <- function(lc, L, na, k, d) {
  Fcoeff^2 * system_config$ri_definition$sigma_n^2 * (
    (2/pi * (L * (k^4 * lc^3) * na^2) /
       (1 + (k * lc)^2 * (4 + na^2)) /
       (1 + (k * lc)^2 * 4)
    ) + 1/4 * (1 - 1/sqrt(1 + (k * lc * na)^2))
  )
}


for idf = 1:length(dfs) %
    df = dfs(idf); % Db from the paper

    eqn2 = @(lmax_r) 6.*(df-3)./df .* (1 - (lmin./lmax_r).^df)./((lmin./lmax_r).^3 - (lmin./lmax_r).^df); %Eqn 2. from the Paper. Array of Nfs for corresponding values of lmax_r % This is equation 2 in the paper.
    lmax = S2D.utility.numericalInversion(eqn2, lmax_r, Nf);  % to get lmax corresponding to Nf. Nf = (lmax/lmin)^D cannot be inverted since we know Db and not D
    sigmaLUT(idf) = sqrt(integral(@(lc) (Pr(lc,lmin,lmax,df).* s2_to_int(lc, thick, system_config.center_lambda, ri_chromatin, system_config.na_c, lmin, lmax, df)), lmin, lmax));
end

    # Loop equivalent to for idf = 1:length(dfs)
    for (idf in 1:length(dfs)) {
      df <- dfs[idf]

      # Anonymous function equivalent to MATLAB's version
      eqn2 <- function(lmax_r) {
        6 * (df - 3) / df * (1 - (lmin / lmax_r)^df) / ((lmin / lmax_r)^3 - (lmin / lmax_r)^df)
      }

      # In R, you need a separate function for numerical inversion, as it's not a built-in function like in MATLAB
      # Assuming S2D.utility.numericalInversion is a custom function in MATLAB, you need an R equivalent
      lmax <- numericalInversion(eqn2, lmax_r, Nf) # Replace with the appropriate R function

      # For the integral, you can use the 'integrate' function in R
      # Assuming Pr and s2_to_int are already defined as functions in R
      integral_result <- integrate(function(lc) Pr(lc, lmin, lmax, df) * s2_to_int(lc, thick, system_config$center_lambda, ri_chromatin, system_config$na_c, lmin, lmax, df), lmin, lmax)$value

      sigmaLUT[idf] <- sqrt(integral_result)
    }

# Assuming sigmaLUT and dfs are already defined
# Fit a third-degree polynomial
fit <- lm(dfs ~ poly(sigmaLUT, 3, raw = TRUE))

# Extract coefficients
VSigma <- coef(fit)


# Assuming sigmaLUT and dfs are already defined
# Fit a third-degree polynomial
fit <- lm(dfs ~ poly(sigmaLUT, 3, raw = TRUE))

# Extract coefficients
VSigma <- coef(fit)

# Assuming the polynomial model is stored in 'fit'
dOut <- predict(fit, newdata = data.frame(sigmaLUT = sigmaIn))

lmin = 1 # minimal fractal distance that is possible

# No need to assign these as local variables, I will address later
# ri_glass = system_config.ri_definition.ri_glass;
# ri_chromatin = system_config.ri_definition.ri_chromatin;
# ri_media = system_config.ri_definition.ri_media;
# ri_immersion = system_config.immersion_ri;










