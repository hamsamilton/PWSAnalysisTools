setwd("/Users/samhamilton/Library/Mobile Documents/com~apple~CloudDocs/Thesis_Aim1/libraries/PWSAnalysisTools")
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
           Image_ROI             = str_c(Image,`ROI Name`,ROI),
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
           title = Cond) +facet.themes
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
