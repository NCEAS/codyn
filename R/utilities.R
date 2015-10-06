#' Convert from a longform abundance dataframe to a time by species dataframe.
#'
#' @param df A dataframe containing time.var, species.var and abundance.var columns
#' @param time.var The name of the time column from df
#' @param species.var The name of the species column from df
#' @param abundance.var The name of the abundance column from df
#' @return comdat A dataframe of species abundances x time
#' @export
transpose_community <- function(df, time.var, species.var, abundance.var) {
    df<-as.data.frame(df)
    df[species.var]<-if(is.factor(df[[species.var]])==TRUE){factor(df[[species.var]])} else {df[species.var]}  
    df<-df[order(df[time.var], df[species.var]),]
    comdat<-tapply(df[[abundance.var]], list(df[[time.var]], as.vector(df[[species.var]])), sum)
    comdat[is.na(comdat)]<-0
    comdat<-as.data.frame(comdat)
    return(comdat)
}

#' check names of data frames
#'
#' @param given Vector of variable names as supplied by user
#' @param data Data frame containing variables
check_names <- function(given, data) {
    for (i in given){
        assertthat::assert_that(assertthat::has_name(data, i))
    }
}

#' Utility function to ensure only a single record exists for a given species within one replicate, for one time point. 
#' @param df A dataframe containing time.var, species.var, and replicate.var columns
#' @param time.var The name of the time column from df
#' @param species.var The name of the species column from df
#' @param replicate The name of the replicate column from df

check_single <- function(df, time.var, species.var, replicate.var){
  X <- split(df, df[replicate.var]) 
  checksingle <- lapply(X, FUN = function(xx) apply(table(xx[[species.var]], xx[[time.var]]), 2, function(x) any(x>1)))    
  reptest <- unlist(lapply(checksingle, any))    
  yrtest <- lapply(checksingle, which)
  
  if(any(unlist(checksingle))){         
    stop(paste("In replicate(s)", names(reptest)[which(reptest)], "there are more than one record(s) for species at the time point(s)", unlist(lapply(yrtest, names))))
  }
}