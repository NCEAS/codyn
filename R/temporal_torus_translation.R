#' @title Temporal Modification of the Torus Translation
#' @description Calculates a user-specified test statistic on a null ecological community using a temporal modification of the Torus translation.
#' 
#' @param df A data frame containing time, species and abundance columns and an optional column of replicates
#' @param time.var The name of the time column 
#' @param species.var The name of the species column 
#' @param abundance.var The name of the abundance column 
#' @param FUN A function to calculate on the null community
#' @return The temporal_torus_translation function returns an integer output of the user-specified function calculated on a null community. 
#' @details The input data frame needs to contain columns for time, species and abundance; time.var, species.var and abundance.var are used to indicate which columns contain those variables.
#' @export
temporal_torus_translation <- function(df, time.var="year", species.var="species",  abundance.var="abundance", FUN){
  if(!is.numeric(df[,abundance.var])) { stop("Abundance variable is not numeric") }
  out<-FUN(genRand(transpose_community(df, time.var,  species.var, abundance.var)))
  bootout<-unlist(out)
  return(bootout)
}

#' @title Confidence Intervals Using a Modification of the Torus Translation
#' @description Calculates a user-specified test statistic on a null ecological community using a temporal modification of the Torus translation.
#' 
#' @param df A data frame containing time, species and abundance columns and an optional column of replicates
#' @param time.var The name of the time column 
#' @param species.var The name of the species column 
#' @param abundance.var The name of the abundance column 
#' @param FUN A function to calculate on the null community
#' @param bootnumber The number of null model iterations used to calculated confidence intervals
#' @param replicate.var The name of the replication column from df
#' @param li The lower confidence interval, defaults to lowest 2.5\% CI
#' @param ui The upper confidence interval, defaults to upper 97.5\% CI  
#' @param average.replicates If true returns the CIs averaged across replicates; if false returns the CI for each replicate
#' @return The temporal_torus_translation_CI function returns a dataframe with the following attributes:
#' \itemize{
#'  \item{lowerCI: }{A numeric column with the lowest confidence interval value.}
#'  \item{upperCI: }{A numeric column with the highest confidence interval value.}
#'  \item{nullMean: }{A numeric column with the average value of the specified function when calculated on a null community.}
#'  \item{replicate.var: }{A column that has same name and type as the replicate.var column, if replication is specified.}
#' }
#' @details
#' The input data frame needs to contain columns for time, species and abundance; time.var, species.var and abundance.var are used to indicate which columns contain those variables.
#' If multiple replicates are included in the data frame, that column should be specified with replicate.var. Each replicate should reflect a single experimental unit - there must be a single abundance value per species within each time point and replicate.
#' Null model confidence intervals default to the standard lowest 2.5\% and upper 97.5\% of the null distribution, typically these do not need to be change, but they can be user-modified to set more stringent CIs.
#' @export
temporal_torus_translation_CI<-function(df,  time.var="year",species.var="species", abundance.var="abundance", FUN, bootnumber, replicate.var=NA, li=0.025, ui=0.975, average.replicates=TRUE){
  if(is.na(replicate.var)){
    check_single_onerep(df, time.var, species.var)
    
    out<-replicate(bootnumber, temporal_torus_translation(df, time.var, species.var,  abundance.var, FUN))
    lowerCI <- quantile(out, li)
    upperCI <-quantile(out, ui)
    nullmean<-mean(out)
    output<-cbind(lowerCI, upperCI, nullmean)
    row.names(output)<-NULL
  } else {
    df[replicate.var]<-if(is.factor(df[[replicate.var]])==TRUE){factor(df[[replicate.var]])} else {df[replicate.var]} 
    
    df<-df[order(df[[replicate.var]]),]  
    check_single(df, time.var, species.var, replicate.var)
    if(average.replicates==TRUE){
      X<-split(df, df[replicate.var])
    out<-replicate(bootnumber, mean(unlist(lapply(X, temporal_torus_translation, time.var, species.var, abundance.var, FUN)))) 
    lowerCI <- quantile(out, li)
    upperCI <-quantile(out, ui)
    nullmean<-mean(out)
    output<-cbind(lowerCI, upperCI, nullmean)
    row.names(output)<-NULL
  } else {
    X <- split(df, df[replicate.var])
    out<-lapply(X, function(x,time.var, species.var,  abundance.var, FUN, bootnumber){replicate(bootnumber, temporal_torus_translation(x, time.var, species.var, abundance.var, FUN))}, time.var, species.var, abundance.var, FUN, bootnumber)
    lowerCI<-do.call("rbind", lapply(out, quantile, li))
    upperCI<-do.call("rbind", lapply(out, quantile, ui))
    nullmean<-do.call("rbind", lapply(out, mean))
    reps<-unique(df[replicate.var])
    output<-cbind(reps, lowerCI, upperCI, nullmean)
    names(output)[2:3]=c("lowerCI", "upperCI")
  }
  }
  return(as.data.frame(output))
}




############################################################################
#
# Private functions: these are internal functions not intended for reuse.  
# Future package releases may change these without notice. External callers
# should not use them.
#
############################################################################



#' A function to generate a community dataframe with a random start time for each species
#'
#' @param comdat A community dataframe
#' @return rand.use A randomized community dataframe
genRand<-function(comdat){
  comdat2<-rbind(comdat, comdat)
  rand.comdat<-matrix(NA, nrow(comdat), ncol(comdat)) 
  for(i in 1:ncol(comdat)){  
    rand.start<-sample(1:nrow(comdat), 1)
    rand.comdat[,i]<-comdat2[rand.start:(rand.start+nrow(comdat)-1), i]
  }
  rand.use<-rand.comdat[1:nrow(rand.comdat), 1:ncol(rand.comdat)]
  return(rand.use)
}

