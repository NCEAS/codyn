#' @title Temporal modification of the Torus translation
#' @description A function to calculate a vector of null test statistics using a temporal modification of the Torus translation
#' 
#' @param df A dataframe containing time, species and abundance columns
#' @param species.var The name of the species column from df
#' @param time.var The name of the time column from df
#' @param abundance.var The name of the abundance column from df
#' @param FUN A function to calculate on the null community
#' @return a vector of null test statistics calculated from a randomized community matrix in which species autocorrelation has been maintained via a Torus translation
#' @export
temporal_torus_translation<-function(df, time.var, species.var,  abundance.var, FUN){
  out<-FUN(genRand(transpose_community(df, time.var,  species.var, abundance.var)))
  bootout<-unlist(out)
  return(bootout)
}


#' A function that returns confidence intervals calculated from a temporal modification of the torus translation
#'
#' @param df A dataframe containing time, species and abundance columns
#' @param time.var The name of the time column from df
#' @param species.var The name of the species column from df
#' @param abundance.var The name of the abundance column from df
#' @param FUN A function to calculate on the null community
#' @param bootnumber The number of null model iterations used to calculated CIs
#' @param replicate.var The name of the replication column from df
#' @param li The lower confidence interval, defaults to lowest 2.5\% CI
#' @param ui The upper confidence interval, defaults to upper 97.5\% CI  
#' @param average.replicates If true returns CI averaged across reps; if false returns CI for each rep
#' @return output A dataframe containing lowerCI, upperCI and nullmean value of the test statistic
#' @export
temporal_torus_translation_CI<-function(df,  time.var="year",species.var="species", abundance.var="abundance", FUN, bootnumber, replicate.var=NA, li=0.025, ui=0.975, average.replicates=TRUE){
  if(is.na(replicate.var)==TRUE){
    out<-replicate(bootnumber, temporal_torus_translation(df, time.var, species.var,  abundance.var, FUN))
    lowerCI <- quantile(out, li)
    upperCI <-quantile(out, ui)
    nullmean<-mean(out)
    output<-cbind(lowerCI, upperCI, nullmean)
    row.names(output)<-NULL
  }else{
    df[replicate.var]<-if(is.factor(df[[replicate.var]])==TRUE){factor(df[[replicate.var]])} else {df[replicate.var]}
    df<-df[order(df[[replicate.var]]),]  
    if(average.replicates==TRUE){
    X<-split(df, df[replicate.var])
    out<-replicate(bootnumber, mean(unlist(lapply(X, temporal_torus_translation, time.var, species.var, abundance.var, FUN)))) 
    lowerCI <- quantile(out, li)
    upperCI <-quantile(out, ui)
    nullmean<-mean(out)
    output<-cbind(lowerCI, upperCI, nullmean)
    row.names(output)<-NULL
  } else{
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

