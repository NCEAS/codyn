#' A function to calculate a vector of null test statistics using a temporal modification of the Torus translation
#'
#' @param data1 A dataframe containing year, species and abundance columns
#' @param species The name of the species column from data1
#' @param year The name of the year column from data1
#' @param abundance The name of the abundance column from data1
#' @param FUN A function to calculate on the null community
#' @return a vector of null test statistics calculated from a randomized community matrix in which species autocorrelation has been maintained via a Torus translation
#' @export
temporal_torus_translation<-function(data1, species, year, abundance, FUN){
  out<-FUN(genRand(calComDat(data1, species, year, abundance)))
  bootout<-unlist(out)
  return(bootout)
}


#' A function that returns confidence intervals calculated from a temporal modification of the torus translation
#'
#' @param data1 A dataframe containing year, species and abundance columns
#' @param replicate The name of the replication column from data1
#' @param species The name of the species column from data1
#' @param year The name of the year column from data1
#' @param abundance The name of the abundance column from data1
#' @param FUN A function to calculate on the null community
#' @param bootnumber The number of null model iterations used to calculated CIs
#' @param li The lower confidence interval, defaults to lowest 2.5\% CI
#' @param ui The upper confidence interval, defaults to upper 97.5\% CI  
#' @param averagereps If true returns CI averaged across reps; if false returns CI for each rep
#' @return output A dataframe containing lowerCI, upperCI and nullmean value of the test statistic
#' @export
temporal_torus_translation_CI<-function(data1, replicate="replicate", species="species", year="year", abundance="abundance", FUN, bootnumber, li=0.025, ui=0.975, averagereps=TRUE){
  if(is.na(replicate)==TRUE){
    out<-replicate(bootnumber, temporal_torus_translation(data1, species, year, abundance, FUN))
    lowerCI <- quantile(out, li)
    upperCI <-quantile(out, ui)
    nullmean<-mean(out)
    output<-cbind(lowerCI, upperCI, nullmean)
    row.names(output)<-NULL
  }else{
  if(averagereps==TRUE){
    X<-split(data1, data1[replicate])
    out<-replicate(bootnumber, mean(unlist(lapply(X, temporal_torus_translation, species, year, abundance, FUN)))) 
    lowerCI <- quantile(out, li)
    upperCI <-quantile(out, ui)
    nullmean<-mean(out)
    output<-cbind(lowerCI, upperCI, nullmean)
    row.names(output)<-NULL
  } else{
    X <- split(data1, data1[replicate])
    out<-lapply(X, function(x, species, year, abundance, FUN, bootnumber){replicate(bootnumber, temporal_torus_translation(x, species, year, abundance, FUN))}, species, year, abundance, FUN, bootnumber)
    lowerCI<-do.call("rbind", lapply(out, quantile, li))
    upperCI<-do.call("rbind", lapply(out, quantile, ui))
    nullmean<-do.call("rbind", lapply(out, mean))
    reps<-unique(data1[replicate])
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

comdat <- calComDat(data1 = dat1, "species", "year", "abundance")


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

