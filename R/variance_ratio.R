#' A function to calculate the variance ratio and null model mean and CIs within multiple replicates
#'
#' @param data1 A dataframe containing rep, species, year and abundance columns
#' @param rep The name of the replicate column from data1
#' @param species The name of the species column from data1
#' @param year The name of the year column from data1
#' @param abundance The name of the abundance column from data1
#' @param bootnumber The number of null model iterations used to calculated CIs
#' @param averagereps If true returns VR and CI averaged across reps; if false returns VR and CI for each rep
#' @param li The lower confidence interval, defaults to lowest 2.5\% CI
#' @param ui The upper confidence interval, defaults to upper 97.5\% CI  
#'          If true, null VR are calculated within each rep, averaged, and the repeated for length of bootnumber
#' @return output A dataframe containing the replicate name, VR  nullVRCIlow, nullVRCIhigh and nullVRmean
#'          VR is the actual variance ratio
#'          nullVRCIow is the 0.025 CI 
#'          nullVRCIhigh is the 0.975 CI 
#'          nullVRmean is the mean variance ratio calculated on null communities
#' @export
varianceratio<-function(data1, rep="rep", species="species", year="year", abundance="abundance", bootnumber, li=0.025, ui=0.975, averagereps=TRUE){
  if(averagereps==TRUE){
    X<-split(data1, data1[rep])
    out<-replicate(bootnumber, mean(unlist(lapply(X, FUN=temporal_torus_translation, species, year, abundance, calVR)))) 
    lowerCI <- quantile(out, li)
    upperCI <-quantile(out, ui)
    nullmean<-mean(out)
    VR<-mean(unlist(lapply(X, FUN=calVR2, species, year, abundance)))
    output<-cbind(VR, lowerCI, upperCI, nullmean)
    row.names(output)<-NULL
  } else{
    X <- split(data1, data1[rep])
    out<-lapply(X, FUN=calVRrealnull, species, year, abundance, bootnumber)
    reps<-unique(data1[rep])
    output<-cbind(reps, do.call("rbind", out))
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

#' A function to calculate the variance ratio 
#'
#' @param comdat A community dataframe
#' @return var.ratio The variance ratio of the community
calVR<-function(comdat){
  all.cov <- cov(comdat, use="pairwise.complete.obs")
  col.var<-apply(comdat, 2, var)
  com.var <-sum(all.cov)
  pop.var <-sum(col.var)
  var.ratio<-com.var/pop.var
  return(var.ratio)
}

#' A function to calculate the variance ratio from a longform dataframe
#'
#' @param data1 A dataframe containing year, rep, species and abundance columns
#' @param year The name of the year column from data1
#' @param species The name of the species column from data1
#' @param abundance The name of the abundance column from data1
#' @return var.ratio The variance ratio of the community
calVR_longformdata<-function(data1, species, year, abundance){
  com.use<-calComDat(data1, species, year, abundance)
  var.ratio<-calVR(com.use)
  return(var.ratio)
}

#' A function to calculate both the real and mean null variance ratio along with lower 2.5\% CI, upper 97.5\% CI using a temporal modification of the Torus translation
#'
#' @param data1 A dataframe containing year, species and abundance columns
#' @param species The name of the species column from data1
#' @param year The name of the year column from data1
#' @param abundance The name of the abundance column from data1
#' @param bootnumber The number of null model iterations used to calculated CIs
#' @param li The lower confidence interval, defaults to lowest 2.5\% CI
#' @param ui The upper confidence interval, defaults to upper 97.5\% CI  
#' @return output A dataframe 
#'          VR is the actual variance ratio
#'          lowerCI defaults to the 0.025 CI 
#'          upperCI defaults to the 0.975 CI 
#'          nullmean is the mean variance ratio calculated on null communities
calVRrealnull<-function(data1, species, year, abundance, bootnumber, li=0.025, ui=0.975){
  VR<-calVR_longformdata(data1, species, year, abundance)
  nullVR<-temporal_torus_translation_CI(data1, species, year, abundance, FUN=calVR, bootnumber=bootnumber, li=li, ui=ui)
  out<-cbind(VR, nullVR)
  return(out)
}

