#' A function to calculate the variance ratio and null model mean and CIs within multiple replicates
#'
#' @param data1 A dataframe containing replicate, species, year and abundance columns
#' @param replicate The name of the replicate column from data1
#' @param species The name of the species column from data1
#' @param year The name of the year column from data1
#' @param abundance The name of the abundance column from data1
#' @param bootnumber The number of null model iterations used to calculated CIs
#' @param averagereps If true returns VR and CI averaged across reps; if false returns VR and CI for each rep
#' @param li The lower confidence interval, defaults to lowest 2.5\% CI
#' @param ui The upper confidence interval, defaults to upper 97.5\% CI  
#' @return output A dataframe containing the replicate name, VR  nullVRCIlow, nullVRCIhigh and nullVRmean
#'          VR is the actual variance ratio
#'          nullVRCIow is the 0.025 CI 
#'          nullVRCIhigh is the 0.975 CI 
#'          nullVRmean is the mean variance ratio calculated on null communities
#' @export
varianceratio<-function(data1, replicate="replicate", species="species", year="year", abundance="abundance",  bootnumber, li=0.025, ui=0.975, averagereps=TRUE){
if(is.na(replicate)==TRUE){
  VR<-calVR_longformdata(data1, species, year, abundance)}else{
    data1[replicate]<-if(is.factor(data1[[replicate]])==TRUE){factor(data1[[replicate]])} else {data1[replicate]}
    if(averagereps==TRUE){
      X<-split(data1, data1[replicate])
      VR<-mean(unlist(lapply(X, FUN=calVR_longformdata, species, year, abundance)))}else{
        X <- split(data1, data1[replicate])
        VR<-do.call("rbind", lapply(X, FUN=calVR_longformdata, species, year, abundance))
      }
  }
nullout<-temporal_torus_translation_CI(data1, replicate, species, year, abundance,calVR, bootnumber, li, ui, averagereps)
output<-(cbind(nullout, VR))
row.names(output)<-NULL
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
  com.use<-transpose_community(data1, species, year, abundance)
  var.ratio<-calVR(com.use)
  return(var.ratio)
}
