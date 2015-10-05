#' A function to calculate the variance ratio and null model mean and CIs within multiple replicates
#'
#' @param df A dataframe containing replicate.var, species.var, time.var and abundance.var columns
#' @param time.var The name of the time column from df
#' @param species.var The name of the species column from df
#' @param abundance.var The name of the abundance column from df
#' @param bootnumber The number of null model iterations used to calculated CIs
#' @param li The lower confidence interval, defaults to lowest 2.5\% CI
#' @param ui The upper confidence interval, defaults to upper 97.5\% CI  
#' @param replicate.var The name of the replicate column from df
#' @param average.replicates If true returns VR and CI averaged across reps; if false returns VR and CI for each rep
#' @return output A dataframe containing the replicate.var name, VR  nullVRCIlow, nullVRCIhigh and nullVRmean
#'          VR is the actual variance ratio
#'          nullVRCIow is the 0.025 CI 
#'          nullVRCIhigh is the 0.975 CI 
#'          nullVRmean is the mean variance ratio calculated on null communities
#' @export
varianceratio<-function(df, time.var="year", species.var="species",  abundance.var="abundance", bootnumber, replicate.var=NA,
                        li=0.025, ui=0.975,  average.replicates=TRUE){
  stopifnot(is.numeric(df[[time.var]]))
  if(is.na(replicate.var)==TRUE){
  VR<-varianceratio_longformdata(df,time.var, species.var, abundance.var)}else{
    df[replicate.var]<-if(is.factor(df[[replicate.var]])==TRUE){factor(df[[replicate.var]])} else {df[replicate.var]}
    if(average.replicates==TRUE){
      X<-split(df, df[replicate.var])
      VR<-mean(unlist(lapply(X, FUN=varianceratio_longformdata, time.var, species.var, abundance.var)))}else{
        X <- split(df, df[replicate.var])
        VR<-do.call("rbind", lapply(X, FUN=varianceratio_longformdata, time.var,  species.var,abundance.var))
      }
  }
nullout<-temporal_torus_translation_CI(df, time.var, species.var, abundance.var, varianceratio_matrixdata, bootnumber, replicate.var, li, ui, average.replicates)
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
varianceratio_matrixdata<-function(comdat){
  all.cov <- cov(comdat, use="pairwise.complete.obs")
  col.var<-apply(comdat, 2, var)
  com.var <-sum(all.cov)
  pop.var <-sum(col.var)
  var.ratio<-com.var/pop.var
  return(var.ratio)
}

#' A function to calculate the variance ratio from a longform dataframe
#'
#' @param df A dataframe containing time.var, replicate.var, species.var and abundance.var columns
#' @param time.var The name of the time.var column from df
#' @param species.var The name of the species.var column from df
#' @param abundance.var The name of the abundance.var column from df
#' @return var.ratio The variance ratio of the community
varianceratio_longformdata<-function(df, time.var, species.var, abundance.var){
  com.use<-transpose_community(df, time.var, species.var, abundance.var)
  var.ratio<-varianceratio_matrixdata(com.use)
  return(var.ratio)
}
