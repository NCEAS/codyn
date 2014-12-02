#' A function to calculate community stability over time
#' @param replicate The name of the replicate column from data1
#' @param year The name of the year column from data1
#' @param species The name of the species column from data1
#' @param abundance The name of the abundance column from data1
#' @return output The stability of community abundance, calculated as mean/standard deviation
#' @export
community_stability<-function(data1, replicate="replicate", year="year", abundance="abundance"){
 if(is.na(replicate)==FALSE){
  #sum abundance within a replicate and year
  aggform<-as.formula(paste(abundance, "~", replicate, "+", year, sep=""))
  data2<-aggregate(aggform, data=data1, sum)
  X<-split(data2, data2[replicate])
  out<-lapply(X, stability_onerep, abundance)
  reps<-unique(data2[replicate])
  output<-cbind(reps, do.call("rbind", out))
  names(output)=c(replicate, "stability")
  output<-subset(output, is.na(stability)==F)} 
  else{
    #sum abundance within a year
    aggform<-as.formula(paste(abundance, "~", year, sep=""))
    data2<-aggregate(aggform, data=data1, sum)
    output<-stability_onerep(data2, abundance)
  }
  row.names(output)<-NULL
 return(output)
}


############################################################################
#
# Private functions: these are internal functions not intended for reuse.  
# Future package releases may change these without notice. External callers
# should not use them.
#
############################################################################



#' A function to calculate species synchrony over time within one replicate
#'
#' @param data1 A dataframe containing x column
#' @param x The column to calculate stability on
#' @return Stability of x, calculated as the mean/sd
stability_onerep<-function(data1, x){
  return(mean(data1[[x]])/sd(data1[[x]]))
}

