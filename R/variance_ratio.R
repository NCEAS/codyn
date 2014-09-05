library(reshape)


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
#' @param rep The name of the replicate column from data1
#' @param species The name of the species column from data1
#' @param abundance The name of the abundance column from data1
#' @return var.ratio The variance ratio of the community
calVR2<-function(data1, species, year, abundance){
  com.use<-calComDat(data1, species, year, abundance)
  var.ratio<-calVR(com.use)
  return(var.ratio)
}

#' A function to calculate the variance ratio within multiple replicates
#'
#' calVRx a dataframe of replicates with their variance ratio values
#' @param data1 A dataframe containing year, rep, species and abundance columns
#' @param year The name of the year column from data1
#' @param rep The name of the replicate column from data1
#' @param species The name of the species column from data1
#' @param abundance The name of the abundance column from data1
#' @return output A dataframe containing rep and VR (the variance ratio)
#' @export
calVRs<-function(data1, rep, species, year, abundance){
  X <- split(data1, data1[rep])
  out<-lapply(X, FUN=calVR2, species, year, abundance)
  output<-cbind((names(out)), as.data.frame(unlist(out)))
  names(output)<-c(rep, "VR")
  row.names(output)<-NULL
  return(output)
}

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
  rand.use<-rand.comdat[1:nrow(rand.comdat), 2:ncol(rand.comdat)]
  return(rand.use)
}
