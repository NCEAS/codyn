#' A function to calculate the variance ratio 
#'
#' @param comdat A community dataframe
#' @return var.ratio The variance ratio of the community
#' @export
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
#' @export
calVR2<-function(data1, species, year, abundance){
  com.use<-calComDat(data1, species, year, abundance)
  var.ratio<-calVR(com.use)
  return(var.ratio)
}

#' A function to generate a community dataframe with a random start time for each species
#'
#' @param comdat A community dataframe
#' @return rand.use A randomized community dataframe
#' @export
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

#' A function to calculate a null variance ratio from longform data using a temporal modification of the Torus translation
#'
#' @param data1 A dataframe containing year, rep, species and abundance columns
#' @param species The name of the species column from data1
#' @param year The name of the year column from data1
#' @param abundance The name of the abundance column from data1
#' @return randVR A variance ratio calculated from a randomized community matrix in which species autocorrelation has been maintained via a Torus translation
#' @export
calnullVR<-function(data1, species, year, abundance){
  comdat<-calComDat(data1, species, year, abundance)
  rand.dat<-genRand(comdat)
  var.ratio<-calVR(rand.dat)
  randVR <-data.frame(VR=var.ratio)
  return(randVR)
}

#' A function to generate lower 2.5\% CI, upper 97.5\% CI and mean null VR values
#'
#' @param data1 A dataframe containing year, rep, species and abundance columns
#' @param species The name of the species column from data1
#' @param year The name of the year column from data1
#' @param abundance The name of the abundance column from data1
#' @param bootnumber The number of null model iterations used to calculated CIs
#' @return output A dataframe containing nullVRCIlow, nullVRCIhigh and nullVRmean
#'          nullVRCIow is the 0.025 CI and nullVRCIhigh is the 0.975 CI
#' @export
nullVRCI<-function(data1, species, year, abundance, bootnumber){
  out<-replicate(bootnumber, calnullVR(data1, species, year, abundance))
  bootout<-(as.data.frame(unlist(out)))
  names(bootout)<-c("nullVR")
  nullVRlow <- quantile(bootout$nullVR, (.025))
  nullVRhigh<-quantile(bootout$nullVR, 0.975)
  nullVRmean<-mean(bootout$nullVR)
  output<-cbind(nullVRlow, nullVRhigh, nullVRmean)
  row.names(output)<-NULL
  return(output)
}

#' A function to calculate both the real and mean null variance ratio along with lower 2.5\% CI, upper 97.5\% CI using a temporal modification of the Torus translation
#'
#' @param data1 A dataframe containing year, rep, species and abundance columns
#' @param species The name of the species column from data1
#' @param year The name of the year column from data1
#' @param abundance The name of the abundance column from data1
#' @param bootnumber The number of null model iterations used to calculated CIs
#' @return output A dataframe containing VR  nullVRCIlow, nullVRCIhigh and nullVRmean
#'          VR is the actual variance ratio
#'          nullVRCIow is the 0.025 CI 
#'          nullVRCIhigh is the 0.975 CI 
#'          nullVRmean is the mean variance ratio calculated on null communities
#' @export
calVRrealnull<-function(data1, species, year, abundance, bootnumber){
  VR<-calVR2(data1, species, year, abundance)
  nullVR<-nullVRCI(data1, species, year, abundance, bootnumber)
  out<-cbind(VR, nullVR)
  return(out)
}

#' A function to calculate the variance ratio and null model mean and CIs within multiple replicates
#'
#' @param data1 A dataframe containing rep, species, year and abundance columns
#' @param rep The name of the replicate column from data1
#' @param species The name of the species column from data1
#' @param year The name of the year column from data1
#' @param abundance The name of the abundance column from data1
#' @param bootnumber The number of null model iterations used to calculated CIs
#' @param averagereps If true returns VR and CI averaged across reps; if false returns VR and CI for each rep
#'          If true, null VR are calculated within each rep, averaged, and the repeated for length of bootnumber
#' @return output A dataframe containing the replicate name, VR  nullVRCIlow, nullVRCIhigh and nullVRmean
#'          VR is the actual variance ratio
#'          nullVRCIow is the 0.025 CI 
#'          nullVRCIhigh is the 0.975 CI 
#'          nullVRmean is the mean variance ratio calculated on null communities
#' @export
VR<-function(data1, rep, species, year, abundance, bootnumber, averagereps=TRUE){
  if(averagereps==TRUE){
    X<-split(data1, data1[rep])
    out<-replicate(bootnumber, mean(unlist(lapply(X, FUN=calnullVR, species, year, abundance)))) 
    nullVRlow <- quantile(out, (.025))
    nullVRhigh<-quantile(out, 0.975)
    nullVRmean<-mean(out)
    VR<-mean(unlist(lapply(X, FUN=calVR2, species, year, abundance)))
    output<-cbind(VR, nullVRlow, nullVRhigh, nullVRmean)
    row.names(output)<-NULL
  } else{
    X <- split(data1, data1[rep])
    out<-lapply(X, FUN=calVRrealnull, species, year, abundance, bootnumber)
    reps<-unique(data1[rep])
    output<-cbind(reps, do.call("rbind", out))
  }
  return(output)
}
