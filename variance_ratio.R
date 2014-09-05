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

#' A function to calculate the variance ratio within multiple replicates
#'
#' calrealVR a dataframe of replicates with their variance ratio values
#' @param data1 A dataframe containing year, sitesubplot, species and abundance columns
#' @param year The name of the year column from data1
#' @param sitesubplot The name of the replicate column from data1
#' @param species The name of the species column from data1
#' @param abundance The name of the abundance column from data1
#' @return output A dataframe containing sitesubplot and VR (the variance ratio)
#' @export


calrealVR <- function(data1, year, sitesubplot, species, abundance){
  plotnames<-as.numeric (unique(data1$sitesubplot))
  site<-data.frame(cbind(VR=NA, sitesubplot=NA))
  for(i in 1:length(plotnames)){
    current.plot<-plotnames[i]
    subber<-data1[as.numeric(data1$sitesubplot)==current.plot, ]
    melt.all <-melt(subber, id=c("species", "year", "abundance"))
    melt.all$value <- melt.all$variable <-NULL
    comdat<-as.data.frame(cast(melt.all, year ~ species, value="abundance", fill=0))
    com.use<-comdat
    com.use$year<-NULL
    var.ratio<-calVR(com.use)
    calvarbind <-data.frame(VR=var.ratio)
    calvarbind$sitesubplot<-unique(subber$sitesubplot)
    site <-rbind(site, calvarbind)
    output<-site[which(is.na(site$VR)==FALSE),]
  }
  return(output)
}