#' A function to calculate species turnover between two years within a rep
#'
#' @param data1 A dataframe containing year, rep, species and abundance columns
#' @param year The name of the year column from data1
#' @param rep The name of the replicate column from data1
#' @param species The name of the species column from data1
#' @param abundance The name of the abundance column from data1
#' @return turnover The species turnover (appearances + disappearances) relative to the total species richness observed across both years
#' @import reshape
#' @export
pairwise_turnover<-function(data1, rep, species, year, abundance){
  d1<-data1[which(data1["abundance"]>0),]
  fstr<-(paste(species, "+", rep, "~", year, sep=""))
  f<-as.formula(fstr)
  d2<-as.data.frame(cast(d1, f, value=abundance, fill=0))
  d2["disapp"]<-ifelse(d2[3]>0 & d2[4]==0, 1, 0)
  d2["app"]<-ifelse(d2[3]==0 & d2[4]>0, 1, 0)
  disapp<-sum(d2["disapp"])
  app<-sum(d2["app"])
  sppchange<-sum(disapp, app)
  totspp<-as.numeric(nrow(d2))
  reldis<-disapp/totspp
  relapp<-app/totspp
  turnover<-sppchange/totspp
  return(turnover)
}

