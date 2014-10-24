library(reshape)
#' A function to covert from a longform abundance dataframe to a year x species dataframe
#'
#' @param data1 A dataframe containing year, species and abundance columns
#' @param year The name of the year column from data1
#' @param species The name of the species column from data1
#' @param abundance The name of the abundance column from data1
#' @return comdat A dataframe of species abundances x year
calComDat<-function(data1, species, year, abundance){
  data1<-data1[order(data1[year]),]
  fstr<-(paste(year, "~", species, sep=""))
  f<-as.formula(fstr)
  comdat<-as.data.frame(cast(data1, f, value=abundance, fill=0))
  comdat[year]<-NULL
  return(comdat)
}
