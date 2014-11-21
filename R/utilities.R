#' A function to covert from a longform abundance dataframe to a year x species dataframe
#'
#' @param data1 A dataframe containing year, species and abundance columns
#' @param year The name of the year column from data1
#' @param species The name of the species column from data1
#' @param abundance The name of the abundance column from data1
#' @return comdat A dataframe of species abundances x year
calComDat<-function(data1, species, year, abundance){
  data1<-data1[order(data1[year], data1[species]),]
  comdat<-tapply(data1[[abundance]], list(data1[[year]], as.vector(data1[[species]])), sum)
  comdat[is.na(comdat)]<-0
  comdat<-as.data.frame(comdat)
  return(comdat)
}


#' A function to covert from a longform abundance dataframe to a species x year time series
#'
#' @param data1 A dataframe containing year, species and abundance columns
#' @param year The name of the year column from data1
#' @param species The name of the species column from data1
#' @param abundance The name of the abundance column from data1
#' @return comdat A dataframe of species abundances x year
calComTS<-function(data1, species, year, abundance){
  data1<-data1[order(data1[year], data1[species]),]
  comdat<-(tapply(data1[[abundance]], list(data1[[species]], as.vector(data1[[year]])), sum))
  comdat[is.na(comdat)]<-0  
  comdat<-comdat[which(rowSums(comdat)>0),]
  start_time <- min(as.numeric(colnames(comdat)))
  end_time <- max(as.numeric(colnames(comdat)))
  comdat <- as.ts(comdat, start = start_time, end = end_time)
  return(comdat)
}
