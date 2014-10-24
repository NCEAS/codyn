#' A function to covert from a longform abundance dataframe to a year x species dataframe
#'
#' @param data1 A dataframe containing year, species and abundance columns
#' @param year The name of the year column from data1
#' @param species The name of the species column from data1
#' @param abundance The name of the abundance column from data1
#' @return comdat A dataframe of species abundances x year
#' @import reshape
calComDat<-function(data1, species, year, abundance){
  data1<-data1[order(data1[year]),]
  fstr<-(paste(year, "~", species, sep=""))
  f<-as.formula(fstr)
  comdat<-as.data.frame(cast(data1, f, value=abundance, fill=0))
  comdat[year]<-NULL
  return(comdat)
}


#' A function to covert from a longform abundance dataframe to a species x year time series
#'
#' @param data1 A dataframe containing year, species and abundance columns
#' @param year The name of the year column from data1
#' @param species The name of the species column from data1
#' @param abundance The name of the abundance column from data1
#' @return comdat A dataframe of species abundances x year
#' @import reshape
calComTS<-function(data1, species, year, abundance){
  data1 <- data1[order(data1[year]), ]
  fstr <- paste(species, "~", year, sep="")
  f <- as.formula(fstr)
  comdat <- as.data.frame(cast(data1, f, value = abundance, fill=0))
  rownames(comdat) <- comdat[[species]]
  comdat[species] <- NULL
  comdat <- as.matrix(comdat)
  start_time <- min(as.numeric(colnames(comdat)))
  end_time <- max(as.numeric(colnames(comdat)))
  comdat <- as.ts(comdat, start = start_time, end = end_time)
  return(comdat)
}
