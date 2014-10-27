#' A function to calculate species turnover between two years 
#'
#' @param data1 A dataframe containing year, rep, species and abundance columns
#' @param year The name of the year column from data1
#' @param species The name of the species column from data1
#' @param abundance The name of the abundance column from data1
#' @param metric The turnover metric to return; the default, turnover, returns summed appearances and disappearances relative to total species richness across both years
#'          appearance returns the number of appearances in the second year relative to total species richness across both years
#'          disappearance returns the number of disappearances in the second year relative to the total species richness across both years
#' @return output The specificed turnover metric
#' @import reshape
#' @export
pairwise_turnover<-function(data1,  species, year, abundance, metric="turnover"){
  d1<-data1[which(data1["abundance"]>0),]
  fstr<-(paste(species, "~", year, sep=""))
  f<-as.formula(fstr)
  d2<-as.data.frame(cast(d1, f, value=abundance, fill=0))
  d2["disapp"]<-ifelse(d2[2]>0 & d2[3]==0, 1, 0)
  d2["app"]<-ifelse(d2[2]==0 & d2[3]>0, 1, 0)
  disapp<-sum(d2["disapp"])
  app<-sum(d2["app"])
  sppchange<-sum(disapp, app)
  totspp<-as.numeric(nrow(d2))
  reldis<-disapp/totspp
  relapp<-app/totspp
  turnover<-sppchange/totspp
  if(metric=="turnover"){
    output<-turnover} else {
      if(metric=="appearance"){
        output<-relapp} else{
          if(metric=="disappearance"){
            output<-reldis
          }
        }
    }
return(output)
}


#' A function to calculate species turnover between years
#'
#' @param data1 A dataframe containing year,species and abundance columns
#' @param year The name of the year column from data1
#' @param species The name of the species column from data1
#' @param abundance The name of the abundance column from data1
#' @param metric The turnover metric to return; the default, turnover, returns summed appearances and disappearances relative to total species richness across both years
#'          appearance returns the number of appearances in the second year relative to total species richness across both years
#'          disappearance returns the number of disappearances in the second year relative to the total species richness across both years
#' @return output A dataframe containing the specificed turnover metric and year
#' @import reshape
#' @export
allyear_turnover<-function(data1, species, year, abundance, metric="turnover"){
  data1<-data1[order(data1[year]),]
  ## split data by year
  yearlist <- split(data1, data1[[year]])
  ## create consecutive pairs of years
  y1 <- yearlist[-length(yearlist)]
  y2 <- yearlist[-1]
  ## rbind consecutive pairs of years
  yearpair <- Map(function(d1, d2){rbind(d1,d2)}, y1, y2)
  ## calculate turnover for across all years
  out<-lapply(yearpair, FUN=pairwise_turnover, species, year, abundance, metric)
  output<-as.data.frame(do.call("rbind", out))
  names(output)[1]=metric
  ## add year column
  allyr<-as.matrix(unique(data1[year]))
  currentyr<-allyr[2:nrow(allyr)]
  output[year]=(currentyr)
  return(output)
}


