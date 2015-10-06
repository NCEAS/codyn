#' A function to calculate species turnover between years
#'
#' @param df A dataframe containing year,species and abundance columns
#' @param time.var The name of the year column from df
#' @param species.var The name of the species column from df
#' @param abundance.var The name of the abundance column from df
#' @param replicate.var The name of the replicate column from df
#' @param metric The turnover metric to return; the default, total, returns summed appearances and disappearances relative to total species richness across both years
#'          appearance returns the number of appearances in the second year relative to total species richness across both years
#'          disappearance returns the number of disappearances in the second year relative to the total species richness across both years
#' @return output A dataframe containing the specificed turnover metric and year
#' @export
#turnover<-function(data1, replicate="replicate", species="species", year="year", abundance="abundance", metric="total") {
turnover<-function(df, time.var="year", species.var="species", abundance.var="abundance", replicate.var=as.character(NA), metric="total") {
  if(is.na(replicate.var)==TRUE){
    output<-turnover_allyears(df, species.var, time.var, abundance.var)}else{
      df[replicate.var]<-if(is.factor(df[[replicate.var]])==TRUE){factor(df[[replicate.var]])} else {df[replicate.var]}
      X <- split(df, df[replicate.var])
      out<-lapply(X, FUN=turnover_allyears, species.var, time.var, abundance.var, metric)
      ID<-unique(names(out))
      out<-mapply(function(x, y) "[<-"(x, replicate.var, value = y) ,
                  out, ID, SIMPLIFY = FALSE)
      output<-do.call("rbind", out)
    }
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


#' A function to calculate species turnover between two years 
#'
#' @param d1 A dataframe containing a species column from one year
#' @param d2 A dataframe containing a species column from the following year
#' @param species The name of the species column in d1 and d2
#' @param metric The turnover metric to return; the default, total, returns summed appearances and disappearances relative to total species richness across both years
#'          appearance returns the number of appearances in the second year relative to total species richness across both years
#'          disappearance returns the number of disappearances in the second year relative to the total species richness across both years
#' @return output The specificed turnover metric
getturnover <- function(d1, d2, species = "species", metric="total"){
  d1spp<-as.character(unique(d1[[species]]))
  d2spp<-as.character(unique(d2[[species]]))
  commspp<-intersect(d1spp, d2spp)
  disappear<-length(d1spp)-length(commspp)
  appear<-length(d2spp)-length(commspp)
  totrich<-sum(disappear, appear, length(commspp))
  if(metric=="total"){
    output<-((appear+disappear)/totrich)} else {
      if(metric=="appearance"){
        output<-appear/totrich} else{
          if(metric=="disappearance"){
            output<-disappear/totrich
          }
        }
    }
  return(output)
}



#' A function to calculate species turnover between years
#'
#' @param df A dataframe containing year,species and abundance columns
#' @param species The name of the species column from df
#' @param year The name of the year column from df
#' @param abundance The name of the abundance column from df
#' @param metric The turnover metric to return; the default, total, returns summed appearances and disappearances relative to total species richness across both years
#'          appearance returns the number of appearances in the second year relative to total species richness across both years
#'          disappearance returns the number of disappearances in the second year relative to the total species richness across both years
#' @return output A dataframe containing the specificed turnover metric and year
turnover_allyears<-function(df, species, year, abundance, metric="total"){
  df<-df[order(df[year]),]
  df<-df[which(df[[abundance]]>0),]
  ## split data by year
  yearlist <- split(df, df[[year]])
  ## create consecutive pairs of years
  y1 <- yearlist[-length(yearlist)]
  y2 <- yearlist[-1]
  ## rbind consecutive pairs of years
  yearpair <- Map(function(d1, d2){rbind(d1,d2)}, y1, y2)
  ## calculate turnover for across all years
  out <- Map(getturnover, y1, y2, species, metric)
  output<-as.data.frame(unlist(out))
  names(output)[1]=metric
  ## add year column
  allyr<-as.matrix(unique(df[year]))
  currentyr<-allyr[2:nrow(allyr)]
  output[year]=(currentyr)
  return(output)
}
