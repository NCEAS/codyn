#' @title Community Diversity
#' @description Calculates Shannon's or Simpson's diversity of a community.  
#' @param df A data frame containing species and abundance columns and optional columns of time and/or replicate. Note that at least time.var or replicate.var must be specified.
#' @param time.var The name of the optional time column 
#' @param species.var The name of the species column 
#' @param abundance.var The name of the abundance column 
#' @param replicate.var The name of the optional replicate column 
#' @param metric The diversity metric to return:
#' \itemize{
#'  \item{"Shannon": }{The default metric, calculates Shannon's diversity.}
#'  \item{"Simpson": }{Calculates Inverse of Simpson's diversity.}
#' }
#' 
#' @return The community_diversity function returns a data frame with the following attributes:
#' \itemize{
#'  \item{time.var: }{A column that has the same name and type as the time.var column, if time.var is specified.}
#'  \item{replicate.var: }{A column that has same name and type as the replicate.var column, if replicate.var is specified.}
#'  \item{Shannon: }{A numeric column of Shannons diversity if metric = "Shannon"}
#'  \item{Simpson: }{A numeric column of Simpsons diversity if metric = "InvSimpson"}
#' }
#' @examples
#' data(pplots)
#' #Example with both time and replicates
#' community_diversity(subset(pplots, plot==25|plot==6), 
#'                     time.var="year", 
#'                     replicate.var = "plot", 
#'                     abundance.var = "relative_cover") # for Shannon's diversity metric
#'
#' community_diversity(subset(pplots, plot==25|plot==6), 
#'                     time.var="year", 
#'                     replicate.var = "plot", 
#'                     abundance.var = "relative_cover", 
#'                     metric = "InvSimpson") # for Inverse of Simpson's diversity metric
#'
#' #Example with no replicates
#' community_diversity(subset(pplots, plot==25), 
#'                     time.var="year", 
#'                     abundance.var = "relative_cover") # for Shannon's diversity metric
#'
#' #Example with only a single time point
#' community_diversity(subset(pplots, year==2002&plot==25|year==2002&plot==6), 
#'                     replicate.var = "plot", 
#'                     abundance.var = "relative_cover") # for Shannon's diversity metric
#'
#'
#' @export
#'
community_diversity <- function(df,  time.var = NULL, 
                                abundance.var, 
                                replicate.var = NULL,  
                                metric = c("Shannon", "InvSimpson")) {
  
  # verify metric choice
  metric <- match.arg(metric)
  
  if(is.null(replicate.var)) {
    myformula <- as.formula(paste(abundance.var, "~", time.var))
  } else if(is.null(time.var)) {
    myformula <- as.formula(paste(abundance.var, "~", replicate.var))
  } else {
    myformula <- as.formula(paste(abundance.var, "~", time.var, "+", replicate.var))
  }
  
  if(metric == "Shannon") {
      comdiv <- aggregate(myformula, data = df, FUN = function(x) diversity = Shannon(x))
  } else if(metric == "InvSimpson") {
      comdiv <- aggregate(myformula, data = df, FUN = function(x) diversity = InvSimpson(x))
  }
  
  names(comdiv)[names(comdiv) == abundance.var] <- metric
  
  return(comdiv)
}

############################################################################
#
# Private functions: these are internal functions not intended for reuse.
# Future package releases may change these without notice. External callers
# should not use them.
#
############################################################################


#' A function to calculate Simpson's Divsersity from Smith and Wilson 1996
#' @param x the vector of abundances of each species
#' @param N the total abundance
#' @param p the vector of relative abundances of each species
InvSimpson <- function(x, N = sum(x[x!=0&!is.na(x)]), ps = x[x!=0&!is.na(x)]/N, p2=ps*ps ){
  D <- sum(p2)
  1/D
}

#' A function to calculate Shannon's Divsersity 
#' @param x the vector of abundances of each species
#' @param N the total abundance
#' @param p the vector of relative abundances of each species
Shannon <- function(x, N = sum(x[x!=0&!is.na(x)]), ps = x[x!=0&!is.na(x)]/N ){
  -sum(ps*log(ps))
}




