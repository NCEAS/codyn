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
#'  \item{Simpson: }{A numeric column of Simpsons diversity if metric = "InverseSimpson"}
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
#'                     metric = "InverseSimpson") # for Inverse of Simpson's diversity metric
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
community_diversity <- function(df,  time.var = NULL, 
                                abundance.var, 
                                replicate.var = NULL,  
                                metric = c("Shannon", "InverseSimpson")) {
  
  # verify metric choice
  metric <- match.arg(metric)

  # specify aggregate formula from arguments
  if(is.null(replicate.var)) {
    by <- time.var
  } else if(is.null(time.var)) {
    by <- replicate.var
  } else {
    by <- c(time.var, replicate.var)
  }

  # get function for chosen metric, and calculate output
  diversity <- get(metric)
  comdiv <- aggregate.data.frame(df[abundance.var], df[by], FUN = diversity)
  names(comdiv) <- c(by, paste(abundance.var, metric, sep = '.'))
  
  return(comdiv)
}

############################################################################
#
# Private functions: these are internal functions not intended for reuse.
# Future package releases may change these without notice. External callers
# should not use them.
#
############################################################################


# A function to calculate Simpson's Divsersity from Smith and Wilson 1996
# @param x the vector of abundances of each species
# @param N the total abundance
# @param p the vector of relative abundances of each species
InverseSimpson <- function(x, N = sum(x[x!=0&!is.na(x)]), ps = x[x!=0&!is.na(x)]/N, p2=ps*ps ){
  D <- sum(p2)
  1/D
}

# A function to calculate Shannon's Divsersity 
# @param x the vector of abundances of each species
# @param N the total abundance
# @param p the vector of relative abundances of each species
Shannon <- function(x, N = sum(x[x!=0&!is.na(x)]), ps = x[x!=0&!is.na(x)]/N ){
  -sum(ps*log(ps))
}




