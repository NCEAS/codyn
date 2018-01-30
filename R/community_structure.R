#' @title Community Structure
#' @description Calculates species richness and evenness of a community. Evenness may be calculated as either Simpson's (1/D/S) or as EQ.
#' @param df A data frame containing species and abundance columns and optional columns of time and/or replicate. Note that at least time.var or replicate.var must be specified.
#' @param time.var The name of the optional time column 
#' @param species.var The name of the species column 
#' @param abundance.var The name of the abundance column 
#' @param replicate.var The name of the optional replicate column 
#' @param evenness The evenness metric to return:
#' \itemize{
#'  \item{"EQ": }{The default metric, calculates EQ evenness from Smith and Wilson 1996}
#'  \item{"SimpEven": }{Calculates Simpsons eveness.}
#' }
#'  
#' @return The community_structure function returns a data frame with the following attributes:
#' \itemize{
#'  \item{time.var: }{A column that has the same name and type as the time.var column, if time.var is specified.}
#'  \item{replicate.var: }{A column that has same name and type as the replicate.var column, if specified.}
#'  \item{richness: }{A numeric column of species richness}
#'  \item{Shannon: }{A numeric column of EQ evenness if evenness = "EQ"}
#'  \item{Simpson: }{A numeric column of Simpsons evenness if evenness = "SimpEven"}
#' }
#' @references Smith, B. and Wilson, J. B. 1996. A consumer's guide to evnness indices. Oikos 76: 70-82.
#' @example
#' data(pplots)
#' #Example with both time and replicates
#' community_structure(subset(pplots, plot==25|plot==6), 
#'                     time.var="year", 
#'                     replicate.var = "plot", 
#'                     abundance.var = "relative_cover") # for EQ evenness metric
#'
#' community_structure(subset(pplots, plot==25|plot==6), 
#'                     time.var="year", 
#'                     replicate.var = "plot", 
#'                     abundance.var = "relative_cover", 
#'                     metric = "SimpEven") # for Simpson's evenness metric
#'
#' #Example with no replicates
#' community_structure(subset(pplots, plot==25), 
#'                     time.var="year", 
#'                     abundance.var = "relative_cover") # for EQ evenness metric
#'
#' #Example with only a single time point
#' community_structure(subset(pplots, year==2002&plot==25|year==2002&plot==6), 
#'                     replicate.var = "plot", 
#'                     abundance.var = "relative_cover")# For EQ evenness metric
#'
#' @export

community_structure <- function(df,  time.var = NULL, 
                                abundance.var, 
                                replicate.var = NULL, 
                                metric = c("EQ", "SimpsonEvenness")) {
                                  
  # verify metric choice
  metric <- match.arg(metric)

  if(is.null(replicate.var)) {
    myformula <- as.formula(paste(abundance.var, "~", time.var))
  } else if(is.null(time.var)) {
    myformula <- as.formula(paste(abundance.var, "~", replicate.var))
  } else {
    myformula <- as.formula(paste(abundance.var, "~", time.var, "+", replicate.var))
  }

  if(metric == "EQ") {
    comstruct <- do.call(
      data.frame,
      aggregate(myformula, data = df,
                FUN = function(x) c(SpR = S(x), evenness = EQ(x))))
  } else if(metric == "SimpsonEvenness"){
    comstruct <- do.call(
      data.frame,
      aggregate(myformula, data = df,
                FUN = function(x) c(SpR = S(x), evenness = SimpEven(x))))
  } 

  idx <- which(names(comstruct) == paste(abundance.var, 'SpR', sep = '.'))
  names(comstruct)[idx] <- "richness"
  idx <- which(names(comstruct) == paste(abundance.var, 'evenness', sep = '.'))
  names(comstruct)[idx] <- metric
  
  return(comstruct)
}


############################################################################
#
# Private functions: these are internal functions not intended for reuse.
# Future package releases may change these without notice. External callers
# should not use them.
#
############################################################################

#'  A function to calculate E1/D (inverse of Simpson's) from Smith and Wilson 1996
#' @param S the number of species in the sample
#' @param x the vector of abundances of each species
#' @param N the total abundance
#' @param p the vector of relative abundances of each species
SimpEven <- function(x, S = length(x[x != 0 & !is.na(x)]), N = sum(x[x != 0 & !is.na(x)]), ps = x[x != 0 & !is.na(x)]/N, p2 = ps*ps ){
  D <- sum(p2)
  (1/D)/S
}
