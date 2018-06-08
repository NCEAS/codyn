#' @title Community Structure
#' 
#' @description Calculates species richness and evenness of a community. Evenness may be calculated as Simpson's (1/D/S), EQ, or Evar, but only one metric of evenness can be calculated at a time and must be specified.
#'   
#' @param df A data frame containing species and abundance columns and optional columns of time and/or replicate.
#' @param time.var The name of the optional time column 
#' @param abundance.var The name of the abundance column 
#' @param replicate.var The name of the optional replicate column. If specified, replicate must be unique within the dataset and cannot be nested within treatments or blocks. 
#' @param metric The evenness metric to return:
#' \itemize{
#'  \item{"Evar": }{The default metric, calculates Evar evenness from Smith and Wilson 1996}
#'  \item{"SimpsonEvenness": }{Calculates Simpsons' evenness}
#'  \item{"EQ": }{Calculates EQ evenness from Smith and Wilson 1996}
#' }
#'  
#' @return The community_structure function returns a data frame with the
#'   following attributes:
#' \itemize{
#'  \item{time.var: }{A column that has the same name and type as the time.var column, if time.var is specified.}
#'  \item{replicate.var: }{A column that has same name and type as the replicate.var column, if specified.}
#'  \item{richness: }{A numeric column of species richness}
#'  \item{Evar: }{A numeric column of Evar evenness if evenness = "Evar"}
#'  \item{EQ: }{A numeric column of EQ evenness if evenness = "EQ"}
#'  \item{SimpsonEvenness: }{A numeric column of Simpsons evenness if evenness =  "SimpsonEveness"}
#' }
#' @references Smith, B. and Wilson, J. B. 1996. A consumer's guide to evenness indices. Oikos 76: 70-82.
#' @examples
#' data(pplots)
#' #Example with both time and replicates
#' df <- subset(pplots, plot == 25 | plot == 6)
#' community_structure(df, 
#'                     time.var="year", 
#'                     replicate.var = "plot", 
#'                     abundance.var = "relative_cover") # for Evar evenness metric
#'
#' df <- subset(pplots, plot == 25 | plot == 6)
#' community_structure(df,
#'                     time.var="year", 
#'                     replicate.var = "plot", 
#'                     abundance.var = "relative_cover", 
#'                     metric = "SimpsonEvenness") # for Simpson's evenness metric
#'
#' #Example with no replicates
#' df <- subset(pplots, plot == 25)
#' community_structure(df, 
#'                     time.var="year", 
#'                     abundance.var = "relative_cover",
#'                     metric = "EQ") # for EQ evenness metric
#'
#' #Example with only a single time point and no replicates
#' df <- subset(pplots, plot == 25 & year == 2002)
#' community_structure(df, 
#'                     abundance.var = "relative_cover") # for Evar evenness metric
#' @importFrom stats aggregate.data.frame
#' @export

community_structure <- function(df,
                                time.var = NULL, 
                                abundance.var, 
                                replicate.var = NULL, 
                                metric = c("Evar", "SimpsonEvenness", "EQ")) {

  # verify metric choice
  metric <- match.arg(metric)
  
  # check no NAs in abundance column
  if(any(is.na(df[[abundance.var]]))) stop("Abundance column contains missing values")

  # specify aggregate formula from arguments
  if(is.null(replicate.var)) {
    by <- time.var
  } else if(is.null(time.var)) {
    by <- replicate.var
  } else {
    by <- c(time.var, replicate.var)
  }

  # get function for chosen metric, and calculate output
  evenness <- get(metric)
  comstruct <- aggregate.data.frame(df[abundance.var], df[by],
                         FUN = function(x) cbind(S(x), evenness(x)))
  comstruct <- do.call(data.frame, comstruct)
  
  if(any(is.na(comstruct[[paste(abundance.var, 2, sep = ".")]]))) warning("Evenness values contain NAs because there are plots with only one species")
  
  names(comstruct) <- c(by, 'richness', metric)

   return(comstruct)
}


############################################################################
#
# Private functions: these are internal functions not intended for reuse.
# Future package releases may change these without notice. External callers
# should not use them.
#
############################################################################

# A function to calculate 1/D (inverse of Simpson's) from Smith and Wilson 1996
# @param S the number of species in the sample
# @param x the vector of abundances of each species
# @param N the total abundance
# @param ps the vector of relative abundances of each species
# @param p2 the vector of the square of relative abundances
SimpsonEvenness <- function(x, S = length(x[x != 0]), N = sum(x[x != 0]), ps = x[x != 0]/N, p2 = ps*ps ){
  D <- sum(p2)
  (1/D)/S
}

#' Utility function to calculate EQ evenness from Smith and Wilson 1996
#' @param x Vector of abundance of each species
#' If all abundances are equal it returns a 1
#' @importFrom stats lm
EQ <- function(x){
  x1 <- x[x != 0]
  if (length(x1) == 1) {
    return(NA)
  }
  if (abs(max(x1) - min(x1)) < .Machine$double.eps^0.5) {
    return(1)
  }
  r <- rank(-x1, ties.method = "average")
  r_scale <- r/max(r)
  x_log <- log(x1)
  fit <- lm(r_scale~x_log)
  b <- fit$coefficients[[2]]
  -2/pi*atan(b)
}
