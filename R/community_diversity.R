#' @title Community Diversity
#' 
#' @description Calculates Shannon's or Inverse Simpson's diversity of a
#'   community, but only one measure of diversity can be calculated at a time and
#'   must be specified.
#'   
#' @param df A data frame containing species and abundance columns and optional
#'   columns of time and/or replicate.
#' @param time.var The name of the optional time column 
#' @param abundance.var The name of the abundance column 
#' @param replicate.var The name of the optional replicate column. If specified,
#'   replicate must be unique within the dataset and cannot be nested within
#'   treatments or blocks.
#' @param metric The diversity measure to return:
#' \itemize{
#'  \item{"Shannon": }{The default measure, calculates Shannon's diversity.}
#'  \item{"InverseSimpson": }{Calculates inverse of Simpson's diversity.}
#' }
#' 
#' @return The community_diversity function returns a data frame with the
#'   following attributes:
#' \itemize{
#'  \item{time.var: }{A column that has the same name and type as the time.var
#'  column, if time.var is specified.}
#'  \item{replicate.var: }{A column that has same name and type as the
#'  replicate.var column, if replicate.var is specified.}
#'  \item{Shannon: }{A numeric column of Shannon's diversity if metric =
#'  "Shannon"}
#'  \item{InverseSimpson: }{A numeric column of the inverse of Simpson's
#'  diversity if metric = "InverseSimpson"}
#' }
#' @references Magurran, A.E. 2004. Measuring Biological Diversity. Blackwell
#'   Publishing, Malden MA, USA.
#' @examples
#' data(pplots)
#' #Example with both time and replicates
#' df <- subset(pplots, plot == 25 | plot == 6)
#' community_diversity(df,
#'                     time.var="year", 
#'                     replicate.var = "plot", 
#'                     abundance.var = "relative_cover") # for Shannon's diversity measure
#'
#'df <- subset(pplots, plot == 25 | plot == 6)
#' community_diversity(df, 
#'                     time.var="year", 
#'                     replicate.var = "plot", 
#'                     abundance.var = "relative_cover", 
#'                     metric = "InverseSimpson") # for Inverse of Simpson's diversity measure
#'
#' #Example with no replicates
#' df <- subset(pplots, plot == 25)
#' community_diversity(df,
#'                     time.var="year", 
#'                     abundance.var = "relative_cover") # for Shannon's diversity measure
#'                     
#' #Example with no time or replicate
#' df <- subset(pplots, plot == 25 & year == 2002)
#' community_diversity(df,
#'                     abundance.var = "relative_cover") # for Shannon's diversity measure
#' @importFrom stats aggregate.data.frame
#' @export
community_diversity <- function(df,
                                time.var = NULL, 
                                abundance.var, 
                                replicate.var = NULL,  
                                metric = c("Shannon", "InverseSimpson")) {
  
  # verify measure choice
  measure <- match.arg(metric)
  
  # check no NAs in abundance column
  if (any(is.na(df[[abundance.var]]))) stop("Abundance column contains missing values")

  # specify aggregate formula from arguments
  if (is.null(replicate.var)) {
    by <- time.var
  } else if (is.null(time.var)) {
    by <- replicate.var
  } else {
    by <- c(time.var, replicate.var)
  }

  # get function for chosen measure, and calculate output
  diversity <- get(measure)
  comdiv <- aggregate.data.frame(df[abundance.var], df[by], FUN = diversity)
  names(comdiv) <- c(by, measure)
  
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
# @param ps the vector of relative abundances of each species
# @param p the vector of the square of relative abundances
InverseSimpson <- function(x, N = sum(x[x != 0]), ps = x[x != 0]/N, p2=ps*ps ){
  D <- sum(p2)
  1/D
}

# A function to calculate Shannon's Diversity 
# @param x the vector of abundances of each species
# @param N the total abundance
# @param ps the vector of relative abundances of each species
Shannon <- function(x, N = sum(x[x != 0]), ps = x[x != 0]/N ){
  -sum(ps*log(ps))
}




