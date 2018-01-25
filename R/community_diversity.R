#' @title Community Diversity
#' @description 
#' @param df A data frame containing species and abundance columns and optional columns of time point and/or replicates. Note that at least time.var or replicate.var must be specified.
#' @param time.var The name of the optional time column 
#' @param species.var The name of the species column 
#' @param abundance.var The name of the abundance column 
#' @param replicate.var The name of the optional replicate column 
#' @param diversity The diversity metric to return:
#' \itemize{
#'  \item{"Shannon": }{The default metric, calculates Shannon diversity.}
#'  \item{"Simpson": }{Calculates Simpson diversity.}
#' }
#' 
#' @return The community_diversity function returns a data frame with the following attributes:
#' \itemize{
#'  \item{time.var: }{A column that has the same name and type as the time.var column, if time.var is specified.}
#'  \item{replicate.var: }{A column that has same name and type as the replicate.var column, if specified.}
#'  \item{Shannon: }{A numeric column of Shannons diversity if diversity = "Shannon"}
#'  \item{Simpson: }{A numeric column of Simpsons diversity if diversity = "Simpson"}
#' }
#' @export
#'
community_diversity <- function(df,  time.var = NULL, 
                                abundance.var, 
                                replicate.var = NULL,  
                                diversity = "Shannon") {
  
  if(is.null(replicate.var)) {
    
    myformula <- as.formula(paste(abundance.var, "~", time.var))
    
    if(diversity == "Shannon") {
      
      comdiv <- aggregate(myformula, data = df, FUN = function(x) diversity = Shannon(x))
      names(comdiv)[2] <- "Shannon"
      
    } else {
      
      comdiv <- aggregate(myformula, data=df, FUN = function(x) diversity = Simpson(x))
      names(comdiv)[2] <- "Simpson"
      
      }
  } else {
    
    if(is.null(time.var)) {
      
      myformula <- as.formula(paste(abundance.var, "~", replicate.var))
      
      if(diversity == "Shannon"){
        
        comdiv <- aggregate(myformula, data = df, FUN = function(x) diversity = Shannon(x))
        names(comdiv)[2] <- "Shannon"
        
      } else {
        
        comdiv <- aggregate(myformula, data=df, FUN = function(x) diversity = Simpson(x))
        names(comdiv)[2] <- "Simpson"
        
      }
      
    } else {
      
    myformula <- as.formula(paste(abundance.var, "~", time.var, "+", replicate.var))
    
    if(diversity == "Shannon") {
      comdiv <- aggregate(myformula, data = df, FUN = function(x) diversity = Shannon(x))
      names(comdiv)[3] <- "Shannon"
      
    } else {
      
      comdiv <- aggregate(myformula, data = df, FUN = function(x) diversity = Simpson(x))
      names(comdiv)[3] <- "Simpson"
    }
    }
  }
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
Simpson <- function(x, N = sum(x[x!=0&!is.na(x)]), ps = x[x!=0&!is.na(x)]/N, p2=ps*ps ){
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




