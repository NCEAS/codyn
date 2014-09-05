#' Apparent rate of immigration and extinction
#'
#' Counts how frequently a species goes extinct
#' @param presvector A vector of presence (1) and absence (0)
#' @return list of length 2, giving perceived frequencies of immigration and extinction

apparent <- function(presvector = c(0,1,1,0,0,0,1,0,1,0,1,0,0,1)){
  list(immigrate = sum(diff(presvector)==1),
       extinct   = sum(diff(presvector)==-1))
}