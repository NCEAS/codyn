#' Function for calculating mean rank shifts
#'
#' This is a function that calculates mean rank shifts
#' @param comm_data Community data
#' @param timevar the Time variable
#' @param speciesvar the Species variable
#' @return the Truth


meanrank <- function(comm_data = dat, timevar = "year", speciesvar = "species"){
  ## split data by year
  yearlist <- split(as.character(dat[["species"]]), dat[["year"]])
  ## Compare consecutive pairs of years
  y2 <- yearlist[-length(yearlist)]
  y1 <- yearlist[-1]

  Map(intersect, y1, y2)
  ## compare two years
}



