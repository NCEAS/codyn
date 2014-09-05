#' Function for calculating mean rank shifts
#'
#' This is a function that calculates mean rank shifts
#' @param comm_data Community data
#' @param timevar the Time variable
#' @param speciesvar the Species variable
#' @return the Truth


getintersect <- function(d1, d2, dataname = "species"){
  commspp <- intersect(d1[[dataname]], d2[[dataname]])
  names(d1)[-1] <- paste0(names(d1)[-1],"1")
  names(d2)[-1] <- paste0(names(d2)[-1],"2")
  merge(x = d1[d1[[dataname]] %in% commspp, ],
        y = d2[d2[[dataname]] %in% commspp, ])
}

meanrank <- function(comm_data = dat, timevar = "year", speciesvar = "species"){
  ## split data by year
  yearlist <- split(dat, dat[["year"]])
  ## Compare consecutive pairs of years
  y1 <- yearlist[-length(yearlist)]
  y2 <- yearlist[-1]

  commonspp <- Map(getintersect, y1, y2)
  combineddata <- do.call(rbind, commonspp)

  ## generate names to be combined

  transform(combineddata, tipepair = paste0("year1"

}



