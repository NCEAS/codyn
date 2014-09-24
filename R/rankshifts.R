#' Create intersected data frames
#'
#' Create intersections.
#' @param d1 A dataframe
#' @param d2 A dataframe
#' @param dataname The name of the column on which the two datasets will be joined and intersected
getintersect <- function(d1, d2, dataname = "species"){
  commspp <- intersect(d1[[dataname]], d2[[dataname]])
  names(d1)[-1] <- paste0(names(d1)[-1],"1")
  names(d2)[-1] <- paste0(names(d2)[-1],"2")
  merge(x = d1[d1[[dataname]] %in% commspp, ],
        y = d2[d2[[dataname]] %in% commspp, ])
}


#' Function for calculating mean rank shifts
#'
#'
#' This is a function that calculates mean rank shifts
#' @param comm_data Community dataset. Must be in 'long' format.
#' @param timevar The time variable
#' @param speciesvar The species variable
#' @param abdvar The abundance variable
#' @return a dataframe, showing years compared
#' @export
meanrank <- function(comm_data = dat, timevar = "year",
                     speciesvar = "species", abdvar = "abundance"){
  ## split data by year
  yearlist <- split(dat, dat[[timevar]])
  ## Compare consecutive pairs of years
  y1 <- yearlist[-length(yearlist)]
  y2 <- yearlist[-1]

  commonspp <- Map(getintersect, y1, y2)

  names(commonspp) <- Map(function(x, y) paste0(x, "-", y), names(y1), names(y2))

  abdname1 <- paste0(abdvar,"1")
  abdname2 <- paste0(abdvar,"2")

  ranknames <- lapply(commonspp, function(x) cbind(x,
                                                   rank1 = rank(x[[abdname1]]),
                                                   rank2 = rank(x[[abdname2]])
  ))

  rankdiff <- lapply(ranknames,
                     function(x) transform(x, abs_ch_rank = abs(rank2 - rank1)))

  MRS <- sapply(rankdiff, function(x) mean(x$abs_ch_rank))

  data.frame(year_pair = names(MRS), MRS, row.names = NULL)
}



