#' Function for calculating mean rank shifts
#'
#' This is a function that calculates mean rank shifts
#' 
#' @param data1 Community dataset. Must be in 'long' format.
#' @param replicate The replication variable
#' @param year The time variable
#' @param species The species variable
#' @param abundance The abundance variable
#' @return a dataframe, showing years compared
#' @export
meanrankshift <- function(data1, replicate="replicate", year = "year",
                     species = "species", abundance = "abundance"){
  if(is.na(replicate)==TRUE){
    output<-meanrank(data1, year, species, abundance)}else{
        data1[replicate]<-if(is.factor(data1[[replicate]])==TRUE){factor(data1[[replicate]])} else {data1[replicate]}
        X<-split(data1, data1[replicate])
        out<-(lapply(X, FUN=meanrank, year, species, abundance))
        ID<-unique(names(out))
        out<-mapply(function(x, y) "[<-"(x, replicate, value = y) ,
                    out, ID, SIMPLIFY = FALSE)
        output<-do.call("rbind", out)
        row.names(output)<-NULL
      }
return(output)
}
  
  

############################################################################
#
# Private functions: these are internal functions not intended for reuse.  
# Future package releases may change these without notice. External callers
# should not use them.
#
############################################################################



#' Function for calculating mean rank shifts
#'
#'
#' This is a function that calculates mean rank shifts
#' @param comm_data Community dataset. Must be in 'long' format.
#' @param year The time variable
#' @param species The species variable
#' @param abundance The abundance variable
#' @return a dataframe, showing years compared
meanrank <- function(comm_data, year = "year",
                     species = "species", abundance = "abundance"){
    ## split data by year
    yearlist <- split(comm_data, comm_data[[year]])
    ## Compare consecutive pairs of years
    y1 <- yearlist[-length(yearlist)]
    y2 <- yearlist[-1]

    commonspp <- Map(getintersect, y1, y2, dataname = species)

    names(commonspp) <- Map(function(x, y) paste0(x, "-", y), names(y1), names(y2))

    abdname1 <- paste0(abundance,"1")
    abdname2 <- paste0(abundance,"2")
    rank1 <- ""   # Note: initialized rank1 and rank2 simply to eliminate R CMD check NOTE
    rank2 <- ""
    ranknames <- lapply(commonspp, function(x) cbind(x,
                                                     rank1 = rank(x[[abdname1]]),
                                                     rank2 = rank(x[[abdname2]])
    ))

    rankdiff <- lapply(ranknames,
                       function(x) transform(x, abs_ch_rank = abs(rank2 - rank1)))

    MRS <- sapply(rankdiff, function(x) mean(x$abs_ch_rank))

    data.frame(year_pair = names(MRS), MRS, row.names = NULL)
}


#' Create intersected data frames
#'
#' Create intersections.
#' @param d1 A dataframe
#' @param d2 A dataframe
#' @param dataname The name of the column on which the two datasets will be joined and intersected
getintersect <- function(d1, d2, dataname = "species"){
  commspp <- intersect(d1[[dataname]], d2[[dataname]])
  ## select out the dataname columsn from d1 and d2
  d1dataname<-data.frame(d1[[dataname]])
  names(d1dataname)=dataname
  d2dataname<-data.frame(d2[[dataname]])
  names(d2dataname)=dataname
  ## rename d1 and d2 columns
  d1[[dataname]]<-NULL
  d2[[dataname]]<-NULL
  names(d1) <- paste0(names(d1),"1")
  names(d2) <- paste0(names(d2),"2")

  d1<-cbind(d1, d1dataname)
  d2<-cbind(d2, d2dataname)
  merge(x = d1[d1[[dataname]] %in% commspp, ],
        y = d2[d2[[dataname]] %in% commspp, ])
}



