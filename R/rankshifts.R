#' Function for calculating mean rank shifts
#'
#' This is a function that calculates mean rank shifts
#' 
#' @param df dataframe of Community dataset. Must be in 'long' format.
#' @param replicate.var The replication variable
#' @param time.var The time variable
#' @param species.var The species variable
#' @param abundance.var The abundance variable
#' @return a dataframe, showing years compared
#' @export
meanrankshift <- function(df, time.var = "year", species.var = "species",
                    abundance.var = "abundance", replicate.var=as.character(NA)){
  if(is.na(replicate.var)==TRUE){
    output<-meanrank(df, time.var, species.var, abundance.var)}else{
        df[replicate.var]<-if(is.factor(df[[replicate.var]])==TRUE){factor(df[[replicate.var]])} else {df[replicate.var]}
        X<-split(df, df[replicate.var])
        out<-(lapply(X, FUN=meanrank, time.var, species.var, abundance.var))
        ID<-unique(names(out))
        out<-mapply(function(x, y) "[<-"(x, replicate.var, value = y) ,
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
#' @param comm_data dataframe of Community dataset. Must be in 'long' format.
#' @param time.var The time variable
#' @param species.var The species variable
#' @param abundance.var The abundance variable
#' @return a dataframe, showing years compared
meanrank <- function(comm_data, time.var = "year",
                     species.var = "species", abundance.var = "abundance"){
    ## split data by year
    yearlist <- split(comm_data, comm_data[[time.var]])
    ## Compare consecutive pairs of years
    y1 <- yearlist[-length(yearlist)]
    y2 <- yearlist[-1]

    commonspp <- Map(getintersect, y1, y2, dataname = species.var)

    names(commonspp) <- Map(function(x, y) paste0(x, "-", y), names(y1), names(y2))

    abdname1 <- paste0(abundance.var,"1")
    abdname2 <- paste0(abundance.var,"2")
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
#' @param df1 A dataframe
#' @param df2 A dataframe
#' @param dataname The name of the column on which the two datasets will be joined and intersected
getintersect <- function(df1, df2, dataname = "species"){
  commspp <- intersect(df1[[dataname]], df2[[dataname]])
  ## select out the dataname columsn from df1 and df2
  df1dataname<-data.frame(df1[[dataname]])
  names(df1dataname)=dataname
  df2dataname<-data.frame(df2[[dataname]])
  names(df2dataname)=dataname
  ## rename df1 and df2 columns
  df1[[dataname]]<-NULL
  df2[[dataname]]<-NULL
  names(df1) <- paste0(names(df1),"1")
  names(df2) <- paste0(names(df2),"2")

  df1<-cbind(df1, df1dataname)
  df2<-cbind(df2, df2dataname)
  merge(x = df1[df1[[dataname]] %in% commspp, ],
        y = df2[df2[[dataname]] %in% commspp, ])
}



