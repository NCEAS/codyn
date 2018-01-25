#' @title Add Ranks for Replicates
#'@description Ranks species by abundance in each replicate
#' @param df A data frame containing species, abundance and replicate columns and an optional column of time points
#' @param time.var The name of the optional time column 
#' @param species.var The name of the species column 
#' @param abundance.var The name of the abundance column 
#' @param replicate.var The name of the replicate column 
#' 
#' @return The add_ranks function returns a data frame with the following attributes:
#' \itemize{
#'  \item{time.var: }{A column that has the same name and type as the time.var column, if time.var is specified.}
#'  \item{abundance.var: }{A column that has same name and type as the abundance.var column.}
#'  \item{species.var: }{A column that has same name and type as the species.var column.}
#'  \item{replicate.var: }{A column that has same name and type as the replicate.var column.}
#'  \item{rank: }{A numeric column with the species rank; a rank of 1 indicates the species was most abundant in that time period. Species that are not present in that time period have the largest rank value.}
#' }
#' @export
#' 
add_ranks_replicate <- function(df, time.var = NULL, 
                                species.var, 
                                abundance.var, 
                                replicate.var) {
  
  df <- as.data.frame(df)
  
  if(is.null(time.var)) {
    
    df <- subset(df, select = c(replicate.var, species.var, abundance.var))
    
  ##add ranks for present species
  rank_pres <- subset(df, df[[abundance.var]] != 0)
  rank_pres$rank <- ave(rank_pres[[abundance.var]], rank_pres[[replicate.var]], FUN = function(x) rank(-x, ties.method = "average"))

  #adding zeros
  
  # sort and fill zeros 
  
  df2 <- subset(df, select = c(replicate.var,species.var,abundance.var))
  wide <- reshape(df2, idvar = replicate.var, timevar = species.var, direction = "wide")
  wide[is.na(wide)] <- 0
  
  long <- reshape(wide, idvar = replicate.var, ids = replicate.var, time = names(wide), timevar = abundance.var, direction = "long")
  colnames(long)[3] <- abundance.var
  allsp <- long

  ###make zero abundant species have the rank S+1 (the size of the species pool plus 1)
  ##pull out zeros
  zeros <- subset(allsp, allsp[[abundance.var]] == 0)
  
  ##get species richness for each year
  myformula <- as.formula(paste(abundance.var, "~", replicate.var))
  SpR <- aggregate(myformula, FUN = S, data = allsp)
  colnames(SpR)[2] <- "S"
  
  ##merge together make zero abundances rank S+1
  zero_rank <- merge(zeros, SpR, by = c(replicate.var))
  zero_rank$rank <- zero_rank$S+1
  zero_rank <- subset(zero_rank, select = -S)
  
  ##combine all
  rank <- rbind(rank_pres, zero_rank)
  
  } else {
    
    df <- subset(df, select = c(time.var, species.var, abundance.var, replicate.var))
    ##add ranks for present species
    rank_pres <- subset(df, df[[abundance.var]]!=0)
    rank_pres$rep_time <- paste(rank_pres[[replicate.var]], rank_pres[[time.var]], sep="_")
    rank_pres$rank <- ave(rank_pres[[abundance.var]], rank_pres$rep_time, FUN = function(x) rank(-x, ties.method = "average"))
    rank_pres <- subset(rank_pres, select = -rep_time)
    
    #adding zeros
    
    # sort and apply fill_zeros to all time steps
    df <- df[order(df[[time.var]]),]
    X <- split(df, df[time.var])
    out <- lapply(X, FUN = fill_zeros_rep, replicate.var, species.var, abundance.var)
    ID <- unique(names(out))
    out <- mapply(function(x, y) "[<-"(x, time.var, value = y) ,
                  out, ID, SIMPLIFY = FALSE)
    allsp <- do.call("rbind", out)
    
    
    ###make zero abundant species have the rank S+1 (the size of the species pool plus 1)
    ##pull out zeros
    zeros <- subset(allsp, allsp[[abundance.var]] == 0)
    
    ##get species richness for each year
    ## Note to Meghan: This uses a function that I thought was only for community_structure
    ## Might make a separate file of shared diversity functions
    myformula <- as.formula(paste(abundance.var, "~", replicate.var, "+", time.var))
    SpR <- aggregate(myformula, FUN = S, data = allsp)
    colnames(SpR)[3] <- "S"
    
    ##merge together make zero abundances rank S+1
    zero_rank <-merge(zeros, SpR, by = c(time.var, replicate.var))
    zero_rank$rank <- zero_rank$S + 1
    zero_rank <- subset(zero_rank, select = -S)
    
    ##combine all
    rank <- rbind(rank_pres, zero_rank)
  }
  
  return(rank)
}



############################################################################
#
# Private functions: these are internal functions not intended for reuse.
# Future package releases may change these without notice. External callers
# should not use them.
#
############################################################################

