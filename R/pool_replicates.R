#' @title pool replicates into treatments and add ranks
#' @description pool replicates into treatments and add ranks
#' @param df A data frame containing an optional time column and requred species, abundance, replicate and treatment columns
#' @param time.var The name of the optional time column 
#' @param species.var The name of the species column 
#' @param abundance.var The name of the abundance column 
#' @param replicate.var The name of the replicate column 
#' @param treatment.var The name of the treatment column


pool_replicates <- function(df, time.var=NULL, species.var, abundance.var, replicate.var, treatment.var) {
  df<-as.data.frame(df)
  
  if(is.null(time.var)){
    ## isolate rep-trts
    rep_trt<-unique(subset(df, select = c(replicate.var, treatment.var)))
    
    #make wide to add zeros and long again and add back in treatments
    df2 <- subset(df, select = c(replicate.var,species.var,abundance.var))
    wide <- reshape(df2, idvar = replicate.var, timevar = species.var, direction = "wide")
    wide[is.na(wide)] <- 0
    
    long<-reshape(wide, idvar = replicate.var, ids = replicate.var, time = names(wide), timevar = abundance.var, direction = "long")
    colnames(long)[3] <- abundance.var
    allsp <- merge(long, rep_trt, by=replicate.var)
    
  #get averages of each species by treatment
    myformula <- as.formula(paste(abundance.var, "~", treatment.var, "+", species.var))
    spave <- aggregate(myformula, FUN=mean, data=allsp)

  ##add ranks for present species
    rank_pres <- subset(spave, spave[[abundance.var]]!=0)
    rank_pres$rank <- ave(rank_pres[[abundance.var]], rank_pres[[treatment.var]], FUN = function(x) rank(-x, ties.method = "average"))
  
  ###make zero abundant species have the rank S+1 (the size of the species pool plus 1)
  ##pull out zeros
  zeros <- subset(spave, spave[[abundance.var]]==0)
  ##get species richness for each year
  myformula <- as.formula(paste(abundance.var, "~", treatment.var))
  SpR <- aggregate(myformula, FUN = S, data = spave)
  colnames(SpR)[2] <- "S"
  
  ##merge together make zero abundances rank S+1
  zero_rank <- merge(zeros, SpR, by = c(treatment.var))
  zero_rank$rank <- zero_rank$S+1
  zero_rank <- subset(zero_rank, select = -S)
  
  ##combine all
  rankdf<-rbind(rank_pres, zero_rank)
  }
  else{
    rep_trt<-unique(subset(df, select = c(replicate.var, treatment.var)))
    
    # sort and apply fill_zeros to all time steps
    df <- df[order(df[[time.var]]),]
    X <- split(df, df[time.var])
    out <- lapply(X, FUN = fill_zeros_rep, replicate.var, species.var, abundance.var)
    ID <- unique(names(out))
    out <- mapply(function(x, y) "[<-"(x, time.var, value = y) ,
                  out, ID, SIMPLIFY = FALSE)
    out2 <- do.call("rbind", out)
    
    allsp <- merge(out2, rep_trt, by=replicate.var)
    
    #get averages of each species by treatment
    myformula <- as.formula(paste(abundance.var, "~", treatment.var, "+", species.var, "+", time.var))
    spave <- aggregate(myformula, FUN=mean, data=allsp)
    
    ##add ranks for present species
    rank_pres <- subset(spave, spave[[abundance.var]]!=0)
    rank_pres$trt_time <- paste(rank_pres[[treatment.var]], rank_pres[[time.var]], sep="##")
    rank_pres$rank <- ave(rank_pres[[abundance.var]], rank_pres$trt_time, FUN = function(x) rank(-x, ties.method = "average"))
    rank_pres <- subset(rank_pres, select = -trt_time)

    ###make zero abundant species have the rank S+1 (the size of the species pool plus 1)
    ##pull out zeros
    zeros <- subset(spave, spave[[abundance.var]]==0)
    ##get species richness for each year
    myformula <- as.formula(paste(abundance.var, "~", treatment.var, "+", time.var))
    SpR <- aggregate(myformula, FUN = S, data = spave)
    colnames(SpR)[3] <- "S"
    
    ##merge together make zero abundances rank S+1
    zero_rank <- merge(zeros, SpR, by = c(treatment.var, time.var))
    zero_rank$rank <- zero_rank$S+1
    zero_rank <- subset(zero_rank, select = -S)
    
    ##combine all
    rankdf<-rbind(rank_pres, zero_rank)
    
  }
  
  return(rankdf)
}
