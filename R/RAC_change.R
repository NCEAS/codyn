##calculating RAC changes
#' @title Rank Abundance Curve Changes
#' @description 
#' @param df A data frame containing time, species and abundance columns and an optional column of replicates
#' @param time.var The name of the time column 
#' @param species.var The name of the species column 
#' @param abundance.var The name of the abundance column 
#' @param replicate.var The name of the optional replicate column 
#' 
#' 
#' TO DO: Add time.var_pair

RAC_change <- function(df, time.var, species.var, abundance.var, replicate.var=NULL) {
  if(is.null(replicate.var)){
  
  rankdf <- add_ranks_time(df, time.var, species.var, abundance.var, replicate.var=NULL)
  
  time1 <- sort(unique(rankdf[[time.var]]))
  time2 <- c(time1[2:length(time1)], NA)

  # current year rankdf
  df2 <- rankdf[which(rankdf[[time.var]]%in%time2),]
  
  # previous year rank df
  mytimes <- data.frame(cbind(time1, time2))
  names(mytimes) = c(time.var, "dummytime")

  df1 <- merge(rankdf, mytimes)
  names(df1)[1] <- 'time1'
  names(df1)[[ncol(df1)]] <- time.var
  df1 <- subset(df1, !is.na(df1[[time.var]]))
  
  # merge: .x is for previous time point, .y for current time point, time.var corresponds to current (i.e., .y)
  df12 <- merge(df1, df2,  by=c(species.var, time.var), all=T)
  df12<-subset(df12, df12[[paste(abundance.var, ".x", sep = "")]]!=0|df12[[paste(abundance.var, ".y", sep = "")]]!=0)
  df12<-subset(df12, !is.na(df12[[paste(abundance.var, ".x", sep = "")]]) & !is.na(df12[[paste(abundance.var, ".y", sep = "")]]))
  
  # sort and apply to all time pairs
  df12 <- df12[order(df12[[time.var]]),]
  X <- split(df12, df12[[time.var]])
  out <- lapply(X, FUN=SERGL, time.var, "rank.x", "rank.y", paste(abundance.var, ".x", sep = ""),paste(abundance.var, ".y", sep = "")) ## need to check how SERGL assigns time names
  output <- do.call("rbind", out)  
}
  else{
    
    rankdf <- add_ranks_time(df,  time.var, species.var, abundance.var, replicate.var)
    
    time1 <- sort(unique(rankdf[[time.var]]))
    time2 <- c(time1[2:length(time1)], NA)
    
    # current year rankdf
    df2 <- rankdf[which(rankdf[[time.var]]%in%time2),]
    
    # previous year rank df
    mytimes <- data.frame(cbind(time1, time2))
    names(mytimes) = c(time.var, "dummytime")
    
    df1 <- merge(rankdf, mytimes)
    names(df1)[1] <- 'time1'
    names(df1)[[ncol(df1)]] <- time.var
    df1 <- subset(df1, !is.na(df1[[time.var]]))
    
    # merge: .x is for previous time point, .y for current time point, time.var corresponds to current (i.e., .y)
    df12 <- merge(df1, df2,  by=c(species.var,replicate.var, time.var), all=T)
    df12<-subset(df12, df12[[paste(abundance.var, ".x", sep = "")]]!=0|df12[[paste(abundance.var, ".y", sep = "")]]!=0)
    
    #what does this step do? Why is it necessary?
    df12<-subset(df12, !is.na(df12[[paste(abundance.var, ".x", sep = "")]]) & !is.na(df12[[paste(abundance.var, ".y", sep = "")]]))
    df12$splitvariable <- paste(df12[[replicate.var]], df12[[time.var]], sep="##") 
    
    # sort and apply to all time and replicate combinations
    df12 <- df12[order(df12$splitvariable),]
    X <- split(df12, df12$splitvariable)
    
    
    out <- lapply(X, FUN=SERGL, time.var, "rank.x", "rank.y", paste(abundance.var, ".x", sep = ""),paste(abundance.var, ".y", sep = "")) 
    ID <- unique(names(out))
    out <- mapply(function(x, y) "[<-"(x, "splitvariable", value = y) ,
                  out, ID, SIMPLIFY = FALSE)
    output <- do.call("rbind", out)  
    
    outnames <- data.frame(do.call('rbind', strsplit(as.character(output$splitvariable),'##',fixed=TRUE)))
    names(outnames) = c(replicate.var, time.var)
    outnames <- outnames[1]
    
    output$splitvariable <- NULL
    output <- cbind(outnames, output)
    
  }
  return(output)
}

### PRIVATE FUNCTIONS ###


## function for the richness and evenness differences, gains and losses, and rankshifts returning a dataframe with those and the MRSc output
#rename this
SERGL <- function(df, time.var, rank.var1, rank.var2, abundance.var1, abundance.var2){
  #ricness and evenness differences
  s_t1 <- S(df[[abundance.var1]])
  e_t1 <- EQ(as.numeric(df[[abundance.var1]]))
  s_t2 <- S(df[[abundance.var2]])
  e_t2 <- EQ(as.numeric(df[[abundance.var2]]))
  
  sdiff <- abs(s_t1-s_t2)/nrow(df)
  ediff <- abs(e_t1-e_t2)/nrow(df)
  
  #gains and losses
  df$gain <- ifelse(df[[abundance.var1]]==0, 1, 0)
  df$loss <- ifelse(df[[abundance.var2]]==0, 1, 0)
  
  gain <- sum(df$gain)/nrow(df)
  loss <- sum(df$loss)/nrow(df)
  
  time1.1<-unique(df$time1)
  time1.2<-unique(df[[time.var]])
  
  mrsc <- mean(abs(df[[rank.var1]]-df[[rank.var2]])/nrow(df))
  
  time_pair <- paste(time1.1, time1.2, sep="-")
  
  metrics <- data.frame(time=time_pair, richness_change=sdiff, evenness_change=ediff, rank_change=mrsc, gains=gain, losses=loss)
  
  colnames(metrics)[1]<-paste(time.var, "pair", sep="_")
  
  return(metrics)
}
