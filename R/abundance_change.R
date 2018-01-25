#' @title Species Abundance Changes
#' @description 
#' @param df A data frame containing time, species and abundance columns and an optional column of replicates
#' @param time.var The name of the time column 
#' @param species.var The name of the species column 
#' @param abundance.var The name of the abundance column 
#' @param replicate.var The name of the optional replicate column 
#' 
#' @return The abundance_change function returns a data frame with the following attributes:
#' \itemize{
#'  \item{replicate.var: }{A column that has same name and type as the replicate.var column, if replicate.var is specified.}
#'  \item{time.var_pair: }{A characteric column that has the time points to be compared, separated by a dash.}
#'  \item{species.var: }{A column that has same name and type as the species.var column.}
#'  \item{abund_change: }{A numeric column of the abundance differences between timepoints}
#' }
#' @details 
#' @references 
#' @example 
#' @export


abundance_change <- function(df, time.var, 
                             species.var, 
                             abundance.var, 
                             replicate.var = NULL) {
  
  if(is.null(replicate.var)) {
    
    rankdf <- add_ranks_time(df, time.var, species.var, abundance.var, replicate.var = NULL)
    
    time1 <- sort(unique(rankdf[[time.var]]))
    time2 <- c(time1[2:length(time1)], NA)
    
    # current year rankdf
    df2 <- rankdf[which(rankdf[[time.var]]%in%time2),]
    
    # previous year rank df
    mytimes <- data.frame(cbind(time1, time2))
    names(mytimes) = c(time.var, "dummytime")
    
    df1 <- merge(rankdf, mytimes)
    names(df1)[1] = 'time1'
    names(df1)[[ncol(df1)]] = time.var
    df1 <- subset(df1, !is.na(df1[[time.var]]))
    
    # merge: .x is for previous time point, .y for current time point, time.var corresponds to current (i.e., .y)
    df12 <- merge(df1, df2,  by = c(species.var, time.var), all = T)
    df12<-subset(df12, df12[[paste(abundance.var, ".x", sep = "")]] != 0 | df12[[paste(abundance.var, ".y", sep = "")]] != 0)
    df12<-subset(df12, !is.na(df12[[paste(abundance.var, ".x", sep = "")]]) & !is.na(df12[[paste(abundance.var, ".y", sep = "")]]))
    
    # sort and apply turnover to all replicates
    df12 <- df12[order(df12[[time.var]]),]
    X <- split(df12, df12[[time.var]])
    
    
    out <- lapply(X, FUN = abundchange, time.var, species.var, paste(abundance.var, ".x", sep = ""),paste(abundance.var, ".y", sep = "")) 
    output <- do.call("rbind", out)  

  } else {
    
    rankdf <- add_ranks_time(df,  time.var, species.var, abundance.var, replicate.var)
    
    time1 <- sort(unique(rankdf[[time.var]]))
    time2 <- c(time1[2:length(time1)], NA)
    
    # current year rankdf
    df2 <- rankdf[which(rankdf[[time.var]]%in%time2),]
    
    # previous year rank df
    mytimes <- data.frame(cbind(time1, time2))
    names(mytimes) = c(time.var, "dummytime")
    
    df1 <- merge(rankdf, mytimes)
    names(df1)[1] = 'time1'
    names(df1)[[ncol(df1)]] = time.var
    df1 <- subset(df1, !is.na(df1[[time.var]]))
    
    # merge: .x is for previous time point, .y for current time point, time.var corresponds to current (i.e., .y)
    df12 <- merge(df1, df2,  by = c(species.var,replicate.var, time.var), all = T)
    df12<-subset(df12, df12[[paste(abundance.var, ".x", sep = "")]] !=0 | df12[[paste(abundance.var, ".y", sep = "")]] !=0)
    df12<-subset(df12, !is.na(df12[[paste(abundance.var, ".x", sep = "")]]) & !is.na(df12[[paste(abundance.var, ".y", sep = "")]]))
    df12$splitvariable <- paste(df12[[replicate.var]], df12[[time.var]], sep = "##") 
    
    # sort and apply turnover to all replicates
    df12 <- df12[order(df12$splitvariable),]
    X <- split(df12, df12$splitvariable)
    
    
    out <- lapply(X, FUN = abundchange, time.var, species.var, paste(abundance.var, ".x", sep = ""),paste(abundance.var, ".y", sep = "")) 
    ID <- unique(names(out))
    out <- mapply(function(x, y) "[<-"(x, "splitvariable", value = y) ,
                  out, ID, SIMPLIFY = FALSE)
    output <- do.call("rbind", out)  
    
    outnames <- data.frame(do.call('rbind', strsplit(as.character(output$splitvariable), '##', fixed = TRUE)))
    names(outnames) = c(replicate.var, time.var)
    outnames <- outnames[1]
    
    output$splitvariable <- NULL
    output <- cbind(outnames, output)
    
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

abundchange <- function(df, time.var, species.var, abundance.var1, abundance.var2){

  time1.1 <- unique(df$time1)
  time1.2 <- unique(df[[time.var]])
  
  df$abund_change <- df[[abundance.var1]] - df[[abundance.var2]]
  df$time_pair <- paste(time1.1, time1.2, sep = "-")
  
  df <- subset(df, select = c("time_pair", species.var, "abund_change"))
  
  colnames(df)[1] <- paste(time.var, "pair", sep="_")
  
  return(df)
}
