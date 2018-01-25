#' @title Curve Shape Changes
#' @description The curve_change function compares the areas of difference between two curves that are consecutive time steps. It is optional to do this for replicates. For this to work, there must be more than 1 species in a plot at both time points, and when included, the replicate must also be measured at both time points.
#' @param df A data frame containing time, species and abundance columns and an optional column of replicates
#' @param time.var The name of the time column 
#' @param species.var The name of the species column 
#' @param abundance.var The name of the abundance column 
#' @param replicate.var The name of the optional replicate column
#'  
#' @return The curve_change function returns a data frame with the following attributes:
#' \itemize{
#'  \item{replicate.var: }{A column that has same name and type as the replicate.var column, if specified.}
#'  \item{time.var_pair: }{A characteric column that has the time points to be compared, separated by a dash.}
#'  \item{curve_change: }{A numeric column of the differences in curves between time points.}
#' }
#' @export
#'
curve_change <- function(df, time.var, 
                         species.var, 
                         abundance.var, 
                         replicate.var = NULL) {
  

  if(is.null(replicate.var)) {
  
  df <- subset(df, select = c(time.var, species.var, abundance.var))
  relrank <- subset(df, df[[abundance.var]] != 0)
  relrank$rank <- ave(relrank[[abundance.var]], relrank[[time.var]], FUN = function(x) rank(-x, ties.method = "average"))
  relrank$maxrank = ave(relrank$rank, relrank[[time.var]], FUN = function(x) max(x))
  relrank$relrank = relrank$rank/relrank$maxrank
  relrank <- relrank[order(relrank[[time.var]], -relrank[[abundance.var]]),]
  relrank$cumabund <- ave(relrank[[abundance.var]], relrank[[time.var]], FUN = function(x) cumsum(x))
 
  timestep<-sort(unique(relrank[[time.var]]))
  cc_out<-data.frame()
  
  for(i in 1:(length(timestep)-1)) {
    subset_t1 <- relrank[relrank[[time.var]] == timestep[i],]
    subset_t2 <- relrank[relrank[[time.var]] == timestep[i+1],]
    subset_t12 <- rbind(subset_t1, subset_t2)
    
    output <- curvechange(subset_t12, time.var, relrank, cumabund)
    cc_out <- rbind(cc_out, output)
    
  }
  
  } else {
    df <- subset(df, select = c(time.var, species.var, abundance.var, replicate.var))
    relrank <- subset(df, df[[abundance.var]] != 0)
    relrank$rep_time <- paste(relrank[[replicate.var]], relrank[[time.var]], sep = "_")
    relrank$rank <- ave(relrank[[abundance.var]], relrank$rep_time, FUN = function(x) rank(-x, ties.method = "average"))
    relrank$maxrank = ave(relrank$rank, relrank$rep_time, FUN = function(x) max(x))
    relrank$relrank = relrank$rank/relrank$maxrank
    relrank <- relrank[order(relrank[[time.var]], relrank[[replicate.var]], -relrank[[abundance.var]]),]
    relrank$cumabund <- ave(relrank[[abundance.var]], relrank$rep_time, FUN = function(x) cumsum(x))
    
    timestep<-sort(unique(relrank[[time.var]]))
    cc_out<-data.frame()
    
    for(i in 1:(length(timestep)-1)) {
      
      subset_t1 <- relrank[relrank[[time.var]] == timestep[i],]
      subset_t2 <- relrank[relrank[[time.var]] == timestep[i+1],]
      subset_t12 <- rbind(subset_t1, subset_t2)
      
      ##dropping plots that were not measured both years
      plots_t1 <- as.data.frame(unique(subset_t1[[replicate.var]]))
      colnames(plots_t1)[1] <- replicate.var

      plots_t2 <- as.data.frame(unique(subset_t2[[replicate.var]]))
      colnames(plots_t2)[1] <- replicate.var
      
      plots_bothyrs <- merge(plots_t1, plots_t2, by=replicate.var)
      
      subset_t12_2 <- merge(plots_bothyrs, subset_t12, by=replicate.var)
      subset_t12_2[[replicate.var]]<-as.character(subset_t12_2[[replicate.var]])
      
      X <- split(subset_t12_2, subset_t12_2[[replicate.var]])
      out <- lapply(X, FUN=curvechange, time.var, relrank, cumabund) 
      ID <- unique(names(out))
      out <- mapply(function(x, y) "[<-"(x, replicate.var, value = y) ,
                    out, ID, SIMPLIFY = FALSE)
       output <- do.call("rbind", out)
      
       cc_out <- rbind(cc_out, output)
    }
  }
    
    return (cc_out)
}
    
############################################################################
#
# Private functions: these are internal functions not intended for reuse.
# Future package releases may change these without notice. External callers
# should not use them.
#
############################################################################
    
curvechange <- function(df, time.var, relrank, cumabund) {
    
    df <- df[order(df[[time.var]], df$cumabund),]
  
    timestep2 <- unique(df[[time.var]]) 

    df1 <- df[df[[time.var]] == timestep2[1],]
    df2 <- df[df[[time.var]] == timestep2[2],]
    sf1 <- stepfun(df1$relrank, c(0, df1$cumabund))
    sf2 <- stepfun(df2$relrank, c(0, df2$cumabund))
    r <- sort(unique(c(0, df1$relrank, df2$relrank)))
    h <- abs(sf1(r) - sf2(r))
    w <- c(diff(r), 0)
    CC=sum(w*h)
    time_pair <- paste(timestep2[1], timestep2[2], sep="-")
    
    output=data.frame(timepair=time_pair, curve_change=CC)
    colnames(output)[1]<-paste(time.var, "pair", sep = "_")
    
    return(output)
    } 
  