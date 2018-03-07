#' @title Curve Changes
#' @description Calculates the area difference between two rank abundance curves between two consecutive time periods. If replicate is specified, it must be measued in both time points, otherwise it will be dropped for that time period comparision. Also, a replicate must have more than a single species in both time periods.
#' @param df A data frame containing time, species, and abundance columns and an optional column of replicates
#' @param time.var The name of the time column 
#' @param species.var The name of the species column 
#' @param abundance.var The name of the abundance column 
#' @param replicate.var The name of the optional replicate column
#'  
#' @return The curve_change function returns a data frame with the following attributes:
#' \itemize{
#'  \item{time.var_pair: }{A characteric column that has the time points to be compared, separated by a dash.}
#'  \item{curve_change: }{A numeric column of the change in curves between time points.}
#'  \item{replicate.var: }{A column that has same name and type as the replicate.var column, if specified.}
#' }
#' @references Avolio et al.OUR PAPER
#' @examples 
#' data(pplots)
#' # Without replicates
#' df <- subset(pplots, plot == 25)
#' curve_change(df = df,
#'            species.var = "species",
#'            abundance.var = "relative_cover",
#'            time.var = "year")
#'
#' # With replicates
#' df <- subset(pplots, year < 2004 & plot %in% c(6, 25, 32))
#' curve_change(df = df,
#'            species.var = "species",
#'            abundance.var = "relative_cover",
#'            replicate.var = "plot",
#'            time.var = "year")
#' @export
curve_change <- function(df, time.var, 
                         species.var, 
                         abundance.var, 
                         replicate.var = NULL) {
  
  # check no NAs in abundance column
  if(any(is.na(df[[abundance.var]]))) stop("Abundance column contains missing values")

  # check unique species x time x replicate combinations
  if(is.null(replicate.var)){
    check_single_onerep(df, time.var, species.var)
  } else {
    check_single(df, time.var, species.var, replicate.var)
  }
  
  df <- subset(df, select = c(time.var, species.var, abundance.var, replicate.var))
  relrank <- subset(df, df[[abundance.var]] != 0)
  
  # add ranks within each time step and optionally replicate
  by <- c(replicate.var, time.var)
  relrank$rank <- ave(relrank[[abundance.var]], relrank[by], FUN = function(x) rank(-x, ties.method = "average"))
  
  relrank$maxrank = ave(relrank$rank, relrank[by], FUN = function(x) max(x))
  relrank$relrank = relrank$rank/relrank$maxrank
  if(is.null(replicate.var)){
    relrank <- relrank[order(relrank[[time.var]], -relrank[[abundance.var]]),]
  } else {
    relrank <- relrank[order(relrank[[time.var]], relrank[[replicate.var]], -relrank[[abundance.var]]),]
  }
  relrank$cumabund <- ave(relrank[[abundance.var]], relrank[by], FUN = function(x) cumsum(x))
 
  timestep<-sort(unique(relrank[[time.var]]))
  cc_out<-data.frame()
  
  for(i in 1:(length(timestep)-1)) {
    subset_t1 <- relrank[relrank[[time.var]] == timestep[i],]
    subset_t2 <- relrank[relrank[[time.var]] == timestep[i+1],]
    subset_t12 <- rbind(subset_t1, subset_t2)
    
  if(is.null(replicate.var)){
      output <- curvechange(subset_t12, time.var, relrank, cumabund)
      cc_out <- rbind(cc_out, output)
    } else {
      
      #dropping plots that were not measured both years
      plots_t1 <- as.data.frame(unique(subset_t1[[replicate.var]]))
      colnames(plots_t1)[1] <- replicate.var
      
      plots_t2 <- as.data.frame(unique(subset_t2[[replicate.var]]))
      colnames(plots_t2)[1] <- replicate.var
      
      plots_bothyrs <- merge(plots_t1, plots_t2, by=replicate.var)
      
      subset_t12_2 <- merge(plots_bothyrs, subset_t12, by=replicate.var)
      subset_t12_2[[replicate.var]]<-as.character(subset_t12_2[[replicate.var]])
      
      # doing curve change for each replicate between consecutive time periods
      splitvars <- c(replicate.var)
      X <- split(subset_t12_2, subset_t12_2[splitvars])
      out <- lapply(X, FUN = curvechange, time.var, relrank, cumabund)
      unsplit <- lapply(out, nrow)
      unsplit <- rep(names(unsplit), unsplit)
      output <- do.call(rbind, c(out, list(make.row.names = FALSE)))
      output[by] <- do.call(rbind, as.list(unsplit))
      cc_out <- rbind(cc_out, output)
    }
  }
    
  output_order <- c(
    paste(time.var,"pair", sep="_"),
    replicate.var,
    'curve_change')
  
  return(cc_out[intersect(output_order, names(cc_out))])
  
}
    
############################################################################
#
# Private functions: these are internal functions not intended for reuse.
# Future package releases may change these without notice. External callers
# should not use them.
#
############################################################################

# A function calculate the curve changes between two time peroids
# @param df a dataframe
# @param time.var the name of the time column
# @param relrank the name of the relative rank of each species in the sample
# @param cumabund the name of the cumulative abundance of each species in the sample  
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
    CC <- sum(w*h)
    time_pair <- paste(timestep2[1], timestep2[2], sep = "-")
    
    output <- data.frame(timepair = time_pair, curve_change = CC)
    colnames(output)[1]<-paste(time.var, "pair", sep = "_")
    
    return(output)
}
