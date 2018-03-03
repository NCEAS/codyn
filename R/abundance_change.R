#' @title Species Abundance Changes
#' @description Calculates the abundance change for species in a replicate between two consecutive time points.
#' @param df A data frame containing time, species, and abundance columns and an optional column of replicates
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
#'  \item{change: }{A numeric column of the change in abundance between consecutive timepoints. A postive value occurs when a species increases in abundnace over time, and a negative value when a species decreases in abundance over time.}
#' }
#' @references Avolio et al. OUR PAPER
#' @examples 
#' data(pplots)
#' # Without replicates
#' df <- subset(pplots, plot == 25)
#' abundance_change(df = df,
#'                  species.var = "species",
#'                  abundance.var = "relative_cover",
#'                  time.var = "year")
#'
#' # With replicates
#' df <- subset(pplots, year < 2004 & plot %in% c(6, 25, 32))
#' abundance_change(df = df,
#'                  species.var = "species",
#'                  abundance.var = "relative_cover",
#'                  replicate.var = "plot",
#'                  time.var = "year")
#' @export
abundance_change <- function(df, time.var, 
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
  
  # add zeros for species absent from a time period within a replicate
  if (is.null(replicate.var)) {
    allsp <- fill_zeros(df, species.var, abundance.var)
  } else {
    by <- c(replicate.var)
    allsp <- do.call(rbind, c(
      lapply(split(df, df[by], drop = TRUE),
        FUN = fill_zeros, species.var, abundance.var),
      list(make.row.names = FALSE)))
  }

  # rank species in each time and optionally replicate
  by <- c(time.var, replicate.var)
  rankdf <- do.call(rbind, c(
    lapply(split(allsp, allsp[by], drop = TRUE),
           FUN = add_ranks, species.var, abundance.var),
    list(make.row.names = FALSE)))

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
  mergevars <- c(species.var, replicate.var, time.var)
  df12 <- merge(df1, df2,  by = mergevars, all = T)
  df12<-subset(df12, df12[[paste(abundance.var, ".x", sep = "")]] !=0 | df12[[paste(abundance.var, ".y", sep = "")]] !=0)
  
  # sort and apply turnover to all replicates for each time point
  splitvars <- c(replicate.var, time.var)
  X <- split(df12,
             df12[splitvars], 
             sep = "##", drop = TRUE)
  out <- lapply(X, 
                FUN = abundchange, time.var, species.var, paste(abundance.var, ".x", sep = ""), paste(abundance.var, ".y", sep = "")) 
  unsplit <- lapply(out, nrow)
  unsplit <- rep(names(unsplit), unsplit)
  output <- do.call(rbind, c(out, list(make.row.names = FALSE)))
  output[splitvars] <- do.call(rbind, strsplit(unsplit, '##'))
  
  output_order <- c(
    paste(time.var, 'pair', sep = '_'), ## FIXEM why not time.var time.var2 cols?
    replicate.var,
    species.var,
    'change')
  
  return(output[intersect(output_order, names(output))])
}

############################################################################
#
# Private functions: these are internal functions not intended for reuse.
# Future package releases may change these without notice. External callers
# should not use them.
#
############################################################################

# A function to calculate abundance changes for a species between two consecutive time points 
# @param df a dataframe
# @param species.var the name of the species column
# @param time.var the name of the time column
# @param abundance.var1 the name of the abundance column for the first time peroid
# @param abundance.var2 the name of the abundance column for the second time period
abundchange <- function(df, time.var, species.var, abundance.var1, abundance.var2){
  
  time1.1 <- unique(df$time1)
  time1.2 <- unique(df[[time.var]])
  
  df$change <- df[[abundance.var1]] - df[[abundance.var2]]
  df$time_pair <- paste(time1.1, time1.2, sep = "-")
  
  df <- subset(df, select = c("time_pair", species.var, "change"))
  
  colnames(df)[1] <- paste(time.var, "pair", sep="_")
  
  return(df)
}
