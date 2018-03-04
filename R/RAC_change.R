#' @title Rank Abundance Curve Changes
#' @description Calculates change of the five aspects of rank abundance curves
#'   (richness, evenness, rank, species gains, and species losses) for a
#'   replicate between two consecutive time points.
#' @param df A data frame containing time, species, and abundance columns and an
#'   optional column of replicates
#' @param time.var The name of the time column 
#' @param species.var The name of the species column 
#' @param abundance.var The name of the abundance column 
#' @param replicate.var The name of the optional replicate column 
#' @return The RAC_change function returns a data frame with the following
#'   attributes:
#' \itemize{
#'  \item{replicate.var: }{A column that has same name and type as the
#'  replicate.var column, if replicate.var is specified.}
#'  \item{time.var_pair: }{A characteric column that has the time points to be
#'  compared, separated by a dash.}
#'  \item{richness_change: }{A numeric column that is the change in richness
#'  between the two consecutive time peroids for a repicate divided by the total
#'  number of species in both time periods.}
#'  \item{evenness_change: }{A numeric column that is the change in evenness
#'  (measured with EQ) between the two consecutive time peroids for a repicate
#'  divided by the total number of species in both time periods.}
#'  \item{rank_change: }{A numeric column that is the average change in rank of
#'  a species between the two consecutive time peroids for a repicate divided by
#'  the total number of species in both time periods. Species that are not
#'  present in both time periods are given the S+1 rank in the sample it is
#'  absent in, where S is the number of species in that sample.}
#'  \item{gains: }{A numeric column of the number of species that are present at
#'  time period 2 that were not present at time period 1 for a replicate divided
#'  by the total number of species in both time periods. This is equivelant to
#'  the turnover function with metric = "appearances".}
#'  \item{losses: }{A numeric column of the number of species that are not
#'  present at time period 2 but were  present at time period 1 for a replicate
#'  divided by the total number of species in both time periods. This is
#'  equivelant to the turnover function with metric = "disappearance".}
#' }
#' @references Avolio et al.OUR PAPER
#' @examples 
#' data(pplots)
#' # Without replicates
#' df <- subset(pplots, plot == 25)
#' RAC_change(df = df,
#'            species.var = "species",
#'            abundance.var = "relative_cover",
#'            time.var = "year")
#'
#' # With replicates
#' df <- subset(pplots, year < 2004 & plot %in% c(6, 25, 32))
#' RAC_change(df = df,
#'            species.var = "species",
#'            abundance.var = "relative_cover",
#'            replicate.var = "plot",
#'            time.var = "year")
#' @export

RAC_change <- function(df, time.var, species.var, abundance.var, replicate.var=NULL) {
  
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
  
  # merge subsets on time difference of one time step
  cross.var <- time.var
  cross.var2 <- paste(cross.var, 2, sep = '')
  split_by <- c(species.var, replicate.var)
  merge_on <- !(names(rankdf) %in% split_by)
  rankdf_split <- lapply(split(rankdf, rankdf[split_by], drop = TRUE),
    function(x) {
      y <- x[merge_on]
      cross <- merge(x, y, by = NULL, suffixes = c('', '2'))
      f <- factor(cross[[cross.var]])
      f2 <- factor(cross[[cross.var2]], levels = levels(f))
      idx <- (as.integer(f2) - as.integer(f)) == 1
      cross[idx, ]
    })
  ranktog <- do.call(rbind, c(rankdf_split, list(make.row.names = FALSE)))

  # remove species not present in either year
  abundance.var2 <- paste(abundance.var, "2", sep = "")
  idx <- ranktog[[abundance.var]] != 0 | ranktog[[abundance.var2]] != 0
  ranktog <- ranktog[idx, ]
  
  # apply turnover calculation to all replicates for each time point
  by <- c(replicate.var, time.var)
  ranktog_split <- lapply(split(ranktog, ranktog[by], drop = TRUE),
                          FUN = SERGL,
                          species.var, abundance.var, abundance.var2)
  output <- do.call(rbind, c(ranktog_split, list(make.row.names = FALSE)))
  
  output_order <- c(
    time.var, paste(time.var, 2, sep = ''),
    replicate.var,
    'richness_change', 'evenness_change', 'rank_change', 'gains', 'losses')
  
  return(output[intersect(output_order, names(output))])
}

############################################################################
#
# Private functions: these are internal functions not intended for reuse.
# Future package releases may change these without notice. External callers
# should not use them.
#
############################################################################

# A function to calculate RAC changes for a replicate between two consecutive time points 
# @param df a dataframe
# @param time.var the name of the time column
# @param rank.var1 the name of the speices rank column for the first time peroid
# @param rank.var2 the name of the species rank column for the second time period
# @param abundance.var1 the name of the abundance column for the first time peroid
# @param abundance.var2 the name of the abundance column for the second time period
SERGL <- function(df, species.var, abundance.var, abundance.var2) {
  
  out <- c(species.var, 'rank', 'rank2', abundance.var, abundance.var2)
  out <- unique(df[!(names(df) %in% out)])
  if (nrow(out) != 1)
    stop('Input df has not been correctly split.')
  
  # ricness and evenness differences
  s_t1 <- S(df[[abundance.var]])
  e_t1 <- EQ(as.numeric(df[[abundance.var]]))
  s_t2 <- S(df[[abundance.var2]])
  e_t2 <- EQ(as.numeric(df[[abundance.var2]]))
  
  sdiff <- abs(s_t1-s_t2) / nrow(df)
  ediff <- abs(e_t1-e_t2) / nrow(df)
  
  # gains and losses
  df$gain <- ifelse(df[[abundance.var]] == 0, 1, 0)
  df$loss <- ifelse(df[[abundance.var2]] == 0, 1, 0)
  
  gain <- sum(df$gain) / nrow(df)
  loss <- sum(df$loss) / nrow(df)
  
  mrsc <- mean(abs(df[['rank']] - df[['rank2']]) / nrow(df))
  
  metrics <- data.frame(richness_change = sdiff, evenness_change = ediff,
                        rank_change = mrsc, gains = gain, losses = loss)
  
  return(cbind(out, metrics))
}
