#' @title Rank Abundance Curve Changes
#' @description Calculates change of the five aspects of rank abundance curves (richness, evenness, rank, species gains, and species losses) for a replicate between two time points.
#' @param df A data frame containing time, species, and abundance columns and an optional columns of replicates.
#' @param time.var The name of the time column 
#' @param species.var The name of the species column 
#' @param abundance.var The name of the abundance column 
#' @param replicate.var The name of the optional replicate column. If specified, replicate must be unique within the dataset and cannot be nested within treatments or blocks.
#' @param reference.time The name of the optional time point that all other time points should be compared to (e.g. the first year of data). If not specified, each comparison is between consecutive time points (the first and second year, second and third year, etc.)
#' @return The RAC_change function returns a data frame with the following attributes:
#' \itemize{
#'  \item{replicate.var: }{A column that has same name and type as the replicate.var column, if replicate.var is specified.}
#'  \item{time.var: }{A column with the specified time.var and a second column, with '2' appended to the name. Time is subtracted from time2.}
#'  \item{richness_change: }{A numeric column that is the change in richness between the two time periods for a replicate divided by the total number of unique species in both time periods. A positive value occurs when a there is an increase in species richness over time, and a negative value when there is a decreases in species richness over  time.}
#'  \item{evenness_change: }{A numeric column that is the change in evenness(measured with Evar) between the two time periods for a replicate. A positive value occurs when evenness  increases over time, and a negative value when evenness decreases in over time.}
#'  \item{rank_change: }{A numeric column that is the absolute value of the average change in species ranks between the two time periods for a replicate divided by the total number of unique species in both time periods. Species that are not present in both time periods are given the S+1 rank in the sample it is absent in, where S is the number of species in that sample.}
#'  \item{gains: }{A numeric column of the number of species that are present at time period 2 that were not present at time period 1 for a replicate divided by the total number of unique species in both time periods. This is equivalent to the turnover function with metric = "appearances".}
#'  \item{losses: }{A numeric column of the number of species that are not present at time period 2 but were  present at time period 1 for a replicate divided by the total number of unique species in both time periods. This is equivalent to the turnover function with metric = "disappearance".}
#' }
#' @references Avolio et al. submitted to MEE
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
#'            
#' # With reference year
#' df <- subset(pplots, year < 2005 & plot %in% c(6, 25, 32))
#' RAC_change(df = df,
#'            species.var = "species",
#'            abundance.var = "relative_cover",
#'            replicate.var = "plot",
#'            time.var = "year",
#'            reference.time = 2002)
#'
#' @export

RAC_change <- function(df,
                       time.var,
                       species.var,
                       abundance.var,
                       replicate.var = NULL,
                       reference.time = NULL) {

  # validate function call and purge extraneous columns
  args <- as.list(match.call()[-1])
  df <- do.call(check_args, args, envir = parent.frame())
  
  # add zeros for species absent from a time period within a replicate
  by <- c(replicate.var)
  allsp <- split_apply_combine(df, by, FUN = fill_species, species.var, abundance.var)

  # rank species in each time and optionally replicate
  by <- c(time.var, replicate.var)
  rankdf <- split_apply_combine(allsp, by, FUN = add_ranks, abundance.var)

  # merge subsets on time difference of one time step
  cross.var <- time.var
  cross.var2 <- paste(cross.var, 2, sep = '')
  split_by <- c(replicate.var)
  merge_to <- !(names(rankdf) %in% split_by)
  if (is.null(reference.time)) {
    ranktog <- split_apply_combine(rankdf, split_by, FUN = function(x) {
      y <- x[merge_to]
      cross <- merge(x, y, by = species.var, suffixes = c('', '2'))
      f <- factor(cross[[cross.var]])
      f2 <- factor(cross[[cross.var2]], levels = levels(f))
      idx <- (as.integer(f2) - as.integer(f)) == 1
      cross[idx, ]
    })
  } else {
    ranktog <- split_apply_combine(rankdf, split_by, FUN = function(x) {
      y <- x[x[[time.var]] != reference.time, merge_to]
      x <- x[x[[time.var]] == reference.time, ]
      merge(x, y, by = species.var, suffixes = c('', '2'))
    })
  }
  # remove rows with NA for both abundances (preferably only when introduced
  # by fill_species)
  idx <- is.na(ranktog[[abundance.var]])
  abundance.var2 <- paste(abundance.var, 2, sep = '')
  idx2 <- is.na(ranktog[[abundance.var2]])
  ranktog[idx, abundance.var] <- 0
  ranktog[idx2, abundance.var2] <- 0
  idx <- ranktog[[abundance.var]] != 0 | ranktog[[abundance.var2]] != 0
  ranktog <- ranktog[idx, ]
  
  # apply turnover calculation to all replicates for each time point
  by <- c(replicate.var, time.var)
  output <- split_apply_combine(ranktog, by, FUN = SERGL,
    species.var, abundance.var, abundance.var2)
  
  if(any(is.na(output$evenness_change)))
    warning(paste0("evenness_change values contain NAs because there are plots",
                   " with only one species"))

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

  # ricness and evenness differences
  s_t1 <- S(df[[abundance.var]])
  e_t1 <- Evar(as.numeric(df[[abundance.var]]))
  s_t2 <- S(df[[abundance.var2]])
  e_t2 <- Evar(as.numeric(df[[abundance.var2]]))
  
  delta_s <- (s_t2-s_t1) / nrow(df)
  delta_e <- (e_t2-e_t1) 
  
  # gains and lqosses
  df$gain <- ifelse(df[[abundance.var]] == 0, 1, 0)
  df$loss <- ifelse(df[[abundance.var2]] == 0, 1, 0)
  
  gain <- sum(df$gain) / nrow(df)
  loss <- sum(df$loss) / nrow(df)
  
  delta_r <- mean(abs(df[['rank']] - df[['rank2']])) / nrow(df)
  
  metrics <- data.frame(richness_change = delta_s, evenness_change = delta_e, rank_change = delta_r, gains = gain, losses = loss)
  
  return(cbind(out, metrics))
}
