#' @title Curve Change
#' @description Calculates the area difference between two rank abundance curves between two time periods. If replicate is specified, it must be measured in both time points, otherwise it will be dropped for that time period comparison.
#' @param df A data frame containing time, species, and abundance columns and an optional column of replicates
#' @param time.var The name of the time column 
#' @param species.var The name of the species column 
#' @param abundance.var The name of the abundance column 
#' @param replicate.var The name of the optional replicate column. If specified, replicate must be unique within the dataset and cannot be nested within treatments or blocks.
#' @param reference.time The name of the optional time point that all other time points should be compared to (e.g. the first year of data). If not specified, each comparison is between consecutive time points (e.g. first to  second year, second to third year, etc.)
#' @return The curve_change function returns a data frame with the following attributes:
#' \itemize{
#'  \item{time.var: }{A column with the specified time.var and a second column, with '2' appended to the name. Time is subtracted from time2.}
#'  \item{curve_change: }{A numeric column of the change in curves between time points.}
#'  \item{replicate.var: }{A column that has same name and type as the replicate.var column, if specified.}
#' }
#' @references Avolio et al. Submitted to MEE
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
#'            
#' # With reference year
#' df <- subset(pplots, year < 2005 & plot %in% c(6, 25, 32))
#' curve_change(df = df,
#'            species.var = "species",
#'            abundance.var = "relative_cover",
#'            replicate.var = "plot",
#'            time.var = "year",
#'            reference.time = 2002)
#' @export
curve_change <- function(df, time.var, 
                         species.var, 
                         abundance.var, 
                         replicate.var = NULL,
                         reference.time = NULL) {

  # validate function call and purge extraneous columns
  args <- as.list(match.call()[-1])
  df <- do.call(check_args, args, envir = parent.frame())
  
  # add rank abundance function within each time step and optionally replicate
  by <- c(replicate.var, time.var)
  rankabunddf <- split_apply_combine(df, by, FUN = add_rank_abundance,
    species.var, abundance.var)

  # merge subsets on time difference of one time step
  cross.var <- time.var
  cross.var2 <- paste(cross.var, 2, sep = '')
  split_by <- c(replicate.var)
  merge_to <- !(names(rankabunddf) %in% split_by)
  if (is.null(reference.time)) {
    output <- split_apply_combine(rankabunddf, split_by, FUN = function(x) {
      y <- x[merge_to]
      cross <- merge(x, y, by = NULL, suffixes = c('', '2'))
      f <- factor(cross[[cross.var]])
      f2 <- factor(cross[[cross.var2]], levels = levels(f))
      idx <- (as.integer(f2) - as.integer(f)) == 1
      cross[idx, ]
    })
  } else {
    output <- split_apply_combine(rankabunddf, split_by, FUN = function(x) {
      y <- x[x[[time.var]] != reference.time, merge_to]
      x <- x[x[[time.var]] == reference.time, ]
      merge(x, y, by = NULL, suffixes = c('', '2'))
    })
  }
  
  # split on treatment pairs (and block if not null)
  output[['curve_change']] <- mapply(curve_dissim,
    output[['rankabund']], output[['rankabund2']])

  output_order <- c(
    time.var, paste(time.var, '2', sep = ''),
    replicate.var,
    'curve_change')
  
  return(output[intersect(output_order, names(output))])
}
