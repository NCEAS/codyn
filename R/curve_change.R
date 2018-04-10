#' @title Curve Changes
#' @description Calculates the area difference between two rank abundance curves
#'   between two consecutive time periods. If replicate is specified, it must be
#'   measued in both time points, otherwise it will be dropped for that time
#'   period comparision. Also, a replicate must have more than a single species
#'   in both time periods.
#' @param df A data frame containing time, species, and abundance columns and an
#'   optional column of replicates
#' @param time.var The name of the time column 
#' @param species.var The name of the species column 
#' @param abundance.var The name of the abundance column 
#' @param replicate.var The name of the optional replicate column
#'  
#' @return The curve_change function returns a data frame with the following
#'   attributes:
#' \itemize{
#'  \item{time.var_pair: }{A characteric column that has the time points to be
#'  compared, separated by a dash.}
#'  \item{curve_change: }{A numeric column of the change in curves between time
#'  points.}
#'  \item{replicate.var: }{A column that has same name and type as the
#'  replicate.var column, if specified.}
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

  # drop extraneous columns
  args <- as.list(match.call())
  df <- as.data.frame(df[as.character(args[grep('\\.var$', names(args))])])
  
  # check no NAs in abundance column
  if(any(is.na(df[[abundance.var]]))) stop("Abundance column contains missing values")
  
  # check no NAs in species column
  if(any(is.na(df[[species.var]]))) stop("Species names are missing")

  # check unique species x time x replicate combinations
  if(is.null(replicate.var)){
    check_single_onerep(df, time.var, species.var)
  } else {
    check_single(df, time.var, species.var, replicate.var)
  }
  
  df <- subset(df, select = c(time.var, species.var, abundance.var, replicate.var))
  df <- subset(df, df[[abundance.var]] != 0)
  
  # add rank abundance function within each time step and optionally replicate
  by <- c(replicate.var, time.var)
  rankabunddf <- split_apply_combine(df, by, FUN = add_rank_abundance,
    species.var, abundance.var)

  # merge subsets on time difference of one time step
  cross.var <- time.var
  cross.var2 <- paste(cross.var, 2, sep = '')
  split_by <- c(replicate.var)
  merge_on <- !(names(rankabunddf) %in% split_by)
  output <- split_apply_combine(rankabunddf, split_by, FUN = function(x) {
    y <- x[merge_on]
    cross <- merge(x, y, by = NULL, suffixes = c('', '2'))
    f <- factor(cross[[cross.var]])
    f2 <- factor(cross[[cross.var2]], levels = levels(f))
    idx <- (as.integer(f2) - as.integer(f)) == 1
    cross[idx, ]
  })
  
  # split on treatment pairs (and block if not null)
  output[['curve_change']] <- mapply(curve_dissim,
    output[['rankabund']], output[['rankabund2']])

  output_order <- c(
    time.var, paste(time.var, '2', sep = ''),
    replicate.var,
    'curve_change')
  
  return(output[intersect(output_order, names(output))])
}
