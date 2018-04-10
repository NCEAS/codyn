#' @title Species Abundance Changes
#' @description Calculates the abundance change for species in a replicate
#'   between two consecutive time points.
#' @param df A data frame containing time, species, and abundance columns and an
#'   optional column of replicates
#' @param time.var The name of the time column 
#' @param species.var The name of the species column 
#' @param abundance.var The name of the abundance column 
#' @param replicate.var The name of the optional replicate column 
#' 
#' @return The abundance_change function returns a data frame with the following
#'   fields:
#' \itemize{
#'  \item{replicate.var: }{A column with the specified replicate.var, if it is specified.}
#'  \item{time.var: }{A column with the specified time.var and a second column, with '2'
#'  appended to the name, giving the time of the subtracted abundance.}
#'  \item{species.var: }{A column with the specified species.var.}
#'  \item{change: }{A numeric column of the change in abundance between consecutive
#'  timepoints. A postive value occurs when a species increases in abundnace
#'  over time, and a negative value when a species decreases in abundance over
#'  time.}
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
  
  # add zeros for species absent from a time period within a replicate
  by <- c(replicate.var)
  allsp <- split_apply_combine(df, by, FUN = fill_species, species.var, abundance.var)

  # merge subsets on time difference of one time step
  cross.var <- time.var
  cross.var2 <- paste(cross.var, 2, sep = '')
  split_by <- c(replicate.var)
  merge_on <- !(names(allsp) %in% split_by)
  ranktog <- split_apply_combine(allsp, split_by, FUN = function(x) {
      y <- x[merge_on]
      cross <- merge(x, y, by = species.var, suffixes = c('', '2'))
      f <- factor(cross[[cross.var]])
      f2 <- factor(cross[[cross.var2]], levels = levels(f))
      idx <- (as.integer(f2) - as.integer(f)) == 1
      cross[idx, ]
  })
  
  # remove rows with NA for both abundances (preferably only when introduced
  # by fill_species)
  idx <- is.na(ranktog[[abundance.var]])
  abundance.var2 <- paste(abundance.var, 2, sep = '')
  idx2 <- is.na(ranktog[[abundance.var2]])
  ranktog[idx, abundance.var] <- 0
  ranktog[idx2, abundance.var2] <- 0
  output <- ranktog[!(idx & idx2), ]

  # append change column
  output[['change']] <- output[[abundance.var2]] - output[[abundance.var]]
  output_order <- c(
    time.var, paste(time.var, '2', sep = ''),
    replicate.var,
    species.var,
    'change')
  
  return(output[intersect(output_order, names(output))])
}
