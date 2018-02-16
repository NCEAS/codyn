#' @title  Abundance Differences
#' @description Calculates the abundnace difference for species between two samples. There are three ways differences can be calculated. 1) Between treatments within a block (note: block.var and treatment.var need to be specified). 2) Between treatments, pooling all replicates into a single species pool (note: pool = TRUE, treatment.var needs to be specified, and block.var = NULL). 3) All pairwise combinations between all replicates (note: block.var = NULL, pool = FALSE and specifying treatment.var is optional. If treatment.var is specified, the treatment that each replicate belongs to will also be listed in the output).
#' @param df A data frame containing species, abundance, replicate columns and optional time, treatment and block columns.
#' @param time.var The name of the optional time column 
#' @param species.var The name of the species column 
#' @param abundance.var The name of the abundance column 
#' @param replicate.var The name of the replicate column 
#' @param treatment.var The name of the optional treatment column
#' @param block.var The name of the optional block column
#' @param pool An argument to allow abundance values to be pooled within a treatment. The default value is "FALSE", a value of "TRUE" averages abundance of each species within a treatment at a given time point.
#' 
#' @return The abundance_difference function returns a data frame with the following attributes:
#' \itemize{
#'  \item{species.var: }{A column that has same name and type as the species.var column.}
#'  \item{difference: }{A numeric column of the abundance differences between the two samples being compared (replicates or treatments).}
#'  \item{replicate.var: }{A column that has same name and type as the replicate.var column, represents the first replicate being compared. Note, a replicate column will be returned only when pool = FALSE or block.var = NULL.}
#'  \item{replicate.var2: }{A column that has the same type as the replicate.var column, and is named replicate.var with a 2 appended to it, represents the second replicate being compared. Note, a replicate.var column will be returned only when pool = FALSE and block.var = NULL.}
#'  \item{time.var: }{A column that has the same name and type as the time.var column, if time.var is specified.}
#'  \item{treatment.var: }{A column that has same name and type as the treatment.var column, represents the first treatment being compared. A treatment.var column will be returned when pool = TRUE, block.var is specified, or treatment.var is specified.}
#'  \item{treatment.var2: }{A column that has the same type as the treatment.var column, and is named treatment.var with a 2 appended to it, represents the second treatment being compared. A treatment.var column will be returned when pool = TRUE, block.var is specified, or treatment.var is specified.}
#'  \item{block.var: }{A column that has same name and type as the block.var column, if block.var is specified.}
#' }
#' @references Avolio et al. OUR PAPER
#' @examples 
#' data(pplots)
#' # With block and no time
#' df <- subset(pplots, year == 2002 & block < 3)
#' abundance_difference(df = df,
#'                      species.var = "species",
#'                      abundance.var = "relative_cover",
#'                      treatment.var = 'treatment',
#'                      block.var = "block",
#'                      replicate.var = "plot")
#' 
#' # With blocks and time
#' df <- subset(pplots, year < 2004 & block < 3)
#' abundance_difference(df = df,
#'                      species.var = "species",
#'                      abundance.var = "relative_cover",
#'                      treatment.var = 'treatment',
#'                      block.var = "block",
#'                      replicate.var = "plot",
#'                      time.var = "year")
#' 
#' # Pooling by treatment with time
#' df <- subset(pplots, year < 2004)
#' abundance_difference(df = df,
#'                      species.var = "species",
#'                      abundance.var = "relative_cover",
#'                      treatment.var = 'treatment',
#'                      pool = TRUE,
#'                      replicate.var = "plot",
#'                      time.var = "year")
#' 
#' # All pairwise replicates with treatment
#' df <- subset(pplots, year < 2004 & plot %in% c(6, 25, 32))
#' abundance_difference(df = df,
#'                      species.var = "species",
#'                      abundance.var = "relative_cover",
#'                      replicate.var = "plot",
#'                      time.var = "year",
#'                      treatment.var = "treatment")
#' 
#' # All pairwise replicates without treatment
#' df <- subset(pplots, year < 2004 & plot %in% c(6, 25, 32))
#' abundance_difference(df = df,
#'                      species.var = "species",
#'                      abundance.var = "relative_cover",
#'                      replicate.var = "plot",
#'                      time.var = "year")
#' @export
abundance_difference <- function(df, time.var = NULL, species.var, 
                                 abundance.var, replicate.var,
                                 treatment.var = NULL, pool = FALSE, 
                                 block.var = NULL) {

  # check no NAs in abundance column
  if(any(is.na(df[[abundance.var]]))) stop("Abundance column contains missing values")
  
  #check no species are repeated
  if (is.null(time.var)){
    # check there unique species x time combinations
    check_single_onerep(df, replicate.var, species.var)
  }
  else {
    # check unique species x time x replicate combinations
    check_single(df, time.var, species.var, replicate.var)
  }
  
  
  if (!is.null(block.var)) {
    reps_exp <- length(unique(df[[block.var]])) * length(unique(df[[treatment.var]]))
    reps_obs <- length(unique(df[[replicate.var]]))
    if (reps_exp != reps_obs)
      stop("There is not one replicate per treatment in a block")
    cross.var <- treatment.var
  } else if (pool) {
    cross.var <- treatment.var
  } else {
    cross.var <- replicate.var
  }

  if (pool) {
    rankdf <- pool_replicates(df, time.var, species.var, abundance.var, replicate.var, treatment.var)
  } else {
    # rank species in each replicate
    rep_trt <- unique(df[c(replicate.var, treatment.var, block.var)])
    rankdf <- add_ranks_replicate(df, time.var, species.var, abundance.var, replicate.var)
    rankdf <- merge(rankdf, rep_trt, by = replicate.var)
  }
  
  # cross join for pairwise comparisons
  splitvars <- c(species.var, block.var, time.var)
  cross.var2 <- paste(cross.var, 2, sep = '')
  rankdf <- lapply(split(rankdf, rankdf[splitvars]),
                   function(x) {
                     y <- x
                     y[splitvars] <- NULL
                     cross <- merge(x, y, by = NULL, suffixes = c('', '2'))
                     idx <- as.integer(cross[[cross.var]])
                     idx <- idx < as.integer(cross[[cross.var2]])
                     cross[idx,]
                   })
  ranktog <- do.call(rbind, c(rankdf, list(make.row.names = FALSE)))
  
  # split on treatment pairs (and block if not null)
  splitvars <- c(block.var, time.var, cross.var, cross.var2)
  ranktog_split <- split(ranktog,
                         ranktog[splitvars], 
                         sep = "##", drop = TRUE)
  ranktog_split <- lapply(ranktog_split,
                          FUN = abund_diff, species.var, abundance.var)
  unsplit <- lapply(ranktog_split, nrow)
  unsplit <- rep(names(unsplit), unsplit)
  output <- do.call(rbind, c(ranktog_split, list(make.row.names = FALSE)))
  output[splitvars] <- do.call(rbind, strsplit(unsplit, '##'))

  if (is.null(block.var) & !pool & !is.null(treatment.var)) {
    # add treatment for reference
    output <- merge(output, merge(rep_trt, rep_trt, by = NULL))
  }
  
  ## FIXME reset column types based on df

return(output)
  
}

############################################################################
#
# Private functions: these are internal functions not intended for reuse.
# Future package releases may change these without notice. External callers
# should not use them.
#
############################################################################

# A function to calculate abundance differences for a species between two samples 
# @param df a dataframe
# @param species.var the name of the species column
# @param abundance.var the name of the abundance column
abund_diff <- function(df, species.var, abundance.var) {

  abundance.var2 <- paste(abundance.var, 2, sep = '')
  df[['difference']] <- df[[abundance.var]] - df[[abundance.var2]]

  return(df[c(species.var, 'difference')])
}
