#' @title  Abundance Differences
#' 
#' @description Calculates the abundnace difference for species between two
#'   samples. There are three ways differences can be calculated. 1) Between
#'   treatments within a block (note: block.var and treatment.var need to be
#'   specified). 2) Between treatments, pooling all replicates into a single
#'   species pool (note: pool = TRUE, treatment.var needs to be specified, and
#'   block.var = NULL). 3) All pairwise combinations between all replicates
#'   (note: block.var = NULL, pool = FALSE and specifying treatment.var is
#'   optional. If treatment.var is specified, the treatment that each replicate
#'   belongs to will also be listed in the output).
#'   
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
#'  \item{species.var: }{A column that has same name and type as the species.var
#'  column.}
#'  \item{difference: }{A numeric column of the abundance differences between
#'  the two samples being compared (replicates or treatments). A numeric column
#'  of the change in abundance between consecutive timepoints. A postive value
#'  occurs when a species has greater abundance in replicate.var2 than in
#'  replicate.var and/or in treatment.var2 than in treatment.var.}
#'  \item{replicate.var: }{A column that has same name and type as the
#'  replicate.var column, represents the first replicate being compared. Note, a
#'  replicate column will be returned only when pool = FALSE or block.var =
#'  NULL.}
#'  \item{replicate.var2: }{A column that has the same type as the replicate.var
#'  column, and is named replicate.var with a 2 appended to it, represents the
#'  second replicate being compared. Note, a replicate.var column will be
#'  returned only when pool = FALSE and block.var = NULL.}
#'  \item{time.var: }{A column that has the same name and type as the time.var
#'  column, if time.var is specified.}
#'  \item{treatment.var: }{A column that has same name and type as the
#'  treatment.var column, represents the first treatment being compared. A
#'  treatment.var column will be returned when pool = TRUE, block.var is
#'  specified, or treatment.var is specified.}
#'  \item{treatment.var2: }{A column that has the same type as the treatment.var
#'  column, and is named treatment.var with a 2 appended to it, represents the
#'  second treatment being compared. A treatment.var column will be returned
#'  when pool = TRUE, block.var is specified, or treatment.var is specified.}
#'  \item{block.var: }{A column that has same name and type as the block.var
#'  column, if block.var is specified.}
#' }
#' 
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
#' df <- subset(pplots, year < 2004 & plot %in% c(21, 25, 32))
#' abundance_difference(df = df,
#'                      species.var = "species",
#'                      abundance.var = "relative_cover",
#'                      replicate.var = "plot",
#'                      time.var = "year",
#'                      treatment.var = "treatment")
#' 
#' # All pairwise replicates without treatment
#' df <- subset(pplots, year < 2004 & plot %in% c(21, 25, 32))
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
  
  # drop extraneous columns
  args <- as.list(match.call())
  df <- as.data.frame(df[as.character(args[grep('\\.var$', names(args))])])

  # check no NAs in abundance column
  if(any(is.na(df[[abundance.var]]))) stop("Abundance column contains missing values")

  # check no NAs in species column
  if(any(is.na(df[[species.var]]))) stop("Species names are missing")
  
  # check no species are repeated
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
    # pool and rank species in each replicate
    allsp <- pool_replicates(df, time.var, species.var, abundance.var,
      replicate.var, treatment.var)
  } else {
    # add NA for species absent from a replicate
    by <- c(block.var, time.var)
    allsp <- split_apply_combine(df, by, FUN = fill_species,
      species.var, abundance.var)
  }

  # order cross.var if unordered factor
  to_ordered = is.factor(df[[cross.var]]) & !is.ordered(df[[cross.var]])
  if (to_ordered) {
    class(allsp[[cross.var]]) <- c('ordered', class(df[[cross.var]]))
  }
    
  # cross join for pairwise comparisons
  split_by <- c(block.var, time.var)
  merge_to <- !(names(allsp) %in% split_by)
  cross.var2 <- paste(cross.var, 2, sep = '')
  ranktog <- split_apply_combine(allsp, split_by, FUN = function(x) {
    y <- x[merge_to]
    cross <- merge(x, y, by = species.var, suffixes = c('', '2'))
    idx <- cross[[cross.var]] < cross[[cross.var2]]
    cross[idx, ]
  })

  # unorder cross.var if orginally unordered factor
  if (to_ordered) {
    class(ranktog[[cross.var]]) <- class(df[[cross.var]])
  }
  
  # remove rows with NA for both abundances (preferably only when introduced
  # by fill_species)
  idx <- is.na(ranktog[[abundance.var]])
  abundance.var2 <- paste(abundance.var, 2, sep = '')
  idx2 <- is.na(ranktog[[abundance.var2]])
  ranktog[idx, abundance.var] <- 0
  ranktog[idx2, abundance.var2] <- 0
  ranktog <- ranktog[!(idx & idx2), ]

  # take abundance difference
  output <- abund_diff(ranktog, species.var, abundance.var, abundance.var2)

  # order and select output columns
  output_order <- c(
    time.var,
    block.var,
    replicate.var, paste(replicate.var, 2, sep = ''),
    treatment.var, paste(treatment.var, 2, sep = ''),
    species.var,
    'difference')

  return(output[intersect(output_order, names(output))])
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
abund_diff <- function(df, species.var, abundance.var, abundance.var2) {

  df[['difference']] <- df[[abundance.var2]] - df[[abundance.var]]
  df[c(abundance.var, abundance.var2)] <- NULL

  return(df)
}
