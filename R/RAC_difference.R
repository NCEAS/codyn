#' @title Rank Abundance Curve Differences
#' @description Calculates differences between two samples for four comparable
#'   aspects of rank abundance curves (richness, evenness, rank, species
#'   composition). There are three ways differences can be calculated. 1)
#'   Between treatments within a block (note: block.var and treatment.var need
#'   to be specified). 2) Between treatments, pooling all replicates into a
#'   single species pool (note: pool = TRUE, treatment.var needs to be
#'   specified, and block.var will be NULL). 3) All pairwise combinations
#'   between all replicates (note: block.var = NULL, pool = FALSE and specifying
#'   treatment.var is optional. If treatment.var is specified, the treatment
#'   that each replicate belongs to will also be listed in the output).
#' @param df A data frame containing a species, abundance, and replicate columns
#'   and optional time, treatment, and block columns
#' @param time.var The name of the optional time column 
#' @param species.var The name of the species column 
#' @param abundance.var The name of the abundance column 
#' @param replicate.var The name of the replicate column 
#' @param treatment.var The name of the optional treatment column
#' @param block.var The name of the optional block column
#' @param pool An argument to allow abundance values to be pooled within a
#'   treatment. The default value is "FALSE", a value of "TRUE" averages
#'   abundance of each species within a treatment at a given time point.
#' @return The RAC_difference function returns a data frame with the following attributes:
#' \itemize{
#'  \item{time.var: }{A column that has the same name and type as the time.var
#'  column, if time.var is specified.}
#'  \item{block.var: }{A column that has same name and type as the block.var
#'  column, if block.var is specified.}
#'  \item{replicate.var: }{A column that has same name and type as the
#'  replicate.var column, represents the first replicate being compared. Note, a
#'  replicate column will be returned only when pool is FALSE or block.var =
#'  NULL.}
#'  \item{replicate.var2: }{A column that has the same type as the replicate.var
#'  column, and is named replicate.var with a 2 appended to it, represents the
#'  second replicate being compared. Note, a replicate.var column will be
#'  returned only when pool is FALSE and block.var = NULL.}
#'  \item{treatment.var: }{A column that has same name and type as the
#'  treatment.var column, represents the first treatment being compared. A
#'  treatment.var column will be returned when pool is TRUE or block.var is
#'  present, or treatment.var is specified.}
#'  \item{treatment.var2: }{A column that has the same type as the treatment.var
#'  column, and is named treatment.var with a 2 appended to it, represents the
#'  second treatment being compared. A treatment.var column will be returned
#'  when pool is TRUE or block.var is present, or treatment.var is specified.}
#'  \item{richness_diff: }{A numeric column that is the difference between the
#'  compared samples (treatments or replicates) in species richness divided by
#'  the total number of species in both samples.}
#'  \item{evenness_diff: }{A numeric column of the difference between the
#'  compared samples (treatments or replicates) in evenness (measured using the
#'  EQ metric) divided by the total number of species in both samples.}
#'  \item{rank_diff: }{A numeric column of the average difference between the
#'  compared samples (treatments or replicates) in species' ranks divided by the
#'  total number of species in both samples. Species that are not present in
#'  both samples are given the S+1 rank in the sample it is absent in, where S
#'  is the number of species in that sample.}
#'  \item{species_diff: }{A numeric column of the number of species that are
#'  different between the compared samples (treatments or replicates) divided by
#'  the total number of species in both samples. This is equivelant to the
#'  Jaccard Index.}
#' }
#' @references Avolio et al. OUR PAPER
#' @examples 
#' data(pplots)
#' # With block and no time
#' df <- subset(pplots, year == 2002 & block < 3)
#' RAC_difference(df = df,
#'                species.var = "species",
#'                abundance.var = "relative_cover",
#'                treatment.var = 'treatment',
#'                block.var = "block",
#'                replicate.var = "plot")
#'
#' # With blocks and time
#' df <- subset(pplots, year < 2004 & block < 3)
#' RAC_difference(df = df,
#'                species.var = "species",
#'                abundance.var = "relative_cover",
#'                treatment.var = 'treatment',
#'                block.var = "block",
#'                replicate.var = "plot",
#'                time.var = "year")
#' 
#' # Pooling by treatment with time
#' df <- subset(pplots, year < 2004)
#' RAC_difference(df = df,
#'                species.var = "species",
#'                abundance.var = "relative_cover",
#'                treatment.var = 'treatment',
#'                pool = TRUE,
#'                replicate.var = "plot",
#'                time.var = "year")
#' 
#' # All pairwise replicates with treatment
#' df <- subset(pplots, year < 2004 & plot %in% c(21, 25, 32))
#' RAC_difference(df = df,
#'                species.var = "species",
#'                abundance.var = "relative_cover",
#'                replicate.var = "plot",
#'                time.var = "year",
#'                treatment.var = "treatment")
#' 
#' # All pairwise replicates without treatment
#' df <- subset(pplots, year < 2004 & plot %in% c(21, 25, 32))
#' RAC_difference(df = df,
#'                species.var = "species",
#'                abundance.var = "relative_cover",
#'                replicate.var = "plot",
#'                time.var = "year")
#' @export
RAC_difference <- function(df, time.var = NULL, species.var, 
                                abundance.var, replicate.var,
                                treatment.var = NULL, pool = FALSE, 
                                block.var = NULL) {
  
  # check no NAs in abundance column
  if(any(is.na(df[[abundance.var]]))) stop("Abundance column contains missing values")
  
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
    rankdf <- pool_replicates(df, time.var, species.var, abundance.var,
                              replicate.var, treatment.var)
  } else {
    # add zeros for species absent from a replicate within a treatment
    if (is.null(time.var)) {
      df <- fill_zeros(df, species.var, abundance.var) ## FIXME quietly keeps time if time.var = NULL, okay if unique by rep, but not otherwise
    } else {
      by <- c(time.var)
      df <- do.call(rbind, c(
        lapply(split(df, df[by], drop = TRUE),
              FUN = fill_zeros, species.var, abundance.var),
        list(make.row.names = FALSE)))
    }
    # rank species in each replicate
    by <- c(replicate.var, treatment.var, block.var)
    rankdf <- do.call(rbind, c(
      lapply(split(df, df[by], drop = TRUE),
             FUN = add_ranks, species.var, abundance.var),
      list(make.row.names = FALSE)))
  }
  
  # cross join for pairwise comparisons
  split_by <- c(species.var, block.var, time.var)
  merge_on <- !(names(rankdf) %in% split_by)
  cross.var2 <- paste(cross.var, 2, sep = '')
  rankdf <- lapply(split(rankdf, rankdf[split_by]),
                   function(x) {
                     y <- x[merge_on]
                     cross <- merge(x, y, by = NULL, suffixes = c('', '2'))
                     idx <- as.integer(cross[[cross.var]])
                     idx <- idx < as.integer(cross[[cross.var2]])
                     cross[idx,]
                   })
  ranktog <- do.call(rbind, c(rankdf, list(make.row.names = FALSE)))
  
  # split on treatment pairs (and block if not null)
  by <- c(block.var, time.var, cross.var, cross.var2)
  abundance.var2 <- paste(abundance.var, 2, sep = '')
  ranktog_split <- lapply(split(ranktog, ranktog[by], drop = TRUE),
                          FUN = SERSp,
                          species.var, abundance.var, abundance.var2)
  output <- do.call(rbind, c(ranktog_split, list(make.row.names = FALSE)))
  
  output_order <- c(
    time.var,
    block.var,
    replicate.var, paste(replicate.var, 2, sep = ''),
    treatment.var, paste(treatment.var, 2, sep = ''),
    'richness_diff', 'evenness_diff', 'rank_diff', 'species_diff')

  return(output[intersect(output_order, names(output))])

}

############################################################################
#
# Private functions: these are internal functions not intended for reuse.
# Future package releases may change these without notice. External callers
# should not use them.
#
############################################################################

# A function to calculate RAC difference between two samples 
# @param df a dataframe
# @param rank.var the name of the rank column at time 1
# @param rank.var2 the name of the rank column at time 2
# @param abundance.var the name of the abundance column at time 1
# @param abundance.var2 the name of the abundance column at time 2
SERSp <- function(df, species.var, abundance.var, abundance.var2) {
  
  out <- c(species.var, 'rank', 'rank2', abundance.var, abundance.var2)
  out <- unique(df[!(names(df) %in% out)])
  if (nrow(out) != 1)
    stop('Input df has not been correctly split.')
  
  df <- subset(df, df[[abundance.var]] != 0 | df[[abundance.var2]] != 0)
  
  #ricness and evenness differences
  s_t1 <- S(df[[abundance.var]])
  e_t1 <- EQ(as.numeric(df[[abundance.var]]))
  s_t2 <- S(df[[abundance.var2]])
  e_t2 <- EQ(as.numeric(df[[abundance.var2]]))
  
  sdiff <- abs(s_t1-s_t2)/nrow(df)
  ediff <- abs(e_t1-e_t2)/nrow(df)
  
  #Jaccard Index or Number of species not shared  
  spdiff <- df[df[[abundance.var]] == 0|df[[abundance.var2]] == 0,]
  spdiffc <- nrow(spdiff)/nrow(df)
  
  #Mean Rank Difference
  mrsc_diff <- mean(abs(df[['rank']]-df[['rank2']]) / nrow(df))
  
  metrics <- data.frame(richness_diff = sdiff, evenness_diff = ediff,
                        rank_diff = mrsc_diff, species_diff = spdiffc)
  
  return(cbind(out, metrics))
}
