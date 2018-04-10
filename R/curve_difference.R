#' @title Curve Differences
#' @description Calculates the area difference between two rank abundance
#'   curves. There are three ways differences can be calculated. 1) Between all
#'   treatments within a block (note: block.var and treatment.var need to be
#'   specified. 2) Between treatments, pooling all replicates into a single
#'   species pool (note: pool = TRUE, treatment.var needs to be specified, and
#'   block.var = NULL. 3) All pairwise combinations between all replicates
#'   (note:block.var = NULL, pool = FALSE and specifying treatment.var is
#'   optional. If treatment.var is specified, the treatment that each replicate
#'   belongs to will also be listed in the output).
#' @param df A data frame containing a species, abundance, and replicate columns
#'   and optional time, treatment, and block columns
#' @param time.var The name of the optional time column 
#' @param species.var The name of the species column 
#' @param abundance.var The name of the abundance column 
#' @param replicate.var The name of the replicate column 
#' @param treatment.var The name of the optional treatment column
#' @param pool An argument to allow values to be pooled within treatment. The
#'   default value is "FALSE", a value of "TRUE" takes the average abundance of
#'   all species within a treatment at a given time point prior to comparisons.
#' @param block.var The name of the optional block column
#'  
#' @return The curve_difference function returns a data frame with the following
#'   attributes:
#' \itemize{
#'  \item{curve_diff: }{A numeric column of the area difference in curves
#'  between the two samples being compared (replicates or treatments).}
#'  \item{replicate.var: }{A column that has same name and type as the
#'  replicate.var column, represents the first replicate being compared. Note, a
#'  replicate column will be returned only when pool is FALSE or block.var =
#'  NULL.}
#'  \item{replicate.var2: }{A column that has the same type as the replicate.var
#'  column, and is named replicate.var with a 2 appended to it, represents the
#'  second replicate being compared. Note, a replicate.var column will be
#'  returned only when pool is FALSE and block.var = NULL.}
#'  \item{time.var: }{A column that has the same name and type as the time.var
#'  column, if time.var is specified.}
#'  \item{treatment.var: }{A column that has same name and type as the
#'  treatment.var column, represents the first treatment being compared. A
#'  treatment.var column will be returned when pool is TRUE or block.var is
#'  specified, or treatment.var is specified.}
#'  \item{treatment.var2: }{A column that has the same type as the treatment.var
#'  column, and is named treatment.var with a 2 appended to it, represents the
#'  second treatment being compared. A treatment.var column will be returned
#'  when pool is TRUE or block.var is specified, or treatment.var is specified.}
#'  \item{block.var: }{A column that has same name and type as the block.var
#'  column, if block.var is specified.}
#' }
#' @references Avolio et al. OUR PAPER.
#' @examples
#' data(pplots)
#' # With block and no time
#' df <- subset(pplots, year == 2002 & block < 3)
#' curve_difference(df = df,
#'                  species.var = "species",
#'                  abundance.var = "relative_cover",
#'                  treatment.var = 'treatment',
#'                  block.var = "block",
#'                  replicate.var = "plot")
#' # With blocks and time
#' df <- subset(pplots, year < 2004 & block < 3)
#' curve_difference(df = df,
#'                  species.var = "species",
#'                  abundance.var = "relative_cover",
#'                  treatment.var = 'treatment',
#'                  block.var = "block",
#'                  replicate.var = "plot",
#'                  time.var = "year")
#' 
#' # Pooling by treatment with time
#' df <- subset(pplots, year < 2004)
#' curve_difference(df = df,
#'                  species.var = "species",
#'                  abundance.var = "relative_cover",
#'                  treatment.var = 'treatment',
#'                  pool = TRUE,
#'                  replicate.var = "plot",
#'                  time.var = "year")
#' 
#' # All pairwise replicates with treatment
#' df <- subset(pplots, year < 2004 & plot %in% c(21, 25, 32))
#' curve_difference(df = df,
#'                  species.var = "species",
#'                  abundance.var = "relative_cover",
#'                  replicate.var = "plot",
#'                  time.var = "year",
#'                  treatment.var = "treatment")
#' 
#' # All pairwise replicates without treatment
#' df <- subset(pplots, year < 2004 & plot %in% c(21, 25, 32))
#' curve_difference(df = df,
#'                  species.var = "species",
#'                  abundance.var = "relative_cover",
#'                  replicate.var = "plot",
#'                  time.var = "year")
#' @importFrom stats aggregate.data.frame
#' @export
curve_difference <- function(df, time.var = NULL, species.var, 
                                abundance.var, replicate.var,
                                treatment.var = NULL, pool = FALSE, 
                                block.var = NULL) {

  df<-as.data.frame(df)
  
  # check no NAs in abundance column
  if(any(is.na(df[[abundance.var]]))) stop("Abundance column contains missing values")
  
  # check no NAs in species column
  if(any(is.na(df[[species.var]]))) stop("Species names are missing")

  #check no species are repeated
  if (is.null(time.var)){
    # check there unique species x time combinations
    check_single_onerep(df, replicate.var, species.var)
  } else {
    # check unique species x time x replicate combinations
    check_single(df, time.var, species.var, replicate.var)
  }

  if (!is.null(block.var)) {
    reps_exp <- length(unique(df[[block.var]])) * length(unique(df[[treatment.var]]))
    reps_obs <- length(unique(df[[replicate.var]]))
    if (reps_exp != reps_obs)
      stop("There is not one replicate per treatment in a block")
  }

  rep_trt<-unique(subset(df, select = c(replicate.var, treatment.var)))
  
  if (pool) {
    #add zero abundnaces for missing species to get averages
    if (!is.null(time.var)) {
      df <- df[order(df[[time.var]]), ]
    }
    by <- c(time.var)
    allsp <- split_apply_combine(df, by, FUN = fill_zeros,
      species.var, abundance.var)

    # specify aggregate formula from arguments
    by <- c(species.var, treatment.var, time.var)
    spave <- aggregate.data.frame(allsp[abundance.var], allsp[by], FUN = mean)
    spave <- spave[spave[[abundance.var]] != 0, ]
  
    # rank each species by treatment and optionally time
    by <- c(treatment.var, time.var)
    rankabunddf <- split_apply_combine(spave, by, FUN = add_rank_abundance,
      species.var, abundance.var)
  } else {
    #for block samples
    by <- c(block.var, replicate.var, time.var)
    rankabunddf <- split_apply_combine(df, by, FUN = add_rank_abundance,
        species.var, abundance.var)
  }
  
  # determine which variable indicates contrast for comparison
  if (pool) {
    cross.var <- treatment.var
  } else if (is.null(time.var) & !is.null(block.var)) {
    cross.var <- treatment.var
  } else {
    cross.var <- replicate.var
  }

  # order cross.var if unordered factor
  to_ordered = is.factor(df[[cross.var]]) & !is.ordered(df[[cross.var]])
  if (to_ordered) {
    class(rankabunddf[[cross.var]]) <- c('ordered', class(df[[cross.var]]))
  }
  
  # cross join for pairwise comparisons
  split_by <- c(block.var, time.var)
  merge_to <- !(names(rankabunddf) %in% split_by)
  cross.var2 <- paste(cross.var, 2, sep = '')
  output <- split_apply_combine(rankabunddf, split_by, FUN = function(x) {
    y <- x[merge_to]
    cross <- merge(x, y, by = NULL, suffixes = c('', '2'))
    idx <- cross[[cross.var]] < cross[[cross.var2]]
    cross[idx, ]
  })
  
  # unorder cross.var if orginally unordered factor
  if (to_ordered) {
    class(output[[cross.var]]) <- class(df[[cross.var]])
  }
  
  # split on treatment pairs (and block if not null)
  output[['curve_diff']] <- mapply(curve_dissim,
    output[['rankabund']], output[['rankabund2']])
  
  output_order <- c(
    time.var,
    block.var,
    replicate.var, paste(replicate.var, 2, sep = ''),
    treatment.var, paste(treatment.var, 2, sep = ''),
    'curve_diff')
  return(output[intersect(output_order, names(output))])

}