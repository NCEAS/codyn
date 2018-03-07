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

  # specify aggregate formula from arguments

  rep_trt<-unique(subset(df, select = c(replicate.var, treatment.var)))
  
  if (pool) {
    #add zero abundnaces for missing species to get averages
    if (!is.null(time.var)) {
      df <- df[order(df[[time.var]]), ]
    }
    by <- c(time.var)
    allsp <- split_apply_combine(df, by, FUN = fill_zeros, species.var, abundance.var)

    # specify aggregate formula from arguments
    by <- c(species.var, treatment.var, time.var)
    spave <- aggregate.data.frame(allsp[abundance.var], allsp[by], FUN = mean)
    spave <- spave[spave[[abundance.var]] != 0, ]
  
    # rank each species by treatment and (time)
    by <- c(time.var)
    relrankdf1 <- split_apply_combine(spave, by, FUN = relrank,
      species.var, abundance.var, treatment.var)
  } else {
    #for block samples
    by <- c(time.var, block.var)
    relrankdf1 <- split_apply_combine(df, by, FUN = relrank,
      species.var, abundance.var, replicate.var)
  }

  # split and run curve_diff
  if (is.null(block.var) & is.null(time.var)) {
    if (pool){
      output <- curve_diff(relrankdf1, treatment.var, relrank, cumabund)
    } else {
      output <- curve_diff(relrankdf1, replicate.var, relrank, cumabund)
    }
  } else {
    splitvars <- c(time.var, block.var)
    relrankdf1_split <- split(relrankdf1,
                              relrankdf1[splitvars],
                              sep = "##", drop = TRUE)
    if (!is.null(block.var) | pool) {
      #feed treatment into curve_diff
      out <- lapply(relrankdf1_split,
                    FUN = curve_diff, treatment.var, relrank, cumabund)
    } else {
      #feed replicate into curve_diff
      out <- lapply(relrankdf1_split,
                    FUN = curve_diff, replicate.var, relrank, cumabund)
    }
    
    unsplit <- lapply(out, nrow)
    unsplit <- rep(names(unsplit), unsplit)
    output <- do.call(rbind, c(out, list(make.row.names = FALSE)))
    output[splitvars] <- do.call(rbind, strsplit(unsplit, '##'))
  }

  if (is.null(block.var) & !pool & !is.null(treatment.var)) {
    # add treatment for reference
    output <- merge(output,
                    merge(rep_trt, rep_trt, by = NULL, suffixes = c('', '2')))
  }
  
  output_order <- c(
    time.var,
    block.var,
    replicate.var, paste(replicate.var, 2, sep = ''),
    treatment.var, paste(treatment.var, 2, sep = ''),
    'curve_diff')
  return(output[intersect(output_order, names(output))])

}

############################################################################
#
# Private functions: these are internal functions not intended for reuse.
# Future package releases may change these without notice. External callers
# should not use them.
#
############################################################################

# A function to rank species in a sample by replicate
# @param df a dataframe
# @param species.var the name of the species column
# @param abundance.var the name of the abundance column
# @param replicate.var the name of the replicate column
# NOTE: when ranks are assigned by treatment and not replicate, treatment is fed into replicate.var
relrank <- function(df, species.var, abundance.var, replicate.var) {

  relrank <- subset(df, df[[abundance.var]]!=0)
  relrank$rank <- ave(relrank[[abundance.var]], relrank[[replicate.var]],
                      FUN = function(x) rank(-x, ties.method = "average"))
  relrank$maxrank <- ave(relrank$rank, relrank[[replicate.var]],
                         FUN = function(x) max(x))
  relrank$relrank  <- relrank$rank/relrank$maxrank
  relrank <- relrank[order(relrank[[replicate.var]], -relrank[[abundance.var]]), ]
  relrank$cumabund <- ave(relrank[[abundance.var]], relrank[[replicate.var]],
                          FUN = function(x) cumsum(x))

  return(relrank)
}

# A function calculate the curve difference between two treatments
# @param df a dataframe
# @param treatment.var the name of the treatment column
# @param relrank the name of the relative rank of each species in the sample
# @param cumabund the name of the cumulative abundance of each species in the sample
# NOTE: when doing on replicates not treatments, replicate is fed in as the treatment variable
curve_diff <- function(df, treatment.var, relrank, cumabund) {
   
  # cross join for pairwise comparisons
  treatment.var2 <- paste(treatment.var, 2, sep = '')
  compare <- unique(df[treatment.var])
  compare <- merge(compare, compare, by = NULL, suffixes = c('', '2'))
  idx <- as.integer(compare[[treatment.var]])
  idx <- idx < as.integer(compare[[treatment.var2]])
  compare <- compare[idx,]
  trt1 <- compare[[treatment.var]]
  trt2 <- compare[[treatment.var2]]

  cc_out<-data.frame() 
   
  for(i in 1:(length(trt1))) {
    subset_t1 <- df[df[[treatment.var]] == trt1[i],]
    subset_t1 <- subset_t1[order(subset_t1$cumabund),]
    subset_t2 <- df[df[[treatment.var]] == trt2[i],]
    subset_t2 <- subset_t2[order(subset_t2$cumabund),]
     
    sf1 <- stepfun(subset_t1$relrank, c(0, subset_t1$cumabund))
    sf2 <- stepfun(subset_t2$relrank, c(0, subset_t2$cumabund))
    r <- sort(unique(c(0, subset_t1$relrank, subset_t2$relrank)))
    h <- abs(sf1(r) - sf2(r))
    w <- c(diff(r), 0)
    CC=sum(w*h)
     
    output <- data.frame(trt1 = trt1[i], trt2 = trt2[i], curve_diff=CC)
    colnames(output)[1] <- treatment.var
    colnames(output)[2] <- paste(treatment.var, 2, sep = "")
     
    cc_out <- rbind(cc_out, output)
  }
  
  return(cc_out)
} 
