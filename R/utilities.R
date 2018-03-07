#' Convert from a longform abundance dataframe to a time by species dataframe.
#'
#' @param df A dataframe containing time.var, species.var and abundance.var columns
#' @param time.var The name of the time column from df
#' @param species.var The name of the species column from df
#' @param abundance.var The name of the abundance column from df
#' @return A dataframe of species abundances x time
transpose_community <- function(df, time.var, 
                                species.var, 
                                abundance.var) {
  df <- as.data.frame(df)
  
  # remove unused levels if species is a factor
  df[species.var] <- if(is.factor(df[[species.var]]) == TRUE){factor(df[[species.var]])} else {df[species.var]}
  
  # sort by time and species
  df <- df[order(df[[time.var]], df[[species.var]]),]
 
  # cast as a species x time dataframe; NAs to 0s
  comdat <- tapply(df[[abundance.var]], list(df[[time.var]], as.vector(df[[species.var]])), sum)
  comdat[is.na(comdat)] <- 0
  comdat <- as.data.frame(comdat)
  
  # results
  return(comdat)
}

#' check names of data frames
#'
#' @param given Vector of variable names as supplied by user
#' @param data Data frame containing variables
check_names <- function(given, data) {
    for (i in given) {
        assertthat::assert_that(assertthat::has_name(data, i))
    }
}

#' Utility function to warn users that either multiple records exist within replicates, or that data may be spanning mutiple replicates but no replicate.var has been specified
#' @param df A dataframe containing time.var, species.var and abundance.var columns
#' @param time.var The name of the time column from df
#' @param species.var The name of the species column from df
check_single_onerep <- function(df, time.var, species.var){
  if(max(table(df[[time.var]], df[[species.var]]))>1) stop("Either data span multiple replicates with no replicate.var specified or multiple records within years for some species") }

#' Utility function to ensure only a single record exists for a given species within one replicate, for one time point.
#' @param df A dataframe containing time.var, species.var, and replicate.var columns
#' @param time.var The name of the time column from df
#' @param species.var The name of the species column from df
#' @param replicate.var The name of the replicate column from df

check_single <- function(df, time.var, species.var, replicate.var){
  X <- split(df, df[[replicate.var]])
  checksingle <- lapply(X, FUN = function(xx) apply(table(xx[[species.var]], xx[[time.var]]), 2, function(x) any(x>1)))
  reptest <- unlist(lapply(checksingle, any))
  yrtest <- lapply(checksingle, which)

  if(any(unlist(checksingle))){
    if(length(names(reptest)[which(reptest)]) == 1){

    stop(paste("In replicate", names(reptest)[which(reptest)], "there is more than one record for species at the time point", unlist(lapply(yrtest, names))))
    }
      else  {
        toprint <- unlist(lapply(yrtest, names))
    stop("For the following replicates in the following time points, there are more than one records for species: \n", paste(names(toprint), collapse = "\t"), "\n", paste(toprint, collapse = "\t"))
      }
  }
}

#' Utility to check for numeric abundance and time variables
#'
#' @param df A dataframe containing time.var, species.var, and replicate.var columns
#' @param time.var The name of the time column from df
#' @param abundance.var The name of the replicate column from df

check_numeric <- function(df, time.var, abundance.var) {
  if(!is.numeric(df[[time.var]])) { stop("Time variable is not numeric") }
  if(!is.numeric(df[[abundance.var]])) { stop("Abundance variable is not numeric") }
  }

#' Utility function to stop calculations if only one species occurs in at least one replicate
#' @param df A dataframe containing time.var, species.var and abundance.var columns
#' @param species.var The name of the species column from df
#' @param replicate.var The name of the replicate column from df
check_multispp <- function(df, species.var, replicate.var){

  df <- droplevels(df)

  spptable <- table(df[[species.var]], df[[replicate.var]])
  if (min(colSums(spptable > 1)) < 2) {
    stop("One or more replicates consists of only a single species;
       please remove these replicates prior to calculations ")
  }
}

#' Utility function to stop calculations if the species never change in a replicate
#' @param comdat A community dataframe
#' @importFrom stats var
check_sppvar <- function(comdat){
  sppvar <- sum(apply(comdat, 2, var))
  if(sppvar == 0)
    stop("One or more replicates consist of species that never vary;
         please remove these replicates before calculation")}

#' Utility function to calculate richness
#' @param x Vector of abundance of each species
S <- function(x){
  x1 <- x[x!=0]
  stopifnot(x1==as.numeric(x1))
  length(x1)
}

#' Utility function to calculate EQ evenness from Smith and Wilson 1996
#' @param x Vector of abundance of each species
#' If all abundances are equal it returns a 1
#' @importFrom stats lm
EQ <- function(x){
  x1 <- x[x != 0]
  if (length(x1) == 1) {
    return(NA)
  }
  if (abs(max(x1) - min(x1)) < .Machine$double.eps^0.5) {
    return(1)
  }
  r <- rank(-x1, ties.method = "average")
  r_scale <- r/max(r)
  x_log <- log(x1)
  fit <- lm(r_scale~x_log)
  b <- fit$coefficients[[2]]
  -2/pi*atan(b)
}

#' Add zero abundances for missing species, on the assumption that any species
#' in the \code{species.var} column should be included for every group defined
#' by all the remaining colums save \code{abundance.var}.
#'
#' @param df A dataframe with species, abundances, and at least one other column
#'   to group by
#' @param species.var The name of the species column from df
#' @param abundance.var The name of the abundance column from df
#' @return A dataframe with the same columns as df, but with zeros added for
#'   species that are present in df, but not in a particular group.
fill_zeros <- function(df, species.var, abundance.var) {

  if(any(is.na(df[[species.var]]))) stop("Species names are missing")
  
  mergevars <- !(names(df) %in% c(species.var, abundance.var))
  full <- merge(unique(df[mergevars]), unique(df[species.var]))
  df <- merge(df, full, all = TRUE)
  df[is.na(df[abundance.var]), abundance.var] <- 0

  return(df)
}

#' @title Add abundance ranks
#' @description Rank species by abundance, by specified groupig. Species with
#'   zero abundance receive rank S+1, where S is the total number of species in
#'   the group.
#' @param df A data frame containing a single record per species with its abundance
#' @param abundance.var The name of the abundance column
#' 
#' @return The add_ranks function returns a data frame with the following
#'   additional column:
#'   \itemize{
#'     \item{rank: }{A numeric column with the species rank; a rank of 1
#'     indicates the species was most abundant in that time period. All species
#'     that are not present in that time period have the rank value S+1 where S
#'     is the number of species in the sample.
#'     }
#'   }
add_ranks <- function(df, species.var, abundance.var) {
  
  if (nrow(df) != nrow(unique(df[names(df) != abundance.var])))
    stop('Input df has not been correctly split.')
  
  # get species richness, note that zero abunance does not contribute to S
  S <- S(df[[abundance.var]])
  # rank from high to low abundance
  df[['rank']] <- rank(-df[[abundance.var]], ties.method = 'average')
  # adjust ranks for species with zero abundance to S + 1
  df[df[['rank']] > S, 'rank'] <- S + 1
  
  return(df)
}

#' @title Faster split-apply-combine for data frames, when the FUN does not change
#' the number, order, data type or levels of columns in \code{df}.
#' 
#' @param df A data frame
#' @param by The column(s) of the data frame that determine splits
#' @param FUN The function applied to each data frame
#' @param ... Additional parameters to FUN
#' 
#' @source \url{https://stackoverflow.com/a/9730292/687112}
split_apply_combine <- function(df, by, FUN, ...) {
  if (is.null(by)) {
    # just apply
    df <- FUN(df, ...)
  } else {
    # split (names get in the way)
    dfs <- split(df, df[by], drop = TRUE)
    dfs <- unname(dfs)
    # apply
    dfs <- lapply(dfs, FUN = FUN, ...)
    # combine (flatten outer list, then flatten again across vectors from same column)
    dfs <- unlist(dfs, recursive = FALSE, use.names = TRUE)
    hdr <- unique(names(dfs))
    idx <- seq_along(hdr)
    df <- lapply(idx, function(i) unlist(dfs[i == idx], FALSE, FALSE))
    names(df) <- hdr
  }
  as.data.frame(df)
}
