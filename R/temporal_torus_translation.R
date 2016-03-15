#' @title Cyclic shift
#' @description Performs a user-specified function on a null ecological community created via cyclic shifts (Harms et al. 2001, Hallett et al. 2014).
#' The null community is formed by randomly selected different starting years for the time series of each species.
#' This generates a null community matrix in which species abundances vary independently but within-species autocorrelation is maintained.
#' The user-specified function must require a species x time matrix input.
#' @param df A data frame containing time, species and abundance columns and an optional column of replicates
#' @param time.var The name of the time column
#' @param species.var The name of the species column
#' @param abundance.var The name of the abundance column
#' @param FUN A function to calculate on the null community
#' @return The cyclic_shift function returns the same output as the user-specified function, as calculated on a null community.
#' @details The input data frame needs to contain columns for time, species and abundance; time.var, species.var and abundance.var are used to indicate which columns contain those variables.
#' @examples
#' # Calculate a covariance matrix on a null community
#' data(knz_001d)
#' a1_cyclic <- cyclic_shift(subset(knz_001d, subplot=="A_1"), time.var="year",
#' species.var="species", abundance.var="abundance", FUN=cov, bootnumber = 10)
#' @references
#' Hallett, Lauren M., Joanna S. Hsu, Elsa E. Cleland, Scott L. Collins, Timothy L. Dickson, Emily C. Farrer, Laureano A. Gherardi, et al. "Biotic Mechanisms of Community Stability Shift along a Precipitation Gradient." Ecology 95, no. 6 (2014): 1693-1700.
#'
#' Harms, Kyle E., Richard Condit, Stephen P. Hubbell, and Robin B. Foster. "Habitat Associations of Trees and Shrubs in a 50-Ha Neotropical Forest Plot." Journal of Ecology 89, no. 6 (2001): 947-59.
#' @export
cyclic_shift <- function(df, time.var="year",
                         species.var="species",
                         abundance.var="abundance",
                         replicate.var,
                         FUN,
                         bootnumber,
                         average.replicates=TRUE){

  assertthat::assert_that(is.numeric(df[[abundance.var]]))

  ## if you give a replicate, it must be a factor. This gives users responsibility for order.
  if (!missing(replicate.var)) {
    assertthat::assert_that(is.factor(df[[replicate.var]]))
  }

  if (missing(replicate.var)) {
    check_single_onerep(df, time.var = time.var, species.var = species.var)

    out <- replicate(bootnumber,
                   FUN(
                     shuffle_community(
                       transpose_community(df,
                                           time.var,
                                           species.var,
                                           abundance.var)
                     )))

  shift <- structure(list(out = out), class = "cyclic_shift")

  return(shift)
}


#' @title Confidence Intervals Using a Cyclic Shift Significance Testing
#' @description Calculates confidence intervals for a user-specified test statistic that derives from a species x time matrix.
#' It does so by employing a cyclic shift that creates a null community by randomly selecting different starting points for each species' time series. This generates a community in which species abundances vary independently but within-species autocorrelation is maintained (Harms et al. 2001, Hallett et al. 2014).
#' This randomization is repeated a user-specific number of times and confidence intervals are reported for the resultant null distribution of the test statistic.
#' If the data frame includes multiple replicates, the test statistics for the null communities are averaged within each iteration unless specified otherwise.
#' @param df A data frame containing time, species and abundance columns and an optional column of replicates
#' @param time.var The name of the time column
#' @param species.var The name of the species column
#' @param abundance.var The name of the abundance column
#' @param FUN A function to calculate on the null community that requires a species x time matrix
#' @param bootnumber The number of null model iterations used to calculated confidence intervals
#' @param replicate.var The name of the replication column from df
#' @param li The lower confidence interval, defaults to lowest 2.5\% CI
#' @param ui The upper confidence interval, defaults to upper 97.5\% CI
#' @param average.replicates If true returns the CIs averaged across replicates; if false returns the CI for each replicate
#' @return The confint.cyclic_shift function returns a dataframe with the following attributes:
#' \itemize{
#'  \item{lowerCI: }{A numeric column with the lowest confidence interval value.}
#'  \item{upperCI: }{A numeric column with the highest confidence interval value.}
#'  \item{nullMean: }{A numeric column with the average value of the specified test statistic when calculated on a null community.}
#'  \item{replicate.var: }{A column that has same name and type as the replicate.var column, if replication is specified.}
#' }
#' @details
#' This function was developed to be applied to the variance ratio; to use it for that purpose see the variance_ratio function.
#' It is included as a general function here but has very specific requirements. Namely, it is only relevant for functions that return single values, and for which varying species abundances independently is an acceptable way to develop a null test statistic value.
#'
#' If applying this function to test statistics other the variance ratio, the input data frame needs to contain columns for time, species and abundance; time.var, species.var and abundance.var are used to indicate which columns contain those variables.
#' If multiple replicates are included in the data frame, that column should be specified with replicate.var. Each replicate should reflect a single experimental unit - there must be a single abundance value per species within each time point and replicate.
#' Null model confidence intervals default to the standard lowest 2.5\% and upper 97.5\% of the null distribution, typically these do not need to be change, but they can be user-modified to set more stringent CIs.
#' @references
#' Hallett, Lauren M., Joanna S. Hsu, Elsa E. Cleland, Scott L. Collins, Timothy L. Dickson, Emily C. Farrer, Laureano A. Gherardi, et al. "Biotic Mechanisms of Community Stability Shift along a Precipitation Gradient." Ecology 95, no. 6 (2014): 1693-1700.
#'
#' Harms, Kyle E., Richard Condit, Stephen P. Hubbell, and Robin B. Foster. "Habitat Associations of Trees and Shrubs in a 50-Ha Neotropical Forest Plot." Journal of Ecology 89, no. 6 (2001): 947-59.
#' @import stats
#' @export
confint.cyclic_shift <- function(df, time.var="year", species.var="species",
                                 abundance.var="abundance", FUN, bootnumber,
                                 li=0.025, ui=0.975, replicate.var=NA, average.replicates=TRUE){
  if(!is.numeric(df[[abundance.var]])) { stop("Abundance variable is not numeric") }



  if(is.na(replicate.var)){

    lowerCI <- quantile(out, li)
    upperCI <- quantile(out, ui)
    nullmean <- mean(out)
    output <- cbind(lowerCI, upperCI, nullmean)
    row.names(output) <- NULL

  } else {


    check_single(df, time.var, species.var, replicate.var)
    df <- df[order(df[[replicate.var]]),]
    X <- split(df, df[replicate.var])
    lout <- lapply(X, cyclic_shift, time.var, species.var, abundance.var, FUN, bootnumber)

    ## simple workaround, a model for how this function will need to work
    lout <- lapply(lout, `[[`, i = "out")

    if (average.replicates  ==  TRUE) {
      out <- do.call("rbind", lout)
      out <- colMeans(out)
      lowerCI <- quantile(out, li)
      upperCI <- quantile(out, ui)
      nullmean <- mean(out)
      output <- as.data.frame(cbind(lowerCI, upperCI, nullmean))
      row.names(output) <- NULL
    } else {
      lowerCI <- do.call("rbind", lapply(lout, quantile, li))
      upperCI <- do.call("rbind", lapply(lout, quantile, ui))
      nullmean <- do.call("rbind", lapply(lout, mean))
      reps <- unique(df[replicate.var])
      output <- cbind(reps, lowerCI, upperCI, nullmean)
      names(output)[2:3] <- c("lowerCI", "upperCI")
    }
  }
  return(output)
}



############################################################################
#
# Private functions: these are internal functions not intended for reuse.
# Future package releases may change these without notice. External callers
# should not use them.
#
############################################################################



#' A function to generate a community dataframe with a random start time for each species
#'
#' @param comdat A community dataframe
#' @return rand.comdat A randomized community dataframe
shuffle_community <- function(comdat){
  nr <- nrow(comdat)
  nc <- ncol(comdat)
  rand.comdat <- matrix(NA, nrow = nr, ncol = nc)
  rand.start <- sample.int(nr, nc, replace = TRUE)
  for (i in seq_len(nc)) {
    rand.comdat[, i] <- permute::shuffleSeries(comdat[,i])
  }
  rand.comdat <- as.data.frame(rand.comdat)
  names(rand.comdat) <- names(comdat)
  row.names(rand.comdat) <- row.names(comdat)
  return(rand.comdat)
}
