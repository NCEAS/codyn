#' @title Variance Ratio
#' @description Computes the ratio of the variance of aggregate species abundances 
#' in a community to the sum of the variances of individual, component species. A 
#' variance ratio = 1 indicates that species do not covary,  a variance ratio > 1 
#' indicates predominately positive covariance among species and a variance 
#' ratio < 1 indicates predominately negative covariance (Schluter 1984).
#' 
#' Includes a null modeling option to test if the variance ratio significantly 
#' differs from 1. The null community is created by randomly selecting different 
#' starting points for each species' time series, which generates a community in 
#' which species abundances vary independently but within-species autocorrelation 
#' is maintained (Hallett et al. 2014). This randomization is repeated a user-specific 
#' number of times and confidence intervals are reported for the resultant null 
#' distribution of variance ratios. If the dataframe includes multiple replicates, 
#' the variance ratios for the actual and null communities are averaged within each 
#' iteration unless specified otherwise.
#' @param df A data frame containing time, species and abundance columns and an optional column of replicates
#' @param time.var The name of the time column 
#' @param species.var The name of the species column 
#' @param abundance.var The name of the abundance column 
#' @param bootnumber The number of null model iterations used to calculated confidence intervals
#' @param li The lower confidence interval, defaults to lowest 2.5\% CI
#' @param ui The upper confidence interval, defaults to upper 97.5\% CI  
#' @param replicate.var The name of the optional replicate column 
#' @param average.replicates If true returns the variance ratio and CIs averaged 
#' across replicates; if false returns the variance ratio and CI for each replicate
#' @return The variance_ratio function returns a data frame with the following attributes:
#' \itemize{
#'  \item{VR: }{A numeric column with the actual variance ratio value.}
#'  \item{lowerCI: }{A numeric column with the lowest confidence interval value.}
#'  \item{upperCI: }{A numeric column with the highest confidence interval value.}
#'  \item{nullmean: }{A numeric column with the average null variance ratio value.}
#'  \item{replicate.var: }{A column that has same name and type as the replicate.var column, if replication is specified.}
#' }
#' @details
#' The input data frame needs to contain columns for time, species and abundance; 
#' time.var, species.var and abundance.var are used to indicate which columns 
#' contain those variables. If multiple replicates are included in the data frame, 
#' that column should be specified with replicate.var. Each replicate should 
#' reflect a single experimental unit - there must be a single abundance value 
#' per species within each time point and replicate.
#' 
#' Null model confidence intervals default to the standard lowest 2.5\% and 
#' upper 97.5\% of the null distribution, typically these do not need to be change, 
#' but they can be user-modified to set more stringent CIs.
#'  @references
#'  Hallett, Lauren M., Joanna S. Hsu, Elsa E. Cleland, Scott L. Collins, 
#'  Timothy L. Dickson, Emily C. Farrer, Laureano A. Gherardi, et al. (2014) 
#'  "Biotic Mechanisms of Community Stability Shift along a Precipitation Gradient." 
#'  Ecology 95, no. 6: 1693-1700. doi: 10.1890/13-0895.1
#'  
#'  Schluter, Dolph. (1984) "A Variance Test for Detecting Species Associations, 
#'  with Some Example Applications." Ecology 65, no. 3: 998-1005. doi:10.2307/1938071.
#' @export
#' @examples 
#'  data(knz_001d)
#'  
#'  # Calculate the variance ratio and CIs averaged within replicates
#'  # Here the null model is repeated once, for final use it is recommended to set a 
#'  # large bootnumber (eg, 10000)
#'  
#'  res_averagedreplicates <- variance_ratio(knz_001d, time.var = "year", species.var = "species", 
#'  abundance.var = "abundance", bootnumber = 1, replicate = "subplot")
#'  
#'  #Calculate the variance ratio and CIs for each replicate
#'  
#'  res_withinreplicates <- variance_ratio(knz_001d, time.var = "year", species.var = "species", 
#'  abundance.var = "abundance", bootnumber = 1, replicate = "subplot", average.replicates = FALSE)
variance_ratio<-function(df, time.var="year", species.var="species",  abundance.var="abundance", bootnumber, replicate.var=NA,
                        li=0.025, ui=0.975,  average.replicates=TRUE) {
  # check to make sure abundance is numeric data
  check_numeric(df, time.var, abundance.var)
    if(is.na(replicate.var)) {
        check_single_onerep(df, time.var, species.var)
        VR<-variance_ratio_longformdata(df,time.var, species.var, abundance.var)
    } else {
        check_single(df, time.var, species.var, replicate.var)
        df[replicate.var]<-if(is.factor(df[[replicate.var]])==TRUE) {
            factor(df[[replicate.var]])
        } else {
            df[replicate.var]
        }
        if(average.replicates==TRUE) {
           check_multispp(df, species.var, replicate.var)
            df<-df[order(df[[replicate.var]]),]
            X<-split(df, df[replicate.var])
            VR<-mean(unlist(lapply(X, FUN=variance_ratio_longformdata, time.var, species.var, abundance.var))) 
        } else {
            check_multispp(df, species.var, replicate.var)
            df<-df[order(df[[replicate.var]]),]
            X <- split(df, df[replicate.var])
            VR<-do.call("rbind", lapply(X, FUN=variance_ratio_longformdata, time.var,  species.var,abundance.var))
        }
    }
    nullout<-temporal_torus_translation_CI(df, time.var, species.var, abundance.var, variance_ratio_matrixdata, bootnumber, replicate.var, li, ui, average.replicates)
    output<-(cbind(nullout, VR))
    row.names(output)<-NULL
    return(as.data.frame(output)) 
}

############################################################################
#
# Private functions: these are internal functions not intended for reuse.  
# Future package releases may change these without notice. External callers
# should not use them.
#
############################################################################

#' A function to calculate the variance ratio 
#'
#' @param comdat A community dataframe
#' @return var.ratio The variance ratio of the community
#' @import stats
variance_ratio_matrixdata<-function(comdat){
    check_sppvar(comdat)
    all.cov <- cov(comdat, use="pairwise.complete.obs")
    col.var<-apply(comdat, 2, var)
    com.var <-sum(all.cov)
    pop.var <-sum(col.var)
    var.ratio<-com.var/pop.var
    return(var.ratio)
}

#' A function to calculate the variance ratio from a longform dataframe
#'
#' @param df A dataframe containing time.var, replicate.var, species.var and abundance.var columns
#' @param time.var The name of the time.var column from df
#' @param species.var The name of the species.var column from df
#' @param abundance.var The name of the abundance.var column from df
#' @return var.ratio The variance ratio of the community
variance_ratio_longformdata<-function(df, time.var, species.var, abundance.var){
    com.use<-transpose_community(df, time.var, species.var, abundance.var)
    var.ratio<-variance_ratio_matrixdata(com.use)
    return(var.ratio)
}
