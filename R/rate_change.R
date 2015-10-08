#'  Rate of community change over successive time intervals
#'
#' Calculates the slope of the differences in species composition within a community over increasing time 
#' intervals, which provides a  measures of the rate of directional change in community composition. 
#' Differences in species composition are characterized by Euclidean distances, 
#' which are calculated on pair-wise communities across the entire time series. 
#' For example, a data set with six time intervals will have distance values for five one-year time lags 
#' (year 1 vs year 2, year 2 vs year 3 ...), 
#' four two-year time lags (year 1 vs year 3, year 2 vs year 4 ...) and so forth. 
#' These distance values are regressed against the time lag interval. 
#' The slope of the regression line is reported as an indication of the rate and direction of compositional change in the community.
#' @param df A data frame containing time, species and abundance columns and an optional column of replicates
#' @param time.var The name of the time column 
#' @param species.var The name of the species column 
#' @param abundance.var The name of the abundance column 
#' @param replicate.var The name of the optional replicate column 
#' @return The \code{rate_change} function returns a numeric rate change value unless a replication column is specified in the input data frame. 
#' If replication is specified, the function returns a data frame with the following attributes:
#' \itemize{
#'  \item{rate_change: }{A numeric column with the synchrony values.}
#'  \item{replicate.var: }{A column that shares the same name and type as the replicate.var column in the input data frame.}
#' }
#' @details
#' The input data frame needs to contain columns for time, species and abundance; time.var, species.var and abundance.var are used to indicate which columns contain those variables.
#' If multiple replicates are included in the data frame, that column should be specified with replicate.var. Each replicate should reflect a single experimental unit - there must be a single abundance value per species within each time point and replicate.
#' 
#' The \code{rate_change} function uses linear regression to relate Euclidean distances to time lag intervals. 
#' It is recommended that fit of this relationship be verified using \code{rate_change_interval},
#' which returns the full set of community distance values and associated time lag intervals.
#' @references 
#' Collins, S. L., Micheli, F. and Hartt, L. 2000. A method to determine rates andpatterns of variability in ecological communities. - Oikos 91: 285-293. 
#' @examples 
#' data(knz_001d)
#' rate_change(knz_001d[knz_001d$subplot=="A_1",]) # for one subplot
#' rate_change(knz_001d, replicate.var = "subplot") # across all subplots
#' @export
rate_change <- function(df, time.var="year", species.var="species", abundance.var="abundance", replicate.var=NA) {
  check_numeric(df, time.var, abundance.var)
  if(is.na(replicate.var)) {
        check_single_onerep(df, time.var, species.var)
        output <- get_slope(df, time.var, species.var, abundance.var)
    } else {
        check_single(df, time.var, species.var, replicate.var)
        df[replicate.var] <- if(is.factor(df[[replicate.var]])) {
            factor(df[[replicate.var]])
        } else {
            df[replicate.var]
        }
        X <- split(df, df[replicate.var])
        out <- lapply(X, FUN=get_slope, time.var, species.var, abundance.var)
        reps <- unique(df[replicate.var])
        output <- cbind(reps, do.call("rbind", out))
        names(output)=c(replicate.var, "rate_change")
    }
    return(output)
}

#' 
#' Differences in community composition over successive time lag intervals
#' 
#' Calculates the differences in species composition within a community over increasing time 
#' intervals. Differences in species composition are characterized by Euclidean distances, 
#' which are calculated on pair-wise communities across the entire time series. 
#' For example, a data set with 6 time intervals will have distance values for five one-year time lags 
#' (year 1 vs year 2, year 2 vs year 3 ...), 
#' 4 two-year time lags (year 1 vs year 3, year 2 vs year 4 ...) and so forth. 
#' Returns the full set of community distance values and associated time lag intervals.
#' @param df A data frame containing time, species and abundance columns and an optional column of replicates
#' @param time.var The name of the time column 
#' @param species.var The name of the species column 
#' @param abundance.var The name of the abundance column 
#' @param replicate.var The name of the optional replicate column
#' @return The \code{rate_change_interval} function returns a data frame with the following attributes:
#' \itemize{
#'  \item{interval: }{A numeric column containing the interval length between time periods.}
#'  \item{distance: }{A numeric column containing the Euclidean distances.}
#'  \item{replicate.var: }{A column that shares the same name and type as the replicate.var column in the input data frame.}
#' }
#' The input data frame needs to contain columns for time, species and abundance; time.var, species.var and abundance.var are used to indicate which columns contain those variables.
#' If multiple replicates are included in the data frame, that column should be specified with replicate.var. Each replicate should reflect a single experimental unit - there must be a single abundance value per species within each time point and replicate.
#' @references 
#' Collins, S. L., Micheli, F. and Hartt, L. 2000. A method to determine rates andpatterns of variability in ecological communities. - Oikos 91: 285-293. 
#' @examples 
#' data(knz_001d)
#' rate_change_interval(knz_001d[knz_001d$subplot=="A_1",]) # for one subplot
#' rate_change_interval(knz_001d, replicate.var = "subplot") # across all subplots
#' @export
rate_change_interval <- function(df, time.var="year", species.var="species", abundance.var="abundance", replicate.var=NA) {
    stopifnot(is.numeric(df[[time.var]]))
    stopifnot(is.numeric(df[[abundance.var]]))
    if(is.na(replicate.var)) {
        check_single_onerep(df, time.var, species.var)
        output <- get_lagged_distances(df, time.var, species.var, abundance.var)
    } else {
        check_single(df, time.var, species.var, replicate.var)
        df[replicate.var] <- if(is.factor(df[[replicate.var]])) {
            factor(df[[replicate.var]])
        } else {
            df[replicate.var]
        }
        X <- split(df, df[replicate.var])
        out <- lapply(X, FUN=get_lagged_distances, time.var, species.var, abundance.var)
        #reps <- unique(df[replicate.var])
        #output <- cbind(reps, do.call("rbind", out))
        ID <- unique(names(out))
        out_rep <- mapply(function(x, y) "[<-"(x, replicate.var, value = y) ,
                      out, ID, SIMPLIFY = FALSE)
        output <- do.call("rbind", out_rep)
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

#' Get lagged distances for a single replicate
#' @description Returns a data frame with two columns, interval and distance. The interval is
#' the number of time steps between two communities, while distance is the 
#' euclidean distance of community change within one replicate lagged across invtervals.
#' @param df data frame to compute the slope of community change for
#' @param time.var The name of the time column from df
#' @param species.var The name of the species column from df
#' @param abundance.var The name of the abundance column from df
#' @return a data frame containing of time lags by species distances
#' @import stats
get_lagged_distances <- function(df, time.var="year", species.var="species", abundance.var="abundance") {
    df <- transpose_community(df, time.var, species.var, abundance.var)
    DM <- dist(df, method="euclidean", diag = FALSE, upper = FALSE)
    DM <- as.matrix(DM)
    
    rownums = row(DM)
    colnums = col(DM)
    lag_list = lapply(1:(nrow(df)-1), get_lag_i, DM, rownums, colnums)
    results <- data.frame(do.call(rbind, lag_list))
    names(results)=c("interval", "distance")
    return(results)
}

#' Get slope
#' Returns the slope of community change within one replicate.
#' @param df data frame to compute the slope of community change for
#' @param time.var The name of the time column from df
#' @param species.var The name of the species column from df
#' @param abundance.var The name of the abundance column from df
#' @return a slope of time lags by species distances
#' @import stats
get_slope <- function(df, time.var="year", species.var="species", abundance.var="abundance") {
    results <- get_lagged_distances(df, time.var, species.var, abundance.var)
    lm_coefficents <- lm(distance ~ interval, data=results)
    slope <- data.frame(lm_coefficents[1][[1]])
    return(slope[2,])
}

#' Get lagged values from a distance matrix
#' Get lagged values from distance matrix at value i
#' @param i the index of the matrix to lag
#' @param DM the distance matrix from which lagged values are drawn
#' @param rownums number of rows in the distance matrix
#' @param colnums number of columns in the distance matrix
#' @return the lagged values
get_lag_i <- function(i, DM, rownums, colnums) {
    cbind(lag = i, value = DM[rownums == (colnums + i)])
}
