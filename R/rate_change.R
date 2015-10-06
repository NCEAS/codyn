#' A function to calculate community rate changes over time within multiple replicates
#'
#'#' This is an analysis of differences in species composition between samples at increasing time lags. It measures the rate of directional change in community composition. First, a triangular dissimilarity matrix is calculated using Euclidean distance. Then, the Euclidean distance values are plotted against time lags. For example, a data set with 6 time intervals will have 5 one-year time lags (year 1 vs year 2, year 2 vs year 3 ...) and 4 two-year time lags (year 1 vs year 3, year 2 vs year 4 ...). Finally, distance values are regressed against time lag. The slope of the regression line indicates the rate and direction of change.
#' @param df A dataframe containing replicate, time, species and abundance columns.
#' @param replicate.var The name of the replicate column from df. Defaults to NA.
#' @param time.var The name of the time column from df
#' @param species.var The name of the species column from df
#' @param abundance.var The name of the abundance column from df
#' @return output The rate of community change
#' @export
rate_change <- function(df, time.var="time", species.var="species", abundance.var="abundance", replicate.var=NA) {
    if(is.na(replicate.var)) {
        output<-get_slope(df, time.var, species.var, abundance.var)
    } else {
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


############################################################################
#
# Private functions: these are internal functions not intended for reuse.
# Future package releases may change these without notice. External callers
# should not use them.
#
############################################################################

#' Get slope
#' Returns the slope of community change within one replicate.
#' @param df data frame to compute the slope of community change for
#' @param time.var The name of the time column from df
#' @param species.var The name of the species column from df
#' @param abundance.var The name of the abundance column from df
#' @return a slope of time lags by species distances
get_slope = function(df, time.var="time", species.var="species", abundance.var="abundance") {
    df <- transpose_community(df, time.var, species.var, abundance.var)
    DM <- dist(df[-1], method="euclidean", diag = FALSE, upper = FALSE)
    DM <- as.matrix(DM)
    
    rownums = row(DM)
    colnums = col(DM)
    lag_list = lapply(1:(nrow(df)-1), get_lag_i, DM, rownums, colnums)
    
    results <- data.frame(do.call(rbind, lag_list))
    lm_coefficents <- lm(value ~ lag, data=results)
    slope <- data.frame(lm_coefficents[1][[1]])
    return(slope[2,])
}

#' Get lagged values from a distance matrix
#' Get lagged values from distance matrix at value i
#' @i the index of the matrix to lag
#' @DM the distance matrix from which lagged values are drawn
#' @rownums number of rows in the distance matrix
#' @colnums number of columns in the distance matrix
#' @return the lagged values
get_lag_i = function(i, DM, rownums, colnums) {
    cbind(lag = i, value = DM[rownums == (colnums + i)])
}
