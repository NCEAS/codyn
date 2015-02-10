#' A function to calculate community rate changes over time within multiple replicates
#'
#'#' This is an analysis of differences in species composition between samples at increasing time lags. It measures the rate of directional change in community composition. First, a triangular dissimilarity matrix is calculated using Euclidean distance. Then, the Euclidean distance values are plotted against time lags. For example, a data set with 6 time intervals will have 5 one-year time lags (year 1 vs year 2, year 2 vs year 3 ...) and 4 two-year time lags (year 1 vs year 3, year 2 vs year 4 ...). Finally, distance values are regressed against time lag. The slope of the regression line indicates the rate and direction of change.
#' @param data1 A dataframe containing replicate, year, species and abundance columns.
#' @param replicate The name of the replicate column from data1
#' @param year The name of the year column from data1
#' @param species The name of the species column from data1
#' @param abundance The name of the abundance column from data1
#' @import vegan
#' @return output The rate of community change
rate_change <- function(data1, replicate="replicate", year="year", species="species", abundance="abundance") {
	if(is.na(replicate)==TRUE){
    output<-get_slope(data1, year, species, abundance)}else{
        data1[replicate]<-if(is.factor(data1[[replicate]])==TRUE){factor(data1[[replicate]])} else {data1[replicate]}
		X <- split(data1, data1[replicate])
		out <- lapply(X, FUN=get_slope, year, species, abundance)
		reps <- unique(data1[replicate])
		output <- cbind(reps, do.call("rbind", out))
		names(output)=c(replicate, "rate_change")
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



#' A function that returns the slope of community change within one replicate.
#' 
#' @param year The name of the year column from data1
#' @param species The name of the species column from data1
#' @param abundance The name of the abundance column from data1
#' @import vegan
#' @return a slope of year lags by species distances
get_slope = function(data1, year="year", species="species", abundance="abundance") {
		data1 <- calComDat(data1, species, year, abundance)
		DM <- vegdist(data1[-1], method="euclidean", diag = FALSE, upper = FALSE, na.rm = TRUE)
    DM <- as.matrix(DM)
	
    rownums = row(DM)
    colnums = col(DM)

    get_lag_i = function(i){
        cbind(lag = i, value = DM[rownums == (colnums + i)])
    }

    lag_list = lapply(
        1:(nrow(data1)-1),
        get_lag_i)

    results <- data.frame(do.call(rbind, lag_list))
		lm_coefficents <- lm(value ~ lag, data=results)
		slope <- data.frame(lm_coefficents[1][[1]])
		return(slope[2,])
}

