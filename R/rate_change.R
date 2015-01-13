#' A function to calculate community rate changes over time within multiple replicates
#'     
#' @param comm_data Community dataset. Must be in 'wide' format
#' @param year Year variable. Must be the first column of comm_data
#' @return linear model coefficients
#' @export
#' @import vegan
rate_change <- function(comm_data, year="year") {
    # creating distance matrix from species matrix
    DM <- vegdist(comm_data[-1], method="euclidean", diag = FALSE, upper = FALSE, na.rm = TRUE)
    DM <- as.matrix(DM)
    
    # final: matrix converted to long-form
    DM.long <- data.frame(get_lags(DM, comm_data, year))
    
    # regression
    lm_coefficents <- lm(value ~ lag, data=DM.long)
		return(summary(lm_coefficents))
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
#' This is an analysis of differences in species composition between samples at increasing time lags. It measures the rate of directional change in community composition. First, a triangular dissimilarity matrix is calculated using Euclidean distance. Then, the Euclidean distance values are plotted against time lags. For example, a data set with 6 time intervals will have 5 one-year time lags (year 1 vs year 2, year 2 vs year 3 ...) and 4 two-year time lags (year 1 vs year 3, year 2 vs year 4 ...). Finally, distance values are regressed against time lag. The slope of the regression line indicates the rate and direction of change, and the regression coefficient is a measure of signal verses noise.
#' @param DM distance matrix to be used for lag calculations
#' @param comm_data community data frame
#' @param year The year column in comm_data 
#' @return matrix of lag values
get_lags = function(DM, comm_data, year) {
    # label each row and each column.
    rownums = row(DM)
    colnums = col(DM)

    # A function to get all the elements that are lagged by i time steps,
    # then put them into a long-form matrix with i in the first column and
    # the matrix value in the second column.
    get_lag_i = function(i){
        cbind(lag = i, value = DM[rownums == (colnums + i)])
    }
    
    # apply get_lag_i to all lags from 1 to n-1
    # replace n with number of columns 
    lag_list = lapply(
        1:(length(unique(comm_data[year]))-1),
        get_lag_i)
    # squash all the lags from the list into one long-form matrix
    result <- do.call(rbind, lag_list)
    return(result)
}

#' @param data A year by species dataframe
#' @return a slope of year lags by species distances
get_slope = function(data1) {
	
		# creating distance matrix
		DM <- vegdist(data[-1], method="euclidean", diag = FALSE, upper = FALSE, na.rm = TRUE)
    DM <- as.matrix(DM)
	
		# label each row and each column.
    rownums = row(DM)
    colnums = col(DM)

		# A function to get all the elements that are lagged by i time steps,
    # then put them into a long-form matrix with i in the first column and
    # the matrix value in the second column.
    get_lag_i = function(i){
        cbind(lag = i, value = DM[rownums == (colnums + i)])
    }

		# apply get_lag_i to all lags from 1 to n-1
    # replace n with number of columns 	
    lag_list = lapply(
        1:(length(unique(comm_data$year))-1),
        get_lag_i)

		# squash all the lags from the list into one long-form matrix
    results <- data.frame(do.call(rbind, lag_list))
	
		# fitting a regression for lag 
		lm_coefficents <- lm(value ~ lag, data=results)
		slope <- data.frame(lm_coefficents[1][[1]])
		return(slope[2,])
}
