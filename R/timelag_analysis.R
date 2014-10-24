require(vegan)

#' Function for time-lag analysis or directional change in temporal Community composition data 
#'
#'
#' This is a function that calculates the slope of directional change 
#' @param comm_data Community dataset. Must be in 'wide' format
#' @return linear model coefficients

# creating distance matrix from species matrix
# DM <- vegdist(comm_data, method="euclidean", diag = FALSE, upper = FALSE, na.rm = TRUE)
# DM <- as.matrix(DM)

# get the lags from distance matrix
get_lags = function(DM){
    # label each row and each column.
    rownums = row(DM)
    colnums = col(DM)
    # A mini-function to get all the elements that are lagged by i time steps,
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
    do.call(rbind, lag_list)
}

# final: matrix converted to long-form
# DM.long <- data.frame(get_lags(DM))
# regression
# lm(value ~ lag, data=DM.long)



