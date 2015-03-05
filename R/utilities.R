#' Convert from a longform abundance dataframe to a year by species dataframe.
#'
#' @param data1 A dataframe containing year, species and abundance columns
#' @param year The name of the year column from data1
#' @param species The name of the species column from data1
#' @param abundance The name of the abundance column from data1
#' @return comdat A dataframe of species abundances x year
calComDat<-function(data1, species, year, abundance) {
    data1<-data1[order(data1[year], data1[species]),]
    comdat<-tapply(data1[[abundance]], list(data1[[year]], as.vector(data1[[species]])), sum)
    comdat[is.na(comdat)]<-0
    comdat<-as.data.frame(comdat)
    return(comdat)
}

#' check names of data frames
#'
#' @param given Vector of variable names as supplied by user
#' @param data Data frame containing variables
check_names<-function(given, data) {
    for (i in given){
        assertthat::assert_that(assertthat::has_name(data, i))
    }
}