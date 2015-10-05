#' Convert from a longform abundance dataframe to a time by species dataframe.
#'
#' @param df A dataframe containing time.var, species and abundance columns
#' @param time.var The name of the time column from df
#' @param species The name of the species column from df
#' @param abundance The name of the abundance column from df
#' @return comdat A dataframe of species abundances x time
#' @export
transpose_community <- function(df, time.var, species, abundance) {
    df<-as.data.frame(df)
    df[species]<-if(is.factor(df[[species]])==TRUE){factor(df[[species]])} else {df[species]}  
    df<-df[order(df[time.var], df[species]),]
    comdat<-tapply(df[[abundance]], list(df[[time.var]], as.vector(df[[species]])), sum)
    comdat[is.na(comdat)]<-0
    comdat<-as.data.frame(comdat)
    return(comdat)
}

#' check names of data frames
#'
#' @param given Vector of variable names as supplied by user
#' @param data Data frame containing variables
check_names <- function(given, data) {
    for (i in given){
        assertthat::assert_that(assertthat::has_name(data, i))
    }
}