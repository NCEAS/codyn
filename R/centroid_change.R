#'@title Using centoids to calucuate changes in composition and dispersion
#' @description Calculates the changes in composition and dispersion based off a Bray-Curtis dissimilarity matrix. Composition change is the euclidean distance between the centroids of time period 1 and time period 2. Dispersion change is the difference of average dispersion of each replicate to its centroid between consecutive time periods.
#' @param df A data frame containing time, species, abundance and replicate columns and an optional column of treatment
#' @param time.var The name of the time column 
#' @param species.var The name of the species column 
#' @param abundance.var The name of the abundance column 
#' @param replicate.var The name of the replicate column 
#' @param treatment.var the neame of the optional treatment column
#' @return The multivariate_change function returns a data frame with the following attributes:
#' \itemize{
#'  \item{time.var_pair: }{A characteric column that has the time points to be compared, separated by a dash.}
#'  \item{centroid_distance_change: }{A numeric column that is the euclidean distance between the centroids of two consecutive time points.}
#'  \item{dispersion_change: }{A numeric column that is the difference in the average dispersion of the replicates around the centriod for the two consecutive time periods. A negative value indicates replicates are converging over time (there is less dispersion at time peroid 2 than time period 1) and a postive value indicates replicates are diverging over time (there is more dispersion at time period 2 than time period 1.}
#'  \item{treatment.var: }{A column that has same name and type as the treatment.var column, if treatment.var is specified.}
#' }
#' @examples 
#' data(pplots)
#' #With treatment
#' centroid_change(pplots,
#'                 time.var="year", 
#'                 replicate.var = "plot", 
#'                 treatment.var = "treatment", 
#'                 species.var = "species", 
#'                 abundance.var = "relative_cover")
#' 
#' #Without treatment
#' df <- subset(pplots, treatment == "N1P0")
#' centroid_change(df, 
#'                 time.var="year", 
#'                 replicate.var = "plot", 
#'                 species.var = "species", 
#'                 abundance.var = "relative_cover")

#' @importFrom vegan vegdist betadisper
#' @importFrom stats aggregate as.formula
#' @references Avolio et al. 2015; Avolio et al. OUR PAPER, Mari Anderson?
#' @export
centroid_change <- function(df, time.var, species.var, abundance.var, replicate.var, treatment.var = NULL){
  
  # check no NAs in abundance column
  if(any(is.na(df[[abundance.var]]))) stop("Abundance column contains missing values")
  
  # check unique species x time x replicate combinations
  check_single(df, time.var, species.var, replicate.var)
  
  df <- as.data.frame(df)


  if (is.null(treatment.var)) {
    output <- mult_change1(df, time.var, species.var, abundance.var, replicate.var)
  } else {
    df[[treatment.var]]<-as.character(df[[treatment.var]])
    # calculate change for each treatment
    splitvars <- treatment.var
    X <- split(df, df[splitvars])
    out <- lapply(X, FUN = mult_change1, time.var, species.var, abundance.var, replicate.var)
    unsplit <- lapply(out, nrow)
    unsplit <- rep(names(unsplit), unsplit)
    output <- do.call(rbind, c(out, list(make.row.names = FALSE)))
    output[splitvars] <- do.call(rbind, as.list(unsplit))
  }
    
  output_order <- c(
    time.var, paste(time.var,"2", sep=""),
    treatment.var,
    'composition_change', 'dispersion_change')

  return(output[intersect(output_order, names(output))])
}

############################################################################
#
# Private functions: these are internal functions not intended for reuse.
# Future package releases may change these without notice. External callers
# should not use them.
#
############################################################################
# A function calculate the community compositon change (the distance between the centriods of two consecutive time periods) and dispersion change (the difference in the average dispersion of the replicates around the centriod for the two consecutive time periods). For dispersion change a negative value indicates replicates are converging over time and a postive value indicates replicates are diverging over time.
# @param df a dataframe
# @param time.var the name of the time column
# @param species.var the name of the species column
# @param replicate.var the name of the replicate column
mult_change1 <- function(df, time.var, species.var, abundance.var, replicate.var) {
  #get years
  timestep <- sort(unique(df[[time.var]]))
  
  #transpose data
  df2<-subset(df, select = c(time.var, species.var, abundance.var, replicate.var))
  df2$id <- paste(df2[[time.var]], df2[[replicate.var]], sep="##")
  species <- transpose_community(df2, 'id', species.var, abundance.var)
  species$id <- row.names(species)
  speciesid <- do.call(rbind.data.frame, strsplit(species$id, split="##"))
  colnames(speciesid)[1] <- "time_forfixxyz"
  colnames(speciesid)[2] <- replicate.var
  speciesid[[time.var]]<-as.numeric(as.character(speciesid$time_forfixxyz))
  speciesid.2 <- subset(speciesid, select = -time_forfixxyz)
  species2 <- cbind(speciesid.2, species)
  species3 <- subset(species2, select = -id)
  
  #calculate bray-curtis dissimilarities
  bc <- vegdist(species3[,3:ncol(species3)], method="bray")
  
  #calculate distances of each plot to year centroid (i.e., dispersion)
  disp <- betadisper(bc, species3[[time.var]], type="centroid")
  
  #getting distances between centroids over years; these centroids are in BC space, so that's why this uses euclidean distances
  cent_dist <- as.matrix(vegdist(disp$centroids, method="euclidean"))
  
  ##extracting only the comparisions, year x to year x+1.
  cent_dist_yrs <- data.frame(
    time1 = timestep[1:length(timestep)-1],
    time2 = timestep[2:length(timestep)],
    composition_change = diag(
      as.matrix(cent_dist[2:nrow(cent_dist), 1:(ncol(cent_dist)-1)])))
  colnames(cent_dist_yrs)[1] <- time.var
  colnames(cent_dist_yrs)[2] <- paste(time.var, 2, sep="")
  
  #collecting and labeling distances to centroid from betadisper to get a measure of dispersion and then take the mean for a year
  disp2 <- data.frame(time=species2[[time.var]],
                   dist = disp$distances)
  colnames(disp2)[1] <-time.var
  
  myformula <- as.formula(paste("dist", "~", time.var))
  disp2.2<-aggregate(myformula, mean, data=disp2)
  
  ##merge together  
  bc_dis1 <- merge(cent_dist_yrs, disp2.2, by = time.var)
  bc_dis <- merge(bc_dis1, disp2.2, by.x = paste(time.var, 2, sep = ""), by.y = time.var)
  
  #calculate absolute difference
  bc_dis$dispersion_change <- bc_dis$dist.y - bc_dis$dist.x
  
  bc_dis$dist.y <- NULL
  bc_dis$dist.x <- NULL
  
  return(bc_dis)
}
