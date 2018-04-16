#'@title Using dissimilarity-based metrics to calucuate changes in composition and dispersion
#' @description Calculates the changes in composition and dispersion based off a Bray-Curtis dissimilarity matrix. Composition change is the euclidean distance between the centroids of time period 1 and time period 2 and ranges from 0-1, where 0 are identical communities, and 1 and completelty different communities. Since, Composition change is based on plotted distance between centroids, it is context dependent and depends on how many centroids are being plotted, thus result of composition change depends on how many time periods are being measured. Dispersion change is the difference of average dispersion of each replicate to its centroid between consecutive time periods.
#' @param df A data frame containing time, species, abundance and replicate columns and an optional column of treatment
#' @param time.var The name of the time column 
#' @param species.var The name of the species column 
#' @param abundance.var The name of the abundance column 
#' @param replicate.var The name of the replicate column 
#' @param treatment.var the neame of the optional treatment column
#' @return The multivariate_change function returns a data frame with the following attributes:
#' \itemize{
#'  \item{time.var_pair: }{A characteric column that has the time points to be compared, separated by a dash.}
#'  \item{composition_change: }{A numeric column that is the euclidean distance between the centroids of two consecutive time points.}
#'  \item{dispersion_change: }{A numeric column that is the difference in the average dispersion of the replicates around the centriod for the two consecutive time periods. A negative value indicates replicates are converging over time (there is less dispersion at time peroid 2 than time period 1) and a postive value indicates replicates are diverging over time (there is more dispersion at time period 2 than time period 1.}
#'  \item{treatment.var: }{A column that has same name and type as the treatment.var column, if treatment.var is specified.}
#' }
#' @examples 
#' data(pplots)
#' #With treatment
#' multivariate_change(pplots,
#'                     time.var="year", 
#'                     replicate.var = "plot", 
#'                     treatment.var = "treatment", 
#'                     species.var = "species", 
#'                     abundance.var = "relative_cover")
#' #In each year there are 6 replicates and there are 4 years of data for 3 time comparitions, thus 24 total observations in each treatment.
#' 
#' #Without treatment
#' df <- subset(pplots, treatment == "N1P0")
#' multivariate_change(df, 
#'                     time.var="year", 
#'                     replicate.var = "plot", 
#'                     species.var = "species", 
#'                     abundance.var = "relative_cover")
#' #In each year there are 6 replicates and there are 4 years of data for 3 time comparitions, thus 24 total observations.
#'
#' @importFrom vegan vegdist betadisper
#' @importFrom stats aggregate as.formula
#' @references Avolio et al. 2015; Avolio et al. OUR PAPER, Mari Anderson?
#' @export
multivariate_change <- function(df, time.var, species.var, abundance.var, replicate.var, treatment.var = NULL){
  
  # drop extraneous columns
  args <- as.list(match.call())
  df <- as.data.frame(df[as.character(args[grep('\\.var$', names(args))])])
  
  # check no NAs in abundance column
  if(any(is.na(df[[abundance.var]]))) stop("Abundance column contains missing values")
  
  # check unique species x time x replicate combinations
  check_single(df, time.var, species.var, replicate.var)

  df <- as.data.frame(df)

  # calculate change for each treatment
  by <- treatment.var
  output <- split_apply_combine(df, by, FUN = mult_change, time.var,
    species.var, abundance.var, replicate.var, treatment.var)
  
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
# A function calculate the community compositon change (the distance between the
# centriods of two consecutive time periods) and dispersion change (the
# difference in the average dispersion of the replicates around the centriod for
# the two consecutive time periods). For dispersion change a negative value
# indicates replicates are converging over time and a postive value indicates
# replicates are diverging over time.
#
# @param df a dataframe
# @param time.var the name of the time column
# @param species.var the name of the species column
# @param replicate.var the name of the replicate column
mult_change <- function(df, time.var, species.var, abundance.var,
                        replicate.var, treatment.var) {
  #get years
  timestep <- sort(unique(df[[time.var]]))
  
  #transpose data
  idvar = c(treatment.var, replicate.var, time.var)
  species <- reshape(df, idvar = idvar, timevar = species.var, direction = 'wide')
  species[is.na(species)] <- 0

  #calculate bray-curtis dissimilarities
  abund <- species[, -(1:length(idvar))]
  bc <- vegdist(abund, method = "bray")
  
  #calculate distances of each plot to year centroid (i.e., dispersion)
  if (is.null(treatment.var)) {
    message('Composition and dispersion change calculation using ',
            nrow(species),' observations.')
  } else {
    message('Composition and dispersion change calculation using ',
            nrow(species), ' observations at treatment.var value ',
            df[[treatment.var]][[1]])
  }
  disp <- betadisper(bc, species[[time.var]], type = "centroid")
  
  #getting distances between centroids over years; these centroids are in BC
  #space, so that's why this uses euclidean distances
  cent_dist <- as.matrix(vegdist(disp$centroids, method = "euclidean"))
  
  ##extracting only the comparisions, year x to year x+1.
  time.var2 <- paste(time.var, 2, sep = '')
  cent_dist_yrs <- data.frame(
    time1 = timestep[1:length(timestep) - 1],
    time2 = timestep[2:length(timestep)],
    composition_change = diag(
      as.matrix(cent_dist[2:nrow(cent_dist), 1:(ncol(cent_dist) - 1)])))
  colnames(cent_dist_yrs)[1] <- time.var
  colnames(cent_dist_yrs)[2] <- time.var2
  
  #collecting and labeling distances to centroid from betadisper to get a
  #measure of dispersion and then take the mean for a year
  species['dist'] <- disp['distances']
  by <- c(treatment.var, time.var)
  disp2 <- aggregate.data.frame(species['dist'], species[by], mean)
  
  #merge together  
  bc_dis1 <- merge(cent_dist_yrs, disp2, by = time.var)
  bc_dis <- merge(bc_dis1, disp2[c(time.var, 'dist')],
    by.x = time.var2, by.y = time.var)
  
  #calculate absolute difference
  bc_dis$dispersion_change <- bc_dis$dist.y - bc_dis$dist.x
  
  return(bc_dis)
}
