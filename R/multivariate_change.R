#'@title Bray-Curtis dissimilarity of replicates between and within time periods
#' @description Calculates the average changes in Bray-Curtis dissimilarity of replicates between two consecutive time periods and within time periods. Between changes the average Bray-Curtis dissimilaritiy of all pairwise comparisions of replicates between two consecutive time periods. This is a measure of how dissimilar the community composition of two time periods is. Bray-Curtis dissimilarity ranges from 0-1, where 0 are identical communities, and 1 and completelty different communiites. Within change is derived in two steps. First, the average Bray-Curtis dissimilarity of all replicates within a time period is calculated. This is a measure of how homogenous a community is within a time step is. Then, these averages are compared between two consecutive time periods.
#' @param df A data frame containing time, species, abundance and replicate columns and an optional column of treatment
#' @param time.var The name of the time column 
#' @param species.var The name of the species column 
#' @param abundance.var The name of the abundance column 
#' @param replicate.var The name of the replicate column 
#' @param treatment.var the neame of the optional treatment column
#' @return The multivariate_change function returns a data frame with the following attributes:
#' \itemize{
#'  \item{time.var: }{A characteric column that has the first of two time periods that are being compared.}
#'  \item{time.var2: }{A characteric column that has the second of two time periods that are being compared.}
#'  \item{BC_between_change: }{A numeric column that is the average pairwise Bray-Curtis dissimilarity of replicates in two consecutive time periods. 0 - The communities are similar over time. 1 - The communities are changing over time.}
#'  \item{BC_within_change: }{A numeric column that is change between two time periods in the average pairwise Bray-Curtis dissimilarity of replicates within a time period. A positive number indicates that time.var2 has greater varaiblity in community composition than time.var}
#'  \item{treatment.var: }{A column that has same name and type as the treatment.var column, if treatment.var is specified.}
#' }
#' @examples 
#' data(pplots)
#' #With treatment
#' dissimilarity_change(pplots,
#'                     time.var="year", 
#'                     replicate.var = "plot", 
#'                     treatment.var = "treatment", 
#'                     species.var = "species", 
#'                     abundance.var = "relative_cover")
#' 
#' #Without treatment
#' df <- subset(pplots, treatment == "N1P0")
#' dissimilarity_change(df, 
#'                     time.var="year", 
#'                     replicate.var = "plot", 
#'                     species.var = "species", 
#'                     abundance.var = "relative_cover")
#' @importFrom vegan vegdist
#' @importFrom stats aggregate as.formula
#' @references Avolio et al. 2015; Avolio et al. OUR PAPER, Mari Anderson?
#' @export
dissimilarity_change <- function(df, time.var, 
                                 species.var, 
                                 abundance.var, 
                                 replicate.var, 
                                 treatment.var = NULL,
                                 comp_metric = c("centroid_distance", "ave_BC_dissim")){
  
  # verify metric choice
  metric <- match.arg(metric)
  
  composition_metric <- get(comp_metric)
  
  # check no NAs in abundance column
  if(any(is.na(df[[abundance.var]]))) stop("Abundance column contains missing values")
  
  # check unique species x time x replicate combinations
  check_single(df, time.var, species.var, replicate.var)
  
  df <- as.data.frame(df)

  if (is.null(treatment.var)) {
    output <- dissim_change(df, time.var, species.var, abundance.var, replicate.var)
  } else {
    # calculate change for each treatment
    splitvars <- treatment.var
    X <- split(df, df[splitvars])
    out <- lapply(X, FUN = mult_change, time.var, species.var, abundance.var, replicate.var, comp_metric = composition_metric)
    unsplit <- lapply(out, nrow)
    unsplit <- rep(names(unsplit), unsplit)
    output <- do.call(rbind, c(out, list(make.row.names = FALSE)))
    output[splitvars] <- do.call(rbind, as.list(unsplit))
  }
    
  output_order <- c(
    time.var,
    paste(time.var,"2", sep=""),
    treatment.var,
    'BC_between_change', 'BC_within_change')

  return(output[intersect(output_order, names(output))])
}

############################################################################
#
# Private functions: these are internal functions not intended for reuse.
# Future package releases may change these without notice. External callers
# should not use them.
#
############################################################################
# A function calculate the Bray-Curtis dissimilarity change between consequtive time periods and dispersion change (the difference in the average dispersion of the replicates around the centriod for the two consecutive time periods). For dispersion change a negative value indicates replicates are converging over time and a postive value indicates replicates are diverging over time.
# @param df a dataframe
# @param time.var the name of the time column
# @param species.var the name of the species column
# @param replicate.var the name of the replicate column
mult_change <- function(df, time.var, species.var, abundance.var, replicate.var, comp_metric) {
  
  df2<-subset(df, select = c(time.var, species.var, abundance.var, replicate.var))
  df2$id <- paste(df2[[time.var]], df2[[replicate.var]], sep="##")
  species <- codyn:::transpose_community(df2, 'id', species.var, abundance.var)
  species$id <- row.names(species)
  speciesid <- do.call(rbind.data.frame, strsplit(species$id, split="##"))
  colnames(speciesid)[1] <- "time_forfixxyz"
  colnames(speciesid)[2] <- replicate.var
  speciesid[[time.var]]<-as.numeric(as.character(speciesid$time_forfixxyz))
  speciesid.2 <- subset(speciesid, select = -time_forfixxyz)
  species2 <- cbind(speciesid.2, species)
  species3 <- subset(species2, select = -id)
  
  #calculate bray-curtis dissmilarity
  bc <- vegdist(species3[,3:ncol(species3)], method="bray")
  
  #calculate distances of each plot to year centroid (i.e., dispersion)
  disp <- betadisper(bc, species[[time.var]], type="centroid")
  
  
  #### doing compositional differences
  if(comp_metric == "ave_BC_dissim"){
    #extracting lower diagonal
    bc2 <- as.data.frame(cbind(rownames(bc)[which(lower.tri(bc), arr.ind=T)[,1]],
                               colnames(bc)[which(lower.tri(bc), arr.ind=T)[,2]],
                               bc[lower.tri(bc1)]))
    c1 <- as.data.frame(do.call('rbind', strsplit(as.character(bc2$V1), "##", fixed = T)))
    c2 <- as.data.frame(do.call('rbind', strsplit(as.character(bc2$V2), "##", fixed = T)))
    
    bc3 <- cbind(bc2, c1, c2)
    bc3$bc_dissim <- as.numeric(as.character(bc3$V3))
    colnames(bc3)[4] <- paste(time.var, 2, sep="")
    colnames(bc3)[6] <- time.var
    bc3$compare <- ifelse(bc3[[time.var]] == bc3[[paste(time.var, 2, sep="")]], 1, 2)
    
    #between time differences
    bc_between <- subset(bc3, compare == 2)
    myformula2 <- as.formula(paste("bc_dissim", "~", time.var, "+", paste(time.var, 2, sep = "")))
    bc_between_ave <- aggregate(myformula2, mean, data=bc_between)
    colnames(bc_between_ave)[3] <- "BC_dissim_change"
    
    #select only consecutive years
    bc_between_ave$yr1 <- as.integer(as.factor(bc_between_ave[[time.var]]))
    bc_between_ave$yr2 <- as.integer(as.factor(bc_between_ave[[paste(time.var, 2, sep = "")]]))
    bc_between_ave$diff <- bc_between_ave$yr2 - bc_between_ave$yr1
    bc_between_ave2 <- subset(bc_between_ave, diff==1)
    bc_between_ave2$yr1 <- NULL
    bc_between_ave2$yr2 <- NULL
    bc_between_ave2$diff <- NULL
    
    comp <- bc_between_ave2
    
  } else {
    #getting distances between centroids over years; these centroids are in BC space, so that's why this uses euclidean distances
    cent_dist <- as.matrix(vegdist(disp$centroids, method="euclidean"))
    
    ##extracting only the comparisions, year x to year x+1.
    cent_dist_yrs <- data.frame(
      time1 = timestep[1:length(timestep)-1],
      time2 = timestep[2:length(timestep)],
      centroid_distance_change = diag(
        as.matrix(cent_dist[2:nrow(cent_dist), 1:(ncol(cent_dist)-1)])))
    
    comp <- cent_dist_yrs
  }
  
  ##doing dispersion
  #collecting and labeling distances to centroid from betadisper to get a measure of dispersion and then take the mean for a year
  disp2 <- data.frame(time=species2[[time.var]],
                      dist = disp$distances)
  
  myformula <- as.formula(paste("dist", "~", "time"))
  disp2.2<-aggregate(myformula, mean, data=disp2)
  
  ##subtract consecutive years (x+1 - x). A positive value indicates greater dispersion in year x+1 and a negative value indicates less dispersion in year x+1
  disp_yrs <- data.frame(time2 = timestep[2:length(timestep)],
                         dispersion_change = diff(disp2.2$dist))
  
  ###merging all together
  distances <- merge(comp, disp_yrs, by="time2")
  
  return(bc_dis)
}
