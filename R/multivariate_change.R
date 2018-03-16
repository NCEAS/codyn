#'@title Multivariate measures of community change between and within time periods
#' @description  Calculates the changes in composition and dispersion based off a Bray-Curtis dissimilarity matrix. Composition change is measured one of two ways. If "centroid_distance" is specified, community composition changes is measured as the the euclidean distance between the centroids of time period 1 and time period 2. If "ave_BC_dissim" is specified compositon change is measured as the average Bray-Curtis dissimilaritiy of all pairwise comparisions of replicates between two consecutive time periods. Both composition change metrics ranges from 0-1, where 0 are identical communities, and 1 and completelty different communities. Since, centroid distance is based on plotted distance between centroids, it is context dependent and depends on how many centroids are being plotted. The centroid distance between consecutive time periods depend on how mnay time periods are being measured. Average Bray-Curtis Dissimilarity is not context dependent. Dispersion change is the difference of average dispersion of each replicate to its centroid between consecutive time periods.
#' @param df A data frame containing time, species, abundance and replicate columns and an optional column of treatment
#' @param time.var The name of the time column 
#' @param species.var The name of the species column 
#' @param abundance.var The name of the abundance column 
#' @param replicate.var The name of the replicate column 
#' @param treatment.var the name of the optional treatment column
#' @param metric the composition change metric to return:
#'  \itemize{
#'  \item{"ave_BC_dissim": }{The default metric, calculates average of pairwise Bray-Curtis dissimilarity bewteen all replicates in consecutive time periods.}
#'  \item{"centroid_distance": }{Calculates the Euclidena distance between centroids of two consecutive time periods.}
#'  }
#' @return The multivariate_change function returns a data frame with the following attributes:
#' \itemize{
#'  \item{time.var: }{A characteric column that has the first of two time periods that are being compared.}
#'  \item{time.var2: }{A characteric column that has the second of two time periods that are being compared.}
#'  \item{centroid_distance_change: }{A numeric column that is the euclidean distance between the centroids of two consecutive time points, if "centroid_distance" is specified as the composition metric. Centroid distance is based on plotted distance between centroids, it is context dependent and depends on how many centroids are being plotted.}
#'  \item{BC_dissim_change: }{A numeric column that is the average pairwise Bray-Curtis dissimilarity of replicates in two consecutive time periods, if "ave_BC_dissim" is specified as the composition change metric.}
#'  \item{dispersion_change: }{A numeric column that is the difference in the average dispersion of the replicates around the centriod for the two consecutive time periods. A negative value indicates replicates are converging over time (there is less dispersion at time peroid 2 than time period 1) and a postive value indicates replicates are diverging over time (there is more dispersion at time period 2 than time period 1).}
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
#'                     abundance.var = "relative_cover",
#'                     metric = "centroid_distance") # centroid distance as composition metric
#' 
#' #Without treatment
#' df <- subset(pplots, treatment == "N1P0")
#' multivariate_change(df, 
#'                     time.var="year", 
#'                     replicate.var = "plot", 
#'                     species.var = "species", 
#'                     abundance.var = "relative_cover",
#'                     metric = "ave_BC_dissim") # Average Bray-Curtis dissimilarity as composisition metric
#' @importFrom vegan vegdist
#' @importFrom stats aggregate as.formula
#' @references Avolio et al. 2015; Avolio et al. OUR PAPER, Mari Anderson?
#' @export
multivariate_change <- function(df, time.var, 
                                 species.var, 
                                 abundance.var, 
                                 replicate.var, 
                                 treatment.var = NULL,
                                 metric = c("centroid_distance", "ave_BC_dissim")){
  
  # verify metric choice
  metric <- match.arg(metric)
  
  comp_metric <- get(metric)
  
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
  
  df1<-df[order(df[[time.var]], df[[replicate.var]]),]
  df2<-subset(df1, select = c(time.var, species.var, abundance.var, replicate.var))
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
    colnames(comp)[1] <- time.var
    colnames(comp)[2] <- paste(time.var, 2, sep="")
  }
  
  ##doing dispersion
  #collecting and labeling distances to centroid from betadisper to get a measure of dispersion and then take the mean for a year
  disp2 <- data.frame(time=species2[[time.var]],
                      dist = disp$distances)
  colnames(disp2)[1] <-time.var
  
  myformula <- as.formula(paste("dist", "~", time.var))
  disp2.2<-aggregate(myformula, mean, data=disp2)
  
  ##merge together  
  bc_dis1 <- merge(comp, disp2.2, by = time.var)
  bc_dis <- merge(bc_dis1, disp2.2, by.x = paste(time.var, 2, sep = ""), by.y = time.var)
  
  #calculate absolute difference
  bc_dis$disp_change <- bc_dis$dist.y - bc_dis$dist.x
  
  bc_dis$dist.y <- NULL
  bc_dis$dist.x <- NULL
  
  return(bc_dis)
}
