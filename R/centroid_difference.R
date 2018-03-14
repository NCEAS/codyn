#'@title Using centriods to calcualte differences in composition and dispersion
#'@description Calculates the difference in composition and dispersion between treatments based off a Bray-Curtis dissimilarity matrix at a single point in time. Composition difference is the euclidean distance between the centroids of different treatments. Dispersion difference is the difference of average dispersion of each replicate to its centroid between two treatments.
#' @param df A data frame containing an optional time column, species, abundance and replicate, and treatment columns
#' @param time.var The name of the optional time column 
#' @param species.var The name of the species column 
#' @param abundance.var The name of the abundance column 
#' @param replicate.var The name of the replicate column 
#' @param treatment.var the name of the treatment column
#' @return The multivariate_difference function returns a data frame with the following attributes:
#' \itemize{
#'  \item{treatment.var: }{A column that has same name and type as the treatment.var column, if treatment.var is specified.}
#'  \item{treatment.var2: }{A column that has the same type as the treatment.var column, and is named treatment.var with a 2 appended to it.}
#'  \item{composition_diff: }{A numeric column that is the euclidean distance between the centroids of two treatments at a single point in time.}
#'  \item{abs_dispersion_diff: }{A numeric column that is the absolute value of the difference in the average dispersion of the replicates around the centriod for the two treatments.}
#'  \item{trt_greater_disp: }{A column that has same type as the treatment.var column, and specifies which of the two  treatments has greater dispersion.}
#'  \item{time.var: }{A characteric column that has the same name and type as the time.var column, if specified.}
#' }
#' @references Our Avolio et al. paper, Avolio et al. 2015, Marti Anderson?
#' @importFrom vegan vegdist betadisper
#' @importFrom stats aggregate as.formula
#' @examples  
#' data(pplots)
#' #Without time
#' df <- subset(pplots, year == 2002)
#' centroid_difference(df, 
#'                     replicate.var = "plot", 
#'                     treatment.var = "treatment", 
#'                     species.var = "species", 
#'                     abundance.var = "relative_cover")
#' #With time
#' centroid_difference(pplots, 
#'                     time.var = "year", 
#'                     replicate.var = "plot", 
#'                     species.var = "species", 
#'                     abundance.var = "relative_cover", 
#'                     treatment.var = "treatment")
#' @export
centroid_difference <- function(df, time.var = NULL, species.var,
                                    abundance.var, replicate.var, treatment.var) {
  
  # check no NAs in abundance column
  if(any(is.na(df[[abundance.var]]))) stop("Abundance column contains missing values")
  
  if(is.null(time.var)){
    
    # check there unique species x time combinations
    check_single_onerep(df, replicate.var, species.var)
    
    output <- mult_diff(df, species.var, abundance.var, replicate.var, treatment.var)
    
  } else {
    
    # check unique species x time x replicate combinations
    check_single(df, time.var, species.var, replicate.var)
    
    splitvars <- time.var
    X <- split(df, 
               df[splitvars])
    out <- lapply(X, FUN = mult_diff, species.var, abundance.var, replicate.var, treatment.var)
    unsplit <- lapply(out, nrow)
    unsplit <- rep(names(unsplit), unsplit)
    output <- do.call(rbind, c(out, list(make.row.names = FALSE)))
    output[splitvars] <- do.call(rbind, as.list(unsplit))
  }
  
  output_order <- c(
    time.var,
    treatment.var, paste(treatment.var, 2, sep = ''),
    'composition_diff', 'abs_dispersion_diff', 'trt_greater_disp')
  
  return(output[intersect(output_order, names(output))])
}

############################################################################
#
# Private functions: these are internal functions not intended for reuse.
# Future package releases may change these without notice. External callers
# should not use them.
#
############################################################################

# A function calculate the community compositon difference (the distance between the centriods of two treatments) and absolute dispersion difference (the absolute difference in the average dispersion of the replicates around the centriod for the two treatments).
# @param df a dataframe
# @param species.var the name of the species column
# @param replicate.var the name of the replicate column
# @param treatment.var the name of the treatment column

mult_diff <- function(df, species.var, abundance.var, replicate.var, treatment.var){

  #transpose data
  df2<-subset(df, select = c(species.var, abundance.var, replicate.var, treatment.var))
  df2$id <- paste(df2[[treatment.var]], df2[[replicate.var]], sep="##")
  species <- transpose_community(df2, 'id', species.var, abundance.var)
  species$id <- row.names(species)
  speciesid <- do.call(rbind.data.frame, strsplit(species$id, split="##"))
  colnames(speciesid)[1] <- treatment.var
  colnames(speciesid)[2] <- replicate.var
  species2 <- cbind(speciesid, species)
  species3 <- subset(species2, select = -id)
  
  #calculate bray-curtis dissimilarities
  bc <- vegdist(species3[,3:ncol(species3)], method="bray")
  
  #calculate distances of each plot to treatment centroid (i.e., dispersion)
  disp <- betadisper(bc, species3[[treatment.var]], type="centroid")
  
  #getting distances between treatments with euclidean distances
  cent_dist <- as.data.frame(as.matrix(vegdist(disp$centroids, method="euclidean")))
  
  #extracting all treatment differences
  cent_dist2 <- as.data.frame(cbind(rownames(cent_dist)[which(lower.tri(cent_dist, diag=T), arr.ind=T)[,1]],
        colnames(cent_dist)[which(lower.tri(cent_dist, diag=T), arr.ind=T)[,2]],
        cent_dist[lower.tri(cent_dist, diag=T)]))
  cent_dist3 <- cent_dist2[cent_dist2$V1 != cent_dist2$V2,]
  cent_dist3[3]<-as.numeric(as.character(cent_dist3[[3]]))
  
  colnames(cent_dist3)[1] <- paste(treatment.var, 2, sep="")
  colnames(cent_dist3)[2] <- treatment.var
  colnames(cent_dist3)[3] <- "composition_diff"
  
  #collecting and labeling distances to centroid from betadisper to get a measure of dispersion and then take the mean for a treatment
  disp2 <- data.frame(treatment=species3[[treatment.var]],
                      dist = disp$distances)
  
  myformula <- as.formula(paste("dist", "~", treatment.var))
  disp2.2 <- aggregate(myformula, mean, data=disp2)
  
  #mege into get dispersion for each treatment
  cent_dist_disp <- merge(cent_dist3, disp2.2, by = treatment.var)
  cent_dist_disp2 <- merge(cent_dist_disp, disp2.2, by.x = paste(treatment.var, 2, sep = ""), by.y = treatment.var)
  
  #calculate absolute difference
  cent_dist_disp2[[treatment.var]] <- as.character(cent_dist_disp2[[treatment.var]])
  cent_dist_disp2[[paste(treatment.var, 2, sep = "")]] <- as.character(cent_dist_disp2[[paste(treatment.var, 2, sep = "")]])
  
   cent_dist_disp2$abs_dispersion_diff <- abs(cent_dist_disp2$dist.x - cent_dist_disp2$dist.y)
  cent_dist_disp2$trt_greater_disp <- as.character(ifelse(cent_dist_disp2$dist.x > cent_dist_disp2$dist.y, cent_dist_disp2[[treatment.var]], cent_dist_disp2[[paste(treatment.var, 2, sep = "")]]))
  
  cent_dist_disp2$dist.x <- NULL
  cent_dist_disp2$dist.y <- NULL
  

  return(cent_dist_disp2)
}