#'@title Multivariate measures of community differences between pairs of treatments at a single time point.
#' @description  Calculates the difference in composition and dispersion based off a Bray-Curtis dissimilarity matrix. Composition difference is measured one of two ways. If "centroid_distance" is specified, community composition difference is measured as the the euclidean distance between the centroids of treatments. If "ave_BC_dissim" is specified compositon difference is measured as the average Bray-Curtis dissimilaritiy of all pairwise comparisions of replicates or treatments. Both composition difference metrics range from 0-1, where 0 are identical communities, and 1 and completelty different communities. Since, centroid distance is based on plotted distance between centroids, it is context dependent and depends on how many centroids are being plotted. The centroid distance between treatments depends on how many treatments are being compared. Average Bray-Curtis dissimilarity is not context dependent. Dispersion difference is the difference of average dispersion of each replicate to its centroid between two treatments.
#' @param df A data frame containing time, species, abundance and replicate columns and an optional column of treatment
#' @param time.var The name of the time column 
#' @param species.var The name of the species column 
#' @param abundance.var The name of the abundance column 
#' @param replicate.var The name of the replicate column 
#' @param treatment.var the name of the optional treatment column
#' @param metric the composition change metric to return:
#'  \itemize{
#'  \item{"ave_BC_dissim": }{The default metric, calculates average of pairwise Bray-Curtis dissimilarity bewteen all replicates in two treatments.}
#'  \item{"centroid_distance": }{Calculates the Euclidean distance between centroids of two treatments.}
#'  }
#' @return The multivariate_change function returns a data frame with the following attributes:
#' \itemize{
##'  \item{treatment.var: }{A column that has same name and type as the treatment.var column, if treatment.var is specified.}
#'  \item{treatment.var2: }{A column that has the same type as the treatment.var column, and is named treatment.var with a 2 appended to it.}
#'  \item{centroid_distance_diff: }{A numeric column that is the euclidean distance between the centroids of two treatments, if "centroid_distance" is specified as the composition metric. Centroid distance is based on plotted distance between centroids, it is context dependent and depends on how many centroids are being plotted.}
#'  \item{BC_dissim_diff: }{A numeric column that is the average pairwise Bray-Curtis dissimilarity of two treatments, if "ave_BC_dissim" is specified as the composition change metric.}
#'  \item{dispersion_diff: }{A numeric column that is the difference in the average dispersion of the replicates around the centriod for the two consecutive time periods. A negative value indicates there is more dispersion in treatment than treatment2 and a postive value indicates there is more dispersion in treatment2 than treatment.}
#'  \item{time.var: }{A column that has same name and type as the time.var column, if time.var is specified.}
#' }
#' @examples 
#' data(pplots)
#' #With treatment
#' multivariate_difference(pplots,
#'                         time.var="year", 
#'                         replicate.var = "plot", 
#'                         treatment.var = "treatment", 
#'                         species.var = "species", 
#'                         abundance.var = "relative_cover",
#'                         metric = "centroid_distance") # centroid distance as composition metric
#'   
#' #Without treatment
#' df <- subset(pplots, treatment == "N1P0")
#' multivariate_difference(df, 
#'                         time.var="year", 
#'                         replicate.var = "plot", 
#'                         species.var = "species", 
#'                         abundance.var = "relative_cover",
#'                         metric = "ave_BC_dissim") # Average Bray-Curtis dissimilarity as composisition metric
#' @importFrom vegan vegdist
#' @importFrom stats aggregate as.formula
#' @references Avolio et al. 2015; Avolio et al. OUR PAPER, Mari Anderson?
#' @export
multivariate_difference <- function(df, time.var, 
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

  if (is.null(time.var)) {
    output <- dissim_change(df, species.var, abundance.var, replicate.var, treatment.var)
  } else {
    # calculate change for each treatment
    splitvars <- time.var
    X <- split(df, df[splitvars])
    out <- lapply(X, FUN = mult_diff, species.var, abundance.var, replicate.var, treatment.var, comp_metric = comp_metric)
    unsplit <- lapply(out, nrow)
    unsplit <- rep(names(unsplit), unsplit)
    output <- do.call(rbind, c(out, list(make.row.names = FALSE)))
    output[splitvars] <- do.call(rbind, as.list(unsplit))
  }
    
  output_order <- c(
    time.var,
    treatment.var,
    paste(treatment.var,"2", sep=""),
    'BC_dissim_diff', 'centroid_distance_diff', 'dispersion_diff')

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
mult_diff <- function(df, time.var, species.var, abundance.var, replicate.var, comp_metric) {
  
  df1<-df[order(df[[treatment.var]], df[[replicate.var]]),]
  df2<-subset(df1, select = c(treatment.var, species.var, abundance.var, replicate.var))
  df2$id <- paste(df2[[treatment.var]], df2[[replicate.var]], sep="##")
  species <- codyn:::transpose_community(df2, 'id', species.var, abundance.var)
  species1 <- species
  species$id <- row.names(species)
  speciesid <- do.call(rbind.data.frame, strsplit(species$id, split="##"))
  colnames(speciesid)[1] <- treatment.var
  colnames(speciesid)[2] <- replicate.var
  species2 <- cbind(speciesid, species)
  species3 <- subset(species2, select = -id)
  
  #calculate bray-curtis dissmilarity
  bc <- vegdist(species3[,3:ncol(species3)], method="bray")
  
  #calculate distances of each plot to year centroid (i.e., dispersion)
  disp <- betadisper(bc, species3[[treatment.var]], type="centroid")
  
  
  #### doing compositional differences
  if(comp_metric == "ave_BC_dissim"){
    bc1<-as.data.frame(as.matrix(vegdist(species1, method = "bray")))
    #extracting lower diagonal
    bc2 <- as.data.frame(cbind(rownames(bc1)[which(lower.tri(bc1), arr.ind=T)[,1]],
                               colnames(bc1)[which(lower.tri(bc1), arr.ind=T)[,2]],
                               bc1[lower.tri(bc1)]))
    c1 <- as.data.frame(do.call('rbind', strsplit(as.character(bc2$V1), "##", fixed = T)))
    c2 <- as.data.frame(do.call('rbind', strsplit(as.character(bc2$V2), "##", fixed = T)))
    
    bc3 <- cbind(bc2, c1, c2)
    bc3$bc_dissim <- as.numeric(as.character(bc3$V3))
    colnames(bc3)[4] <- paste(treatment.var, 2, sep="")
    colnames(bc3)[6] <- treatment.var
    bc3$compare <- ifelse(bc3[[treatment.var]] == bc3[[paste(treatment.var, 2, sep="")]], 1, 2)
    
    #between time differences
    bc_between <- subset(bc3, compare == 2)
    myformula2 <- as.formula(paste("bc_dissim", "~", treatment.var, "+", paste(treatment.var, 2, sep = "")))
    bc_between_ave <- aggregate(myformula2, mean, data=bc_between)
    colnames(bc_between_ave)[3] <- "BC_dissim_diff"
    
    comp <- bc_between_ave
    
  } else {
    #getting distances between centroids over years; these centroids are in BC space, so that's why this uses euclidean distances
    cent_dist <- as.matrix(vegdist(disp$centroids, method="euclidean"))
    
    cent_dist2 <- as.data.frame(cbind(rownames(cent_dist)[which(lower.tri(cent_dist, diag=T), arr.ind=T)[,1]],
                                      colnames(cent_dist)[which(lower.tri(cent_dist, diag=T), arr.ind=T)[,2]],
                                      cent_dist[lower.tri(cent_dist, diag=T)]))
    cent_dist3 <- cent_dist2[cent_dist2$V1 != cent_dist2$V2,]
    cent_dist3[3]<-as.numeric(as.character(cent_dist3[[3]]))
    
    colnames(cent_dist3)[1] <- paste(treatment.var, 2, sep="")
    colnames(cent_dist3)[2] <- treatment.var
    colnames(cent_dist3)[3] <- "centroid_distance_diff"
    
    comp <- cent_dist3
  }
  
  ##doing dispersion
  disp2 <- data.frame(treatment=species3[[treatment.var]],
                      dist = disp$distances)
  
  myformula <- as.formula(paste("dist", "~", treatment.var))
  disp2.2 <- aggregate(myformula, mean, data=disp2)
  
  #mege into get dispersion for each treatment
  cent_dist_disp <- merge(comp, disp2.2, by = treatment.var)
  cent_dist_disp2 <- merge(cent_dist_disp, disp2.2, by.x = paste(treatment.var, 2, sep = ""), by.y = treatment.var)
  
  cent_dist_disp2$dispersion_diff <- cent_dist_disp2$dist.y - cent_dist_disp2$dist.x
  
  cent_dist_disp2$dist.x <- NULL
  cent_dist_disp2$dist.y <- NULL
  
  return(cent_dist_disp2)
}
