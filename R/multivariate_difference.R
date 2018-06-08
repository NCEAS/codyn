#'@title Using dissimilarity-based metrics to calculate differences in
#'  composition and dispersion between pairs of treatments at a single time
#'  point
#'  
#'@description Calculates the difference in composition and dispersion between  treatments based off a Bray-Curtis dissimilarity matrix at a single point in time. Composition difference is the euclidean distance between the centroids of different treatments. Since centroid distance is based on plotted distance between centroids, it is context dependent and depends on how many centroids are being plotted. The centroid distance between treatments depends on how many treatments are being compared. Dispersion difference is the difference of average dispersion of each replicate to its centroid between two treatments.
#' @param df A data frame containing an optional time column, species, abundance and replicate, and treatment columns
#' @param time.var The name of the optional time column 
#' @param species.var The name of the species column 
#' @param abundance.var The name of the abundance column 
#' @param replicate.var The name of the replicate column. Replicate must be unique within the dataset and cannot be nested within treatments or blocks. 
#' @param treatment.var the name of the treatment column
#' @param reference.treatment The name of the optional treatment that all other treatments will be compared to (e.g. only controls will be compared to all other treatments). If not specified all pairwise treatment comparisons will be made.
#' @return The multivariate_difference function returns a data frame with the
#'   following attributes:
#' \itemize{
#'  \item{treatment.var: }{A column that has same name and type as the treatment.var column, if treatment.var is specified.}
#'  \item{treatment.var2: }{A column that has the same type as the treatment.var  column, and is named treatment.var with a 2 appended to it.}
#'  \item{composition_diff: }{A numeric column that is the euclidean distance  between the centroids of two treatments at a single point in time.}
#'  \item{abs_dispersion_diff: }{A numeric column that is the absolute value of the difference in the average dispersion of the replicates around the centroid for the two treatments.}
#'  \item{trt_greater_disp: }{A column that has same type as the treatment.var  column, and specifies which of the two  treatments has greater dispersion.}
#'  \item{time.var: }{A characteristic column that has the same name and type as the time.var column, if specified.}
#' }
#' @references Avolio et al. Submitted to MEE, Avolio et al. 2015, Marti Anderson et al. 2006
#' @importFrom vegan vegdist betadisper
#' @importFrom stats aggregate reshape
#' @examples  
#' data(pplots)
#' # Without time
#' df <- subset(pplots, year == 2002)
#' multivariate_difference(df, 
#'                         replicate.var = "plot", 
#'                         treatment.var = "treatment", 
#'                         species.var = "species", 
#'                         abundance.var = "relative_cover")
#' # There are 6 replicates for each of three treatments, thus 18 total
#' # observations.
#' 
#' # Without time and with reference treatment
#' df <- subset(pplots, year == 2002)
#' multivariate_difference(df, 
#'                         replicate.var = "plot", 
#'                         treatment.var = "treatment", 
#'                         species.var = "species", 
#'                         abundance.var = "relative_cover",
#'                         reference.treatment = "N1P0")
#' # There are 6 replicates for each of three treatments, thus 18 total
#' # observations.
#' 
#' # With time
#' multivariate_difference(pplots, 
#'                         time.var = "year", 
#'                         replicate.var = "plot", 
#'                         species.var = "species", 
#'                         abundance.var = "relative_cover", 
#'                         treatment.var = "treatment")
#' # In each year there are 6 replicates for each of three treatments, for a
#' # total of 18 observations.
#' @export
multivariate_difference <- function(df,
                                    time.var = NULL,
                                    species.var,
                                    abundance.var,
                                    replicate.var,
                                    treatment.var,
                                    reference.treatment = NULL) {
  
  # validate function call and purge extraneous columns
  args <- as.list(match.call()[-1])
  df <- do.call(check_args, args, envir = parent.frame())
  
  by <- time.var
  output <- split_apply_combine(df, by, FUN = mult_diff, time.var, species.var,
    abundance.var, replicate.var, treatment.var, reference.treatment)
  
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

# A function calculate the community compositon difference (the distance between
# the centroids of two treatments) and absolute dispersion difference (the
# absolute difference in the average dispersion of the replicates around the
# centroid for the two treatments).
# @param df a dataframe
# @param species.var the name of the species column
# @param replicate.var the name of the replicate column
# @param treatment.var the name of the treatment column

mult_diff <- function(df, time.var, species.var, abundance.var,
                      replicate.var, treatment.var, reference.treatment) {

  #transpose data
  idvar <- c(time.var, replicate.var, treatment.var)
  species <- reshape(df, idvar = idvar, timevar = species.var,
    direction = 'wide')
  species[is.na(species)] <- 0
  
  #calculate bray-curtis dissimilarities
  abund <- species[, -(1:length(idvar))]
  bc <- vegdist(abund, method = "bray")
  
  #calculate distances of each plot to treatment centroid (i.e., dispersion)
  if (is.null(time.var)) {
    message('Composition and dispersion difference calculation using ',
            nrow(species),' observations.')
  } else {
    message('Composition and dispersion difference calculation using ',
            nrow(species), ' observations at time.var value ',
            df[[time.var]][[1]])
  }
  disp <- betadisper(bc, species[[treatment.var]], type = "centroid")
  
  #getting distances between treatments with euclidean distances
  cent_dist <- as.matrix(vegdist(disp$centroids, method = "euclidean"))
  
  #extracting all treatment differences
  treatment.var2 <- paste(treatment.var, 2, sep = '')
  if (is.null(reference.treatment)) {
    lt <- lower.tri(cent_dist, diag = T)
    idx <- which(lt, arr.ind = T)
    cent_dist2 <- as.data.frame(cbind(
      rownames(cent_dist)[idx[,1]],
      colnames(cent_dist)[idx[,2]],
      cent_dist[lt]))
    cent_dist3 <- cent_dist2[cent_dist2$V1 != cent_dist2$V2,]
    cent_dist3[3] <- as.numeric(as.character(cent_dist3[[3]]))
  } else {
    treatments <- colnames(cent_dist)
    idx <- treatments == reference.treatment
    cent_dist3 <- data.frame(
      V1 = treatments[!idx],
      V2 = treatments[idx],
      V3 = cent_dist[which(!idx), idx])
  }
  
  colnames(cent_dist3)[1] <- treatment.var2
  colnames(cent_dist3)[2] <- treatment.var
  colnames(cent_dist3)[3] <- "composition_diff"
  
  # collecting and labeling distances to centroid from betadisper to get a
  # measure of dispersion and then take the mean for a treatment
  species['dist'] <- disp['distances']
  by <- c(time.var, treatment.var)
  disp2 <- aggregate.data.frame(species['dist'], species[by], mean)
  
  # merge into get dispersion for each treatment
  cent_dist_disp <- merge(cent_dist3, disp2, by = treatment.var)
  cent_dist_disp2 <- merge(cent_dist_disp, disp2[c(treatment.var, 'dist')],
    by.x = treatment.var2, by.y = treatment.var)
  
  # calculate absolute difference
  cent_dist_disp2[[treatment.var]] <- as.character(
    cent_dist_disp2[[treatment.var]])
  cent_dist_disp2[[paste(treatment.var, 2, sep = "")]] <- as.character(
    cent_dist_disp2[[paste(treatment.var, 2, sep = "")]])
  
  cent_dist_disp2$abs_dispersion_diff <- abs(
    cent_dist_disp2$dist.x - cent_dist_disp2$dist.y)
  cent_dist_disp2$trt_greater_disp <- as.character(ifelse(
    cent_dist_disp2$dist.x > cent_dist_disp2$dist.y,
    cent_dist_disp2[[treatment.var]],
    cent_dist_disp2[[treatment.var2]]))
  
  return(cent_dist_disp2)
}
