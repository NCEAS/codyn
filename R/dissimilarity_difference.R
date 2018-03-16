#'@title Bray curtis dissimilarity differences of replicates between and within treatments
#'@description Calculates the average difference in Bray-Curtis dissimilarity of replicates treatments and within treatments at a single point in time. Between treatment differences is the average Bray-Curtis dissimilaritiy of all pairwise comparisions of replicates between two treatments. This is a measure of how dissimilar the community composition of two treatments is. Bray-Curtis dissimilarity ranges from 0-1, where 0 are identical communities, and 1 and completelty different communiites. Within difference is derived in two steps. First, the average Bray-Curtis dissimilarity of all replicates within a treatment is calculated. This is a measure of how homogenous a community is within a treatment. Then, these averages are compared between two treatments.
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
#'  \item{BC_between_diff: }{A numeric column that is the average pairwise Bray-Curtis dissimilarity of replicates in two treatments. 0 - The communities of treatments are similar. 1 - The communities of treatments are very different.}
#'  \item{BC_within_diff: }{A numeric column that is difference between two treatments in the average pairwise Bray-Curtis dissimilarity of replicates within each treatment. A positive number indicates that treatment.var2 has varaiblity in community composition than treatment.var.}
#'  \item{time.var: }{A characteric column that has the same name and type as the time.var column, if specified.}
#' }
#' @references Our Avolio et al. paper, Avolio et al. 2015, Marti Anderson?
#' @importFrom vegan vegdist
#' @importFrom stats aggregate as.formula
#' @examples  
#' data(pplots)
#' #Without time
#' df <- subset(pplots, year == 2002)
#' dissimilarity_difference(df, 
#'                         replicate.var = "plot", 
#'                         treatment.var = "treatment", 
#'                         species.var = "species", 
#'                         abundance.var = "relative_cover")
#' #With time
#' dissimilarity_difference(pplots, 
#'                         time.var = "year", 
#'                         replicate.var = "plot", 
#'                         species.var = "species", 
#'                         abundance.var = "relative_cover", 
#'                         treatment.var = "treatment")
#' @export
dissimilarity_difference <- function(df, time.var = NULL, species.var,
                                    abundance.var, replicate.var, treatment.var) {
  
  # check no NAs in abundance column
  if(any(is.na(df[[abundance.var]]))) stop("Abundance column contains missing values")
  
  if(is.null(time.var)){
    
    # check there unique species x time combinations
    check_single_onerep(df, replicate.var, species.var)
    
    output <- dissim_diff(df, species.var, abundance.var, replicate.var, treatment.var)
    
  } else {
    
    # check unique species x time x replicate combinations
    check_single(df, time.var, species.var, replicate.var)
    
    splitvars <- time.var
    X <- split(df, 
               df[splitvars])
    out <- lapply(X, FUN = dissim_diff, species.var, abundance.var, replicate.var, treatment.var)
    unsplit <- lapply(out, nrow)
    unsplit <- rep(names(unsplit), unsplit)
    output <- do.call(rbind, c(out, list(make.row.names = FALSE)))
    output[splitvars] <- do.call(rbind, as.list(unsplit))
  }
  
  output_order <- c(
    time.var,
    treatment.var, paste(treatment.var, 2, sep = ''),
    'BC_between_diff', 'BC_within_diff')
  
  return(output[intersect(output_order, names(output))])
}

############################################################################
#
# Private functions: these are internal functions not intended for reuse.
# Future package releases may change these without notice. External callers
# should not use them.
#
############################################################################

# A function calculate the average Bray-Curtis dissimilarity of all pairwise comparisions of replicates either between two treatments or within a treatment. Once this is calcualted, for within_diff, the differences of within_dissimliarity is calculated betweent two treatments.
# @param df a dataframe
# @param species.var the name of the species column
# @param replicate.var the name of the replicate column
# @param treatment.var the name of the treatment column

dissim_diff <- function(df, species.var, abundance.var, replicate.var, treatment.var){
  df1<-df[order(df[[treatment.var]], df[[replicate.var]]),]
  #transpose data
  df2<-subset(df1, select = c(species.var, abundance.var, replicate.var, treatment.var))
  df2$id <- paste(df2[[treatment.var]], df2[[replicate.var]], sep="##")
  species <- codyn:::transpose_community(df2, 'id', species.var, abundance.var)
 
  #calculate BC dissimilarity
  bc <- as.data.frame(as.matrix(vegdist(species, method="bray")))
  
  #extracting lower diagonal
  bc2 <- as.data.frame(cbind(rownames(bc)[which(lower.tri(bc), arr.ind=T)[,1]],
                             colnames(bc)[which(lower.tri(bc), arr.ind=T)[,2]],
                             bc[lower.tri(bc)]))
  c1 <- as.data.frame(do.call('rbind', strsplit(as.character(bc2$V1), "##", fixed = T)))
  c2 <- as.data.frame(do.call('rbind', strsplit(as.character(bc2$V2), "##", fixed = T)))
  
  bc3 <- cbind(bc2, c1, c2)
  bc3$bc_dissim <- as.numeric(as.character(bc3$V3))
  colnames(bc3)[4] <- paste(treatment.var, 2, sep="")
  colnames(bc3)[6] <- treatment.var
  bc3$compare <- ifelse(bc3[[treatment.var]] == bc3[[paste(treatment.var, 2, sep="")]], 1, 2)
  
  #within treatment differences
  bc_within <- subset(bc3, compare == 1)
  myformula <- as.formula(paste("bc_dissim", "~", treatment.var))
  bc_within_ave <- aggregate(myformula, mean, data=bc_within)
  colnames(bc_within_ave)[2] <- "BC_dissim_within"
  
  #between treatment differences
  bc_between <- subset(bc3, compare == 2)
  myformula2 <- as.formula(paste("bc_dissim", "~", treatment.var, "+", paste(treatment.var, 2, sep = "")))
  bc_between_ave <- aggregate(myformula2, mean, data=bc_between)
  colnames(bc_between_ave)[3] <- "BC_between_diff"
  
  #mege into get bc_within differences for each treatment
  bc_dis1 <- merge(bc_between_ave, bc_within_ave, by = treatment.var)
  bc_dis <- merge(bc_dis1, bc_within_ave, by.x = paste(treatment.var, 2, sep = ""), by.y = treatment.var)
  bc_dis$BC_within_diff <- bc_dis$BC_dissim_within.y - bc_dis$BC_dissim_within.x
  bc_dis$BC_dissim_within.x <- NULL
  bc_dis$BC_dissim_within.y <- NULL

  return(bc_dis)
}