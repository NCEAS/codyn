#'@title Multivariate differences in composition and dispersion
#' @description 
#' @param df A data frame containing an optional time column, species, abundance and replicate, and treatment columns
#' @param time.var The name of the optional time column 
#' @param species.var The name of the species column 
#' @param abundance.var The name of the abundance column 
#' @param replicate.var The name of the replicate column 
#' @param treatment.var the name of the treatment column
#' 
multivariate_difference <- function(df, time.var=NULL, species.var, abundance.var, replicate.var, treatment.var){
  
  if(is.null(time.var)){
    
    output <- mult_diff(df, species.var, abundance.var, replicate.var, treatment.var)
    
  }
  
  else{
    
    df <- df[order(df[[time.var]]),]
    X <- split(df, df[time.var])
    out <- lapply(X, FUN = mult_diff, species.var, abundance.var, replicate.var, treatment.var)
    ID <- unique(names(out))
    out <- mapply(function(x, y) "[<-"(x, time.var, value = y) ,
                  out, ID, SIMPLIFY = FALSE)
    output <- do.call("rbind", out)
    
  }
  
  return(output)
}

###private functions

mult_diff <- function(df, species.var, abundance.var, replicate.var, treatment.var){
  require(vegan)
  
  #transpose data
  df2<-subset(df, select = c(species.var, abundance.var, replicate.var, treatment.var))
  df2$id <- paste(df2[[treatment.var]], df2[[replicate.var]], sep="##")
  species<-codyn:::transpose_community(df2, 'id', species.var, abundance.var)
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
  
  #getting distances between treatments; these centroids are in BC space, so that's why this uses euclidean distances
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