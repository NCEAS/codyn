# @title Pool replicates into treatments and add ranks
# @description This fucntion first takes the average abundnace of each species
#   in a treatment across all replicates to create a single pooled
#   replicate for each treatment. Then it add ranks to each species following
#   the add_ranks function.
# @param df A data frame containing an optional time column and requred species,
#   abundance, replicate and treatment columns
# @param time.var The name of the optional time column 
# @param species.var The name of the species column 
# @param abundance.var The name of the abundance column 
# @param replicate.var The name of the replicate column 
# @param treatment.var The name of the treatment column
# @return The add_ranks function returns a data frame with the following attributes:
# \itemize{
#  \item{treatment.var: }{A column that as the same name and type as the treatment.var column.}
#  \item{species.var: }{A column that has same name and type as the species.var column.}
#  \item{time.var: }{{A column that has same name and type as the time.var column.}
#  \item{abundance.var: }{A column that has same name and type as the abundance.var column.}
#  \item{replicate.var: }{A column that has same name and type as the replicate.var column.}
#  \item{rank: }{A numeric column with the species rank; a rank of 1 indicates the species was most abundant in that time period. All species that are not present in that time period have the rank value S+1 where S is the number of species in the sample.}
# }
pool_replicates <- function(df, time.var=NULL, species.var, abundance.var,
                            replicate.var, treatment.var) {
  
  df <- as.data.frame(df)

  # isolate rep-trts
  rep_trt <- unique(subset(df, select = c(replicate.var, treatment.var)))
  
  # add zeros for species absent from a replicate within a treatment
  if (is.null(time.var)) {
    allsp <- fill_zeros(df, species.var, abundance.var)
  } else {
    by <- c(time.var)
    allsp <- split_apply_combine(df, by, FUN = fill_zeros, species.var, abundance.var)
  }
  
  # get averages of each species by treatment and (optionally) time
  by <- c(time.var, treatment.var, species.var)
  spave <- aggregate.data.frame(allsp[abundance.var], allsp[by], FUN = mean)
  
  # add ranks after splitting into (pooled) communities
  by <- c(treatment.var, time.var)
  rankdf <- split_apply_combine(spave, by, FUN = add_ranks, abundance.var)
  
  return(rankdf)
  
}
