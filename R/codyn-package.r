#' Ecological Community Dynamics
#' @description Temporal diversity indices and community stability metrics for ecologists.
#' @details The functions in codyn implement metrics that are explicitly temporal, and include the option to  calculate them over multiple replicates. Functions fall into two categories: temporal diversity indices and community stability metrics. The diversity indices in codyn are temporal analogs to traditional diversity indices such as richness and rank-abundance curves. Specifically, codyn includes functions to calculate species turnover, mean rank shifts and lags in community similarity between time points. The community stability metrics in codyn calculate overall stability and patterns of species covariance and synchrony over time. Finally, codyn contains vignettes that describe methods and reproduce figures from published papers to help users contextualize and apply functions to their own data.
#' 
#' 
#' @author Lauren M. Hallett
#' @docType package
#' @name codyn
#' @aliases codyn
NULL

#' Konza data from Collins et. all 2008
#'
#' A dataset of tallgrass prairie plant composition at one annually burned and one unburned 
#' site over time at the Konza Prairie LTER, Manhattan Kansas (Collins et al. 2008).
#' 
#' \itemize{
#'   \item replicate The replicate type, i.e. "burned", "unburned"
#'   \item year The sampling frequency
#'   \item species The sampled species
#'   \item abundance The measure of species abundance
#' }
#'
#' @source
#' Scott L. Collins 1987. Interaction of Disturbances in Tallgrass Prairie: A Field Experiment. 
#' Ecology 68:1243–1250. http://dx.doi.org/10.2307/1939208
#' @docType data
#' @keywords datasets
#' @name collins08
#' @usage data(collins08)
#' @format A data frame with 2058 rows and 4 variables
NULL

#' Data from Konza Prairie, watershed 001d
#'
#' Plant composition within multiple replicates at an annually burned tallgrass 
#' prairie site in the Konza Prairie LTER, Manhattan KS (Watershed 001D). 
#' 
#' species,year,subplot,abundance
#' \itemize{
#'   \item species The sampled species
#'   \item year The sampling frequency
#'   \item subplot The sampled subplot
#'   \item abundance The measure of species abundance
#' }
#'
#' @docType data
#' @keywords datasets
#' @name knz_001d
#' @usage data(knz_001d)
#' @format A data frame with 8768 rows and 4 variables
NULL

