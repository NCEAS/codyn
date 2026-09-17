#' Konza data from Collins et al. 2008
#'
#' A dataset of tallgrass prairie plant composition at one annually burned and one unburned
#' site over time at the Konza Prairie LTER, Manhattan Kansas (Collins et al. 2008).
#'
#' A data frame containing a column for replicate, year, species and abundance:
#' \itemize{
#'   \item replicate: A factor column of spatial replicates with two levels ("annually burned" and "unburned")
#'   \item year: An integer column of sampling time points
#'   \item species: A factor column of species sampled
#'   \item abundance: A numeric column of abundance values
#' }
#'
#' @source
#' Collins, Scott L., Katharine N. Suding, Elsa E. Cleland, Michael Batty, Steven C. Pennings, Katherine L. Gross, James B. Grace, Laura Gough, Joe E. Fargione, and Christopher M. Clark. (2008) "Rank clocks and plant community dynamics." Ecology 89, no. 12: 3534-41.
#' @docType data
#' @keywords datasets
#' @name collins08
#' @usage data(collins08)
#' @format A data frame with 2058 rows and 4 variables
"collins08"

#' Data from Konza Prairie, watershed 001d
#'
#' Plant composition within multiple replicates at an annually burned tallgrass
#' prairie site in the Konza Prairie LTER, Manhattan KS (Watershed 001d).
#'
#' A data frame containing a column for species, year, subplot and abundance:
#' \itemize{
#'   \item species: A factor column of species sampled
#'   \item year: An integer column of sampling time points
#'   \item subplot: A factor column of spatial replicates with 20 levels
#'   \item abundance: A numeric column of abundance values
#' }
#'
#' @source
#' Konza Prairie LTER Dataset ID: PVC02, watershed 1D
#'
#' Collins, S. L. (2000) Disturbance frequency and community stability in native tallgrass prairie. American Naturalist 155:311-325.
#' @docType data
#' @keywords datasets
#' @name knz_001d
#' @usage data(knz_001d)
#' @format A data frame with 8768 rows and 4 variables
"knz_001d"

#' Phosphorus plots data from Avolio et al. 2014
#' 
#' A dataset of tallgrass prairie plant composition in a nitrogen and phosphorus addition experiment at Konza Prairie, Manhattan Kansas (Avolio et al. 2014). This dataset is a subset of the full dataset.
#'
#' A data frame containing a column for replicate, year, species, abundance, block and treatment :
#' \itemize{
#'   \item plot: An integer column of spatial replicates with 18 levels (6-48)
#'   \item year: An integer column of sampling time points
#'   \item species: A factor column of species sampled
#'   \item relative_cover: A numeric column of relative cover values
#'   \item block: An integer column of dummy blocking variable, grouping treatment plots into blocks
#'   \item treatment: A factor column of nitrogen and phosphorus treatments applied to the plots
#' }
#'
#' @source
#' Avolio, ML, Koerner, S, La Pierre, K, Wilcox, K, Wilson, GTW, Smith, MD, Collins, S. 2014. Changes in plant community composition, not diversity, to a decade of nitrogen and phosphorus additions drive changes in aboveground productivity in a tallgrass prairie. Journal of Ecology. 102: 1649-1660.  
#' @docType data
#' @keywords datasets
#' @name pplots
#' @usage data(pplots)
#' @format A data frame with 1232 rows and 6 variables
"pplots"
