#' Community Dynamics Metrics
#' @description Univariate and multivariate temporal and spatial diversity indices, rank abundance curves, and community stability metrics for ecologists.
#' @details The functions in \code{codyn} implement metrics that are either explicitly temporal and include the option to  calculate them over multiple replicates, or spatial and include the option to calculate them over multiple time points.
#' Functions fall into five categories: static diversity indices, temporal diversity indices, spatial diversity indices, rank abundance curves, and community stability metrics.
#' The diversity indices in \code{codyn} are temporal and spatial analogs to traditional diversity indices.
#' Specifically, \code{codyn} includes functions to calculate community richness, evenness and diversity at a given point in space and time. In addition, \code{codyn} contains functions to calculate species turnover, mean rank shifts, and lags in community similarity between two time points. For the components of rank abundance curves, the shape of the rank abundance curve, species abundances and multivariate metrics of community composition, \code{codyn} contains functions to calculate these metrics either between time points or a single time point between two paired replicates.
#' The community stability metrics in \code{codyn} calculate overall stability and patterns of species covariance and synchrony over time.
#' Finally, \code{codyn} contains vignettes that describe methods and reproduce figures from published papers to help users contextualize and apply functions to their own data.
#' Work on this package was supported by National Science Foundation grant #1262458 to the National Center for Ecological Analysis and Synthesis (NCEAS), the University of Wisconsin, and the University of New Mexico, and a SESYNC Synthesis Postdoctoral Fellowship to MLA. 
#'
#' @author
#' - Lauren Hallett <lauren.m.hallett@gmail.com>
#' - Meghan L. Avolio <meghan.avolio@jhu.edu>
#' - Ian T. Carroll <icarroll@sesync.org>
#' - Sydney K. Jones <syd@sevilleta.unm.edu>
#' - A. Andrew A. MacDonald <aammacdonald@gmail.com>
#' - Dan F. B. Flynn <flynn@fas.harvard.edu>
#' - Peter Slaughter <slaughter@nceas.ucsb.edu>
#' - Julie Ripplinger <julie.ripplinger@asu.edu>
#' - Scott L. Collins <scollins@sevilleta.unm.edu>
#' - Corinna Gries <cgries@wisc.edu>
#' - Matthew B. Jones <jones@nceas.ucsb.edu>
#'
#' @section Functions:
#' \itemize{
#'  - [`turnover()`]: Calculates species turnover between time periods 
#'  - [`rank_shift()`]: Calculates the mean relative change in species rank abundances
#'  - [`rate_change()`]: Calculates the rate change in a community over time
#'  - [`rate_change_interval()`]: Produces a data frame containing differences in species composition between samples at increasing time intervals
#'  - [`community_stability()`]: Calculates community stability over time
#'  - [`variance_ratio()`]: Computes the ratio of the variance of aggregate species abundances in a community
#'  - [`synchrony()`]: Calculates the degree synchrony in species abundances
#'  - [`temporal_torus_translation()`]: Calculates a test statistic on a null ecological community created via cyclic shifts. \code{confint} provides mean and confidence intervals of this null distribution
#'  - [`community_structure()`]: Calculates richness and evenness (using specified metric) for a replicate
#'  - [`community_diversity()`]: Calculates diversity (using specified metric) for a replicate
#'  - [`RAC_change()`]: Calculates changes in species richness, evenness, species' ranks, gain, and losses for a replicate over time
#'  - [`abundance_change()`]: For each species in a replicate, calculates changes in abundance over time
#'  - [`curve_change()`]: Calculates changes in the shape of the RAC curve for each replicate over time
#'  - [`multivariate_change()`]: Calculates changes in community composition and dispersion over time
#'  - [`RAC_difference()`]: Calculates differences in species richness, evenness, species' ranks, shared species between paired samples at a single point in time
#'  - [`abundance_difference()`]: Calculates differences in abundance for each species in paired samples at a single point in time
#'  - [`curve_difference()`]: Calculates differences in the shape of the RAC between paired samples at a single point in time
#'  - [`multivariate_difference()`]: Calculates differences in community composition and dispersion of all replicates between treatments at a single point in time
#' }
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end

## mockable bindings: start
## mockable bindings: end
NULL