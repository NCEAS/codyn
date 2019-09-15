#' @title Using dissimilarity-based measures to calculate changes in composition
#'   and dispersion
#'
#' @description Calculates the changes in composition and dispersion based off a
#'   Bray-Curtis dissimilarity matrix. Composition change is the pairwise
#'   distance between centroids of compared time periods and ranges from 0-1,
#'   where identical communities give 0 and completely different
#'   communities give 1. Dispersion change is the difference between time periods in
#'   the dispersion of replicates, i.e. the average distance between a replicate
#'   and its centroid.
#'
#' @inheritParams RAC_change
#' @param df A data frame containing time, species, abundance and replicate
#'   columns and an optional column of treatment.
#' @param replicate.var The name of the replicate column. Replicate identifiers
#'   must be unique within the dataset and cannot be nested within treatments or
#'   blocks.
#' @param treatment.var The name of the optional treatment column.
#'
#' @return The multivariate_change function returns a data frame with the
#'   following attributes:
#' \itemize{
#'  \item{time.var: }{A column with the specified time.var and a second column,
#'  with '2' appended to the name. Time is subtracted from time2 for dispersion
#'  change.}
#'  \item{composition_change: }{A numeric column that is the distance
#'  between the centroids of two time points, or NA if a real distance
#'  could not be calculated.}
#'  \item{dispersion_change: }{A numeric column that is the difference in the
#'  average dispersion of the replicates around the centroid for the two time
#'  periods. A negative value indicates replicates are converging over time
#'  (there is less dispersion at time period 2 than time period 1) and a
#'  positive value indicates replicates are diverging over time (there is more
#'  dispersion at time period 2 than time period 1.}
#'  \item{treatment.var: }{A column that has same name and type as the
#'  treatment.var column, if treatment.var is specified.}
#' }
#' @examples
#' data(pplots)
#' # With treatment
#' multivariate_change(pplots,
#'                     time.var="year",
#'                     replicate.var = "plot",
#'                     treatment.var = "treatment",
#'                     species.var = "species",
#'                     abundance.var = "relative_cover")
#' # In each year there are 6 replicates and there are 4 years of data for 3
#' # time comparisons, thus 24 total observations in each treatment.
#'
#' # With treatment and reference year
#' multivariate_change(pplots,
#'                     time.var="year",
#'                     replicate.var = "plot",
#'                     treatment.var = "treatment",
#'                     species.var = "species",
#'                     abundance.var = "relative_cover",
#'                     reference.time = 2002)
#' # In each year there are 6 replicates and there are 4 years of data for 3
#' # time comparisons, thus 24 total observations in each treatment.
#'
#' # Without treatment
#' df <- subset(pplots, treatment == "N1P0")
#' multivariate_change(df,
#'                     time.var="year",
#'                     replicate.var = "plot",
#'                     species.var = "species",
#'                     abundance.var = "relative_cover")
#' # In each year there are 6 replicates and there are 4 years of data for 3
#' # time comparisons, thus 24 total observations.
#'
#' @importFrom vegan vegdist
#' @importFrom stats aggregate reshape
#' @references Avolio et al. 2015; Avolio et al. Submitted, Marti Anderson et al. 2006.
#' @export
multivariate_change <- function(df,
                                time.var,
                                species.var,
                                abundance.var,
                                replicate.var,
                                treatment.var = NULL,
                                reference.time = NULL){

  # validate function call and purge extraneous columns
  args <- as.list(match.call()[-1])
  df <- do.call(check_args, args, envir = parent.frame())

  # calculate replicate centers and dispersion [by treatment]
  by <- c(treatment.var)
  centers <- split_apply_combine(df, by, FUN = pca_centers,
      time.var, species.var, treatment.var, replicate.var)

  # merge subsets on time difference of one time step
  time.var2 <- paste(time.var, 2, sep = '')
  split_by <- c(treatment.var)
  merge_to <- !(names(centers) %in% split_by)
  if (is.null(reference.time)) {
      output <- split_apply_combine(centers, split_by, FUN = function(x) {
          y <- x[merge_to]
          cross <- merge(x, y, by = NULL, suffixes = c('', '2'))
          f <- factor(cross[[time.var]])
          f2 <- factor(cross[[time.var2]], levels = levels(f))
          idx <- (as.integer(f2) - as.integer(f)) == 1
          cross[idx, ]
      })
  } else {
      output <- split_apply_combine(centers, split_by, FUN = function(x) {
          y <- x[x[[time.var]] != reference.time, merge_to]
          x <- x[x[[time.var]] == reference.time, ]
          merge(x, y, by = NULL, suffixes = c('', '2'))
      })
  }

  # compute time.var2 change from time.var
  output$dispersion_change <- output$dispersion2 - output$dispersion
  idx <- abs(Im(output$dispersion_change)) > sqrt(.Machine$double.eps)
  output$dispersion_change <- Re(output$dispersion_change)
  if (any(idx)) {
    output$dispersion_change[idx] <- NA
    warning('NA(s) produced during dispersion change calculation.')
  }
  output$composition_change <- mapply(
      function(x, y) {z <- x - y; sqrt(sum(z*z))},
      output$center2, output$center
  )
  idx <- abs(Im(output$composition_change)) > sqrt(.Machine$double.eps)
  output$composition_change <- Re(output$composition_change)
  if (any(idx)) {
    output$composition_change[idx] <- NA
    warning('NA(s) produced during centroid change calculation.')
  }


  output_order <- c(
    time.var, time.var2, treatment.var,
    'composition_change', 'dispersion_change')

  return(output[intersect(output_order, names(output))])
}

#' @title Using dissimilarity-based measures to calculate differences in
#'   composition and dispersion between pairs of treatments at a single time
#'   point
#'
#' @description Calculates the difference in composition and dispersion between
#'   treatments based off a Bray-Curtis dissimilarity matrix at a single point
#'   in time. Composition difference is the pairwise distance between centroids
#'   of compared treatments and ranges from 0-1, where identical communities
#'   give 0 and completely different communities give 1. Dispersion difference
#'   is the difference between treatments in the dispersion of replicates, i.e.
#'   the average distance between a replicate and its centroid.
#'
#' @inheritParams RAC_difference
#' @param df A data frame containing a species, abundance, replicate, and
#'   treatment columns and optional time column.
#' @param treatment.var The name of the treatment column.
#'
#' @return The multivariate_difference function returns a data frame with the
#'   following attributes:
#' \itemize{
#'  \item{treatment.var: }{A column that has same name and type as the
#'  treatment.var column, if treatment.var is specified.}
#'  \item{treatment.var2: }{A column that has the same type as the treatment.var
#'  column, and is named treatment.var with a 2 appended to it.}
#'  \item{composition_diff: }{A numeric column that is the euclidean distance
#'  between the centroids of two treatments at a single point in time.}
#'  \item{abs_dispersion_diff: }{A numeric column that is the absolute value of
#'  the difference in the average dispersion of the replicates around the
#'  centroid for the two treatments.}
#'  \item{trt_greater_disp: }{A column that has same type as the treatment.var
#'  column, and specifies which of the two  treatments has greater dispersion.}
#'  \item{time.var: }{A characteristic column that has the same name and type as
#'  the time.var column, if specified.}
#' }
#' @references Avolio et al. Submitted, Avolio et al. 2015, Marti Anderson et al. 2006
#' @importFrom vegan vegdist
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

  # calculate replicate centers and dispersion [by treatment]
  by <- c(time.var)
  centers <- split_apply_combine(df, by, FUN = pca_centers,
                                 treatment.var, species.var, time.var, replicate.var)

  # order treatment.var if unordered factor
  to_ordered <- is.factor(centers[[treatment.var]]) &
    !is.ordered(centers[[treatment.var]]) &
    is.null(reference.treatment)
  if (to_ordered) {
    class(centers[[treatment.var]]) <- c('ordered', class(centers[[treatment.var]]))
  }

  # merge subsets on treatment differences [split by time]
  treatment.var2 <- paste(treatment.var, 2, sep = '')
  split_by <- c(time.var)
  merge_to <- !(names(centers) %in% split_by)
  if (is.null(reference.treatment)) {
    output <- split_apply_combine(centers, split_by, FUN = function(x) {
      y <- x[merge_to]
      cross <- merge(x, y, by = NULL, suffixes = c('', '2'))
      idx <- cross[[treatment.var]] < cross[[treatment.var2]]
      cross[idx, ]
    })
  } else {
    output <- split_apply_combine(centers, split_by, FUN = function(x) {
      y <- x[x[[treatment.var]] != reference.treatment, merge_to]
      x <- x[x[[treatment.var]] == reference.treatment, ]
      merge(x, y, by = NULL, suffixes = c('', '2'))
    })
  }

  # unorder treatment.var if orginally unordered factor
  if (to_ordered) {
    x <- class(output[[treatment.var]])
    class(output[[treatment.var]]) <- x[x != 'ordered']
    class(output[[treatment.var2]]) <- x[x != 'ordered']
  }

  # compute treatment.var2 differences from treatment.var
  output$dispersion_diff <- output$dispersion2 - output$dispersion
  output$trt_greater_disp <- output[[treatment.var]]
  idx <- abs(Im(output$dispersion_diff)) > sqrt(.Machine$double.eps)
  output$dispersion_diff <- Re(output$dispersion_diff)
  if (any(idx)) {
    output$dispersion_diff[idx] <- NA
    output$trt_greater_disp[idx] <- NA
    warning('NA(s) produced during dispersion difference calculation.')
  }
  idx <- !is.na(output$dispersion_diff) & output$dispersion_diff > 0
  output$trt_greater_disp[idx] <- output[[treatment.var2]][idx]
  output$abs_dispersion_diff <- abs(output$dispersion_diff)
  output$composition_diff <- mapply(
    function(x, y) {z <- x - y; sqrt(sum(z*z))}, # not dsp*Conj(dsp), not Euclidean!
    output$center2, output$center
  )
  idx <- abs(Im(output$composition_diff)) > sqrt(.Machine$double.eps)
  output$composition_diff <- Re(output$composition_diff)
  if (any(idx)) {
    output$composition_diff[idx] <- NA
    warning('NA(s) produced during centroid difference calculation.')
  }

  output_order <- c(
    time.var, treatment.var, treatment.var2,
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
# A function calculate the community compositon change (distance between the
# centriods) and dispersion change (
# difference in the average dispersion of the replicates around the centriod for
# the two consecutive time periods). For dispersion change a negative value
# indicates replicates are converging over time and a postive value indicates
# replicates are diverging over time.
#
pca_centers <- function(df, cluster.var, species.var, split.var, replicate.var) {

    # spread species out as columns
    idvar <- c(cluster.var, split.var, replicate.var)
    species <- reshape(df, idvar = idvar, timevar = species.var, direction = 'wide')
    species[is.na(species)] <- 0

    # update user
    if (is.null(split.var)) {
        message('Composition and dispersion change calculation using ',
            nrow(species),' observations.')
    } else {
        message('Composition and dispersion change calculation using ',
            nrow(species), ' observations for ',
            df[[split.var]][[1]])
    }

    # the Bray-Curtis dissimilarity matrix
    a <- species[, -(1:length(idvar))]
    d <- as.matrix(vegdist(a, 'bray'))

    # perform PCoA aka Torgerson-Gower Scaling, while
    # violating the assumption that d results from a metric
    b <- -1/2 * dblctr(d^2)
    eig <- eigen(b, symmetric = TRUE)
    prc <- eig$vectors %*% diag(sqrt(as.complex(eig$values)), nrow=nrow(b))

    # centroids and dispersion (Marti Anderson, doi:10.1111/j.1541-0420.2005.00440.x)
    by <- c(split.var, cluster.var)
    ctr <- aggregate(prc, species[by], mean)
    idx <- match(interaction(species[by]), interaction(ctr[by]))
    dsp <- prc - ctr[idx, -(1:length(by))]
    dsp <- sqrt(rowSums(dsp*dsp)) # not dsp*Conj(dsp), not Euclidean
    dsp <- aggregate(dsp, species[by], mean)

    # return a data frame with a list column containing centroids and
    # a numeric column containing dispersion
    result <- ctr[by]
    result$center <- unname(split(unname(
        as.matrix(ctr[, -(1:length(by))])
        ), 1:nrow(ctr)))
    result$dispersion <- dsp$x
    return(result)
}

# Double centering of a dissimilarity matrix
# @param d dissimilarity matrix
dblctr <- function(d) {
    n <- nrow(d)
    c <- diag(1, n) - matrix(1, n, n)/n
    return(c %*% d %*% c)
}
