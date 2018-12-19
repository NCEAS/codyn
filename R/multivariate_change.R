#'@title Using dissimilarity-based metrics to calculate changes in composition and dispersion
#' @description Calculates the changes in composition and dispersion based off a Bray-Curtis dissimilarity matrix. Composition change is the euclidean distance between the centroids of compared time periods and ranges from 0-1, where 0 are identical communities, and 1 and completely different communities. Since composition change is based on plotted distance between centroids, it is context dependent and depends on how many centroids are being plotted, thus result of composition change depends on how many time periods are being measured. Dispersion change is the difference of average dispersion of each replicate to its centroid between time periods.
#' @param df A data frame containing time, species, abundance and replicate columns and an optional column of treatment
#' @param time.var The name of the time column
#' @param species.var The name of the species column
#' @param abundance.var The name of the abundance column
#' @param replicate.var The name of the replicate column. Replicate must be unique within the dataset and cannot be nested within treatments or blocks.
#' @param treatment.var the name of the optional treatment column
#' @param reference.time The name of the optional time point that all other time points should be compared to (e.g. the first year of data). If not specified, each comparison is between consecutive time points (e.g. first to  second year, second to third year, etc.)
#' @return The multivariate_change function returns a data frame with the
#'   following attributes:
#' \itemize{
#'  \item{time.var: }{A column with the specified time.var and a second column, with '2' appended to the name. Time is subtracted from time2 for dispersion change.}
#'  \item{composition_change: }{A numeric column that is the euclidean distance between the centroids of two time points.}
#'  \item{dispersion_change: }{A numeric column that is the difference in the   average dispersion of the replicates around the centroid for the two time periods. A negative value indicates replicates are converging over time (there is less dispersion at time period 2 than time period 1) and a positive value indicates replicates are diverging over time (there is more dispersion at time period 2 than time period 1.}
#'  \item{treatment.var: }{A column that has same name and type as the treatment.var column, if treatment.var is specified.}
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
#' @importFrom vegan vegdist betadisper
#' @importFrom stats aggregate reshape
#' @references Avolio et al. 2015; Avolio et al. Submitted MEE, Mari Anderson et al. 2006.
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

  # notify user ## FIXME, too much
  if (is.null(treatment.var)) {
    message('Composition and dispersion change calculation using ',
    length(unique(df[[time.var]])),' observations.')
  } else {
    lapply(split(df, df[[treatment.var]], drop = TRUE), FUN = function(df) {
      message('Composition and dispersion change calculation using ',
        length(unique(df[[time.var]])), ' observations at treatment.var value ',
        df[[treatment.var]][[1]])
    })
  }

  # calculate replicate centers and dispersion for each time [and treatment]
  by <- c(treatment.var, time.var)
  centers <- split_apply_combine(df, by, FUN = find_centers,
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
  output$composition_change <- mapply(FUN = function(x, y) {
    x <- data.frame(var = names(x), center2 = x)
    y <- data.frame(var = names(y), center = y)
    xy <- merge(x, y, all = TRUE)
    xy[is.na(xy)] <- 0
    bcdism(xy$center2, t(as.matrix(xy$center)))
  }, output$center2, output$center)

  output_order <- c(
    time.var, time.var2, treatment.var,
    'composition_change', 'dispersion_change')

  return(output[intersect(output_order, names(output))])
}

############################################################################
#
# Private functions: these are internal functions not intended for reuse.
# Future package releases may change these without notice. External callers
# should not use them.
#
############################################################################
# A function calculate the community compositon change (the distance between the
# centriods of two consecutive time periods) and dispersion change (the
# difference in the average dispersion of the replicates around the centriod for
# the two consecutive time periods). For dispersion change a negative value
# indicates replicates are converging over time and a postive value indicates
# replicates are diverging over time.
#
find_centers <- function(df, time.var, species.var, treatment.var, replicate.var) {

    # spread species out as columns
    idvar <- c(treatment.var, replicate.var, time.var)
    species <- reshape(df, idvar = idvar, timevar = species.var, direction = 'wide')
    species[is.na(species)] <- 0

    # the abundance.var matrix
    a <- species[, -(1:length(idvar))]
    # idx <- colMeans(a) == 0
    # a <- a[, !idx]

    # the sum of squared dissimilarities
    f <- function(x) {
        d <- bcdism(x, a)
        return(sum(d^2))
    }
    # the gradient of f wrt each abundance.var, for method = 'bray' only
    grad <- function(x) {
        d <- bcdism(x, a)
        s <- sign(t(x - t(a)))
        t(s - d) %*% (1/(sum(x) + rowSums(a))*d)
    }
    # numerical optimization constrained for non-negative abundances
    ## FIXME probably will error with a zero column
    copt <- constrOptim(
        colMeans(a), f, grad,
        ui = diag(ncol(a)), ci = rep(0, ncol(a))
    )

    # return a data frame with one row with a list column containing center and
    # a numeric column containing dispersion
    centers <- species[1, c(treatment.var, time.var), drop = FALSE]
    centers$dispersion <- mean(bcdism(a, copt$par)) ## FIXME mean of dissimilarities, right?
    centers$center <- list(copt$par)
    return(centers)
}

# bray curtis dissimilarity between each row of a and the vector x
bcdism <- function(x, a) {
    colSums(abs(x - t(a))) / colSums(x + t(a))
}
