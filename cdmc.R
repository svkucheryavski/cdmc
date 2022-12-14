################################################################################
# Implementation of algorithm for calibration design with multiple components  #
# described in: http://dx.doi.org/10.1039/c4an00227j with some extra options   #
################################################################################


##############################################################
# Helper functions                                           #
##############################################################


#' Check the design matrix
#'
#' @param X
#' Matrix or data frame with numeric columns with values between 0 and 1
#'
#' @return
#' Same matrix as input
#'
#' @details
#' Stops execution and show error message if X is not a matrix or a data frame
#' with numeric columns, if number of columns in X is smaller than number of rows
#' or if values of X are not in range [0, 1]
#'
checkX <- function(X) {

   if (is.null(dim(X))) {
      stop("Parameter 'X' must be a matrix or a data frame with numeric columns.")
   }

   if (nrow(X) <= ncol(X)) {
      stop("Number of rows in 'X' must be larger than the number of columns.")
   }

   if (min(X) < 0 || max(X) > 1) {
      stop("Values in 'X' must be between 0 and 1.")
   }

   return (X)
}


#' Multivariate binning
#'
#' @param X
#' Matrix or data frame with numeric columns.
#' @param breaks
#' A vector with breaks for binning, can be also a list of vectors for each column of X
#'
#' @return
#' A distribution table with counts
#'
binning <- function(X, breaks) {

   if (!is.list(breaks)) {
      breaks <- rep(list(breaks), ncol(X))
   }

   if (length(breaks) != ncol(X)) {
      stop("Parameter 'break' must be a vector or a list of vectors - one for each column of X.")
   }

  return (table(lapply(1:ncol(X), function(i) cut(X[, i], breaks[[i]]))))
}


#' Maps values of matrix columns from [0, 1] to [xmin, xmax] range
#'
#' @param X
#' Matrix or data frame with numeric columns with values in range [0, 1].
#' @param xmin
#' Vector with lower boundaries for each column of X
#' @param xmax
#' vector with upper boundaries for each column of X
#'
#' @export
mapValues <- function(X, xmin, xmax) {
   X <- checkX(X)
   sapply(1:ncol(X), function(i) X[, i] * (xmax[i] - xmin[i]) + xmin[i])
}


#' Quantize values in range [0, 1] to n levels
#'
#' @param X
#' Matrix or data frame with numeric columns with values in range [0, 1].
#' @param n
#' Number of levels or a vector with the numbers for each column of X
#'
#' @return
#' Matrix X with quantized values
#'
#' @export
quantizeLevels <- function(X, n) {

   X <- checkX(X)

   if (length(n) == 1 && ncol(X) > 1) {
      n <- rep(n, ncol(X))
   }

   return (sapply(1:ncol(X), function(i) (round(X[, i] * (n[i] - 1)) / (n[i] - 1))))
}



##############################################################
# Functions which compute various quality parameters for X   #
##############################################################

#' Compute maximum absolute correlation between all pairs of columns
#'
#' @param X
#' Matrix or data frame with numeric columns.
#'
#' @return
#' A positive value from [0, 1] - the maximum absolute correlation
#'
#' @details
#' The method computes correlation matrix first, then sets diagonal elements to zeros
#' and return maximum absolute value of the matrix.
#'
getMaxAbsCor <- function(X) {
   return (max(abs(cor(X) - diag(1, ncol(X)))))
}


#' Compute smallest distance between all pairs of points in X (rows of X)
#'
#' @param X
#' Matrix or data frame with numeric columns.
#'
#' @return
#' A number - the smallest distance.
#'
#'
getMinDist <- function(X) {
   return (min(dist(X)))
}


#' Compute maximum absolute deviation between the cumulative distribution of values in X
#' and theoretical expectations
#'
#' @param X
#' Matrix or data frame with numeric columns with values in range [0, 1]
#'
#' @details
#' First the variable space is split into smallest part of subvolumes for given dimension of X.
#' Then cumulative distribution of values from X among the subvolumes is computed and is
#' compared with theoretical expectations (assuming that distribution is uniform). The maximum
#' absolute difference between the empirical and expected values is returned.
#' See full description in [1].
#'
#' @references
#' 1. Kirsanov et. al. Analyst, 2014, 139, 4303–4309 (doi:10.1039/C4AN00227J)
#'
#' @return
#' A number - maximum absolute deviation
#'
#' @export
getDmax <- function(X) {

   # check that X meets necessary requirements
   X <- checkX(X)

   # number of components and number of mixtures
   n <- ncol(X)
   N <- nrow(X)

   # number of segments for each component
   k <- floor(N ^ (1/n) + 0.0000001)

   # number of subvolumes
   m <- k^n

   # compute expected values for cumulative sum of frequencies (i * N / m)
   volInd <- seq_len(m)
   expected <- volInd * N / m

   # find distribution of the values among the subvolumes
   breaks <- seq(-0.00000001, 1, length = k + 1)

   # compute cumulative sum of the frequencies in subvolumes
   S <- cumsum(binning(X, breaks))

   # compute and return maximum absolute deviation
   return (max(abs(S - expected)))
}

#' Compute three quality parameters for a design matrix
#'
#' @param X
#' Design matrix
#'
#' @return
#' vector with three values, Dmax, rmin, and largest absolute correlation
#'
#' @export
getDMQuality <- function(X) {
   X <- sapply(1:ncol(X), function(i) (X[, i] - min(X[, i])) / (max(X[, i]) - min(X[, i])))
   return (c(Dmax = getDmax(X), rmin = getMinDist(X), maxCor = getMaxAbsCor(X)))
}

######################################################################
# Functions which check if X and Xp do not satisfy various criteria  #
######################################################################

#' Check if Dmax of Xp is not better (larger or equal) than Dmax of X
#'
#' @param X
#' Design matrix under optimization
#' @param Xp
#' The optimized version of X
#'
#' @details
#' Compute Dmax value for Xp and X and return TRUE is this value for Xp is larger than this
#' value for X (so the optimization loop will continue)
#'
#' @export
critMaxDeviation <- function(X, Xp = NULL) {
   if (is.null(Xp)) return (TRUE)
   return ((getDmax(X) - getDmax(Xp)) <= 0)
}

#' Check if rmin of Xp is not better (smaller or equal) than rmin of X
#'
#' @param X
#' Design matrix under optimization
#' @param Xp
#' The optimized version of X
#'
#' @details
#' Compute rmin value (minimum distance between any two points) for Xp and X and return TRUE if
#' this value for Xp is smaller than this value for X (so the optimization loop will continue)
#'
#' @export
critMinDistance <- function(X, Xp = NULL) {
   if (is.null(Xp)) return (TRUE)
   return ((getMinDist(Xp) - getMinDist(X)) <= 0)
}


#' Check if maximum absolute correlation between columns for Xp is not better than for X
#'
#' @param X
#' Design matrix under optimization
#' @param Xp
#' The optimized version of X
#'
#' @details
#' Compute maximum absolute correlation value for all column pairs of Xp and X and return TRUE if
#' this value for Xp is larger than this value for X (so the optimization loop will continue)
#'
#' @export
critMaxCorrelation <- function(X, Xp = NULL) {
   if (is.null(Xp)) return (TRUE)
   I <- diag(1, ncol(X))
   return ((getMaxAbsCor(X) - getMaxAbsCor(Xp)) <= 0)
}


######################################################################
# Main optimization procedure and wrappers                           #
######################################################################

#' Optimize X values by random shifts with optimization constraints
#'
#' @param X
#' Design matrix with current values to be optimized
#' @param critFun
#' Handler to a function which computes optimization criterion for continue the while loop
#' @param maxIter
#' Maximum number of iterations
#'
#' @return
#' Matrix of the same size as X with optimized values
#'
#' @details
#'
#' The function takes current design matrix X (either initialized with random numbers or after
#' previous optimization step) and optimizes it by adding/subtracting random values to each
#' column. The algorithm is described in [1].
#'
#' Function for `critFun` should take the current optimized design matrix `Xp`, and the original `X`,
#' and return a logical value. If the returned value is TRUE, the `Xp` is not good enough and the
#' optimization must continue. If the returned value is FALSE, the optimization loop stops and next
#' iteration is taken.
#'
#' Both X and the result of optimization must have values ranging from 0 to 1.
#'
#' @references
#' 1. Kirsanov et. al. Analyst, 2014, 139, 4303–4309 (doi:10.1039/C4AN00227J)
#'
#' @export
getOptimizedX <- function(X, critFun, constrFun = NULL, maxIter = 100) {


   # check that X meets requirements
   checkX(X)

   # number of components and number of mixtures
   n <- ncol(X)
   N <- nrow(X)

   # number of segments for each component
   k <- floor(N ^ (1/n))

   # number of subvolumes
   m <- k^n

   # set dxmax parameter (maximum shift for each component) as a half of segment size
   dxmax <- 1 / k

   # apply constraints
   if (!is.null(constrFun)) {
      X <- constrFun(X)
   }

   # outer loop
   for (i in seq_len(maxIter)) {

      # initialize alpha and dX
      alpha <- 1
      dX <- sapply(1:n, function(i) runif(N, -dxmax/2, dxmax/2))

      Xp <- NULL
      niter <- 1
      while (critFun(X, Xp)) {

         # compute optimized X by shifting the values
         Xp <- X + alpha * dX

         # check if values are outside boundaries set them to boundary values
         Xp[Xp < 0] <- 0
         Xp[Xp > 1] <- 1

         # apply constraints
         if (!is.null(constrFun)) {
            Xp <- constrFun(Xp)
         }

         # decrease alpha for next iteration
         alpha <- alpha / 2

         # if after maxIter iterations no improvement is observed - stop internal loop
         niter <- niter + 1
         if (niter > maxIter) {
             Xp <- X
             break
         }
      }

      # reinitialize the X matrix with optimized one and repeat the loop
      X <- Xp
   }
   return(Xp)
}


#' Generates design matrix with N mixtures
#'
#' @param N
#' Number of mixtures to generate
#' @param xmin
#' Vector with lowest concentration boundaries for each component
#' @param xmax
#' Vector with largest concentration boundaries for each component
#' @param varNames
#' Component names (will be used as column names of generated design matrix)
#' @param maxIter
#' Maximum number of iterations
#'
#' @return
#' A data frame with generated design matrix
#'
#' @detail
#' The function implements the whole procedure for generating of design matrix X with N
#' mixtures described in [1].
#'#'
#' @references
#' 1. Kirsanov et. al. Analyst, 2014, 139, 4303–4309 (doi:10.1039/C4AN00227J)
#'
#' @examples
#'
#' # Generate design matrix with 30 mixtures and 3 components in the following concentration
#' # ranges: C1 [0, 10], C2 [10, 100], C3 [50, 300]
#'
#' dm <- cdmc(30, c(0, 10, 50), c(10, 100, 300))
#' show(dm)
#'
#' # Do the same as in previous example but quantize concentration levels so
#' # there are 11 levels for C1 and C2 and 7 levels for C3
#'
#' dm <- cdmc(30, c(0, 10, 50), c(10, 100, 300), nlevels = c(11, 11, 7))
#' show(dm)
#'
#' @export
cdmc <- function(N, xmin, xmax, varNames = NULL, nlevels = NULL, maxIter = 30) {

   # check inputs

   if (length(xmin) < 2 || length(xmin) != length(xmax)) {
      stop("Parameters 'xmin' and 'xmax' must be vectors of the same size.")
   }

   # set number of components
   n <- length(xmin)

   if (N <= n) {
      stop("Number of mixtures (parameters 'N') must be larger than number of components.")
   }

   if (!is.null(varNames) && length(varNames) != n) {
      stop("Number of elements in vector 'varNames' must be the same as number of components.")
   }

   # define constraints function
   constrFun <- if (is.null(nlevels)) NULL else (function(X) quantizeLevels(X, nlevels))

   # initialize X
   X <- sapply(1:n, function(i) runif(N, 0, 1))

   # Stage 2: optimize X to get DMax smaller
   X <- getOptimizedX(X, constrFun = constrFun, maxIter = maxIter,
      critFun = critMaxDeviation)

   # Stage 3: optimize X to get Dmax smaller and rmin larger
   X <- getOptimizedX(X, constrFun = constrFun, maxIter = maxIter,
      critFun = function(X, Xp) critMaxDeviation(X, Xp) || critMinDistance(X, Xp))

   # map the values to real concentration, assign colnames and return
   X <- mapValues(X, xmin, xmax)
   colnames(X) <- if (is.null(varNames)) paste0("C", 1:n) else varNames

   return (X)
}