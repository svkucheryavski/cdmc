################################################################################
# Implementation of algorithm for calibration design with multiple components  #
# described in: http://dx.doi.org/10.1039/c4an00227j                           #
################################################################################

#' Compute maximum absolute deviation between the distribution of values in X
#' and theoretical expectations
#'
#' @param X
#' Design matrix with current values to be optimized
#' @param xmin
#' Vector with lower boundaries for each component
#' @param xmax
#' Vector with upper boundaries for each component
#'
#' @details
#' See full description in [1].
#'
#' @references
#' 1. Kirsanov et. al. Analyst, 2014, 139, 4303–4309 (doi:10.1039/C4AN00227J)
#'
#' @export
getDmax <- function(X, xmin, xmax) {

   # number of components and number of mixtures
   n <- ncol(X)
   N <- nrow(X)

   # number of segments for each component
   k <- floor(N ^ (1/n))

   # number of subvolumes
   m <- k^n

   # compute expected values for cumulative sum of frequencies (i * N / m)
   volInd <- seq_len(m)
   expected <- volInd * N / m

   # find distribution of the values among the subvolumes
   Xints <- sapply(1:n, function(i) cut(X[, i], breaks = seq(xmin[i], xmax[i], length = k + 1)))

   # compute cumulative sum of the frequencies in subvolumes
   S <- cumsum(xtabs(~., as.data.frame(Xints)))

   # compute and return maximum absolute deviation
   return (max(abs(S - expected)))
}

#' Optimize X values by random shifts with optimization constraints
#'
#' @param X
#' Design matrix with current values to be optimized
#' @param xmin
#' Vector with lower boundaries for each component
#' @param xmax
#' Vector with upper boundaries for each component
#' @param critFun
#' Handler to a function which computes optimization criterion for continue the while loop
#' @param xmax
#' Vector with upper boundaries for each component
#' @param maxIter
#' Maximum number of iterations
#'
#' @return
#' Matrix of the same size as X with optimized values
#'
#' @details
#'
#' The function takes current design matrix X (either initialized with random numbers or after
#' previous optimization step) and optimize it by adding/subtracting random values to each
#' column. The algorithm is described in [1].
#'
#' Function for `critFun` should take `DXpmax`, `DXmax`, `rXpmin`, `rXmin` values for current
#' state of optimized X (Xp) as arguments and returns a logical result based on the values, e.g.
#' `function(d1, d2, r1, r2) (d1 > d2)` will continue the loop if `DXpmax > DXpmax`.
#'
#' @references
#' 1. Kirsanov et. al. Analyst, 2014, 139, 4303–4309 (doi:10.1039/C4AN00227J)
#'
#' @export
getOptimizedX <- function(X, xmin, xmax, critFun, maxIter = 10) {

   # number of components and number of mixtures
   n <- ncol(X)
   N <- nrow(X)

   # number of segments for each component
   k <- floor(N ^ (1/n))

   # number of subvolumes
   m <- k^n

   # set dxmax parameter (maximum shift for each component) as a half of segment size
   dxmax <- (xmax - xmin) / (2 * k)

   # multiply xmin/xmax values for easier comparison
   Xmin <- matrix(xmin, N, n, byrow = TRUE)
   Xmax <- matrix(xmax, N, n, byrow = TRUE)

   # compute initial values for DXmax and rXmin
   DXmax <- getDmax(X, xmin, xmax)
   rXmin <- min(dist(X))

   # outer loop
   for (i in seq_len(maxIter)) {

      # initialize alpha and dX
      alpha <- 1
      dX <- sapply(1:n, function(i) runif(N, -dxmax[i]/2, dxmax[i]/2))

      # initialize DXpmax and rXpmin for optimization loop
      # so we have at least one iteration
      DXpmax <- DXmax + 1
      rXpmin <- rXmin * 0.1

      while (critFun(DXpmax, DXmax, rXpmin, rXmin)) {

         # compute optimized X by shifting the values
         Xp <- X + alpha * dX

         # check if values are outside boundaries set them to boundary values
         Xp[Xp < Xmin] <- Xmin[Xp < Xmin]
         Xp[Xp > Xmax] <- Xmax[Xp > Xmax]

         # compute DXpmax and rXpmax
         DXpmax <- getDmax(Xp, xmin, xmax)
         rXpmin <- min(dist(Xp))

         # decrease alpha for next iteration
         alpha <- alpha / 2
      }

      # reinitialize the X matrix with optimized one and repeat the loop
      X <- Xp
      DXmax <- DXpmax
      rXmin <- rXpmin
   }

   return(Xp)
}
