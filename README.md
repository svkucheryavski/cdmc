# Calibration design with multiple components

This repository contains R code implementing algorithm for generating a design matrix for mixtures. For example for training multivariate models to predict concentration of individual components in mixtures of chemical components (e.g. based on spectroscopic data). The algorithm is described in:

D. Kirsanov, V. Panchuk, M. Agafonova-Moroz, M. Khaydukova, A. Lumpov, V. Semenov, A. Legin. *A sample-effective calibration design for multiple components* Analyst, 2014, 139, 4303. [doi: 10.1039/C4AN00227J](https://10.1039/C4AN00227J)

The algorithm generates a design matrix of size $N \times n$, where $N$ is a number of mixtures and $n$ is a number of components. The algorithm finds/optimizes the values of the matrix to meet two criteria:

1. The mixtures are evenly distributed in the $n$-dimensional component space.
2. The individual mixtures are spread in this space as much as possible.

The implementation consists of two R functions. Function `getDmax()` computes maximum absolute deviation (MAD) between the distribution frequencies of any mixture matrix, `X`, and the theoretical expectations. Function `getOptimizedX()` runs optimization loop with some user defined criteria (for example, minimization of the MAD or/and maximization of minimum distance between two points in the component space).

Function `getOptimized()` requires a handler to anonymous function which will compute the logical value based on four parameters: MAD for previous X and new optimized X (`d1`, `d2`), minimum distance between mixtures in previous X and new optimized X (`r1`, `r2`). The function should take these four values as arguments and compute a logical value based on them. If logical value is `TRUE` it means that the optimized X is not good enough and loop must continue, if the logical value is `FALSE`, it means that the optimized X meets the optimization criteria and next iteration can be taken.

Here are some examples:

```r
# continue if MAD for optimized X is larger than MAD for previous X
f1  = function(d1, d2, r1, r2) (d1 > d2)

# continue if MAD for optimized X is larger than MAD for previous X
# or if minimum distance for optimized X is smaller than minimum distance for previous X
f2  = function(d1, d2, r1, r2) (d1 > d2 || r1 < r2)
```

The code chunk below shows how algorithm works using two stages. This code can be also found in file `demo.R`.

```r
# load the two functions from "cdmc.R"
source("cdmc.R")

# here we will create design matrix for 4 components
# first we define boundaries for each component
# (e.g. concentration of C1 is between 1 and 20mM, etc)
xmin <- c(1, 5, 10, 100)
xmax <- c(20, 100, 300, 300)

# then we define how many mixtured we need
N <- 30

# and maximum number of iterations
maxIter <- 100


# Stage 1: Initialize X with random values bounded by xmin and xmax
n <- length(xmin)
X <- sapply(1:n, function(i) runif(N, xmin[i], xmax[i]))

# Compute MAD and mimimal distance between any two mixtures in X
cat("\nInitial X:\n")
print(c(DXmax = getDmax(X, xmin, xmax), rXmin = min(dist(X))))

# Stage 2. Optimize X to get MAD smaller
X <- getOptimizedX(X, xmin, xmax, critFun = function(d1, d2, r1, r2) (d1 > d2))

# Compute MAD and mimimal distance between any two mixtures in the optimized X
cat("\nAfter stage 2:\n")
print(c(DXmax = getDmax(X, xmin, xmax), rXmin = min(dist(X))))

# Stage 3: Optimize X to get MAD smaller and minimal distance between two mixtures larger
X <- getOptimizedX(X, xmin, xmax, critFun = function(d1, d2, r1, r2) (d1 > d2 || r1 < r2))

# Compute MAD and mimimal distance between any two mixtures in the optimized X
cat("\nAfter stage 3:\n")
print(c(DXmax = getDmax(X, xmin, xmax), rXmin = min(dist(X))))

# show the boxplot and scatter plots for the design matrix
boxplot(X)
plot(as.data.frame(X))
```
