rm(list = ls())
source("cdmc.R")

# boundaries for each component
xmin <- c(1, 5, 10, 100)
xmax <- c(20, 100, 300, 300)

# number of mixtures
N <- 30

# maximum number of iterations
maxIter <- 100


# initialize X with random values bounded by xmin and xmax
n <- length(xmin)
X <- sapply(1:n, function(i) runif(N, xmin[i], xmax[i]))
cat("\nInitial X:\n")
print(c(DXmax = getDmax(X, xmin, xmax), rXmin = min(dist(X))))

# optimize X to get DXmax smaller
X <- getOptimizedX(X, xmin, xmax, critFun = function(d1, d2, r1, r2) (d1 > d2))
cat("\nAfter stage 2:\n")
print(c(DXmax = getDmax(X, xmin, xmax), rXmin = min(dist(X))))

# optimize X to get DXmax smaller and rXmin larger
X <- getOptimizedX(X, xmin, xmax, critFun = function(d1, d2, r1, r2) (d1 > d2 || r1 < r2))
cat("\nAfter stage 3:\n")
print(c(DXmax = getDmax(X, xmin, xmax), rXmin = min(dist(X))))

# show the boxplot and scatter plots
boxplot(X)
plot(as.data.frame(X))
