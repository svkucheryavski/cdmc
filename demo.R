rm(list = ls())
source("cdmc.R")

# common parameters for all examples

## number of mixtures
N <- 30

## lower and upper boundaries for concentration levels
xmin <- c(0, 10, 50)
xmax <- c(10, 100, 300)

# Part 1. Use a built in wrapper ----

# Generate design matrix with 30 mixtures and 3 components in the following concentration
# ranges: C1 [0, 10], C2 [10, 100], C3 [50, 300]
dm1 <- cdmc(N, xmin, xmax)
getDMQuality(dm1)

show(dm1)
boxplot(dm1)
plot(as.data.frame(dm1))

# Do the same as in previous example but quantize concentration levels so
# there are 11 levels for C1, 10 levels for C2, and 6 levels for C3
dm2 <- cdmc(N, xmin, xmax, nlevels = c(11, 10, 6))
getDMQuality(dm2)

show(dm2)
boxplot(dm2)
plot(as.data.frame(dm2))


# Part 2. Construct your own optimization routines ----

# initialize X with random values bounded by xmin and xmax
n <- length(xmin)
X <- sapply(1:n, function(i) runif(N, 0, 1))
cat("\nInitial X:\n")
print(getDMQuality(X))

# Criterion 1: optimize X to get DXmax smaller
X1 <- getOptimizedX(X, critFun = critMaxDeviation)
cat("\nMinimization of DXmax:\n")
print(getDMQuality(X1))

# Criterion 2: optimize X to get rmin larger
X2 <- getOptimizedX(X, critFun = critMinDistance)
cat("\nMaximization of rXmin:\n")
print(getDMQuality(X2))

# Criterion 3: optimize X to get maximum correlation smaller
X3 <- getOptimizedX(X, critFun = critMaxCorrelation)
cat("\nMinimization of correlation between components:\n")
print(getDMQuality(X3))

# Criteria 1 and 2:
X12 <- getOptimizedX(X, critFun = function(X, Xp) critMaxDeviation(X, Xp) || critMinDistance(X, Xp))
cat("\nCriteria 1 and 2:\n")
print(getDMQuality(X12))

# Criteria 2 and 3:
X23 <- getOptimizedX(X, critFun = function(X, Xp) critMaxCorrelation(X, Xp) || critMinDistance(X, Xp))
cat("\nCriteria 2 and 3:\n")
print(getDMQuality(X23))


# Criteria 2 and 3 and quantization of values to 10 levels:
Xq <- quantizeLevels(X, 10)

X23q <- getOptimizedX(Xq,
   maxIter = 100,
   constrFun = function(X) quantizeLevels(X, 10),
   critFun = function(X, Xp) critMinDistance(X, Xp) || critMaxDeviation(X, Xp)
)

cat("\nCriteria 2 and 3 + quantization of values to 10 levels:\n")
q <- rbind(getDMQuality(Xq), getDMQuality(X23q))
rownames(q) <- c("Initial", "After optimization")
print(q)

# show  scatter plots for the last case
plot(as.data.frame(X23q))
