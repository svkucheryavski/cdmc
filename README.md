# Calibration design with multiple components

This repository contains R code implementing algorithm for generating design matrices for mixtures. For example, in case if you need to train multivariate models for predicting concentration of individual components in chemical mixtures (e.g. based on spectroscopic data). The original algorithm is described in:

D. Kirsanov, V. Panchuk, M. Agafonova-Moroz, M. Khaydukova, A. Lumpov, V. Semenov, A. Legin. *A sample-effective calibration design for multiple components* Analyst, 2014, 139, 4303. [doi: 10.1039/C4AN00227J](https://10.1039/C4AN00227J)

The algorithm generates a design matrix (DM) of size $N \times n$, where $N$ is a number of mixtures and $n$ is a number of components. The algorithm finds/optimizes the values of the matrix using two criteria:

1. The mixtures are evenly distributed in the $n$-dimensional component space.
2. The individual mixtures are spread in this space as much as possible.

You can use this implementation in two ways. The simplest way is to use a wrapper, function `cdmc()`, which does everything. The more complicated but also more flexible way is to build your own optimization routine by combining different criteria and constrains. Both ways are described below. All examples can be also found in file `demo.R`.

## Generate design matrix by using a wrapper

Function `cdmc()` is a wrapper, which does everything automatically with a small set of parameters that can be tuned. The function implements the algorithm described in [1] as is and has only one additional option — quantization of concentration levels.

In order to generate a design matrix you need to provide three main parameters: `N` is a number of mixtures to generate, `xmin` is a vector with lower concentration boundaries for each component and `xmax`, which is a vector with the upper boundaries.

The code below generates design matrix with 30 mixtures with the following concentration
ranges: $C_1 \in [0, 10]$, $C_2 \in [10, 100]$, $C_3 \in [50, 300]$ (just remember to load the functions from the main file by running `source("cdmc.R")`):

```r
dm1 <- cdmc(30, xmin = c(0, 10, 50), xmax = c(10, 100, 300))
```

By default, the function `cdmc()` uses 30 iterations for optimization, you can increase this number by providing additional parameter `maxIter`:

```r
dm1 <- cdmc(30, xmin = c(0, 10, 50), xmax = c(10, 100, 300), maxIter = 100)
```

With more iterations it could take a bit longer time, but chances to get better DM are also higher.
You can check the design matrix visually, for example by making scatter plots for all pairs of components:

```r
plot(as.data.frame(dm1))
```

Or by computing different quality characteristics. There are three main characteristics which can be used to assess the DM quality:

* `Dmax` — the maximum absolute deviation (MAD) between empirical cumulative distribution of the generated mixtures and the theoretical expected cumulative distribution. This characteristic shows how uniform the distribution of DM values is. The ideal `Dmax` is $0$.
* `rmin` — the smallest distance between all pairs of points, corresponding to the mixtures (rows of DM). The larger the value — the better. This characteristic is computed for normalized concentrations, where every component values are in range $[0, 1]$. If this characteristic is equal to $0$ it means that the design matrix contains two identical mixtures.
* `maxCor` — shows the largest absolute correlation between all pairs of components (columns of DM). The smaller — the better, however a value $<0.3$ can be considered as acceptable.

You can compute and show the three quality characteristics by using function `getDMQuality()`:

```r
getDMQuality(dm1)
```

You can also quantize concentration values to a given number of levels. The number can be set individual for each component or the same for all.


The code below shows how to generate DM for the same number of mixtures and concentration limits, as in the example above, but the concentration values will be quantized so there are 11 levels for $C_1 \in {0, 1, 2, ..., 10}$, 10 levels for $C_2 \in {10, 20, 30, ..., 100}$, and 6 levels for $C_3 \in {50, 100, 150, 200, 250, 300}$.

```r
dm2 <- cdmc(30, xmin = c(0, 10, 50), xmax = c(10, 100, 300), nlevels = c(11, 10, 6))
```

Note, that in this case DM can contain identical mixtures, especially if number of levels is small. It is recommended to check the `rmin` parameter and repeat the procedure if it is equal to 0.

Because in the method the initial values for DM as well as increments inside the optimization loop are selected randomly, different runs of the procedure provide different results. To get best results, it is recommended to use large number of iteration and/or run the procedure several times until the generated DM has satisfactory quality.


## Build your own optimization routine

You can also build your own optimization routine by selecting optimization criteria (or even add a manual one) and constraints.

The main function for optimization is `getOptimizedX()`. It takes initially generated DM (`X`) or DM from previous optimization run and try to improve it by shifting randomly its values. When it shifts the values it creates a new matrix, with optimized values, `Xp`. After that, it compares the two matrices in order to see if `Xp` is "better" than `X`. If yes, it replaces `X` with `Xp` and takes another iteration. If not, it tries to reduce the magnitude of random shifts and check if it makes an improvement.

There are two ways to change the optimization result.

1. Introduce constraints. Constraints are transformations which can be applied to the computed values (`Xp`) before comparing its quality with the previous version of DM. For example, you can quantize the DM values or do any other things. This can be done by providing a constraint function, which takes a matrix and returns a matrix of the same size, after applying the constraints.

2. Select which criteria to use when compare the quality of `Xp` and its previous version (`X`). This can be done by using a criterion function which takes both matrices as arguments and returns a logical value. If this value is `TRUE` the `Xp` is *not* good enough and the loop continues. If this value is `FALSE` the `Xp` is better than the previous version of DM and must replace the previous version.

### Changing the optimization criteria

There are three built in function for checking the optimization criteria:

* `critMaxDeviation(X, Xp)` — returns `FALSE` only if maximum absolute deviation for `Xp` is smaller than for `X`, so values in `Xp` are more evenly distributed.
* `critMinDistance(X, Xp)` — returns `FALSE` only if smallest distance between all points in `Xp` is larger than in `X`, so values in `Xp` are more spread in the component space.
* `critMaxCorrelation(X, Xp)` — returns `FALSE` only if maximum absolute correlation between all columns of `Xp` is smaller than for `X`, so values in `Xp` are less correlated.

Here are some examples.

Optimization which is based only on minimization of *Dmax*:

```r
# number of components
n <- 3

# number of mixtures
N <- 30

# Stage 1: Initialize X with random values from [0, 1]
X <- sapply(1:n, function(i) runif(N, 0, 1))

# Stage 2. Optimize X to get MAD smaller
X <- getOptimizedX(X, critFun = critMaxDeviation)
```

Optimization which is based only on maximization of smallest distance between two points (in this and all other examples we skip definition of `n` and `N`).

```r
# Stage 1: Initialize X with random values from [0, 1]
X <- sapply(1:n, function(i) runif(N, 0, 1))

# Stage 2: Optimize X to get rmin larger
X <- getOptimizedX(X, critFun = critMinDistance)
```

You can also combine the criteria or create your own. The code below shows how to implement all stages described in [1], where on stage 1 an initial version of X is generated, on stage 2 it is optimized to get larger MAD and on stage 3 it is further optimized using both criteria: minimization of MAD and maximization of smallest distance between two mixtures.

```r
# Stage 1: Initialize X with random values from [0, 1]
X <- sapply(1:n, function(i) runif(N, 0, 1))

# Stage 2. Optimize X to get MAD smaller
X <- getOptimizedX(X, critFun = critMaxDeviation)

# Stage 3. Run another optimization loop with both criteria
X <- getOptimizedX(X,
   critFun = function(X, Xp) critMaxDeviation(X, Xp) || critMinDistance(X, Xp))
```

In all examples above the concentration of components in DM are in range [0, 1]. To map them to the real ranges function `mapValues()` can be used:

```r
DM <- mapValues(X, xmin = c(0, 10, 100), xmax = c(10, 100, 300))
```

More examples can be found in `demo.R`.

