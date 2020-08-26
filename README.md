# **BDgraph** 
  
![](https://www.r-pkg.org/badges/version/BDgraph) ![](https://www.r-pkg.org/badges/last-release/BDgraph) ![](https://cranlogs.r-pkg.org/badges/BDgraph) 
![](https://cranlogs.r-pkg.org/badges/grand-total/BDgraph) 


The `R` package **BDgraph** provides statistical tools for Bayesian structure learning in undirected graphical models for *continuous*, *count*, *binary*, and *mixed data*. The package is implemented the recent improvements in the Bayesian graphical models literature, including [Mohammadi and Wit (2015)](https://projecteuclid.org/euclid.ba/1422468425), [Mohammadi et al. (2017)](https://rss.onlinelibrary.wiley.com/doi/full/10.1111/rssc.12171), [Dobra and Mohammadi (2018)](https://projecteuclid.org/euclid.aoas/1532743478), and [Letac et al. (2018)](https://arxiv.org/abs/1706.04416). Besides, the package contains several functions for simulation and visualization, as well as several multivariate datasets taken from the literature. To speed up the computations, the computationally intensive tasks of the package are implemented in `C++` in parallel using **OpenMP**.

## Installation

You can install the latest version from CRAN using:

``` r
install.packages( "BDgraph" )
```

``` r
require( "BDgraph" )
```

## Example 1: Gaussian Graphical Models

Here is a simple example to see the preformance of the package for the Gaussian graphical models. Frist, by using the function `bdgraph.sim` to simulate 100 observations (n = 100) from a multivariate
Gaussian distribution with 8 variables (p = 8) and “scale-free” graph structure, as follows:

``` r
data.sim = bdgraph.sim( n = 100, p = 8, graph = "scale-free", vis = TRUE )
round( head( data.sim $ data, 4 ), 2 )
```

Since the generated data are Gaussian, we run the `bdgraph` function by choosing `method = "ggm"`, as follows:

``` r
bdgraph.obj <- bdgraph( data = data.sim, method = "ggm", iter = 5000 )

summary( ssgraph.obj )
```

To compare the result with the true graph

``` r
compare( data.sim, bdgraph.obj, main = c( "Target", "BDgraph" ), vis = TRUE )
```

Now, as an alternative, we run the `bdgraph.mpl` function which is based on the GGMs and marginal pseudo-likelihood, as follows:

``` r
bdgraph.obj_mpl <- bdgraph.mpl( data = data.sim, method = "ggm", iter = 5000 )

summary( bdgraph.obj_mpl )
```

We could compare the results of both algorithms with the true graph as follows:

``` r
compare( data.sim, bdgraph.obj, bdgraph.obj_mpl, 
         main = c( "Target", "BDgraph", "BDgraph_mpl" ), vis = TRUE )
```

## Example 2: Gaussian Copula Graphical Models

Here is a simple example to see the preformance of the package for the mixed data using Gaussian copula graphical models. Frist, by using the function `bdgraph.sim` to simulate 100 observations (n = 100) from mixed data (`type = "mixed"`) with 7 variables (p = 7) and “random” graph structure, as follows:

``` r
data.sim = bdgraph.sim( n = 100, p = 7, type = "mixed", graph = "random", vis = TRUE )
round( head( data.sim $ data, 4 ), 2 )
```

Since the generated data are mixed data, we are using  run the `bdgraph` function by choosing `method = "gcgm"`, as follows:

``` r
bdgraph.obj <- bdgraph( data = data.sim, method = "gcgm", iter = 5000 )

summary( ssgraph.obj )
```

To compare the result with the true graph, we could run

``` r
compare( data.sim, bdgraph.obj, main = c( "Target", "BDgraph" ), vis = TRUE )
```

For more examples see [Mohammadi and Wit (2019)](https://www.jstatsoft.org/article/view/v089i03).



