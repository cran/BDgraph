<!-- README.md is generated from NEWS.Rmd. Please edit that file -->

## **BDgraph** NEWS <a href='https://CRAN.R-project.org/package=BDgraph'><img src='man/figures/logo.png' align="right" height="170" /></a>

### **BDgraph** Version 2.73

* Function `mse()` is added to the package.

### **BDgraph** Version 2.72

* Function `predict.bdgraph()` is added to the package.
* Function `posterior.predict()` is added to the package.

* In functions `bdgraph()`, `bdgraph.mpl()`, and `bdgraph.dw()`, option `g.prior` is changed from `0.5` to `0.2`.

### **BDgraph** Version 2.71

* Function `auc()` is added to the package.
* In function `plotroc()`, option `calibrate` is added.

* Data set `churn` is removed.

### **BDgraph** Version 2.70

* Bug fixed for `plotroc()` function.

* In function `plotroc()`, options `color`, `typeline`, `size`, `xlab`, and `ylab`, and some more options are added.

### **BDgraph** Version 2.69

* In function `bdgraph.sim()`, options `type="t"` and `type="alternative-t"` are added.
* In function `bdgraph.sim()`, option `nu` is added for the cases `type="t"` and `type="alternative-t"`.

* Bug fixed for `compare()` function.
* In function `compare()`, option `est` is replaced.

* In function `plotroc()`, option `est` is replaced.

### **BDgraph** Version 2.68

* In functions `bdgraph()`, `bdgraph.mpl()`, and `bdgraph.dw()`, option `verbose` is added.
* In function `plotcoda()`, option `verbose` is added.

* A new file is added to the vignette of package.
* Bug fixed for `bdgraph()` function with option `method = gcgm` and `algorithm = rj-dmh`.

### **BDgraph** Version 2.65

* Function `bdgraph.dw()` is added to the package.
* Function `bdw.reg()`    is added to the package.

* In functions `graph.sim()` and `bdgraph.sim()`, option `graph = "smallworld"` is added.
* In functions `plot.graph()` and `plot.bdgraph()`, options `layout`, `vertex.size`, `vertex.color`, and `vertex.label.dist` are added.

* In function `bdgraph.sim()`, options `range.mu` and `range.dispersion` are added.

### **BDgraph** Version 2.64

* Function `roc()`      is added to the package.
* Function `conf.mat()` is added to the package.
* Function `conf.mat.plot()` is added to the package.
* Function `sparsity()` is added to the package.

* Functions `ddweibull()`, `pdweibull()`, `qdweibull()`, `rdweibull()`  are added to the package.

### **BDgraph** Version 2.63

* In functions `bdgraph()`, `bdgraph.mpl()`, and `bdgraph.ts()`, option `print` is removed.
* In function `bdgraph.sim()`, option `type="dw"` is added.
* In function `bdgraph.sim()`, option `type="discrete"` is changed to `type="count"`.

* Bug fixed for `gnorm()` function.

* `README` is added to the package.

### **BDgraph** Version 2.62

* Bug fixed related to the class function; class(.) == *

### **BDgraph** Version 2.61

* Bug fixed related to the C stack by Fortran functions.
* Function `bf()`         is added to the package.
* Function `adj2link()`   is added to the package.
* Function `link2adj()`   is added to the package.
* Function `bdgraph.ts()` is removed.

### **BDgraph** Version 2.60

* Function `precision()`  is added to the package.
* Function `covariance()` is added to the package.

### **BDgraph** Version 2.58

* Bug fixed for `compare()` function.

* Function `get_graph()` is added to the package and used inside function `compare()`.
* Functions `get_g_prior()`, `get_g_start()`, `get_S_n_p()`, `get_K_start()`, and `get_cores()` are added to the package and used inside function `bdgraph()`.

### **BDgraph** Version 2.57

* In function `compare()`, option `est4` is added.

### **BDgraph** Version 2.56

* In functions `rgwish()`, `rgcwish()` and `gnorm()`, option `adj.g` is changed to `adj`.
	
* In function `compare()`, option `sim.obj`      is changed to `target`.
* In function `compare()`, option `bdgraph.obj`  is changed to `est`.
* In function `compare()`, option `bdgraph.obj2` is changed to `est2`.
* In function `compare()`, option `bdgraph.obj3` is changed to `est3`.

* In function `plotroc()`, option `sim.obj`      is changed to `target`.
* In function `plotroc()`, option `bdgraph.obj`  is changed to `est`.
* In function `plotroc()`, option `bdgraph.obj2` is changed to `est2`.
* In function `plotroc()`, option `bdgraph.obj3` is changed to `est3`.
* In function `plotroc()`, option `bdgraph.obj4` is changed to `est4`.
	
* Bug fixed for `bdgraph()` function with option `method = gcgm`.
* Functions `summary.bdgraph`, `plot.bdgraph`, and `print.bdgraph` are modified.

### **BDgraph** Version 2.53

* Files `configure` and `configure.ac` are removed and files `Makevars` and `Makevars.win` are modified accordingly.

* In functions `bdgraph()`, `bdgraph.mpl()`, and `bdgraph.ts()`, option `save.all` is changed to `save`.
* In functions `bdgraph()` and `bdgraph.mpl()`, option `multi.update` is changed to `jump`.
* In functions `bdgraph()` and `bdgraph.ts()`, option `prior.df` is changed to `df.prior`.
	
* In function `pgraph()`, option `adj_g` is changed to `adj`.

### **BDgraph** Version 2.52

* Option `threshold` is added to function `rgwish()`.
* Option `threshold` is added to function `bdgraph()`.
* Option `not.cont`  is added to function `bdgraph()`.
* In function `compare()`, option `colnames` is changed to `main`.

* In functions `bdgraph()` and `bdgraph.mpl()`, option `g.space` is removed. Instead user can use option `g.prior`.

* Bug fixed for `bdgraph.ts()` function.
* Data set `churn` is added.

### **BDgraph** Version 2.51

* Bug fixed In functions `rgwish()` and `rcgwish()`.
* vignette is added to the package.

### **BDgraph** Version 2.47

* Function `graph.sim`  is added to the package.
* Function `plot.graph` is added to the package.

### **BDgraph** Version 2.46

* Function `rmvnorm` is added to the package.
	 
### **BDgraph** Version 2.44

* Functions `bdgraph.ts`, `rcwish`, and `rcgwish` are added to the package.

### **BDgraph** Version 2.42

* Function `local_rates_ggm_mpl` in cpp is added for bdgraph.mpl function.
* `configure` and `configure.ac` are added to the package for Makevars file.
* Functions `bdgraph.ts`, `rcwish`, and `rcgwish` are removed from the package.
	 
### **BDgraph** Version 2.41

* Option `dgm-binary` is added to function `bdgraph()` which is for running Hill-Climbing algorihtm for binary data.
* `configure` and `configure.ac` are added to configure checks of C++ codes.
	
### **BDgraph** Version 2.38

* Function `bdgraph.mpl()` is added to the package, which is based on Marginal Pseudo-Likelihood estimation.
* In the algorithm `gcgm`, step `copula` are implimented in parallel using OpenMP in C++. 
* Function `bdgraph.ts` is implimented in parallel using OpenMP in C++.
* New reference related to the ratio of normalizing constant is added to manual.	

### **BDgraph** Version 2.40

* Option `g.prior` is added to function `bdgraph()`, for prior distribution of the graph.
* Option `cores`   is added to function `bdgraph()`, for determining the number of cores to use for parallel execution.
* Function `transfer()` is added to the package.
   
### **BDgraph** Version 2.36

* The BDMCMC algorithms are implimented in parallel using OpenMP in C++. 

### **BDgraph** Version 2.28

* The Title in Description is changed.
* Functions `bdgraph.ts()`, `rgcwish()`, and `rcwish()` are added to the package.
	
### **BDgraph** Version 2.24

* Function `phat()`   is changed to `plinks()`.
* Function `prob()`   is changed to `pgraph()`.
* Function `log_Ig()` is changed to `gnorm()`.

### **BDgraph** Version 2.23

* Function `I.g()` is chenged to `log_Ig()` and it is implemented in C++.

### **BDgraph** Version 2.20

* Reversible jump MCMC algorithm is added to the `bdgraph()` fonction.

### **BDgraph** Version 2.19

* The Title in Description is changed.
* Function `I.g` is added.
   

