# **BDgraph** 

![](https://www.r-pkg.org/badges/version/BDFgraph) ![](https://www.r-pkg.org/badges/last-release/BDgraph) ![](https://cranlogs.r-pkg.org/badges/BDgraph) 
![](https://cranlogs.r-pkg.org/badges/grand-total/BDgraph) 

## **BDgraph** NEWS

### **BDgraph** Version 2.63

* In the functions `bdgraph()`, `bdgraph.mpl()`, and `bdgraph.ts()`, option `print` is removed.
* In the function `bdgraph.sim()` option `type="dw"` is added.
* In the function `bdgraph.sim()` option `type="discrete"` is changed to `type="count"`.

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

* In the function `compare()` option `est4` is added.

### **BDgraph** Version 2.56

* In the functions `rgwish()`, `rgcwish()` and `gnorm()`, option `adj.g` is changed to `adj`.
	
* In the function `compare()`, option `sim.obj`      is changed to `target`.
* In the function `compare()`, option `bdgraph.obj`  is changed to `est`.
* In the function `compare()`, option `bdgraph.obj2` is changed to `est2`.
* In the function `compare()`, option `bdgraph.obj3` is changed to `est3`.

* In the function `plotroc()`, option `sim.obj`      is changed to `target`.
* In the function `plotroc()`, option `bdgraph.obj`  is changed to `est`.
* In the function `plotroc()`, option `bdgraph.obj2` is changed to `est2`.
* In the function `plotroc()`, option `bdgraph.obj3` is changed to `est3`.
* In the function `plotroc()`, option `bdgraph.obj4` is changed to `est4`.
	
* Bug fixed for `bdgraph()` function with option `method = gcgm`.
* Functions `summary.bdgraph`, `plot.bdgraph`, and `print.bdgraph` are modified.

### **BDgraph** Version 2.53

* Files `configure` and `configure.ac` are removed and files `Makevars` and `Makevars.win` are modified accordingly.

* In the functions `bdgraph()`, `bdgraph.mpl()`, and `bdgraph.ts()` option `save.all` is changed to `save`.
* In the functions `bdgraph()` and `bdgraph.mpl()`, option `multi.update` is changed to `jump`.
* In the functions `bdgraph()` and `bdgraph.ts()`, option `prior.df` is changed to `df.prior`.
	
* In the function `pgraph()`, option `adj_g` is changed to `adj`.

### **BDgraph** Version 2.52

* Option `threshold` is added to the function `rgwish()`.
* Option `threshold` is added to the function `bdgraph()`.
* Option `not.cont`  is added to the function `bdgraph()`.
* In the function `compare()`, option `colnames` is changed to `main`.

* In functions `bdgraph()` and `bdgraph.mpl()`, option `g.space` is removed. Instead user can use option `g.prior`.

* Bug fixed for `bdgraph.ts()` function.
* Data set `churn` is added.

### **BDgraph** Version 2.51

* Bug fixed for functions `rgwish()` and `rcgwish()`.
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
   

