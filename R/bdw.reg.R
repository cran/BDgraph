## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
#     Copyright (C) 2012 - 2022  Reza Mohammadi                                |
#                                                                              |
#     This file is part of BDgraph package.                                    |
#                                                                              |
#     BDgraph is free software: you can redistribute it and/or modify it under |
#     the terms of the GNU General Public License as published by the Free     |
#     Software Foundation; see <https://cran.r-project.org/web/licenses/GPL-3>.|
#                                                                              |
#     Maintainer: Reza Mohammadi <a.mohammadi@uva.nl>                          |
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
# Bayesian parameter estimation for discrete Weibull regression
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |

bdw.reg = function( data, formula = NA, iter = 5000, burnin = NULL, 
                    dist.q = dnorm, dist.beta = dnorm, 
                    par.q = c( 0, 1 ), par.beta = c( 0, 1 ), par.pi = c( 1, 1 ), 
                    initial.q = NULL, initial.beta = NULL, initial.pi = NULL, 
                    ZI = FALSE, scale.proposal = NULL, adapt = TRUE, print = TRUE )
{
  if( is.null( burnin ) ) burnin = round( iter * 0.75 )
  
  if( !is.vector( data ) )
  {
    mf = stats::model.frame( formula, data = data )
    
    y  = stats::model.response( mf, "numeric" )
    # x  = as.matrix( mf )
    # x[ , 1 ] = 1
    x  = stats::model.matrix(formula, data = data)
  }else{
    y = data
    x = matrix( 1, ncol = 1, nrow = length( y ) )
  }
  
  n      = length( y )
  ncol_x = ncol( x )
  ind_y0 = ifelse( y == 0, 1, 0 )
  
  #setting initial parameters
  
  if( is.null( initial.q ) ) 
  {
    theta.q = rep( 0, times = ncol_x )
    q = 1 - sum( y == 0 ) / n
    if( q == 1 ) q = 0.999
    
    theta.q[ 1 ] = log( q / ( 1 - q ) )
  }else{
    theta.q = initial.q
  }
  
  if( is.null( initial.beta ) ) 
  {
    theta.beta = rep( 0, times = ncol_x )
  }else{
    theta.beta = initial.beta
  }
  
  if( is.null( scale.proposal ) ) scale.proposal = 2.38 ^ 2 / ( 2 * ncol_x )
  
  if( !ZI ){
    pii = 1
    z   = rep( 1, n )
  }else{
    pii    = ifelse( is.null( initial.pi ), sum( y != 0 ) / n + sum( y == 0 ) / ( n * 2 ), initial.pi   )
    # z = ifelse( y == 0, 0, 1 ) 
    prob_z = ifelse( y != 0, 1, pii )
    z      = stats::rbinom( n = n, size = 1, prob = prob_z )
  }
  
  fit = stats::optim( par = c( theta.q, theta.beta ),
                      fn = BDgraph::log_post_cond_dw,
                      x = x, y = y, z = z, dist.q = dist.q, par.q = par.q, dist.beta = dist.beta, par.beta = par.beta,
                      control = list( "fnscale" = -1, maxit = 10000 ), hessian = TRUE )
  
  #Initial values: initial posterior mode
  theta.q    = fit $ par[ 1 : ncol_x ]
  theta.beta = fit $ par[ ( ncol_x + 1 ) : ( 2 * ncol_x ) ]
  
  para = c( theta.q, theta.beta )
  if( ZI ) para = c( para, pii )
  
  #Setting the covariance matrix of the proposal distribution for the regression coefficients
  fisher_info_inv = BDgraph::near_positive_definite( -fit $ hessian[ 1 : ( 2 * ncol_x ), 1 : ( 2 * ncol_x ) ] )
  
  if( !isSymmetric( fisher_info_inv ) || isTRUE( class( try( solve( fisher_info_inv ), silent = TRUE ) ) == "try-error" ) )
  {
    # sigma_proposal = diag( 0.5, length( para ) )
    sigma_proposal = diag( 0.5, 2 * ncol_x )
  }else{
    fisher_info     = solve( fisher_info_inv )
    sigma_proposal = scale.proposal * fisher_info
  }
  
  sample       = matrix( nrow = iter + 1, ncol = length( para ) )
  sample[ 1, ] = para
  
  count_accept = 0
  
  for( i_mcmc in 1:iter ) 
  {
    if( ( print == TRUE ) && ( i_mcmc %% round( iter / 100 ) == 0 ) )
    {
      info = paste( round( i_mcmc / iter * 100 ), '% done, Acceptance rate = ', 
                    round( count_accept / iter * 100, 2 ),'%' )
      cat( '\r ', info, rep( ' ', 20 ) )
    }
    
    theta.q    = sample[ i_mcmc, 1 : ncol_x ]
    theta.beta = sample[ i_mcmc, ( ncol_x + 1 ) : ( 2 * ncol_x ) ]
    
    if( ZI )
    {
      pii  = sample[ i_mcmc, ncol( sample ) ]
      dens = BDgraph::ddweibull_reg( x, y, theta.q, theta.beta )
      
      prob_z0 = ind_y0 * ( 1 - pii )
      prob_z1 = dens   * pii
      prob_z  = prob_z1 / ( prob_z1 + prob_z0 )
      prob_z[ is.na( prob_z ) ] = 1
      
      #Gibbs sampling for z: sample from posterior distribution given current values of theta.q and theta.beta
      z = stats::rbinom( n = n, size = 1, prob = prob_z )
      
      #Gibbs sampling pi: sample from z
      n_z0 = sum( z == 0 )
      n_z1 = sum( z == 1 )
      ### Sample new pi from posterior distribution
      pii = stats::rbeta( n = 1, shape1 = par.pi[ 1 ] + n_z1, shape2 = par.pi[ 2 ] + n_z0 )
    }
    
    #Random Metropolis-Hasting for theta.q and theta.beta (given Z and pi)
    
    #Log-posterior at current parameters
    log_posterior_t = BDgraph::log_post_cond_dw( par = c( theta.q, theta.beta ), par.q = par.q, par.beta = par.beta,
                                                 dist.q = dist.q, dist.beta = dist.beta, x = x, y = y, z = z )
    #Adaptive MH
    if( adapt )
      if( i_mcmc >= 100 && i_mcmc %% 100 == 0 )
        sigma_proposal = BDgraph::near_positive_definite( stats::cov( sample[ ( i_mcmc - 99 ) : i_mcmc, 1 : ( 2 * ncol_x ) ] ) )
    
    if( !isSymmetric( sigma_proposal ) || isTRUE( class( try( chol.default( sigma_proposal ), silent = TRUE ) ) == "try-error" ) )
      sigma_proposal = diag( 0.5, ncol( fisher_info ) )
    
    b1 = BDgraph::rmvnorm( n = 1, mean = rep( 0, 2 * ncol_x ), sigma = sigma_proposal )
    
    #New proposed values of theta.q and theta.beta
    value_proposal = sample[ i_mcmc, 1 : ( 2 * ncol_x ) ] + as.vector( b1 ) 
    
    #Log-posterior at proposed parameters (given the same Z and pi)
    log_posterior_proposal = BDgraph::log_post_cond_dw( par = value_proposal,
                                                        par.q  = par.q, par.beta = par.beta,
                                                        dist.q = dist.q, dist.beta = dist.beta,
                                                        x = x, y = y, z = z )
    
    log_alpha = min( 0, log_posterior_proposal - log_posterior_t ) 
    
    if( log( stats::runif( 1 ) ) <= log_alpha ) 
    {
      sample[ i_mcmc + 1, 1 : ( 2 * ncol_x ) ] = value_proposal
      count_accept = count_accept + 1
    }else
      sample[ i_mcmc + 1, 1 : ( 2 * ncol_x ) ] = sample[ i_mcmc, 1 : ( 2 * ncol_x ) ]
    
    if( ZI )
      sample[ i_mcmc + 1, ncol( sample ) ] = pii
  }
  
  cat( "\n" )
  
  sample = sample[ ( burnin + 1 ) : ( iter + 1 ), ]
  
  sample.mean = apply( sample, 2, mean )
  theta.q     = sample.mean[ 1 : ncol_x ]
  theta.beta  = sample.mean[ ( ncol_x + 1 ) : ( 2 * ncol_x ) ]
  
  if( ZI )
    pii = sample.mean[ 2 * ncol_x + 1 ]
  
  q_reg    = 1 / ( 1 + exp( - x %*% theta.q ) )
  beta_reg = exp( x %*% theta.beta )
  
  if( ncol_x == 1 )
  {
    q_reg    = q_reg[ 1 ]
    beta_reg = beta_reg[ 1 ]
  }
  
  return( list( sample = sample, q.est = q_reg, beta.est = beta_reg, 
                pi.est = pii, accept.rate = count_accept / iter ) )
}
    
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
# Discrete Weibull density based on regression 
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
ddweibull_reg = function( x, y, theta.q, theta.beta )
{
    q_reg    = 1 / ( 1 + exp( - x %*% theta.q ) )
    beta_reg = exp( x %*% theta.beta )
    
    density = q_reg ^ ( y ^ beta_reg ) - q_reg ^ ( ( y + 1 ) ^ beta_reg )
    
    return( density )
}

## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
## Posterior of q and beta given Y & Z: P( X, Y | Z ) P( Z ) P( q ) P( beta )
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
log_post_cond_dw = function( par, x, y, z, dist.q, par.q, dist.beta, par.beta )
{
    ncol_x = ncol( x )
    
    theta.q    = par[ 1 : ncol_x ]   
    theta.beta = par[ ( ncol_x + 1 ) : ( 2 * ncol( x ) ) ]
    
    dens_dw_z1 = BDgraph::ddweibull_reg( x = as.matrix( x[ z == 1, ] ), y = y[ z == 1 ], 
                                theta.q = theta.q, theta.beta = theta.beta )
    
    # log_lik = sum( log( dens_dw_z1[ dens_dw_z1 != 0 ] ) ) 
    dens_dw_z1[ dens_dw_z1 == 0 ] = .Machine $ double.xmin
    log_lik = sum( log( dens_dw_z1 ) )
    
    log_prior_q = sum( dist.q( theta.q, par.q[ 1 ], par.q[ 2 ], log = TRUE ) )
    log_prior_b = sum( dist.beta( theta.beta, par.beta[ 1 ], par.beta[ 2 ], log = TRUE ) )
  
    log_post = log_lik + log_prior_q + log_prior_b 
  
    return( log_post )
}

## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
near_positive_definite = function( m )
{
    # Copyright 2003-05 Korbinian Strimmer
    # Rank, condition, and positive definiteness of a matrix; Method by Higham 1988
    
    d = ncol( m )
    eigen_m = eigen( m )
    
    eigen_vectors = eigen_m $ vectors
    eigen_values  = eigen_m $ values
    
    delta = 2 * d * max( abs( eigen_values ) ) * .Machine $ double.eps
    
    # factor two is just to make sure the resulting
    # matrix passes all numerical tests of positive definiteness
    
    tau = pmax( 0, delta - eigen_values )
    dm  = eigen_vectors %*% diag( tau, d ) %*% t( eigen_vectors )
    
    return( m + dm )
}

## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
