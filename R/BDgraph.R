## All functions for the "BDgraph" package
###########################################################
## function for creating matrix H from data
Hmatrix=function(Ti,p)
{
   H=diag(p)
   for (i in 1:(p-1)){
      for (j in (i+1):p) H[i,j]=Ti[i,j]/Ti[j,j]
      }
   return(H)
}
## function for only obtaining elements of Psi matrix
Psi=function(A,b,H,p)
{
  psi=0*A
  for (i in 1:p){
    for (j in i:p){
        if (i==j) psi[i,j]=sqrt(rchisq(1,b+sum(A[i,])))
        if (A[i,j]==1) psi[i,j]=rnorm(1)
        if (A[i,j]==0 & i!=j){
            psi[i,j]=-sum(psi[i,c(i:(j-1))]*H[c(i:(j-1)),j])
            if (i>1){
               for (r in 1:(i-1)){
                 psi[i,j]=psi[i,j]-((sum(psi[r,c(r:i)]*H[c(r:i),i]))*
                 (sum(psi[r,c(r:j)]*H[c(r:j),j])))/(psi[i,i])
               }
            }
        }
    }
  }
  psi
}
# Sampling from G-Wishart distribution; algorithm 2 of our paper.
sample.gwishart=function(A,b,D,round=3)
{
   if (b<=2){
      stop("parameter 'b' in G-Wishart distribution has to be more than 2")
   }
   p=nrow(A)
   Ti=chol(solve(D)) 
   H=Hmatrix(Ti,p)
   psi=Psi(A,b,H,p)
   round(t(psi%*%Ti)%*%(psi%*%Ti),round)
}
#Algorithm 3.2: function for Monte Carlo approxiation for expection in normalizing constant
Exp.MC=function(A,b,H,MCiter,p)
{
   f_T=vector()
   for (i in 1:MCiter){
       psi=Psi(A,b,H,p)
       f_T[i]=exp(-sum((1-(diag(p)+A))*psi*psi)/2)
   }
   mean(f_T)
}
# function for computing Normalizing constans of G-Wishart distribution according to ATAY-KAYIS AND mASSAM (2005)
I.g=function(A,b,D,MCiter=500)
{
   if (b<=2){
      stop("parameter 'b' in G-Wishart distribution has to be more than 2")
   }
   p=nrow(A)
   Ti=chol(solve(D))
   H=Hmatrix(Ti,p)
   Exp.f_T=Exp.MC(A,b,H,MCiter,p)
   c_dT=0
   for (i in 1:p){
      c_dT=c_dT+((sum(A[i,])/2)*log(pi)+((b+2*sum(A[i,]))/2)*log(2)+
      lgamma((b+sum(A[i,]))/2)+(b+sum(A[i,])+sum(A[,i]))*log(Ti[i,i]))
   }
   c_dT=exp(c_dT)
   paste(c("Normalizing constant=",c_dT*Exp.f_T),collapse="")
}
# sampling from precistion matrix K for our BDMCMC algorithm
sampleK=function(A,b,H,Ts,p)
{
   psi=0*A
   for (i in 1:p){
      for (j in i:p){
        if (i==j) psi[i,j]=sqrt(rchisq(1,b+sum(A[i,])))
        if (A[i,j]==1) psi[i,j]=rnorm(1)
        if (A[i,j]==0 & i!=j){
            psi[i,j]=-sum(psi[i,c(i:(j-1))]*H[c(i:(j-1)),j])
            if (i>1){
               for (r in 1:(i-1)){
               psi[i,j]=psi[i,j]-((sum(psi[r,c(r:i)]*H[c(r:i),i]))*
               (sum(psi[r,c(r:j)]*H[c(r:j),j])))/(psi[i,i])
               }
            }
        }
      }
   }
   t(psi%*%Ts)%*%(psi%*%Ts)
}
#auxiliary function to get covariance matrix
get.S=function(data,n,tol=1e-5,meanzero)
{
  if (ncol(data)!=nrow(data)){
     n = nrow(data)
	 if (meanzero==TRUE) S = t(data)%*%data
     if (meanzero==FALSE) S = n*cov(data)
  } else {
     if (sum(abs(data-t(data)))>tol){
        n = nrow(data)
	    if (meanzero==TRUE) S = t(data)%*%data
        if (meanzero==FALSE) S = n*cov(data)
    }
  }
  return(list(S=S,n=n))
}
# Algorithm 2.1: BD-MCMC algorithm for low-dimentional problem (graphs with more less 8 nodes)
bdmcmc.low=function(data,n=NULL,meanzero=FALSE,iter=5000,burn=2000,distance=1,
gamma.b=1,Gprior="Uniform",b=3,D=NULL,A=NULL,MCiter=10,print=FALSE)
{
  if (iter<=burn){
    stop("Number of iterations have to be more than the number of burn-in iterations")
  }
  if (gamma.b<=0){
    stop("birth rate 'gamma.b' has to be positive value")
  }
  Sn = get.S(data,n,meanzero=meanzero)
  if (is.null(Sn$n)&is.null(n)){
    stop("You have to specify the number of observations 'n'")
  }
  S = Sn$S
  n = Sn$n
  p = ncol(S) 
  if (is.null(A)) {
     A = 0 * S
     A[upper.tri(A)] = 1
     } 
  if (is.null(D)) D = diag(p)
  bs=b+n; Ds=D+S
  Dsinv=solve(Ds)
  Ts=chol(Dsinv)
  H=diag(p)
  Hs=Hmatrix(Ts,p)
  K=sampleK(A,bs,Hs,Ts,p)
  Ksample=Asample=rates_sample=allA=list()
  lambda.s=vector() ## waiting time in every state
  for (g in 1:iter){
    if (print==T){
	  cat(paste("time =", format(Sys.time(), "%X")),
      paste(c("Sum.links = ",sum(A)),collapse=""), fill = TRUE,
      labels = paste("{",paste(c("iter=",g),collapse=""),"}:",sep=""))
	}
    rates=0*K
    for (i in 1:(p-1)){
       for (j in (i+1):p){
          if (A[i,j]==0) rates[i,j]=gamma.b
          if (A[i,j]==1){
             sigma=1
             mu=bs*Dsinv[i,j]
             k_xi=rnorm(1,mu,sigma)
             b_xi=dnorm(k_xi,mu,sigma)
             nustar=sum(A[i,])
             Epsi=Exp.MC(A,b,H,MCiter,p)
             Aminus=A
             Aminus[i,j]=0
             Epsiminus=Exp.MC(Aminus,b,H,MCiter,p)
             Kminus=sampleK(Aminus,bs,Hs,Ts,p)
             if (sum(A)==0 & pr==0) pr=1
             rates[i,j]=(((sum(A))^pr)*((gamma.b)^(1-pr)))*2*sqrt(pi)*gamma.b*(b_xi)*
             exp(lgamma((bs+nustar)/2)-lgamma((bs+nustar-1)/2)+
             log(Epsi)-log(Epsiminus)+((bs-2)/2)*(log(abs(det(Kminus)))-log(abs(det(K))))-
             (sum(diag(Ds%*%(Kminus-K))))/2 )
             if (is.infinite(rates[i,j])==TRUE) rates[i,j]=gamma(170)
          }
       }
    }
    allA[[g]]=A
    if (g > burn){
        Ksample[[g-burn]]=K
        Asample[[g-burn]]=A
        lambda.s[g-burn]=sum(rates)
    }
    pmatrix <- melt(rates)
    rows <- which(rmultinom(1,1,pmatrix$value)==1, arr.ind=TRUE)[,'row']
    indices <- pmatrix[rows, c('X1','X2')]
    A[indices$X1,indices$X2]=A[indices$X1,indices$X2]+(-1)^(A[indices$X1,indices$X2])
    K=sampleK(A,bs,Hs,Ts,p)
  }
  return(list(Ksample=Ksample,Asample=Asample,lambda.s=lambda.s,allA=allA))
}
## Algorithm 2.1: BD-MCMC algorithm for high-dimentional problem (graphs with more than 8 nodes)
bdmcmc.high=function(data,n=NULL,meanzero=FALSE,iter=5000,burn=2000,
distance=1,gamma.b=1,Gprior="Uniform",b=3,D=NULL,A=NULL,print=FALSE)
{
  if (iter<=burn){
    stop("Number of iterations have to be more than the number of burn-in iterations")
  }
  if (gamma.b<=0){
    stop("birth rate 'gamma.b' has to be positive value")
  }
  id <- pmatch(Gprior,c("Uniform","Poisson"))[1]
  if(!is.na(id)){
     if(id==1) pr=0
     if(id==2) pr=1
    }
  Sn = get.S(data,n,meanzero=meanzero)
  if (is.null(Sn$n)&is.null(n)){
    stop("If you provide the covariance matrix, you have to specify the number of observations")
  }
  S = Sn$S
  n = Sn$n
  p = ncol(S) 
  if (is.null(A)) {
     A = 0 * S
     A[upper.tri(A)] = 1
     } 
  if (is.null(D)) D = diag(p)
  bs=b+n; Ds=D+S
  Ts=chol(solve(Ds))
  H=diag(p)
  Hs=Hmatrix(Ts,p)
  K=sampleK(A,bs,Hs,Ts,p)
  Ksample=Asample=allA=list()
  lambda.s=vector() ## waiting time for every state
  cont=0
  for (g in 1:iter){
    if (print==T){
	  cat(paste("time =", format(Sys.time(), "%X")),
      paste(c("Sum.links = ",sum(A)),collapse=""), fill = TRUE,
      labels = paste("{",paste(c("iter=",g),collapse=""),"}:",sep=""))
	}
    rates=0*K
    for (i in 1:(p-1)){
       for (j in (i+1):p){
          if (A[i,j]==0) rates[i,j]=gamma.b
          if (A[i,j]==1){
             Aminus=A
             Aminus[i,j]=0
             Kminus=sampleK(Aminus,bs,Hs,Ts,p)
             if (sum(A)==0 & pr==0) pr=1
             rates[i,j]=(((sum(A))^pr)*((gamma.b)^(1-pr)))*exp((n/2)*(log(abs(det(Kminus)))-
             log(abs(det(K))))+(sum(diag(S%*%(K-Kminus))))/2 )
             if (is.infinite(rates[i,j])==TRUE) rates[i,j]=gamma(170)
          }
       }
    }
    allA[[g]]=A
    if (g > burn && g%%distance==0){
        cont=cont+1
        Ksample[[cont]]=K
        Asample[[cont]]=A
        lambda.s[cont]=sum(rates)
    }
    pmatrix <- melt(rates)
    rows <- which(rmultinom(1,1,pmatrix$value)==1, arr.ind=TRUE)[,'row']
    indices <- pmatrix[rows, c('X1','X2')]
    A[indices$X1,indices$X2]=A[indices$X1,indices$X2]+(-1)^(A[indices$X1,indices$X2])
    K=sampleK(A,bs,Hs,Ts,p)
  }
  return(list(Ksample=Ksample,Asample=Asample,lambda.s=lambda.s,allA=allA))
}
## Main function: BDMCMC algorithm for selecting the best graphical model
bdmcmc=function(data,n=NULL,meanzero=FALSE,iter=5000,burn=2000,distance=1,
gamma.b=1,Gprior="Uniform",b=3,D=NULL,A=NULL,MCiter=10,print=FALSE)
{
  p = ncol(data)
  if (p<8){
    return(bdmcmc.low(data,n,meanzero,iter,burn,distance,gamma.b,Gprior,b,D,A,MCiter,print))
  } else {
    return(bdmcmc.high(data,n,meanzero,iter,burn,distance,gamma.b,Gprior,b,D,A,print))
  }
}
# function for comuting probability of all links in graph
phat=function(output,round=3)
{
   As=output$Asample
   lambda.s=output$lambda.s
   p=nrow(As[[1]])
   phat=0*As[[1]]
   for (i in 1:(p-1)){
      for (j in (i+1):p){
         for (k in 1:length(As)){
             phat[i,j]=phat[i,j]+As[[k]][i,j]/lambda.s[k]
         }
         phat[i,j]=phat[i,j]/(sum(1/lambda.s))
      }
   }
   print("Posterior edge inclusion probabilities")
   return(round(phat,round))
}
# function for computing the probability of one especial graph
prob.g=function(A,output)
{
   As=output$Asample
   Ks=output$Ksample
   lambda.s=output$lambda.s
   tlambda=vector()
   for (i in 1:length(As)){
       if (identical(A,As[[i]])==TRUE) tlambda=c(tlambda,lambda.s[i])
   }
   sum(1/tlambda)/sum(1/lambda.s)
}
# plot for probability of graphs according to number of their links
plotLinks=function(output,xlim=c(0,(nrow(output$As[[1]]))*(nrow(output$As[[1]])-1)/2))
{
   p=nrow(output$As[[1]])
   As=output$Asample
   lambda.s=output$lambda.s
   numinator=c(rep(0,xlim[2]-xlim[1]+1))
   for (i in 1:length(As)){
      for (j in xlim[1]:xlim[2]){
          if (sum(As[[i]])==j) {numinator[j-xlim[1]+1]=numinator[j-xlim[1]+1]+1/lambda.s[i]}
      }
   }
   plot(x=xlim[1]:xlim[2],y=numinator/sum(1/lambda.s),type="h",main="",ylab="Pr(number of links in graph|data)",xlab="number of links in graph")
}
# function for checking the convergency of the BDMCMC algorithm
plotConvergency=function(output,distance=1)
{
  allA=output$allA
  p=nrow(allA[[1]])
  allA.new=list()
  for (i in 1:length(allA)){
     g=i*distance
     if (g>length(allA)) break
     allA.new[[i]]=allA[[g]]
  }   
  length.allA=length(allA.new)
  ff=matrix(0,p*(p-1)/2,length.allA)
  for (g in 1:length.allA){
     for (k in 1:g){
        con=0
        for (i in 1:(p-1)){
           for (j in (i+1):p){
              con=con+1
              ff[con,g]=ff[con,g]+allA.new[[k]][i,j]/g
           }
        }
     }
  }
  matplot(x=distance*(1:length.allA),y=t(ff),type="l",lty=1,col=1,xlab="number of iterations",ylab="cumulative occupancy fractions for each links")
}
# Program for computing the probability of all possible graphical models by using the result of BDMCMC algorithm
prob.allg=function(output)
{
   As.lim=output$Asample
   lambda.s=output$lambda.s
   lambda.lim=lambda.s
   list.A=list.lambda=list()
   fre=vector()
   i=0
   while (length(As.lim)>1){
      i=i+1
      list.A[[i]]=As.lim[[1]]
      list.lambda[[i]]=vector()
      free=0
      As.lim2=As.lim
      lambda.lim2=lambda.lim
      for (g in 1:length(As.lim)){
         if (identical(list.A[[i]],As.lim[[g]])==TRUE){
            list.lambda[[i]]=c(list.lambda[[i]],lambda.lim[g])
            As.lim2=As.lim2[-(g-free)]
		        lambda.lim2=lambda.lim2[-(g-free)]
		        free=free+1
         }
      }
      As.lim=As.lim2
      lambda.lim=lambda.lim2
      fre[i]=free
   }
   prob.A=vector()
   for (i in 1:length(list.A)){
       prob.A[i]=sum(1/list.lambda[[i]])/sum(1/lambda.s)
   }
   return(list(list.A=list.A,prob.A=prob.A))
}                                              
# function summary of the result
select.g=function (output,g=2)
{
  output.allg=prob.allg(output)
  list.A=output.allg$list.A
  prob.A=output.allg$prob.A
  for (i in 1:g){
    dev.new()
    G=network(list.A[[which(prob.A==sort(prob.A,decreasing=T)[i])]],directed=F)
    if (i==1){
         main = "BEST GRAPH: graph with highest probability"
    } else {
      if (i==2){
        main = "Graph with 2nd highest probability"
      } else {
        if (i==3){
          main = "Graph with 3rd highest probability"
        } else {
          main = paste(c("Graph with ", i, "th highest probability"),collapse="")
        }
      }
    }
    plot.network(G, label=network.vertex.names(G), main=main,
    sub=paste(c("Posterior probability of graph=",round(sort(prob.A,decreasing=TRUE)[i],3)),collapse=""))
  }
}



