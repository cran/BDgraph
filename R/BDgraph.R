## All functions for the "BDgraph" package
###########################################################
## function for creating matrix H from data
Hmatrix=function(D){
   p=nrow(D)
   Ti=chol(solve(D))  
   H=diag(p)
   for (i in 1:(p-1)){
      for (j in (i+1):p) H[i,j]=Ti[i,j]/Ti[j,j]
      }
   return(H)
}
## function for only obtaining elements of Psi matrix
Psi=function(A,b,H){
  p=nrow(A)
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
sample.G.Wishart=function(A,b,D){
   p=nrow(A)
   Ti=chol(solve(D)) 
   H=Hmatrix(D)
   psi=Psi(A,b,H)
   round(t(psi%*%Ti)%*%(psi%*%Ti),3)
}
#Algorithm 3.2: function for Monte Carlo approxiation for expection in normalizing constant
Exp.MC=function(A,b,H,MC.iter,p){
   f_T=vector()
   for (i in 1:MC.iter){
       psi=Psi(A,b,H)
       f_T[i]=exp(-sum((1-(diag(p)+A))*psi*psi)/2)
   }
   mean(f_T)
}
# function for computing Normalizing constans of G-Wishart distribution according to ATAY-KAYIS AND mASSAM (2005)
I_G=function(A,b,D,MC.iter=300){
   p=nrow(A)
   Ti=chol(solve(D))
   H=Hmatrix(D)
   Exp.f_T=Exp.MC(A,b,H,MC.iter,p)
   c_dT=0
   for (i in 1:p){
      c_dT=c_dT+((sum(A[i,])/2)*log(pi)+((b+2*sum(A[i,]))/2)*log(2)+
      lgamma((b+sum(A[i,]))/2)+(b+sum(A[i,])+sum(A[,i]))*log(Ti[i,i]))
   }
   c_dT=exp(c_dT)
   paste(c("Normalizing constant=",c_dT*Exp.f_T),collapse="")
}
# sampling from precistion matrix K for our BDMCMC algorithm
sampleK=function(A,b,H,Ts,p){
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
get.S=function(S,n,tol=1e-5){
  if (ncol(S)!=nrow(S)){
    n = nrow(S)
    S = n*cov(S)
  } else {
    if (sum(abs(S-t(S)))>tol){
      n = nrow(S)
      S = n*cov(S)
    }
  }
  return(list(S=S,n=n))
}
# Algorithm 2.1: BD-MCMC algorithm for low-dimentional problem (graphs with more less 8 nodes)
BDMCMC.low=function(S,n=NULL,iter=5000,burn=2000,distance=1,lambda_b=1,b=3,D=NULL,A=NULL,MC.iter=10){
  Sn = get.S(S,n)
  if (is.null(Sn$n)&is.null(n)){
    stop("If you provide the covariance matrix, you have to specify the number of observations")
  }
  S = Sn$S
  n = Sn$n
  p = ncol(S) 
  if (is.null(A)) A = rbind(cbind(0,diag(rep(1,p-1))),0) 
  if (is.null(D)) D = diag(p)
  bs=b+n; Ds=D+S
  Dsinv=solve(Ds)
  Ts=chol(Dsinv)
  H=diag(p)
  Hs=Hmatrix(Ds)
  K=sampleK(A,bs,Hs,Ts,p)
  Ksample=Asample=rates_sample=allA=list()
  lambda.s=vector() ## waiting time in every state
  for (g in 1:iter){
    rates=0*K
    for (i in 1:(p-1)){
       for (j in (i+1):p){
          if (A[i,j]==0) rates[i,j]=lambda_b
          if (A[i,j]==1){
             sigma=1
             mu=bs*Dsinv[i,j]
             k_xi=rnorm(1,mu,sigma)
             b_xi=dnorm(k_xi,mu,sigma)
             nustar=sum(A[i,])
             Epsi=Exp.MC(A,b,H,MC.iter,p)
             Aminus=A
             Aminus[i,j]=0
             Epsiminus=Exp.MC(Aminus,b,H,MC.iter,p)
             Kminus=sampleK(Aminus,bs,Hs,Ts,p)
             rates[i,j]=2*sqrt(pi)*lambda_b*(b_xi)*
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
BDMCMC.high=function(S,n=NULL,iter=5000,burn=2000,distance=1,lambda_b=1,b=3,D=NULL,A=NULL){
  Sn = get.S(S,n)
  if (is.null(Sn$n)&is.null(n)){
    stop("If you provide the covariance matrix, you have to specify the number of observations")
  }
  S = Sn$S
  n = Sn$n
  p = ncol(S) 
  if (is.null(A)) A = rbind(cbind(0,diag(rep(1,p-1))),0) 
  if (is.null(D)) D = diag(p)
  bs=b+n; Ds=D+S
  Ts=chol(solve(Ds))
  H=diag(p)
  Hs=Hmatrix(Ds)
  K=sampleK(A,bs,Hs,Ts,p)
  Ksample=Asample=allA=list()
  lambda.s=vector() ## waiting time for every state
  cont=0
  for (g in 1:iter){
    rates=0*K
    for (i in 1:(p-1)){
       for (j in (i+1):p){
          if (A[i,j]==0) rates[i,j]=lambda_b
          if (A[i,j]==1){
             Aminus=A
             Aminus[i,j]=0
             Kminus=sampleK(Aminus,bs,Hs,Ts,p)
             rates[i,j]=lambda_b*exp((n/2)*(log(abs(det(Kminus)))-
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
BDMCMC=function(S,n=NULL,iter=5000,burn=2000,distance=1,lambda_b=1,b=3,D=NULL,A=NULL,MC.iter=10){
  p = ncol(S)
  if (p<8){
    return(BDMCMC.low(S,n,iter,burn,distance,lambda_b,b,D,A,MC.iter))
  } else {
    return(BDMCMC.high(S,n,iter,burn,distance,lambda_b,b,D,A))
  }
}
# function for comuting probability of all links in graph
Phat=function(output){
   As=output$Asample
   lambda.s=output$lambda.s
   p=nrow(As[[1]])
   Phat=0*As[[1]]
   for (i in 1:(p-1)){
      for (j in (i+1):p){
         for (k in 1:length(As)){
             Phat[i,j]=Phat[i,j]+As[[k]][i,j]/lambda.s[k]
         }
         Phat[i,j]=Phat[i,j]/(sum(1/lambda.s))
      }
   }
   print("Posterior edge inclusion probabilities")
   return(Phat)
}
# function for computing the probability of one especial graph
prob.graph=function(A,output){
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
plot_links=function(output,min=0,max=(nrow(output$As[[1]]))*(nrow(output$As[[1]])-1)/2){
   p=nrow(output$As[[1]])
   As=output$Asample
   lambda.s=output$lambda.s
   numinator=c(rep(0,max-min+1))
   for (i in 1:length(As)){
      for (j in min:max){
          if (sum(As[[i]])==j) {numinator[j-min+1]=numinator[j-min+1]+1/lambda.s[i]}
      }
   }
   plot(x=min:max,y=numinator/sum(1/lambda.s),type="h",main="",ylab="Pr(number of links in graph|data)",xlab="number of links in graph")
}
# function for checking the convergency of the BDMCMC algorithm
plot_cumulative=function(output,distance=1){
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
prob.allgraphs=function(output){
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
selection.result=function (output,g=2){
  output.allgraphs=prob.allgraphs(output)
  list.A=output.allgraphs$list.A
  prob.A=output.allgraphs$prob.A
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



