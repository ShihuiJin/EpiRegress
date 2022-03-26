#EpiRegress
library(MASS)
library(LaplacesDemon)
library(coda)
#setwd('~/OneDrive - National University of Singapore/reimagine/coding')
#Rcpp::sourceCpp('global_rcpp.cpp')
Rcpp::sourceCpp('~/OneDrive - National University of Singapore/reimagine/coding/global_rcpp.cpp')
#Rcpp::sourceCpp("C:/Users/shihui.j/OneDrive - National University of Singapore/reimagine/coding/global_rcpp.cpp")
logposterior_c=function(current){
  if(current$tau<=1){
    current$logposterior=-999999
    return(current)
  }
  Rt=Rt_c(rep(current$a,end), current$b, X_sd, end, length(current$b),delta)
  current$Rt=Rt
  current$loglikelihood=loglikelihood_c(start, end, current$phi, current$tau, Rt, w, I, I_imp)
  current$logposterior=current$loglikelihood+sum(dlaplace(current$b,0,1/current$lambda,log=TRUE))
  current
}

mh=function(current,old)
{
  if(current$phi<0|current$tau<=1) return(old)
  current=logposterior_c(current)
  logaccprob=current$logposterior-old$logposterior
  lu=-rexp(1)
  if(lu>logaccprob) current=old
  current
}


#X_sd: standardized covariate matrix
#I,I_imp: local and imported case counts
#delta: time delay in reporting
MCMC=function(X_sd, I, I_imp,delta,MCMCiterations,thin,var,start,end,select){
  logpos=rep(0,MCMCiterations)
  R=matrix(0,MCMCiterations,end)
  current=list(a=0,b=rep(0,length(select)),phi=0.1,tau=5,lambda=5)
  store=matrix(0,MCMCiterations,ncol(X_sd)+3)
  for (iteration in 1:(MCMCiterations*thin))
  {
    if ((iteration-1)%%1e5==0&iteration!=1) print(paste0('iteration: ',iteration-1,";  time: ",Sys.time()))
    old=logposterior_c(current)
    new=mvrnorm(1,c(current$a,current$b,current$phi, current$tau),var)
    current$a=new[1]
    current$b=new[2:(ncol(store)-2)]
    current$phi=new[ncol(store)-1]
    current$tau=new[ncol(store)]
    current=mh(current,old)
    if(iteration%%thin==0){
      store[iteration/thin,1]=current$a
      store[iteration/thin,2:(ncol(store)-2)]=current$b
      store[iteration/thin,ncol(store)-1]=current$phi
      store[iteration/thin,ncol(store)]=current$tau
      R[iteration/thin,]=current$Rt
      #logpos[iteration/thin]=current$logposterior
      logpos[iteration/thin]=current$loglikelihood
    }
  }
  list(R, store,logpos)
}
#remove(logpos,R, current, store)

EpiRegress=function(X_sd, I, I_imp,start,end,w,delta=0,select=c(1:ncol(X_sd)), draws=1e4,thin=1,raw=FALSE){
  X_sd_new=X_sd[,select]
  #current=list(a=0,b=rep(0,length(select)),phi=0.1,tau=5,lambda=5)
  var=diag(c(rep(1e-6,length(select)+2),5e-6),length(select)+3,length(select)+3)
  store=MCMC(X_sd_new, I, I_imp,delta,draws,thin,var,start,end,select)[[2]]
  burn_in=draws/5
  var=cov(store[burn_in:draws,])/300
  #var=as.matrix(read.csv('var_p.csv'))
  result=MCMC(X_sd_new, I, I_imp,delta,draws,thin,var, start,end,select)
  R=result[[1]]
  store=result[[2]]
  logpos=result[[3]]
  R2=matrix(0,end,6)
  for(t in 1:nrow(R2))
  {
    R2[t,1]=mean(R[burn_in:draws,t])
    R2[t,2:6]=quantile(R[burn_in:draws,t],c(0.025,0.25,0.5,0.75,0.975))
  }
  colnames(R2)=c('mean','95CrI_l', '1st_quantile', 'median','3rd_quantile','95CrI_u')
  R2=as.data.frame(R2)
  DIC=-4*mean(logpos)+2*loglikelihood_c(start, end, mean(store[burn_in:draws,ncol(store)-1]), mean(store[burn_in:draws,ncol(store)]), R2[,1], w, I, I_imp)
  #mean=R2[,1]
  #median=R2[,4]
  #CrI=R2[,c(2,6)]
  beta_quantile=matrix(0,ncol(store)-3,6)
  for(i in 1:nrow(beta_quantile))
  {
    beta_quantile[i,1]=mean(store[,i+1])
    beta_quantile[i,2:6]=quantile(store[,i+1],c(0.025,0.25,0.5,0.075,0.975))
  }
  colnames(beta_quantile)=c('mean','95CrI_l', '1st_quantile', 'median','3rd_quantile','95CrI_u')
  beta_quantile=as.data.frame(beta_quantile)
  if(raw==TRUE){
    return(list(Rt_estimates=R2[start:(end-1),], beta_estimates=beta_quantile, DIC=DIC))
  }else{
    return(list(Rt_estimates=R2[start:(end-1),], beta_estimates=beta_quantile, pars=store,DIC=DIC))
  }
  #remove(start,end,select,delta,MCMCiterations,w,var,result)
}
