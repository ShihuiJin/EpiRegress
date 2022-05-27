#rjags version of EpiRegress

#input
#w: serial interval distribution
#X: the raw covariate matrix starting from at least |w| days before the start time
#I: incidence curve for local case counts
#I_imp: incidence curve for imported case counts
#start, end: start and end time for Rt's to be estimated, the last day ("end" day)'s Rt estimate is obtained through extrapolate rather than MCMC inference
#select: list of variables selected to include in the model

# #e.g. 
# start=2; end=nrow(X)
# select=1:ncol(X)
# #serial interval distribution
# w=dlnorm(seq(1,50),1.132,0.742)
# w=w/sum(w)

EpiRegress=function(X, w, start, end, select=c(1:ncol(X)), burn_in=50000, draw=5000, thin=10){
  #preparation of the data used
  estim_period=max(start-length(w),0): end
  factor=rep(0,ncol(X))
  for(i in 1:length(factor))
  {
    factor[i]=sd(X[estim_period,i])
  }
  #normalizing X covariates, those with constant values are removed from inference
  X_sd=as.matrix(do.call('cbind', lapply(which(factor!=0), function(i) (X[,i]-mean(X[estim_period, i]))/sd(X[estim_period, i]))))
  select1=c(1:ncol(X_sd))[which(factor!=0)%in%select]
  
  data=list(X=X_sd[,select1], I=I, I.imp=I_imp, inv.w=rev(w), n=length(select1), m=length(w), start=start, end=end, lambda=5)
  init = list(a=0, b = rep(0,length(select1))) 
  setwd(paste0(address,'EpiRegress_code_0527' ))
  library(rjags)
  jagmod = jags.model("MCMCscript_EpiRegress.txt", inits = init, data = data, n.chains = 1)
  update(jagmod, n.iter = burn_in, progress.bar = "text")
  post = coda.samples(jagmod, c("a",'b','phi','tau','R','mu'), n.iter = draw*thin, progress.bar = "text", thin = thin)
  mcmc_output = as.data.frame(as.matrix(post))
  
  #analyzing the output values
  R2=as.matrix(do.call('rbind',lapply(start:end, function(t) summary_stat(mcmc_output[,paste0('R[',t,']')]))))
  beta_summary=as.matrix(do.call('rbind',lapply(1:length(select), function(i) summary_stat(mcmc_output[,paste0('b[',i,']')]))))
  colnames(R2)[1]=colnames(beta_summary)[1]='mean'
  case_interval=as.matrix(do.call('rbind',lapply(1:(end-start+1), function(t) case_interval_t(t))))
  mean_case_summary=as.matrix(do.call('rbind',lapply(1:(end-start+1), function(t) summary_stat(mcmc_output[,paste0('mu[',t,']')]))))
  #dic
  loglik_raw=as.matrix(do.call('cbind',lapply(1:(end-start+1), function(t) case_prob(t))))
  loglik_sum_raw=rowSums(loglik_raw)
  dic= -2*mean(loglik_sum_raw)+1/2*mean(var(-2*loglik_sum_raw))
  
  #return beta and Rt estimates, expected case counts and 95% prediction interval of case counts)
  list(beta=beta_summary,Rt=R2,dic=dic, mean_case=mean_case_summary[,1], case_interval=case_interval)
}

#functions for analysis
{
  summary_stat=function(data){
    c(mean(data), quantile(data, c(0.025,0.25,0.5,0.75,0.975)))
  }
  case_interval_t=function(t){
    sim=as.vector(do.call('cbind', lapply(1:nrow(mcmc_output), function(i) rnbinom(20, mcmc_output[i,paste0('mu[',t,']')]/(mcmc_output[i,'tau']-1),1/mcmc_output[i,'tau']))))
    quantile(sim, c(0.025,0.975))
  }
  case_prob=function(t){
    as.matrix(do.call('rbind', lapply(1:nrow(mcmc_output), function(i) dnbinom(I[t+start-1], mcmc_output[i,paste0('mu[',t,']')]/(mcmc_output[i,'tau']-1),1/mcmc_output[i,'tau'],log=TRUE))))
  }
}


plot_Rt=function(R2){
  n=nrow(R2)
  library(grid)
  library(scales)
  ym=ceiling(R2)
  png('Rt_plot.png',height=8,width=8,units='cm', res=300, pointsize=10) 
  pushViewport(plotViewport(c(4,4,1,1), xscale=c(0,n),yscale=c(0,ym), clip=TRUE))
  grid.polygon(x=c(1:n,rev(1:n)), y=c(R2[,2],rev(R2[,6])), default.units='native', gp = gpar(fill = scales::alpha("rosybrown3",0.3),col=scales::alpha("rosybrown3",0.05) , alpha = 1))
  grid.lines(1:n, R2[,4], default.units='native', gp=gpar(col=scales::alpha('salmon3',0.6), cex=3))
  grid.lines(0:n, c(1,1), default.units='native', gp=gpar(col=scales::alpha('palevioletred4',0.3), cex=3))
  popViewport()
  pushViewport(plotViewport(c(4,4,1,1), xscale=c(0,n),yscale=c(0,ym)))
  grid.yaxis()
  xtk=seq(0,n,30);xlb=1:(length(xtk)-1)
  grid.xaxis(at=xtk,label=FALSE)
  grid.text(xlb,x=unit(0.5*(xtk[-1]+xtk[-length(xtk)]),'native'),y=unit(-1,'lines'))
  grid.text('Month',y=unit(-2.5,'lines'))
  grid.text(expression(R[t]),x=unit(-2.5,'lines'),rot=90)
  popViewport()
  dev.off()
}
