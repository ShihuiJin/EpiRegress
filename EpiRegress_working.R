#working file of EpiRegress
#select a region
region='SG'
#input X covariate matrix, standardized version as X_sd
source('read_X.r')
#input case counts (autochthonous & imported) as I & I_imp 
source('read_cases.r')

#load EpiRegress functions
source('EpiRegress.r')


#For function EpiRegress
#input
#X_sd: X covariate matrix
#I, I_imp: local and imported cases
#start, end:start (>2) and end time for inference (Rt estimates that are involved in the likelihood calculation is actually from start-length(w) to end-1)
#w: serial interval distribution
#delta: time lag
#select: covariates that are used in inference
#draws: number of draws used for point estimates and credible intervals
#thin: thin in MCMC iterations
#raw: if exact MCMC draws are provided
#output
#Rt_estimates: Rt point estimates and credible intervals for each day from the start day to one day before the end day
#beta_estimates: point estimates and credible intervals for each coefficient 
#DIC: DIC value of the model, used to assess the fit
#pars: exact MCMC draws for the (hyper-)parameters in the model

#E.g.
start=322; end=nrow(X_sd)
#serial interval distribution
w=dlnorm(seq(1,50),1.132,0.742); w=w/sum(w)
#other parameters
delta=0; select=c(1:ncol(X_sd)); draws=1e4; thin=1; raw=FALSE
out=EpiRegress(X_sd, I, I_imp,start,end,w,delta,select, draws,thin,raw)
out$DIC #DIC of the model
R2=matrix(0,end,6)
R2[start:(end-1),]=out$Rt_estimates #Rt estimates from the first day with record (not the 'start' day) to the end day (the last day is extrapolated)
store=out$pars

