#working file of EpiRegress

#For function EpiRegress
#input
#input
#w: serial interval distribution
#X: the raw covariate matrix starting from at least |w| days before the start time
#I: incidence curve for local case counts
#I_imp: incidence curve for imported case counts
#start, end: start(>=2) and end time for Rt's to be estimated, the last day ("end" day)'s Rt estimate is obtained through extrapolate rather than MCMC inference
#select: list of variables selected to include in the model
#draws: number of draws used for point estimates and credible intervals
#thin: thin in MCMC iterations

#output
#Rt: Rt point estimates (posterior mean and median) and credible intervals (50% and 95%, in the form of quantiles) for each day from the start day to the end day
#beta_estimates: point estimates (posterior mean and median) and credible intervals (50% and 95%, in the form of quantiles)
#dic: DIC value of the model, used to assess the fit

#e.g.
#select a region
region='SG'
#input X covariate matrix, standardized version as X_sd
source('read_X.r')
#input case counts (autochthonous & imported) as I & I_imp 
source('read_cases.r')

#load EpiRegress functions
source('EpiRegress_rjags.r')

#set parameters
start=300; end=nrow(X)
select=1:ncol(X)
#serial interval distribution
w=dlnorm(seq(1,50),1.132,0.742)
w=w/sum(w)

out=EpiRegress(X, w, start, end)

#plot of Rt estimates
#input: Rt estimates by EpiRegress
#output: a 8cm*8cm png file of R_t estimates with posterior medians as point estimates and 95% CrI
plot(out$Rt)

