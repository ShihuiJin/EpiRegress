#read X matrix, including reading cases for the dataset till the end of 2021
#mobility data
if (region=='SG'){
  mobility_2020=read.csv("google_mobility_2020.csv")
  mobility_2021=read.csv("google_mobility_2021.csv")
}else{
  mobility_2020=read.csv(paste0('2020_',region,'_Region_Mobility_Report.csv'))
  mobility_2021=read.csv(paste0('2021_',region,'_Region_Mobility_Report.csv'))
}

#days in each month,2021
m=c(31,28,31,30,31,30,31,31,30,31,30,31)
m=cumsum(m)

#vaccination
vax=read.csv(paste0('vaccination_',region,'.csv'))
#2020 part
X_2020=mobility_2020[(ncol(mobility_2020)-5):ncol(mobility_2020)]
X_2020$vaccination=rep(0,nrow(X_2020))
X_2020$delta=rep(0,nrow(X_2020))
#2021 part
X_2021=mobility_2021[(ncol(mobility_2021)-5):ncol(mobility_2021)]
X_2021$vaccination=rep(0,nrow(X_2021))
if(region=='SG'){
  p_vax=rev(taRifx::destring(vax$per.100.people))
  X_2021$vaccination=p_vax[1:nrow(X_2021)]
}
if(region=="TW"){
  p_vax=rev(taRifx::destring(vax$per.100.people)) #vaccination doses per 100 people
  X_2021$vaccination[(21+m[2]):nrow(X_2021)]=p_vax[1:(nrow(X_2021)-m[2]-20)] #for TW, vaccination data from 21/03,2021
}
if(region=="NZ"){
  p_vax=rev(taRifx::destring(vax$per.100.people)) #vaccination doses per 100 people
  X_2021$vaccination[(17+m[1]):nrow(X_2021)]=p_vax[1:(nrow(X_2021)-m[1]-16)] #for NZ, vaccination data from 17/02,2021
}
if(region=="NSW"){
  p_vax=rev(taRifx::destring(vax$total.doses))
  X_2021$vaccination[(16+m[1]):nrow(X_2021)]=p_vax[1:(nrow(X_2021)-m[1]-15)]/8.2e4 #for AU_NSW, vaccination data from 16/02,2021, total population 8.2 million
} 

#delta variant proportion
delta_count=as.matrix(read.csv(paste0('delta_count_',region,'.csv'))) #week data, starting from week 0
if(region=="NSW") {
  delta_day=rep(0,nrow(X_2021))
  for(i in 1:nrow(X_2021))
  {
    k=ceiling((i-2)/7)+1
    if(k<=nrow(delta_count)){
      delta_day[i]=taRifx::destring(delta_count[k,1])/taRifx::destring(delta_count[k,2])*100
      if(is.na(delta_day[i])==1)delta_day[i]=0
    }else{
      delta_day[i]=100
    }
  }
}else if(region=="TW"){   #change delta_variant proportion manually for Taiwan due to unknown data. 
  delta_day=rep(0,nrow(X_2021))
  for(i in 1:nrow(X_2021))
  {
    k=ceiling((i-2)/7)+1
    if(k<=27){
      delta_day[i]=delta_count[k,1]/delta_count[k,2]*100
      if(is.na(delta_day[i])==1)delta_day[i]=0
    }else{
      delta_day[i]=100
    }
  }
}else{
  delta_day=rep(0,nrow(X_2021))
  for(i in 1:nrow(X_2021))
  {
    k=ceiling((i-2)/7)+1
    if(k<=nrow(delta_count)){
      delta_day[i]=delta_count[k,1]/delta_count[k,2]*100
      if(is.na(delta_day[i])==1)delta_day[i]=0
    }else{
      delta_day[i]=100
    }
  }
}
X_2021$delta=delta_day

#phase data
if(region=="SG"){
  #different phases in Singapore
  #2020
  m=c(31,29,31,30,31,30,31,31,30,31,30,31)
  m=cumsum(m)
  #phase 2: June 18 - Dec 27, 2020, May 8 - 15, 2021
  X_2020$P2=rep(0,nrow(X_2020))
  X_2020$P2[(m[5]+18):(m[11]+27)-(366-nrow(X_2020))]=1
  #phase 2 heightend alert: May 16 - June 13, July 22 - Aug 9,2021
  X_2020$P2h=rep(0,nrow(X_2020))
  #phase 3: Dec 28, 2020 - May 7,2021
  X_2020$P3=rep(0,nrow(X_2020))
  X_2020$P3[(nrow(X_2020)-3):nrow(X_2020)]=1
  #phase 3 heightend alert: June 14 - July 21, 2021
  X_2020$P3h=rep(0,nrow(X_2020))
  #phase Preparatory Stage, Aug 10 - Sept 26, 2021
  X_2020$PP=rep(0,nrow(X_2020))
  #phase Stationary, Sept 27 - Sept 30, 2021
  X_2020$SP=rep(0,nrow(X_2020))
  #2021
  m=c(31,28,31,30,31,30,31,31,30,31,30,31)
  m=cumsum(m)
  #phase 2: May 8 - 15
  X_2021$P2=rep(0,nrow(X_2021))
  X_2021$P2[(m[4]+8):(m[4]+15)]=1
  #phase 2 heightend alert: May 16 - June 13, July 22 - Aug 9, 2021
  X_2021$P2h=rep(0,nrow(X_2021))
  X_2021$P2h[(m[4]+16):(m[5]+13)]=1
  X_2021$P2h[(m[6]+22):(m[7]+9)]=1
  #phase 3: Dec 28, 2020 - May 7, 2021
  X_2021$P3=rep(0,nrow(X_2021))
  X_2021$P3[1:(m[4]+7)]=1
  #phase 3 heightend alert: June 14 - July 21, 2021
  X_2021$P3h=rep(0,nrow(X_2021))
  X_2021$P3h[(m[5]+14):(m[6]+21)]=1
  #phase Preparatory Stage, Aug 10 - Sept 26, 2021
  X_2021$PP=rep(0,nrow(X_2021))
  X_2021$PP[(m[7]+10):(m[8]+26)]=1
  #phase Stationary, Sept 27 - Sept 30, 2021
  X_2021$SP=rep(0,nrow(X_2021))
  X_2021$SP[(m[8]+27):nrow(X_2021)]=1
}

X=rbind(as.matrix(X_2020),as.matrix(X_2021[,(ncol(X_2021)-ncol(X_2020)+1):ncol(X_2021)]))

#exclude one dummy variable for Singapore
if (region=='SG') X=X[,-ncol(X_2020)]

#Oxford policy data
oxford_policy=as.matrix(read.csv(paste0('Oxford_new_',region,'.csv'))) #20 variables
for(i in 2:ncol(oxford_policy))
{
  X=cbind(X,taRifx::destring(oxford_policy[,i]))
}
X_sd=X

#normalizing X
for(i in 1:ncol(X))
{
  if(var(X[,i])>0){
    X_sd[,i]=(X[,i]-mean(X[,i]))/sqrt(var(X[,i]))
  }else{
    X_sd[,i]=0
  }
}

#variance
factor=rep(0,ncol(X))
for(i in 1:length(factor))
{
  factor[i]=sqrt(var(X[,i]))
}
write.csv(factor,paste0('factor_',region,'.csv'),row.names=FALSE)

X_sd=X_sd[,which(factor!=0)]

rm(delta_count,oxford_policy,vax,X_2020,X_2021)
rm(delta_day,i,k,m,p_vax)
#rm(mobility_2020,mobility_2021)
