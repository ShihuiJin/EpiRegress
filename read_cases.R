#read number of cases
if (region=='SG'){
  covid_2020=read.csv('covid_cases_breakdown_2020.csv')
  covid_2021=read.csv('covid_cases_breakdown_2021.csv')
  I=c(covid_2020$Total.Community[(nrow(covid_2020)-nrow(mobility_2020)+1):nrow(covid_2020)], taRifx::destring(covid_2021$Community))
  I_imp=c(covid_2020$Imported[(nrow(covid_2020)-nrow(mobility_2020)+1):nrow(covid_2020)], covid_2021$Imported)
  #I_imp=c(covid_2020$Imported[(nrow(covid_2020)-nrow(mobility_2020)+1):nrow(covid_2020)]+taRifx::destring(covid_2020$Dorm.Residents[(nrow(covid_2020)-nrow(mobility_2020)+1):nrow(covid_2020)]), covid_2021$Imported+covid_2021$Dorm.Residents)
}else{
  covid=read.csv(paste0('case_',region,'.csv'))
  #for TW: breakdown case count starting from 01/01/2020; total case count starting from 22/01/2020
  if (region=="TW"){
    I=covid$local[(366-nrow(mobility_2020)+1): (366+nrow(mobility_2021))]
    I_imp=covid$imported[(366-nrow(mobility_2020)+1): (366+nrow(mobility_2021))]
  }
  #for NSW: case count starting from 26/01/2020
  if (region=='NSW'){
    I=rev(covid$LOCAL)[(366-25-nrow(mobility_2020)+1): (366-25+nrow(mobility_2021))]
    I=taRifx::destring(I)
    I_imp=rev(covid$OVERSEAS)[(366-25-nrow(mobility_2020)+1): (366-25+nrow(mobility_2021))]
    I_imp=taRifx::destring(I_imp)
  }
  #for NZ: case count starting from 28/02/2020
  if (region=='NZ'){
    case=rep(0,nrow(covid)+27+31)
    case[(28+31):(length(case))]=covid$local
    I=case[(366-nrow(mobility_2020)+1): (366+nrow(mobility_2021))]
    case[(28+31):(length(case))]=covid$imported
    I_imp=case[(366-nrow(mobility_2020)+1): (366+nrow(mobility_2021))]
    remove(case)
  }
}