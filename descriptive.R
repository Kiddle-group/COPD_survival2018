print('Perform descriptive analysis')

patient_char <- data.frame(all=NA,none=NA,onetwo=NA,threemore=NA)

cols <- c(31,59:93) # column numbers for Barnett comorbidities

# number of comorbidities per patient
surv_copd$comorbid <- apply(surv_copd[,cols],1,'sumBool')
  
# bin patients based on that number
none <- surv_copd$comorbid==0
onetwo <- surv_copd$comorbid>0 & surv_copd$comorbid<=2
threemore <- surv_copd$comorbid>=3
  
# populate demographics table with information
patient_char['n',] <- c(dim(surv_copd)[1],length(which(none)),length(which(onetwo)),length(which(threemore)))
patient_char['percent',] <- signif(patient_char['n',] / dim(surv_copd)[1]) *100
patient_char['Female',] <- c(length(which(surv_copd[,'gender']==2)),length(which(surv_copd[none,'gender']==2)),
                             length(which(surv_copd[onetwo,'gender']==2)),length(which(surv_copd[threemore,'gender']==2)))
patient_char['Female percent',] <- signif(patient_char['Female',] / patient_char['n',]) *100
patient_char['age median',] <- c(median(surv_copd[,'age_copd']),median(surv_copd[none,'age_copd']),
                                 median(surv_copd[onetwo,'age_copd']),median(surv_copd[threemore,'age_copd']))
patient_char['age LQ',] <- c(quantile(surv_copd[,'age_copd'],0.25),quantile(surv_copd[none,'age_copd'],0.25),
                              quantile(surv_copd[onetwo,'age_copd'],0.25),quantile(surv_copd[threemore,'age_copd'],0.25))
patient_char['age UQ',] <- c(quantile(surv_copd[,'age_copd'],0.75),quantile(surv_copd[none,'age_copd'],0.75),
                             quantile(surv_copd[onetwo,'age_copd'],0.75),quantile(surv_copd[threemore,'age_copd'],0.75))

patient_char['bmi median',] <- c(median(surv_copd[!surv_copd$bmiMiss,'bmi']),median(surv_copd[!surv_copd$bmiMiss & none,'bmi']),
                                 median(surv_copd[!surv_copd$bmiMiss & onetwo,'bmi']),median(surv_copd[!surv_copd$bmiMiss & threemore,'bmi']))
patient_char['bmi LQ',] <- c(quantile(surv_copd[!surv_copd$bmiMiss,'bmi'],0.25),quantile(surv_copd[!surv_copd$bmiMiss & none,'bmi'],0.25),
                             quantile(surv_copd[!surv_copd$bmiMiss & onetwo,'bmi'],0.25),quantile(surv_copd[!surv_copd$bmiMiss & threemore,'bmi'],0.25))
patient_char['bmi UQ',] <- c(quantile(surv_copd[!surv_copd$bmiMiss,'bmi'],0.75),quantile(surv_copd[!surv_copd$bmiMiss & none,'bmi'],0.75),
                             quantile(surv_copd[!surv_copd$bmiMiss & onetwo,'bmi'],0.75),quantile(surv_copd[!surv_copd$bmiMiss & threemore,'bmi'],0.75))
patient_char['bmi missing',] <- c(length(which(surv_copd[,'bmiMiss'])),length(which(surv_copd[none,'bmiMiss'])),
                             length(which(surv_copd[onetwo,'bmiMiss'])),length(which(surv_copd[threemore,'bmiMiss'])))
patient_char['bmi missing percent',] <- signif(patient_char['bmi missing',] / patient_char['n',]) * 100

patient_char['Never smoker',] <- c(length(which(surv_copd[,'smokstatus']==2)),length(which(surv_copd[none,'smokstatus']==2)),
                                length(which(surv_copd[onetwo,'smokstatus']==2)),length(which(surv_copd[threemore,'smokstatus']==2)))
patient_char['Never smoker percent',] <- signif(patient_char['Never smoker',] / patient_char['n',]) *100
patient_char['Ex smoker',] <- c(length(which(surv_copd[,'smokstatus']==3)),length(which(surv_copd[none,'smokstatus']==3)),
                                length(which(surv_copd[onetwo,'smokstatus']==3)),length(which(surv_copd[threemore,'smokstatus']==3)))
patient_char['Ex smoker percent',] <- signif(patient_char['Ex smoker',] / patient_char['n',]) *100
patient_char['Current smoker',] <- c(length(which(surv_copd[,'smokstatus']==1)),length(which(surv_copd[none,'smokstatus']==1)),
                             length(which(surv_copd[onetwo,'smokstatus']==1)),length(which(surv_copd[threemore,'smokstatus']==1)))

patient_char['Current smoker percent',] <- signif(patient_char['Current smoker',] / patient_char['n',]) *100

patient_char['Smoking status missing',] <- c(numNA(surv_copd[,'smokstatus']),numNA(surv_copd[none,'smokstatus']),
                                  numNA(surv_copd[onetwo,'smokstatus']),numNA(surv_copd[threemore,'smokstatus']))
patient_char['Smoking status missing percent',] <- signif(patient_char['Smoking status missing',] / patient_char['n',]) * 100

patient_char['airflow median',] <- c(median(surv_copd[!surv_copd$airflowMiss,'airflow']),median(surv_copd[!surv_copd$airflowMiss & none,'airflow']),
                                     median(surv_copd[!surv_copd$airflowMiss & onetwo,'airflow']),median(surv_copd[!surv_copd$airflowMiss & threemore,'airflow']))
patient_char['airflow LQ',] <- c(quantile(surv_copd[!surv_copd$airflowMiss,'airflow'],0.25),quantile(surv_copd[!surv_copd$airflowMiss & none,'airflow'],0.25),
                                 quantile(surv_copd[!surv_copd$airflowMiss & onetwo,'airflow'],0.25),quantile(surv_copd[!surv_copd$airflowMiss & threemore,'airflow'],0.25))
patient_char['airflow UQ',] <- c(quantile(surv_copd[!surv_copd$airflowMiss,'airflow'],0.75),quantile(surv_copd[!surv_copd$airflowMiss & none,'airflow'],0.75),
                                 quantile(surv_copd[!surv_copd$airflowMiss & onetwo,'airflow'],0.75),quantile(surv_copd[!surv_copd$airflowMiss & threemore,'airflow'],0.75))
patient_char['airflow missing',] <- c(length(which(surv_copd[,'airflowMiss'])),length(which(surv_copd[none,'airflowMiss'])),
                                      length(which(surv_copd[onetwo,'airflowMiss'])),length(which(surv_copd[threemore,'airflowMiss'])))
patient_char['airflow missing percent',] <- signif(patient_char['airflow missing',] / patient_char['n',]) * 100



patient_char['IMD median',] <- c(median(surv_copd[,'imd2010_20'],na.rm=T),median(surv_copd[ none,'imd2010_20'],na.rm=T),
                                 median(surv_copd[ onetwo,'imd2010_20'],na.rm=T),median(surv_copd[threemore,'imd2010_20'],na.rm=T))
patient_char['IMD LQ',] <- c(quantile(surv_copd[,'imd2010_20'],0.25,na.rm=T),quantile(surv_copd[ none,'imd2010_20'],0.25,na.rm=T),
                             quantile(surv_copd[ onetwo,'imd2010_20'],0.25,na.rm=T),quantile(surv_copd[threemore,'imd2010_20'],0.25,na.rm=T))
patient_char['IMD UQ',] <- c(quantile(surv_copd[,'imd2010_20'],0.75,na.rm=T),quantile(surv_copd[ none,'imd2010_20'],0.75,na.rm=T),
                             quantile(surv_copd[onetwo,'imd2010_20'],0.75,na.rm=T),quantile(surv_copd[threemore,'imd2010_20'],0.75,na.rm=T))
patient_char['IMD missing',] <- c(numNA(surv_copd[,'imd2010_20']),numNA(surv_copd[none,'imd2010_20']),
                                  numNA(surv_copd[onetwo,'imd2010_20']),numNA(surv_copd[threemore,'imd2010_20']))
patient_char['IMD missing percent',] <- signif(patient_char['IMD missing',] / patient_char['n',]) * 100

ind_all <- 1:dim(surv_copd)[1]

patient_char['Num deaths 5 ',] <- c(num_dead(ind_all,5),num_dead(which(none),5),
                                                      num_dead(which(onetwo),5),num_dead(which(threemore),5))
patient_char['Percent deaths 5',] <- (patient_char['Num deaths 5',]*100) / patient_char['n',]

write.csv(patient_char,file='../output/demographics.csv') 



# now create table given number with and prevalence of each co-morbidity

comorbid_prevalence <- data.frame(colnames(surv_copd)[comorbidities =cols],num=NA,prev=NA)

# loop over barnett comorbidities
for (i in 1:length(cols)){
  
  print(i)
  
  # count number with
  comorbid_prevalence[i,2] <- length(which(surv_copd[,cols[i]])) 
  
}

# turn into prevalence
comorbid_prevalence$prev <- (comorbid_prevalence$num / dim(surv_copd)[1])*100

# order in terms of prevalence
comorbid_prevalence <- comorbid_prevalence[order(comorbid_prevalence$prev,decreasing = T),]

write.csv(comorbid_prevalence,file = '../output/comorbid_prevalence.csv')
