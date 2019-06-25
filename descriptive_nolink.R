print('Perform descriptive analysis')



cols <- c(28,56:90) # column numbers for Barnett comorbidities

tod_before <- which(copd_eligible$tod < copd_eligible$first_copd)

copd_eligible <- copd_eligible[-tod_before,]


copd_eligible$fiveYears <- (copd_eligible$first_copd +(5*365))

copd_eligible$dead5 <- copd_eligible$fiveYears > copd_eligible$deathdate

copd_eligible$dead5[is.na(copd_eligible$dead5)] <- FALSE


copd_eligible$censor <- copd_eligible$tod < copd_eligible$fiveYears

copd_eligible$censor[is.na(copd_eligible$censor)] <- FALSE

copd_eligible$censoredAlive <- copd_eligible$censor & !copd_eligible$dead5

surv_copd <- copd_eligible

# number of comorbidities per patient

surv_copd$comorbid <- apply(surv_copd[,cols],1,'sumBool')
  

patient_char <- data.frame(all=NA,none=NA,onetwo=NA,threemore=NA)

# bin patients based on that number
none <- surv_copd$comorbid==0
onetwo <- surv_copd$comorbid>0 & surv_copd$comorbid<=2
threemore <- surv_copd$comorbid>=3

surv_copd <- data.table(surv_copd)
  
# populate demographics table with information
patient_char['n',] <- c(dim(surv_copd)[1],length(which(none)),length(which(onetwo)),length(which(threemore)))
patient_char['percent',] <- signif(patient_char['n',] / dim(surv_copd)[1]) *100
patient_char['Female',] <- c(length(which(surv_copd[,'gender']==2)),length(which(surv_copd[none,'gender']==2)),
                             length(which(surv_copd[onetwo,'gender']==2)),length(which(surv_copd[threemore,'gender']==2)))
patient_char['Female percent',] <- signif(patient_char['Female',] / patient_char['n',]) *100
patient_char['age median',] <- c(median(surv_copd[,age_copd]),median(surv_copd[none,age_copd]),
                                 median(surv_copd[onetwo,age_copd]),median(surv_copd[threemore,age_copd]))
patient_char['age LQ',] <- c(quantile(surv_copd[,age_copd],0.25),quantile(surv_copd[none,age_copd],0.25),
                              quantile(surv_copd[onetwo,age_copd],0.25),quantile(surv_copd[threemore,age_copd],0.25))
patient_char['age UQ',] <- c(quantile(surv_copd[,age_copd],0.75),quantile(surv_copd[none,age_copd],0.75),
                             quantile(surv_copd[onetwo,age_copd],0.75),quantile(surv_copd[threemore,age_copd],0.75))

patient_char['bmi median',] <- c(median(surv_copd[!surv_copd$bmiMiss,bmi]),median(surv_copd[!surv_copd$bmiMiss & none,bmi]),
                                 median(surv_copd[!surv_copd$bmiMiss & onetwo,bmi]),median(surv_copd[!surv_copd$bmiMiss & threemore,bmi]))
patient_char['bmi LQ',] <- c(quantile(surv_copd[!surv_copd$bmiMiss,bmi],0.25),quantile(surv_copd[!surv_copd$bmiMiss & none,bmi],0.25),
                             quantile(surv_copd[!surv_copd$bmiMiss & onetwo,bmi],0.25),quantile(surv_copd[!surv_copd$bmiMiss & threemore,bmi],0.25))
patient_char['bmi UQ',] <- c(quantile(surv_copd[!surv_copd$bmiMiss,bmi],0.75),quantile(surv_copd[!surv_copd$bmiMiss & none,bmi],0.75),
                             quantile(surv_copd[!surv_copd$bmiMiss & onetwo,bmi],0.75),quantile(surv_copd[!surv_copd$bmiMiss & threemore,bmi],0.75))
patient_char['bmi missing',] <- c(length(which(surv_copd[,bmiMiss])),length(which(surv_copd[none,bmiMiss])),
                             length(which(surv_copd[onetwo,bmiMiss])),length(which(surv_copd[threemore,bmiMiss])))
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

patient_char['airflow median',] <- c(median(surv_copd[!surv_copd$airflowMiss,airflow]),median(surv_copd[!surv_copd$airflowMiss & none,airflow]),
                                     median(surv_copd[!surv_copd$airflowMiss & onetwo,airflow]),median(surv_copd[!surv_copd$airflowMiss & threemore,airflow]))
patient_char['airflow LQ',] <- c(quantile(surv_copd[!surv_copd$airflowMiss,airflow],0.25),quantile(surv_copd[!surv_copd$airflowMiss & none,airflow],0.25),
                                 quantile(surv_copd[!surv_copd$airflowMiss & onetwo,airflow],0.25),quantile(surv_copd[!surv_copd$airflowMiss & threemore,airflow],0.25))
patient_char['airflow UQ',] <- c(quantile(surv_copd[!surv_copd$airflowMiss,airflow],0.75),quantile(surv_copd[!surv_copd$airflowMiss & none,airflow],0.75),
                                 quantile(surv_copd[!surv_copd$airflowMiss & onetwo,airflow],0.75),quantile(surv_copd[!surv_copd$airflowMiss & threemore,airflow],0.75))
patient_char['airflow missing',] <- c(length(which(surv_copd[,airflowMiss])),length(which(surv_copd[none,airflowMiss])),
                                      length(which(surv_copd[onetwo,airflowMiss])),length(which(surv_copd[threemore,airflowMiss])))
patient_char['airflow missing percent',] <- signif(patient_char['airflow missing',] / patient_char['n',]) * 100

patient_char['England',] <- c(length(which(surv_copd[,'region']<=10)),length(which(surv_copd[none,'region']<=10)),
                                   length(which(surv_copd[onetwo,'region']<=10)),length(which(surv_copd[threemore,'region']<=10)))

patient_char['England percent',] <- signif(patient_char['England',] / patient_char['n',]) * 100

patient_char['NI',] <- c(length(which(surv_copd[,'region']==11)),length(which(surv_copd[none,'region']==11)),
                              length(which(surv_copd[onetwo,'region']==11)),length(which(surv_copd[threemore,'region']==11)))

patient_char['NI percent',] <- signif(patient_char['NI',] / patient_char['n',]) * 100

patient_char['Scotland',] <- c(length(which(surv_copd[,'region']==12)),length(which(surv_copd[none,'region']==12)),
                         length(which(surv_copd[onetwo,'region']==12)),length(which(surv_copd[threemore,'region']==12)))

patient_char['Scotland percent',] <- signif(patient_char['Scotland',] / patient_char['n',]) * 100

patient_char['Wales',] <- c(length(which(surv_copd[,'region']==13)),length(which(surv_copd[none,'region']==13)),
                         length(which(surv_copd[onetwo,'region']==13)),length(which(surv_copd[threemore,'region']==13)))

patient_char['Wales percent',] <- signif(patient_char['Wales',] / patient_char['n',]) * 100

patient_char['Num censor 5 ',] <- c(length(which(surv_copd$censoredAlive)),length(which(surv_copd$censoredAlive[none])),
                                    length(which(surv_copd$censoredAlive[onetwo])),length(which(surv_copd$censoredAlive[threemore])))
patient_char['Percent censor 5',] <- (patient_char['Num censor 5',]*100) / patient_char['n',]

patient_char['Num deaths 5 ',] <- c(length(which(surv_copd$dead5)),length(which(surv_copd$dead5[none])),
                                    length(which(surv_copd$dead5[onetwo])),length(which(surv_copd$dead5[threemore])))
patient_char['Percent deaths 5',] <- (patient_char['Num deaths 5',]*100) / patient_char['n',]

write.csv(patient_char,file='../output/demographics_nolink.csv') 






surv_copd2 <- surv_copd[!surv_copd$censoredAlive,]

patient_char <- data.frame(all=NA,none=NA,onetwo=NA,threemore=NA)

# bin patients based on that number
none <- surv_copd2$comorbid==0
onetwo <- surv_copd2$comorbid>0 & surv_copd2$comorbid<=2
threemore <- surv_copd2$comorbid>=3

surv_copd2 <- data.table(surv_copd2)

# populate demographics table with information
patient_char['n',] <- c(dim(surv_copd2)[1],length(which(none)),length(which(onetwo)),length(which(threemore)))
patient_char['percent',] <- signif(patient_char['n',] / dim(surv_copd2)[1]) *100
patient_char['Female',] <- c(length(which(surv_copd2[,'gender']==2)),length(which(surv_copd2[none,'gender']==2)),
                             length(which(surv_copd2[onetwo,'gender']==2)),length(which(surv_copd2[threemore,'gender']==2)))
patient_char['Female percent',] <- signif(patient_char['Female',] / patient_char['n',]) *100
patient_char['age median',] <- c(median(surv_copd2[,age_copd]),median(surv_copd2[none,age_copd]),
                                 median(surv_copd2[onetwo,age_copd]),median(surv_copd2[threemore,age_copd]))
patient_char['age LQ',] <- c(quantile(surv_copd2[,age_copd],0.25),quantile(surv_copd2[none,age_copd],0.25),
                             quantile(surv_copd2[onetwo,age_copd],0.25),quantile(surv_copd2[threemore,age_copd],0.25))
patient_char['age UQ',] <- c(quantile(surv_copd2[,age_copd],0.75),quantile(surv_copd2[none,age_copd],0.75),
                             quantile(surv_copd2[onetwo,age_copd],0.75),quantile(surv_copd2[threemore,age_copd],0.75))

patient_char['bmi median',] <- c(median(surv_copd2[!surv_copd2$bmiMiss,bmi]),median(surv_copd2[!surv_copd2$bmiMiss & none,bmi]),
                                 median(surv_copd2[!surv_copd2$bmiMiss & onetwo,bmi]),median(surv_copd2[!surv_copd2$bmiMiss & threemore,bmi]))
patient_char['bmi LQ',] <- c(quantile(surv_copd2[!surv_copd2$bmiMiss,bmi],0.25),quantile(surv_copd2[!surv_copd2$bmiMiss & none,bmi],0.25),
                             quantile(surv_copd2[!surv_copd2$bmiMiss & onetwo,bmi],0.25),quantile(surv_copd2[!surv_copd2$bmiMiss & threemore,bmi],0.25))
patient_char['bmi UQ',] <- c(quantile(surv_copd2[!surv_copd2$bmiMiss,bmi],0.75),quantile(surv_copd2[!surv_copd2$bmiMiss & none,bmi],0.75),
                             quantile(surv_copd2[!surv_copd2$bmiMiss & onetwo,bmi],0.75),quantile(surv_copd2[!surv_copd2$bmiMiss & threemore,bmi],0.75))
patient_char['bmi missing',] <- c(length(which(surv_copd2[,bmiMiss])),length(which(surv_copd2[none,bmiMiss])),
                                  length(which(surv_copd2[onetwo,bmiMiss])),length(which(surv_copd2[threemore,bmiMiss])))
patient_char['bmi missing percent',] <- signif(patient_char['bmi missing',] / patient_char['n',]) * 100

patient_char['Never smoker',] <- c(length(which(surv_copd2[,'smokstatus']==2)),length(which(surv_copd2[none,'smokstatus']==2)),
                                   length(which(surv_copd2[onetwo,'smokstatus']==2)),length(which(surv_copd2[threemore,'smokstatus']==2)))
patient_char['Never smoker percent',] <- signif(patient_char['Never smoker',] / patient_char['n',]) *100
patient_char['Ex smoker',] <- c(length(which(surv_copd2[,'smokstatus']==3)),length(which(surv_copd2[none,'smokstatus']==3)),
                                length(which(surv_copd2[onetwo,'smokstatus']==3)),length(which(surv_copd2[threemore,'smokstatus']==3)))
patient_char['Ex smoker percent',] <- signif(patient_char['Ex smoker',] / patient_char['n',]) *100
patient_char['Current smoker',] <- c(length(which(surv_copd2[,'smokstatus']==1)),length(which(surv_copd2[none,'smokstatus']==1)),
                                     length(which(surv_copd2[onetwo,'smokstatus']==1)),length(which(surv_copd2[threemore,'smokstatus']==1)))

patient_char['Current smoker percent',] <- signif(patient_char['Current smoker',] / patient_char['n',]) *100

patient_char['Smoking status missing',] <- c(numNA(surv_copd2[,'smokstatus']),numNA(surv_copd2[none,'smokstatus']),
                                             numNA(surv_copd2[onetwo,'smokstatus']),numNA(surv_copd2[threemore,'smokstatus']))
patient_char['Smoking status missing percent',] <- signif(patient_char['Smoking status missing',] / patient_char['n',]) * 100

patient_char['airflow median',] <- c(median(surv_copd2[!surv_copd2$airflowMiss,airflow]),median(surv_copd2[!surv_copd2$airflowMiss & none,airflow]),
                                     median(surv_copd2[!surv_copd2$airflowMiss & onetwo,airflow]),median(surv_copd2[!surv_copd2$airflowMiss & threemore,airflow]))
patient_char['airflow LQ',] <- c(quantile(surv_copd2[!surv_copd2$airflowMiss,airflow],0.25),quantile(surv_copd2[!surv_copd2$airflowMiss & none,airflow],0.25),
                                 quantile(surv_copd2[!surv_copd2$airflowMiss & onetwo,airflow],0.25),quantile(surv_copd2[!surv_copd2$airflowMiss & threemore,airflow],0.25))
patient_char['airflow UQ',] <- c(quantile(surv_copd2[!surv_copd2$airflowMiss,airflow],0.75),quantile(surv_copd2[!surv_copd2$airflowMiss & none,airflow],0.75),
                                 quantile(surv_copd2[!surv_copd2$airflowMiss & onetwo,airflow],0.75),quantile(surv_copd2[!surv_copd2$airflowMiss & threemore,airflow],0.75))
patient_char['airflow missing',] <- c(length(which(surv_copd2[,airflowMiss])),length(which(surv_copd2[none,airflowMiss])),
                                      length(which(surv_copd2[onetwo,airflowMiss])),length(which(surv_copd2[threemore,airflowMiss])))
patient_char['airflow missing percent',] <- signif(patient_char['airflow missing',] / patient_char['n',]) * 100

patient_char['England',] <- c(length(which(surv_copd2[,'region']<=10)),length(which(surv_copd2[none,'region']<=10)),
                              length(which(surv_copd2[onetwo,'region']<=10)),length(which(surv_copd2[threemore,'region']<=10)))

patient_char['England percent',] <- signif(patient_char['England',] / patient_char['n',]) * 100

patient_char['NI',] <- c(length(which(surv_copd2[,'region']==11)),length(which(surv_copd2[none,'region']==11)),
                         length(which(surv_copd2[onetwo,'region']==11)),length(which(surv_copd2[threemore,'region']==11)))

patient_char['NI percent',] <- signif(patient_char['NI',] / patient_char['n',]) * 100

patient_char['Scotland',] <- c(length(which(surv_copd2[,'region']==12)),length(which(surv_copd2[none,'region']==12)),
                               length(which(surv_copd2[onetwo,'region']==12)),length(which(surv_copd2[threemore,'region']==12)))

patient_char['Scotland percent',] <- signif(patient_char['Scotland',] / patient_char['n',]) * 100

patient_char['Wales',] <- c(length(which(surv_copd2[,'region']==13)),length(which(surv_copd2[none,'region']==13)),
                            length(which(surv_copd2[onetwo,'region']==13)),length(which(surv_copd2[threemore,'region']==13)))

patient_char['Wales percent',] <- signif(patient_char['Wales',] / patient_char['n',]) * 100

patient_char['Num censor 5 ',] <- c(length(which(surv_copd2$censoredAlive)),length(which(surv_copd2$censoredAlive[none])),
                                    length(which(surv_copd2$censoredAlive[onetwo])),length(which(surv_copd2$censoredAlive[threemore])))
patient_char['Percent censor 5',] <- (patient_char['Num censor 5',]*100) / patient_char['n',]

patient_char['Num deaths 5 ',] <- c(length(which(surv_copd2$dead5)),length(which(surv_copd2$dead5[none])),
                                    length(which(surv_copd2$dead5[onetwo])),length(which(surv_copd2$dead5[threemore])))
patient_char['Percent deaths 5',] <- (patient_char['Num deaths 5',]*100) / patient_char['n',]

write.csv(patient_char,file='../output/demographics_nocensor.csv') 




# now create table given number with and prevalence of each co-morbidity

comorbid_prevalence <- data.frame(colnames(surv_copd)[comorbidities =cols],num=NA,prev=NA)

# loop over barnett comorbidities
for (i in 1:length(cols)){
  
  print(i)
  
  # count number with
  comorbid_prevalence[i,2] <- length(which(copd_eligible[,cols[i]])) 
  
}

# turn into prevalence
comorbid_prevalence$prev <- (comorbid_prevalence$num / dim(surv_copd)[1])*100

# order in terms of prevalence
comorbid_prevalence <- comorbid_prevalence[order(comorbid_prevalence$prev,decreasing = T),]

write.csv(comorbid_prevalence,file = '../output/comorbid_prevalence_nolink.csv')
