# ALBUMIN

# identify latest instances from test table
test_relevant <- test_before[test_before$enttype == 152,c('patid','eventdate','data1','data2','first_copd','enttype')]
albumin <- find_latest(test_relevant)

#QC and add
albumin <- albumin[albumin$data2<70,] # remove outlier
albumin <- albumin[albumin$data1==3 & !is.na(albumin$data2),]
albumin$albumin <- albumin$data2
surv_copd <- merge(surv_copd,albumin[,c('patid','albumin')],by='patid',all.x=T)


# C-REACTIVE PROTEIN

# identify latest instances from test table
test_relevant <- test_before[test_before$enttype == 280,c('patid','eventdate','data1','data2','first_copd','enttype')]
crp <- find_latest(test_relevant)

# QC and add
crp <- crp[crp$data1==3 & !is.na(crp$data2),]
crp <- crp[crp$data2<370,] # remove outliers, visual, jump in plot of ordered values
crp$crp <- crp$data2
surv_copd <- merge(surv_copd,crp[,c('patid','crp')],by='patid',all.x=T)


# PLATELETS

# identify latest instances from test table
test_relevant <- test_before[test_before$enttype == 189,c('patid','eventdate','data1','data2','first_copd','enttype')]
plate <- find_latest(test_relevant)

# QC and add
plate <- plate[plate$data1 == 3 & !is.na(plate$data2),]
plate <- plate[plate$data2<1280,] # remove outliers, visual, jump in plot of ordered values
plate$plate <- plate$data2
surv_copd <- merge(surv_copd,plate[,c('patid','plate')],by='patid',all.x=T)
