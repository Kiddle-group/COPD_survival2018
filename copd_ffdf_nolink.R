# first need to convert files to csv in linux
# cat copd_cohort_march2017_Extract_Therapy_01.txt | tr "\\t" ","  > Therapy_01.csv

print('Making ffdfs, takes a long time but is worth it, only first run ')

copd_eligible <- data.table(copd_eligible)

# turn copd_eligible into ffdf, memory efficient object
eligible <- as.ffdf(copd_eligible)

print('CLINICAL')

# CLINICAL
#
# I have commented out an attempt to go from raw data as i can't correctly format the dates.
# Could be solved by splitting into multiple smaller files, but i went with a Quint group
# pre-existing file instead.

# have trouble reformatting the dates, so just go from dta file

load.ffdf('../output/proc_data/clinical')

clinical_eligible <- merge(clinical,eligible,by='patid')

# does this work?
clinical_eligible$eventdate <- as.character(clinical_eligible$eventdate) #works fine
clinical_eligible$eventdate <- as.Date(clinical_eligible$eventdate,format='%d/%m/%Y')

clinical_before_nolink <- clinical_eligible[clinical_eligible$eventdate < clinical_eligible$first_copd,]
save.ffdf(clinical_before_nolink,dir = '../output/proc_data/clinical_before_nolink')










#save.ffdf(clinical,dir = '../output/proc_data/clinical')

#copd_eligible$first_copd <- as.Date(copd_eligible$first_copd,format='%Y-%m-%d')



# load Quint group version of clinical with nice date formats
#clinical <- read.dta13('../data/Final + new linkage/Final_cohort/clinical_base.dta')
# clinical <- as.ffdf(clinical)

# extract only data for eligible and save
#clinical_eligible <- merge(clinical[,c('patid','eventdate','medcode','enttype','adid')],copd_eligible,by='patid')
#tmp <- data.table(clinical_eligible)
#tmp_ff <- as.ffdf(tmp)
#clinical_eligible <- tmp_ff
#save.ffdf(clinical_eligible,dir = '../output/proc_data/clinical_eligible')
#rm(tmp,tmp_ff)


# post-hoc turn dates into as.Date (does this work?)
#clinical_eligible$eventdate <- as.character(clinical_eligible$eventdate)
#clinical_eligible$eventdate <- as.Date(clinical_eligible$eventdate,format='%d/%m/%Y')

# extract just the entries that occur before COPD diagnosis and save
#clinical_before <- clinical_eligible[clinical_eligible$eventdate < clinical_eligible$first_copd,]
#save.ffdf(clinical_before,dir = '../output/proc_data/clinical_before')





print('TEST')

# Test
#

# read in raw files to ffdf objet
# dates are in d/m/Y, so not the right format for dates here, so have to read as character: ''
load.ffdf('../output/proc_data/test')

# extract only data for eligible and save
test_eligible <- merge(test,eligible,by='patid')

# post-hoc turn dates into as.Date

test_eligible$eventdate <- as.character(test_eligible$eventdate)
test_eligible$eventdate <- as.Date(test_eligible$eventdate,format='%d/%m/%Y')

# extract just the entries that occur before COPD diagnosis and save
test_before_nolink <- test_eligible[test_eligible$eventdate < test_eligible$first_copd,]
save.ffdf(test_before_nolink,dir = '../output/proc_data/test_before_nolink')

print('REFERRAL')

# Referral
#

# dates are in d/m/Y, so not the right format for dates here, so have to read as character: ''
load.ffdf('../output/proc_data/referral')

# extract only data for eligible and save
referral_eligible <- merge(referral,eligible,by='patid')

# post-hoc turn dates into as.Date

referral_eligible$eventdate <- as.character(referral_eligible$eventdate)
referral_eligible$eventdate <- as.Date(referral_eligible$eventdate,format='%d/%m/%Y')

# extract just the entries that occur before COPD diagnosis and save

referral_before_nolink <- referral_eligible[referral_eligible$eventdate < referral_eligible$first_copd,]
save.ffdf(referral_before_nolink,dir = '../output/proc_data/referral_before_nolink')

print('IMMUNISATION')

# Immunisation
#

# dates are in d/m/Y, so not the right format for dates here, so have to read as character: ''
load.ffdf('../output/proc_data/immunisation')

# extract only data for eligible and save
immunisation_eligible <- merge(immunisation,eligible,by='patid')

# post-hoc turn dates into as.Date

immunisation_eligible$eventdate <- as.character(immunisation_eligible$eventdate)
immunisation_eligible$eventdate <- as.Date(immunisation_eligible$eventdate,format='%d/%m/%Y')

# extract just the entries that occur before COPD diagnosis and save
immunisation_before_nolink <- immunisation_eligible[immunisation_eligible$eventdate < immunisation_eligible$first_copd,]
save.ffdf(immunisation_before_nolink,dir = '../output/proc_data/immunisation_before_nolink')

print('ADDITIONAL')

# Additional
# no dates, so no need to filter

# dates are in d/m/Y, so not the right format for dates here, so have to read as character: ''
load.ffdf('../output/proc_data/additional')  

# extract only data for eligible and save
additional_eligible <- merge(additional,eligible,by='patid')
save.ffdf(additional_eligible_nolink,dir = '../output/proc_data/additional_eligible_nolink')

print('THERAPY')

# Therapy
#

# dates are in d/m/Y, so not the right format for dates here, so have to read as character: ''

load.ffdf('../output/proc_data/therapy')

# extract only data for eligible and save
therapy_eligible <- merge(therapy,eligible,by='patid')


# post-hoc turn dates into as.Date

therapy_eligible$eventdate <- as.character(therapy_eligible$eventdate)
therapy_eligible$eventdate <- as.Date(therapy_eligible$eventdate,format='%d/%m/%Y')

# extract just the entries that occur before COPD diagnosis and save
therapy_before_nolink <- therapy_eligible[therapy_eligible$eventdate < therapy_eligible$first_copd,]
save.ffdf(therapy_before_nolink,dir = '../output/proc_data/therapy_before_nolink')

