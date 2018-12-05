 load('../../newscratch/steven/CPRD_COPD/output/final.RData')

coef_df <- data.frame(logOddsRatio=fit$coefficients,confint.default(fit))

coef_df <- coef_df[order(coef_df$logOddsRatio),]

coef_df$label <- rownames(coef_df)

colnames(coef_df)[2:3] <- c('log_lower','log_upper')

coef_df$oddsRatio <- exp(coef_df$logOddsRatio)
coef_df$upper <- exp(coef_df$log_upper)
coef_df$lower <- exp(coef_df$log_lower)



wtf <- c('Five_Year_Mortality','Chronic Kidney Disease','Alcohol problems','Anorexia/bulimia','Atrial fibrillation','Blindness/low vision','Bronciectasis','Chronic liver disease','Sinusitis','Coronary heart disease','Dementia','Diabetes','Diverticular disease','Hearing loss','Heart failure','Hypertension',
'Inflammatory bowel disease','Learning disability','Multiple sclerosis','Peripheral vascular disorder','Parkinsons','Prostate disorders','Substance abuse','Connective tissue disorders',
'Stroke','Thyroid disorders','Anxiety','Asthma','Cancer','Constipation','Depression','Epilepsy','Painful condition','Irritable bowel syndrome','Migrane','Psoriasis/eczema','Pyschosis/bipolar','Normalised_MM_burden')

wtf <- sort(wtf)

coef_df$label <- rev(c('Intercept','Never smoker','Ex smoker','eGFR NR*','Female','IBS',wtf[32],wtf[4],wtf[17],wtf[35],'CRP NR',wtf[20],wtf[38],'Albumin','BMI',wtf[22],
'Airflow obstruction','eGFR','Platelets','CRP','IMD','GORD','Platelets NR','Albumin NR',wtf[33],'Age',wtf[c(11,7,13,29,23,6,3,15,12,36,16,18,34)],'BMI NR',wtf[31],'FEV1 NR',
wtf[c(5,37,1,8,21)]))

coef_df$label <- factor(coef_df$label,levels=unique(coef_df$label))


ggplot(data=coef_df,aes(x=label,y=oddsRatio,ymin=lower,ymax=upper)) + geom_pointrange() + geom_hline(yintercept=1,lty=2) + coord_flip() + xlab('Variable') + ylab('Survival odds ratio (95% CI)') + theme_bw() +ylim(0.25,2)

write.csv(coef_df,file='finalModel.csv')
