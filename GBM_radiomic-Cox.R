##### Cox regression with radiomics features.
rm(list=ls())

library(survival)
load('./GBM_dataset.RData')


radiomics_df1 = read.csv('./TCGA_GBM_radiomicFeatures_Train.csv',header=T)
radiomics_df2 = read.csv('./TCGA_GBM_radiomicFeatures_Test.csv',header=T)
#dim(radiomics_df1);dim(radiomics_df2)
radiomics_df = rbind(radiomics_df1,radiomics_df2)

radiomics_df = dplyr::inner_join(gbm_dat[,1:7], radiomics_df, by = c('id'='ID'))
radiomics_df = radiomics_df[,-which(colnames(radiomics_df)== "Date")]
radiomics_df = radiomics_df[,!apply(is.na(radiomics_df),2,any)]
dim(radiomics_df)

glmnet.fcox.cv <- glmnet::cv.glmnet(
  x=as.matrix(radiomics_df[,-(1:3)]),
  y=Surv(time = radiomics_df$os, event = radiomics_df$status),
  penalty.factor = c(0,0,0,0,rep(1,ncol(radiomics_df)-7)),
  family = 'cox', nfolds = 10, #default=10
  alpha=1)  ##Lasso
#plot(glmnet.fcox.cv)

glmnet.fcox.fit = glmnet::glmnet(x=as.matrix(radiomics_df[,-(1:3)]),
                                 y=Surv(time = radiomics_df$os, event = radiomics_df$status),
                                 penalty.factor = c(0,0,0,0,rep(1,ncol(radiomics_df)-7)),
                                 family = 'cox', lambda = glmnet.fcox.cv$lambda[which.min(glmnet.fcox.cv$cvm)],
                                 alpha=1)  ##Lasso
#coef(glmnet.fcox.fit)



#Selected features.
sig.index = which(abs(as.vector(coef(glmnet.fcox.fit)))>0)
radiomics_df_cox = radiomics_df[,c(2:3,sig.index+3)]

###LOOCV predictive risk
loocv.risk = 0
for(ll in 1:nrow(radiomics_df_cox)){
  loocv.cox = survival::coxph(Surv(time = os, event = status)~.,data = radiomics_df_cox[-ll,])
  loocv.risk[ll] = predict(loocv.cox, radiomics_df_cox[ll,-(1:2)], type='risk')
}

radiomics_df$group = NA
radiomics_df$group[loocv.risk>=median(loocv.risk)] = 'high'
radiomics_df$group[loocv.risk<median(loocv.risk)] = 'low'


logrank = survival::survdiff( survival::Surv(time=os,event=status)~group, data = radiomics_df)
km_plot_radiomics = survminer::ggsurvplot(survfit(Surv(os, status) ~ group, data = radiomics_df), 
                                          conf.int = TRUE, pval = sprintf("p-value = %.0e", logrank$pvalue), pval.coord = c(1150, 0.9), 
                                          xlim = c(0, max(radiomics_df$os)),
                                          legend.labs = c("High","Low"),title=NULL)
#km_plot_radiomics

#base::save.image('./GBM_radiomic-Cox.RData')