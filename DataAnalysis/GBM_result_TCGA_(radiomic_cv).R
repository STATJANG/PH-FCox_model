##### Cox regression with radiomics features.
rm(list=ls())
library(survival)

gbm_dat = readRDS('./GBM_R/Revision(SIM)/gbm_tcga_clinical.rds')
gbm_df = gbm_dat[,c('id','os','status')]
gbm_df$Xc = I(gbm_dat[,c('age', 'gender',  'ntumor','kps','roi_F')])

radio_df = read.csv('./data_BraTS/TCGA_GBM_radiomicFeatures_Train.csv',header=T)
radio_df = rbind(radio_df, read.csv('./data_BraTS/TCGA_GBM_radiomicFeatures_Test.csv',header=T))
radio_df = radio_df[,-which(colnames(radio_df)== "Date")]
radio_df = radio_df[,!apply(is.na(radio_df),2,any)]
colnames(radio_df)[1] = 'id'
radio_df = radio_df %>%
  filter(id %in% gbm_df$id) %>%
  slice(match(gbm_df$id, id))
dim(radio_df)
all(radio_df$id == gbm_dat$id)
gbm_df$Xr = I(radio_df[,-1])


lam_vec = exp(seq(-6,-2,length.out=30)) 

radio_cox = function(ll){
  set.seed(777)
  cv_sample = sample(x = 1:nrow(gbm_dat), size = nrow(gbm_dat), replace = FALSE)
  cv_idx = list(0)
  for(ii in 1:9){ cv_idx[[ii]] = cv_sample[(1:13)+13*(ii-1)]}
  cv_idx[[10]] = cv_sample[(13*9+1):length(cv_sample)]
  
  tmp_risk = NA
  for(cc in 1:length(cv_idx)){
    cat('\n cc is ',cc,'\n')
    
    df_tr = gbm_df[-cv_idx[[cc]],]
    colnames(df_tr$Xc) = c('age', 'gender', 'ntumor', 'kps','roi_F')
    df_te = gbm_df[(cv_idx[[cc]]),]
    colnames(df_te$Xc) = c('age', 'gender', 'ntumor', 'kps','roi_F')
    
    df_te$Xc[,c(1,3,4)] = I(sweep(df_te$Xc[,c(1,3,4)],2,colMeans(df_tr$Xc[,c(1,3,4)]),'-'))
    df_te$Xc[,c(1,3,4)] = I(sweep(df_te$Xc[,c(1,3,4)],2,apply(df_tr$Xc[,c(1,3,4)],2,sd),'/'))
    df_tr$Xc[,c(1,3,4)] = I(scale(df_tr$Xc[,c(1,3,4)], scale = TRUE, center = TRUE))
    
    df_te$Xr = I(sweep(df_te$Xr,2,colMeans(df_tr$Xr),'-'))
    df_te$Xr = I(sweep(df_te$Xr,2,apply(df_tr$Xr,2,sd),'/'))
    df_tr$Xr = I(scale(df_tr$Xr, scale = TRUE, center = TRUE))
    
    
    penal_vec = c(rep(0,ncol(gbm_df$Xc)),rep(1,ncol(gbm_df$Xr)))
    radio.fit = glmnet::glmnet(x=as.matrix(df_tr[,c('Xc','Xr')]),
                               y=Surv(time = df_tr$os, event = df_tr$status),
                               penalty.factor = penal_vec, 
                               standardize = FALSE,
                               family = 'cox', lambda = lam_vec[ll],
                               alpha=1)  ##Lasso
    
    tmp_risk[cv_idx[[cc]]] = as.matrix(df_te[,c('Xc','Xr'),drop=F])%*%coef(radio.fit)
  }
  
  km_df = data.frame(gbm_df[,c('id','os','status')],risk_group = tmp_risk<median(tmp_risk))
  pval_tmp = survdiff(Surv(time=os,event=status)~risk_group, data = km_df)$pvalue
  c_index_tmp = concordance(coxph(Surv(time=os,event=status)~risk_group, data = km_df))$concordance
  
  return(list(pred_pval = pval_tmp, pred_c_index = c_index_tmp,
              pred_df = km_df))
  
}

pval_vec = c_index_vec = NA
for (ll in 1:length(lam_vec)){
  tryCatch(
    expr = {
      pred_tmp = radio_cox(ll)
      pval_vec[ll] = pred_tmp$pred_pval
      c_index_vec[ll] = pred_tmp$pred_c_index
      cat('\n ll is ',ll,'\n')
    },
    error = function(e) {
      message("Here is the official R error message: ", e$message)
    })
}

which.min(pval_vec)
radio_best = radio_cox(which.min(pval_vec))

km_plot_radiomics = survminer::ggsurvplot(survfit(Surv(os, status) ~ risk_group, data = radio_best$pred_df), 
                                          conf.int = TRUE, pval = sprintf("p-value = %.0e", radio_best$pred_pval), pval.coord = c(1150, 0.9), 
                                          xlim = c(0, max(radio_best$pred_df$os)),
                                          legend.labs = c("High","Low"),title=NULL)

km_plot_radiomics
#base::save.image('./GBM_R/Revision(SIM)/GBM_result_TCGA_(radiomic_cv).RData')