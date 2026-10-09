##### Cox regression with SECT features.
rm(list=ls())
library(survival)

gbm_dat = readRDS('./GBM_R/Revision(SIM)/gbm_tcga_clinical.rds')
gbm_df = gbm_dat[,c('id','os','status')]
gbm_df$Xc = I(gbm_dat[,c('age', 'gender',  'ntumor','kps','roi_F')])

####Load ec_curve
tmp_env = new.env()
#ec_curve_list = list(NULL,NULL,NULL) #1, 2, 4 
ec_curve1 = ec_curve2 = ec_curve4 = NULL
for(ee in 1:nrow(gbm_dat)){
  load(paste0('./GBM_ec_curve/',gbm_dat$id[ee],'_mni152_tumor.nii/ec_curve_data_SECT_1.RData'), envir = tmp_env)
  ec_curve1 = rbind(ec_curve1, tmp_env$ec_curve_data)
  
  load(paste0('./GBM_ec_curve/',gbm_dat$id[ee],'_mni152_tumor.nii/ec_curve_data_SECT_2.RData'), envir = tmp_env)
  ec_curve2 = rbind(ec_curve2, tmp_env$ec_curve_data)
  
  load(paste0('./GBM_ec_curve/',gbm_dat$id[ee],'_mni152_tumor.nii/ec_curve_data_SECT_4.RData'), envir = tmp_env)
  ec_curve4 = rbind(ec_curve4, tmp_env$ec_curve_data)
  if(ee%%10==0) print(ee)
}
rownames(ec_curve1) = rownames(ec_curve2) = rownames(ec_curve4) = gbm_dat$id

ps_list = vector(mode='list', length=3)
ps_list[[1]] = ec_curve1[match(gbm_dat$id, rownames(ec_curve1)), ]
ps_list[[2]] = ec_curve2[match(gbm_dat$id, rownames(ec_curve2)), ]
ps_list[[3]] = ec_curve4[match(gbm_dat$id, rownames(ec_curve4)), ]

sect_cox = function(ll){
  
  set.seed(777)
  cv_sample = sample(x = 1:nrow(gbm_dat), size = nrow(gbm_dat), replace = FALSE)
  cv_idx = list(0)
  for(ii in 1:9){ cv_idx[[ii]] = cv_sample[(1:13)+13*(ii-1)]}
  cv_idx[[10]] = cv_sample[(13*9+1):length(cv_sample)]
  lam_vec = exp(seq(-6,-2,length.out=30))
  
  
  cat('\n ll = ',ll , '\n')
  
  tmp_risk = NA
  
  ##Cross-validation
  for(cc in 1:length(cv_idx)){
    df_tr = gbm_df[-cv_idx[[cc]],]
    colnames(df_tr$Xc) = c('age', 'gender', 'ntumor', 'kps', 'roi_F')
    df_te = gbm_df[(cv_idx[[cc]]),]
    colnames(df_te$Xc) = c('age', 'gender', 'ntumor', 'kps', 'roi_F')
    
    df_te$Xc[,c(1,3,4)] = I(sweep(df_te$Xc[,c(1,3,4)],2,colMeans(df_tr$Xc[,c(1,3,4)]),'-'))
    df_te$Xc[,c(1,3,4)] = I(sweep(df_te$Xc[,c(1,3,4)],2,apply(df_tr$Xc[,c(1,3,4)],2,sd),'/'))
    df_tr$Xc[,c(1,3,4)] = I(scale(df_tr$Xc[,c(1,3,4)], scale = TRUE, center = TRUE))
    
    df_tr$X0 = I(ps_list[[1]][-cv_idx[[cc]],])
    df_tr$X1 = I(ps_list[[2]][-cv_idx[[cc]],])
    df_tr$X2 = I(ps_list[[3]][-cv_idx[[cc]],])
    df_te$X0 = I(ps_list[[1]][(cv_idx[[cc]]),,drop=F])
    df_te$X1 = I(ps_list[[2]][(cv_idx[[cc]]),,drop=F])
    df_te$X2 = I(ps_list[[3]][(cv_idx[[cc]]),,drop=F])
    
    #Centering for FPC
    df_te$X0 = I(sweep(df_te$X0,2,colMeans(df_tr$X0),'-'))
    df_te$X1 = I(sweep(df_te$X1,2,colMeans(df_tr$X1),'-'))
    df_te$X2 = I(sweep(df_te$X2,2,colMeans(df_tr$X2),'-'))
    
    df_tr$X0 = I(scale(df_tr$X0, center=TRUE, scale=FALSE ))
    df_tr$X1 = I(scale(df_tr$X1, center=TRUE, scale=FALSE ))
    df_tr$X2 = I(scale(df_tr$X2, center=TRUE, scale=FALSE ))
    
    #FPC
    num_pc=20
    svd0 = svd(as.matrix(df_tr$X0/sqrt(nrow(df_tr$X0))), nv=num_pc, nu=num_pc)  
    svd1 = svd(as.matrix(df_tr$X1/sqrt(nrow(df_tr$X1))), nv=num_pc, nu=num_pc)
    svd2 = svd(as.matrix(df_tr$X2/sqrt(nrow(df_tr$X2))), nv=num_pc, nu=num_pc)
    eigenval0 = cumsum(svd0$d^2)/sum(svd0$d^2)  
    eigenval1 = cumsum(svd1$d^2)/sum(svd1$d^2) 
    eigenval2 = cumsum(svd2$d^2)/sum(svd2$d^2) 
    hh0 = min(min(which(eigenval0>0.9)),num_pc)
    hh1 = min(min(which(eigenval1>0.9)),num_pc)
    hh2 = min(min(which(eigenval2>0.9)),num_pc)
    eigenfun0 = svd0$v[,1:hh0]
    eigenfun1 = svd1$v[,1:hh1]
    eigenfun2 = svd2$v[,1:hh2]
    
    df_tr$X0 =  I(as.matrix(df_tr$X0)%*%eigenfun0)
    df_tr$X1 =  I(as.matrix(df_tr$X1)%*%eigenfun1)
    df_tr$X2 =  I(as.matrix(df_tr$X2)%*%eigenfun2)
    colnames(df_tr$X0) = paste0('dim0.',1:hh0)
    colnames(df_tr$X1) = paste0('dim1.',1:hh1)
    colnames(df_tr$X2) = paste0('dim2.',1:hh2)
    
    df_tr$X0i =  I(df_tr$X0*matrix(rep(gbm_dat$roi_F[-cv_idx[[cc]]],hh0),ncol=hh0))
    df_tr$X1i =  I(df_tr$X1*matrix(rep(gbm_dat$roi_F[-cv_idx[[cc]]],hh1),ncol=hh1))
    df_tr$X2i =  I(df_tr$X2*matrix(rep(gbm_dat$roi_F[-cv_idx[[cc]]],hh2),ncol=hh2))
    colnames(df_tr$X0i) = paste0('dim0.',1:hh0,'i')
    colnames(df_tr$X1i) = paste0('dim1.',1:hh1,'i')
    colnames(df_tr$X2i) = paste0('dim2.',1:hh2,'i')
    
    df_te$X0 =  I(as.matrix(df_te$X0)%*%eigenfun0)
    df_te$X1 =  I(as.matrix(df_te$X1)%*%eigenfun1)
    df_te$X2 =  I(as.matrix(df_te$X2)%*%eigenfun2)
    colnames(df_te$X0) = paste0('dim0.',1:hh0)
    colnames(df_te$X1) = paste0('dim1.',1:hh1)
    colnames(df_te$X2) = paste0('dim2.',1:hh2)
    
    df_te$X0i =  I(df_te$X0*matrix(rep(gbm_dat$roi_F[(cv_idx[[cc]])],hh0),ncol=hh0))
    df_te$X1i =  I(df_te$X1*matrix(rep(gbm_dat$roi_F[(cv_idx[[cc]])],hh1),ncol=hh1))
    df_te$X2i =  I(df_te$X2*matrix(rep(gbm_dat$roi_F[(cv_idx[[cc]])],hh2),ncol=hh2))
    colnames(df_te$X0i) = paste0('dim0.',1:hh0,'i')
    colnames(df_te$X1i) = paste0('dim1.',1:hh1,'i')
    colnames(df_te$X2i) = paste0('dim2.',1:hh2,'i')
    
    #Standardization
    df_te$X0 = I(sweep(df_te$X0,2,colMeans(df_tr$X0),'-'))
    df_te$X1 = I(sweep(df_te$X1,2,colMeans(df_tr$X1),'-'))
    df_te$X2 = I(sweep(df_te$X2,2,colMeans(df_tr$X2),'-'))
    df_te$X0 = I(sweep(df_te$X0,2,apply(df_tr$X0,2,sd),'/'))
    df_te$X1 = I(sweep(df_te$X1,2,apply(df_tr$X1,2,sd),'/'))
    df_te$X2 = I(sweep(df_te$X2,2,apply(df_tr$X2,2,sd),'/'))
    
    df_te$X0i = I(sweep(df_te$X0i,2,colMeans(df_tr$X0i),'-'))
    df_te$X1i = I(sweep(df_te$X1i,2,colMeans(df_tr$X1i),'-'))
    df_te$X2i = I(sweep(df_te$X2i,2,colMeans(df_tr$X2i),'-'))
    df_te$X0i = I(sweep(df_te$X0i,2,apply(df_tr$X0i,2,sd),'/'))
    df_te$X1i = I(sweep(df_te$X1i,2,apply(df_tr$X1i,2,sd),'/'))
    df_te$X2i = I(sweep(df_te$X2i,2,apply(df_tr$X2i,2,sd),'/'))
    
    
    df_tr$X0 = I(scale(df_tr$X0, center=TRUE, scale=TRUE ))
    df_tr$X1 = I(scale(df_tr$X1, center=TRUE, scale=TRUE ))
    df_tr$X2 = I(scale(df_tr$X2, center=TRUE, scale=TRUE ))
    df_tr$X0i = I(scale(df_tr$X0i, center=TRUE, scale=TRUE ))
    df_tr$X1i = I(scale(df_tr$X1i, center=TRUE, scale=TRUE ))
    df_tr$X2i = I(scale(df_tr$X2i, center=TRUE, scale=TRUE ))
    
    penal_vec = c(rep(0,ncol(df_tr$Xc)-1),rep(1,hh0+hh1+hh2))
    lam.fit1 = glmnet(y=Surv(df_tr$os, df_tr$status),
                      x=as.matrix(df_tr[,c('Xc','X0','X1','X2')])[,-c(5)], #age,gender,ntumor
                      family='cox',alpha=1,standardize = FALSE,
                      penalty.factor = penal_vec, #type.measure = 'C',
                      lambda=lam_vec[ll])
    tmp_risk[(cv_idx[[cc]])] = as.matrix(df_te[,c('Xc','X0','X1','X2'),drop=F])[,-c(4),drop=F]%*%coef(lam.fit1)
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
      pred_tmp = sect_cox(ll)
      pval_vec[ll] = pred_tmp$pred_pval
      c_index_vec[ll] = pred_tmp$pred_c_index
      cat('\n ll is ',ll,'\n')
    },
    error = function(e) {
      message("Here is the official R error message: ", e$message)
    })
}


which.min(pval_vec)
sect_best = sect_cox(which.min(pval_vec))
km_plot_sect = survminer::ggsurvplot(survfit(Surv(os, status) ~ risk_group, data = sect_best$pred_df), 
                                          conf.int = TRUE, pval = sprintf("p-value = %.0e", sect_best$pred_pval), pval.coord = c(1150, 0.9), 
                                          xlim = c(0, max(sect_best$pred_df$os)),
                                          legend.labs = c("High","Low"),title=NULL)

km_plot_sect
#base::save.image('./GBM_R/Revision(SIM)/GBM_result_TCGA_(sect_cv).RData')