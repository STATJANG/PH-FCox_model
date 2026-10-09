rm(list=ls())
library(spatstat)
library(glmnet)
library(survival)
library(irlba)
library(ggplot2)
library(ggthemes)

zip_path = paste0('./GBM_R/Revision(SIM)/loocv_gbm_upenn.zip') #training set
loocv_files <- unzip(zip_path, list = TRUE)[-1,]
p_df = NA
for(ll in 1:nrow(loocv_files)){
  tmp = read.table(unz(zip_path,paste0('loocv_',ll,'.txt')),header=T)
  p_df = rbind(p_df , c(apply(tmp[,1:3],2,min),apply(tmp[,-(1:3)],2,max)))
}
p_df = p_df[-1,]

apply(p_df[,1:3],2,min,na.rm=T)
apply(p_df[,-(1:3)],2,max,na.rm=T)

optimal_sm_vec = c( apply(p_df[,1:3],2,which.min), apply(p_df[,-(1:3)],2,which.max) )

optimal_ll_vec = c(
  which.min(read.table(unz(zip_path,paste0('loocv_',optimal_sm_vec[1],'.txt')),header=T)[,1]),
  which.min(read.table(unz(zip_path,paste0('loocv_',optimal_sm_vec[2],'.txt')),header=T)[,2]),
  which.min(read.table(unz(zip_path,paste0('loocv_',optimal_sm_vec[3],'.txt')),header=T)[,3]),
  which.max(read.table(unz(zip_path,paste0('loocv_',optimal_sm_vec[4],'.txt')),header=T)[,4]),
  which.max(read.table(unz(zip_path,paste0('loocv_',optimal_sm_vec[5],'.txt')),header=T)[,5]),
  which.max(read.table(unz(zip_path,paste0('loocv_',optimal_sm_vec[6],'.txt')),header=T)[,6]))




## Validation 
fcox_validation = function(sm, ll, train_clinical, test_clinical, train_pdg, test_pdg, train_pdg_range){
  
  sm_vector = unname(unlist(sm_par_mat[sm,]))
  
  ps_list = list(0)
  for(dd in 0:2){
    pixnum = train_pdg_range[dd+1,2] - train_pdg_range[dd+1,1]
    ps_grd = seq(train_pdg_range[dd+1,1]+0.5, train_pdg_range[dd+1,2]-0.5, length.out = pixnum)
    ps_grd = cbind( rep(ps_grd,each=pixnum), rep(ps_grd,pixnum) )
    ps_tmp = 0
    
    for(ii in 1:nrow(train_clinical)){
      PD_bd = train_pdg[[ii]]
      PD_bd = PD_bd[PD_bd$dim == dd, c("birth","death")]
      
      weight_tmp = PD_bd$death - PD_bd$birth
      weight_tmp = abs(cbind(weight_tmp, PD_bd$birth, PD_bd$death))
      weight_tmp = apply(weight_tmp,1,max) 
      
      surface_z = density(ppp(PD_bd[,1], PD_bd[,2], train_pdg_range[dd+1,], train_pdg_range[dd+1,]),
                          weights = weight_tmp, sigma = sm_par_mat[sm,dd+1], 
                          dimyx = c(pixnum,pixnum))$v #kernel = "gaussian" (default)
      surface_z = surface_z[ps_grd[,1]<=ps_grd[,2]]
      surface_z[surface_z<0] = 0
      ps_tmp = rbind(ps_tmp,surface_z)
    }
    ps_list[[dd+1]] = ps_tmp[-1,]
    rownames(ps_list[[dd+1]]) = names(train_pdg)
  }
  ps_list[[1]] = ps_list[[1]][match(train_clinical$id, rownames(ps_list[[1]])), ]
  ps_list[[2]] = ps_list[[2]][match(train_clinical$id, rownames(ps_list[[2]])), ]
  ps_list[[3]] = ps_list[[3]][match(train_clinical$id, rownames(ps_list[[3]])), ]
  
  train_df = train_clinical[,c('id','os','status')]
  train_df$Xc = I(train_clinical[,c('age', 'gender', 'ntumor', 'roi_F')])
  
  train_df$Xc[,c('age','ntumor')] = I(as.matrix(scale(train_df$Xc[,c('age','ntumor')], scale = TRUE, center = TRUE)))
  train_df$X0m = I(as.matrix(scale(ps_list[[1]], center=TRUE, scale=FALSE )))
  train_df$X1m = I(as.matrix(scale(ps_list[[2]], center=TRUE, scale=FALSE )))
  train_df$X2m = I(as.matrix(scale(ps_list[[3]], center=TRUE, scale=FALSE )))
  
  #FPC
  num_pc=20
  svd0 = svd(as.matrix(train_df$X0m/sqrt(nrow(train_df$X0m))), nv=num_pc, nu=num_pc)  
  svd1 = svd(as.matrix(train_df$X1m/sqrt(nrow(train_df$X1m))), nv=num_pc, nu=num_pc)
  svd2 = svd(as.matrix(train_df$X2m/sqrt(nrow(train_df$X2m))), nv=num_pc, nu=num_pc)
  eigenval0 = cumsum(svd0$d^2)/sum(svd0$d^2)  
  eigenval1 = cumsum(svd1$d^2)/sum(svd1$d^2) 
  eigenval2 = cumsum(svd2$d^2)/sum(svd2$d^2) 
  hh0 = min(which(eigenval0>0.9))
  hh1 = min(which(eigenval1>0.9))
  hh2 = min(which(eigenval2>0.9))
  eigenfun0 = svd0$v[,1:hh0,drop=F]
  eigenfun1 = svd1$v[,1:hh1,drop=F]
  eigenfun2 = svd2$v[,1:hh2,drop=F]
  
  #'m'=main, 'i'=interaction
  train_df$X0m =  I(as.matrix(train_df$X0m)%*%eigenfun0)
  train_df$X1m =  I(as.matrix(train_df$X1m)%*%eigenfun1)
  train_df$X2m =  I(as.matrix(train_df$X2m)%*%eigenfun2)
  colnames(train_df$X0m) = paste0('dim0.',1:hh0)
  colnames(train_df$X1m) = paste0('dim1.',1:hh1)
  colnames(train_df$X2m) = paste0('dim2.',1:hh2)
  
  train_df$X0i =  I(train_df$X0*matrix(rep(train_clinical$roi_F,hh0),ncol=hh0))
  train_df$X1i =  I(train_df$X1*matrix(rep(train_clinical$roi_F,hh1),ncol=hh1))
  train_df$X2i =  I(train_df$X2*matrix(rep(train_clinical$roi_F,hh2),ncol=hh2))
  colnames(train_df$X0i) = paste0('dim0.',1:hh0,'i')
  colnames(train_df$X1i) = paste0('dim1.',1:hh1,'i')
  colnames(train_df$X2i) = paste0('dim2.',1:hh2,'i')
  
  train_df$X0m = I(as.matrix(scale(train_df$X0m, center=TRUE, scale=TRUE )))
  train_df$X1m = I(as.matrix(scale(train_df$X1m, center=TRUE, scale=TRUE )))
  train_df$X2m = I(as.matrix(scale(train_df$X2m, center=TRUE, scale=TRUE )))
  train_df$X0i = I(as.matrix(scale(train_df$X0i, center=TRUE, scale=TRUE )))
  train_df$X1i = I(as.matrix(scale(train_df$X1i, center=TRUE, scale=TRUE )))
  train_df$X2i = I(as.matrix(scale(train_df$X2i, center=TRUE, scale=TRUE )))
  
  penal_vec = c(rep(0,ncol(train_df$Xc)-1),rep(1,hh0+hh1+hh2))
  lam.fit1 = glmnet(y=Surv(train_df$os, train_df$status),
                    x=as.matrix(train_df[,c('Xc','X0m','X1m','X2m')])[,-4], #remove frontal-lobe indicator
                    family='cox',alpha=1,standardize = FALSE,
                    penalty.factor = penal_vec, #type.measure = 'C',
                    lambda=lam_vec[ll])
  
  penal_vec = c(rep(0,ncol(train_df$Xc)),rep(1,hh0+hh1+hh2))
  lam.fit2 = glmnet(y=Surv(train_df$os, train_df$status),
                    x=as.matrix(train_df[,c('Xc','X0m','X1m','X2m')]),
                    family='cox',alpha=1,standardize = FALSE,
                    penalty.factor = penal_vec, #type.measure = 'C',
                    lambda=lam_vec[ll])
  
  penal_vec = c(rep(0,ncol(train_df$Xc)),rep(1,(hh0+hh1+hh2)*2))
  lam.fit3 = glmnet(y=Surv(train_df$os, train_df$status),
                    x=as.matrix(train_df[,c('Xc','X0m','X0i','X1m','X1i','X2m','X2i')]),
                    family='cox',alpha=1,standardize = FALSE,
                    penalty.factor = penal_vec, #type.measure = 'C',
                    lambda=lam_vec[ll])
  
  
  # Functional coefficients.
  X0m_coef = coef(lam.fit3)[grep("X0m", rownames(coef(lam.fit3)))]      
  X1m_coef = coef(lam.fit3)[grep("X1m", rownames(coef(lam.fit3)))]      
  X2m_coef = coef(lam.fit3)[grep("X2m", rownames(coef(lam.fit3)))] 
  X0i_coef = coef(lam.fit3)[grep("X0i", rownames(coef(lam.fit3)))]      
  X1i_coef = coef(lam.fit3)[grep("X1i", rownames(coef(lam.fit3)))]      
  X2i_coef = coef(lam.fit3)[grep("X2i", rownames(coef(lam.fit3)))]
  coef_list_m = list(X0m=eigenfun0%*%(X0m_coef/svd0$d[1:hh0]),
                     X1m=eigenfun1%*%(X1m_coef/svd1$d[1:hh1]),
                     X2m=eigenfun2%*%(X2m_coef/svd2$d[1:hh2]))
  coef_list_i = list(X0i=eigenfun0%*%(X0i_coef/svd0$d[1:hh0]),
                     X1i=eigenfun1%*%(X1i_coef/svd1$d[1:hh1]),
                     X2i=eigenfun2%*%(X2i_coef/svd2$d[1:hh2]))      
  
  fcoef_list = list(0)      
  for(dd in 0:2){
    fcoef_list[[dd+1]] = vector(mode='list',length = 2)
    names(fcoef_list[[dd+1]]) = c('non_frontal','frontal')
    
    pixnum = train_pdg_range[dd+1,2] - train_pdg_range[dd+1,1]
    ps_grd = seq(train_pdg_range[dd+1,1]+0.5, train_pdg_range[dd+1,2]-0.5, length.out = pixnum)
    ps_grd = cbind( rep(ps_grd,each=pixnum), rep(ps_grd,pixnum) )
    
    coef_df = data.frame(ps_grd[ps_grd[,1]<=ps_grd[,2],],
                         non_frontal = coef_list_m[[dd+1]], 
                         frontal = coef_list_m[[dd+1]]+coef_list_i[[dd+1]] )
    colnames(coef_df) = c('x','y','non_frontal','frontal')
    coef_lim = max(abs(coef_df[,3:4]))
    coef_lim = c(-1,1)*coef_lim  
    
    fcoef_list[[dd+1]]$'non_frontal' = ggplot(coef_df,aes(x,y,fill = non_frontal))+
      geom_raster()+
      theme_tufte()+
      labs(x="birth",y="death",title = NULL,fill=NULL) +
      geom_hline(yintercept=0)+
      geom_vline(xintercept=0)+
      geom_abline(slope=1,intercept=0)+
      scale_fill_gradient2(limits = coef_lim ) +
      theme(legend.position = c(0.9, 0.25))+
      theme(plot.title = element_text(hjust = 0.4))
    
    fcoef_list[[dd+1]]$'frontal' = ggplot(coef_df,aes(x,y,fill = frontal))+
      geom_raster()+
      theme_tufte()+
      labs(x="birth",y="death",title = NULL, fill=NULL) +
      geom_hline(yintercept=0)+
      geom_vline(xintercept=0)+
      geom_abline(slope=1,intercept=0)+
      scale_fill_gradient2(limits = coef_lim ) +
      theme(legend.position = c(0.9, 0.25))+
      theme(plot.title = element_text(hjust = 0.4))
  }      
  
  
  
  
  
  
  # Fit the model to test_clinical
  ps_list = list(0)
  for(dd in 0:2){
    pixnum = train_pdg_range[dd+1,2] - train_pdg_range[dd+1,1]
    ps_grd = seq(train_pdg_range[dd+1,1]+0.5, train_pdg_range[dd+1,2]-0.5, length.out = pixnum)
    ps_grd = cbind( rep(ps_grd,each=pixnum), rep(ps_grd,pixnum) )
    
    ps_tmp = 0
    
    for(ii in 1:nrow(test_clinical)){
      PD_bd = test_pdg[[ii]]
      PD_bd = PD_bd[PD_bd$dim == dd, c("birth","death")]
      PD_bd = PD_bd[PD_bd$birth>= train_pdg_range[dd+1,1],]
      PD_bd = PD_bd[PD_bd$death<= train_pdg_range[dd+1,2],]
      
      weight_tmp = PD_bd$death - PD_bd$birth
      weight_tmp = abs(cbind(weight_tmp, PD_bd$birth, PD_bd$death))
      weight_tmp = apply(weight_tmp,1,max) 
      
      surface_z = density(ppp(PD_bd[,1], PD_bd[,2], train_pdg_range[dd+1,], train_pdg_range[dd+1,]),
                          weights = weight_tmp, sigma = sm_par_mat[sm,dd+1], 
                          dimyx = c(pixnum,pixnum))$v #kernel = "gaussian" (default)
      surface_z = surface_z[ps_grd[,1]<=ps_grd[,2]]
      surface_z[surface_z<0] = 0
      ps_tmp = rbind(ps_tmp,surface_z)
    }
    ps_list[[dd+1]] = ps_tmp[-1,]
    rownames(ps_list[[dd+1]]) = names(test_pdg)
  }
  ps_list[[1]] = ps_list[[1]][match(test_clinical$id, rownames(ps_list[[1]])), ]
  ps_list[[2]] = ps_list[[2]][match(test_clinical$id, rownames(ps_list[[2]])), ]
  ps_list[[3]] = ps_list[[3]][match(test_clinical$id, rownames(ps_list[[3]])), ]
  
  test_df = test_clinical[,c('id','os','status')]
  test_df$Xc = I(test_clinical[,c('age', 'gender', 'ntumor', 'roi_F')])
  test_df$Xc[,c('age','ntumor')] = I(as.matrix(scale(test_df$Xc[,c('age','ntumor')], scale = TRUE, center = TRUE)))
  test_df$X0m = I(as.matrix(scale(ps_list[[1]], center=TRUE, scale=FALSE )))
  test_df$X1m = I(as.matrix(scale(ps_list[[2]], center=TRUE, scale=FALSE )))
  test_df$X2m = I(as.matrix(scale(ps_list[[3]], center=TRUE, scale=FALSE )))
  
  test_df$X0m =  I(as.matrix(test_df$X0)%*%eigenfun0)
  test_df$X1m =  I(as.matrix(test_df$X1)%*%eigenfun1)
  test_df$X2m =  I(as.matrix(test_df$X2)%*%eigenfun2)
  colnames(test_df$X0m) = paste0('dim0.',1:hh0)
  colnames(test_df$X1m) = paste0('dim1.',1:hh1)
  colnames(test_df$X2m) = paste0('dim2.',1:hh2)
  
  test_df$X0i =  I(test_df$X0*matrix(rep(test_clinical$roi_F,hh0),ncol=hh0))
  test_df$X1i =  I(test_df$X1*matrix(rep(test_clinical$roi_F,hh1),ncol=hh1))
  test_df$X2i =  I(test_df$X2*matrix(rep(test_clinical$roi_F,hh2),ncol=hh2))
  colnames(test_df$X0i) = paste0('dim0.',1:hh0,'i')
  colnames(test_df$X1i) = paste0('dim1.',1:hh1,'i')
  colnames(test_df$X2i) = paste0('dim2.',1:hh2,'i')
  
  test_df$X0m = I(as.matrix(scale(test_df$X0m, center=TRUE, scale=TRUE )))
  test_df$X1m = I(as.matrix(scale(test_df$X1m, center=TRUE, scale=TRUE )))
  test_df$X2m = I(as.matrix(scale(test_df$X2m, center=TRUE, scale=TRUE )))
  test_df$X0i = I(as.matrix(scale(test_df$X0i, center=TRUE, scale=TRUE )))
  test_df$X1i = I(as.matrix(scale(test_df$X1i, center=TRUE, scale=TRUE )))
  test_df$X2i = I(as.matrix(scale(test_df$X2i, center=TRUE, scale=TRUE )))
  
  tmp_risk1 = as.matrix(test_df[,c('Xc','X0m','X1m','X2m'),drop=F])[,-4]%*%coef(lam.fit1)
  tmp_risk2 = as.matrix(test_df[,c('Xc','X0m','X1m','X2m'),drop=F])%*%coef(lam.fit2)
  tmp_risk3 = as.matrix(test_df[,c('Xc','X0m','X0i','X1m','X1i','X2m','X2i'),drop=F])%*%coef(lam.fit3)
  
  km_df1 = data.frame(test_df[,c('id','os','status')], est_risk = as.vector(tmp_risk1))
  km_df2 = data.frame(test_df[,c('id','os','status')], est_risk = as.vector(tmp_risk2))
  km_df3 = data.frame(test_df[,c('id','os','status')], est_risk = as.vector(tmp_risk3))
  
  return(list(fit1 = lam.fit1, fit2 = lam.fit2, fit3 = lam.fit3,
              df1 = km_df1, df2 = km_df2, df3 = km_df3,
              fcoef=fcoef_list))
}


sm_par_mat = expand.grid( seq(.3, 3, length.out=10), seq(.3, 3, length.out=10), seq(.3, 3, length.out=10))
lam_vec = exp(seq(-6,-2,length.out=30)) #guessted from cv.glmnet

gbm_upenn_clinical = readRDS('./GBM_R/Revision(SIM)/gbm_upenn_clinical.rds') #train
gbm_tcga_clinical = readRDS('./GBM_R/Revision(SIM)/gbm_tcga_clinical.rds')  #test

env_upenn =  new.env()
env_tcga = new.env()
load('./GBM_R/Revision(SIM)/GBM_pdg_UPENN.RData', envir = env_upenn)
load('./GBM_R/Revision(SIM)/GBM_pdg_TCGA.RData', envir = env_tcga)
pdg_range_common = cbind(
  apply(cbind(env_tcga$pdg_range[,1], env_upenn$pdg_range[,1]),1,min),
  apply(cbind(env_tcga$pdg_range[,2], env_upenn$pdg_range[,2]),1,max))



km3 = fcox_validation(sm = optimal_sm_vec[3], 
                          ll=optimal_ll_vec[3],
                          train_clinical = gbm_upenn_clinical,
                          test_clinical = gbm_tcga_clinical,
                          train_pdg = env_upenn$pdg,
                          test_pdg = env_tcga$pdg,
                          train_pdg_range = pdg_range_common)
km_df3 = km3$df3
km_df3$risk_group = km_df3$est_risk < median(km_df3$est_risk)
cindex3 = concordance(coxph(Surv(time=os,event=status)~est_risk, data = km_df3))
pval3 = survdiff(Surv(time=os,event=status)~risk_group, data = km_df3)$pvalue
km_surv3 = survminer::ggsurvplot(survfit(Surv(time=os,event=status)~risk_group, data = km_df3),
                                    conf.int = TRUE,  pval = sprintf("p-value = %.0e", pval3), 
                                    pval.coord = c(1150, 0.9), 
                                    xlim = c(0, 3000),
                                    #legend.labs = c("Low","High"),title="ROI-specific PS")
                                    legend.labs = c("High","Low"),title=NULL)


##Clinical model (Upenn = train set)
fit.upenn.cox = coxph(Surv(os,status)~age+gender+ntumor+roi_F,data=gbm_upenn_clinical)
risk_vec = predict(fit.upenn.cox, type='lp', newdata = gbm_tcga_clinical)
km_df_tcga = data.frame(gbm_tcga_clinical[,c('id','os','status')],risk_group = risk_vec<median(risk_vec)+0)
survminer::ggsurvplot(survfit(Surv(os,status)~risk_group,data=km_df_tcga)) 
pval_clin = survdiff(Surv(os,status)~risk_group,data=km_df_tcga)$pvalue
cindex_clin = concordance(coxph(Surv(os,status)~risk_group,data=km_df_tcga))

#Save results.
km_surv_clin = survminer::ggsurvplot(survfit(Surv(time=os,event=status)~risk_group, data = km_df_tcga),
                                             conf.int = TRUE,  pval = round(pval_clin,4), 
                                             pval.coord = c(1150, 0.9), 
                                             xlim = c(0, 3000),
                                             legend.labs = c("High","Low"),
                                              title=NULL)
gg_surv <- ggpubr::ggarrange(
  km_surv_clin$plot, km_surv3$plot,
  labels = c("(a)", "(b)"),
  nrow = 1,ncol = 2)
#ggsave('./Revision(SIM)/SIM_GBM_surv(upenn_to_tcga).eps',units='in', height=3.5, width=7.5)

#Flipped cases
print(
km_df_tcga[(km_df_tcga$risk_group==TRUE & km_df3$risk_group==FALSE),'id',drop=F],
row.names=F)

print(
  km_df_tcga[(km_df_tcga$risk_group==FALSE & km_df3$risk_group==TRUE),'id',drop=F],
  row.names=F)



gg_fcoef <- ggpubr::ggarrange(
  km3$fcoef[[1]]$non_frontal, km3$fcoef[[2]]$non_frontal, km3$fcoef[[3]]$non_frontal, 
  km3$fcoef[[1]]$frontal, km3$fcoef[[2]]$frontal, km3$fcoef[[3]]$frontal,
  labels = c("(a)", "(c)", "(e)","(b)", "(d)", "(f)"),
  nrow = 2,ncol = 3)
ggsave('./GBM_R/Revision(SIM)/SIM_GBM_coef(upenn_to_tcga).eps',units='in', 
device = cairo_ps, height=6, width=10,
font.label = list(size = 20, face = "bold"))

#base::save.image('./Revision(SIM)/GBM_result_TCGA(upenn_to_tcga).RData')