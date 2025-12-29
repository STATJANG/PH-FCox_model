#rm(list=ls())


#Load clinical variables
load(paste0('./GBM/GBM_dataset.RData'))



##Clinical variables-only model
loocv.ahr.clinical.only = 0
for(loocv in 1:nrow(gbm_dat)){
  loocv.model = coxph(Surv(time = os, event = status)~gender+age+     ntumor+ks,data = gbm_dat[-loocv,])
  loocv.ahr.clinical.only[loocv] = 
    predict(loocv.model, data.frame(gbm_dat[loocv,c('gender','age','ntumor','ks'),drop=F]),type = 'risk') 
}
gbm_tmp_dat = data.frame( gbm_dat[,c('status','os')], group = (loocv.ahr.clinical.only  < median(loocv.ahr.clinical.only )) )
logrank = survival::survdiff( survival::Surv(time=os,event=status)~group, data = gbm_tmp_dat)
km_plot_clinical = survminer::ggsurvplot(survfit(Surv(os, status) ~ group, data = gbm_tmp_dat), 
                                         conf.int = TRUE, pval = sprintf("p-value = %.0e", logrank$pvalue), pval.coord = c(1150, 0.9), 
                                         xlim = c(0, max(gbm_tmp_dat$os)),
                                         legend.labs = c("High","Low"),title=NULL)







#### SVC setup
svc_roi = c('Frontal')
sm_par_mat = expand.grid( seq(.3, 3, length.out=10), seq(.3, 3, length.out=10), seq(.3, 3, length.out=10))
fpc.max.order = 10

source(paste0('./GBM/GBM_functions.R'))

EVratio_tab_list = vector(mode='list', length = length(svc_roi))  
coef_list = vector(mode='list', length = length(svc_roi))  
lambda_list = vector(mode='list', length = length(svc_roi))
sigma_list = vector(mode='list', length = length(svc_roi))
pred_risk_list = vector(mode='list', length = length(svc_roi))  
fpredictor_list = vector(mode='list', length = length(svc_roi))

names(EVratio_tab_list) = svc_roi
names(coef_list) = svc_roi
names(lambda_list) = svc_roi
names(sigma_list) = svc_roi
names(pred_risk_list) = svc_roi
names(fpredictor_list) = svc_roi


# Load optimized hyper parameters
load('./GBM/GBM_PH_FCox.RData')
# Persistent Surface
load('./GBM/GBM_pd.RData') 

index_F = which.min(prescv[,1])
lambda_F = prescv[index_F,2]
sm_vector_F = unname(unlist(sm_par_mat[index_F,]))
sigma_list$'Frontal' = sm_vector_F


# Persistence surface
gbm_ps_F  = persistent_surface_list(PD_list = gbm_pd, minmax_table = persistent_range(gbm_pd), weight_type='mdw', sd_vector = sm_vector_F )
gbm_fp_F  = E_score(gbm_ps_F$surface , max_dim = 2) 
hh.F  = sapply(gbm_fp_F$propEV  , function(x){ min(min(which(x>0.9)),fpc.max.order) })  


# Functional features
fp_mat_F  = data.frame(gbm_fp_F$Eigenscore[[1]][ ,1:hh.F[1] ] , gbm_fp_F$Eigenscore[[2]][,1:hh.F[2]], gbm_fp_F$Eigenscore[[3]][,1:hh.F[3]])

# Interaction term
fp_mat_F  = data.frame(fp_mat_F , fp_mat_F *as.data.frame( matrix( rep(gbm_dat$'roi_F' ,ncol(fp_mat_F )), ncol = ncol(fp_mat_F ), nrow = nrow(fp_mat_F ))) )

#w = baseline, i=interaction.
colname_F1  = c(apply(cbind( rep(apply(cbind(rep('ps.w',3),0:2),1,paste,collapse=".dim"), hh.F) , c(1:hh.F[1] ,1:hh.F[2] ,1:hh.F[3])  ),1,paste,collapse=".") )
colname_F2  = c(apply(cbind( rep(apply(cbind(rep('ps.i',3),0:2),1,paste,collapse=".dim"), hh.F) , c(1:hh.F[1] ,1:hh.F[2] ,1:hh.F[3])  ),1,paste,collapse=".") )
colname_F  = c(colname_F1 , colname_F2 )
#colname_FT = as.vector(apply(matrix(colname_FT,ncol=1),2,function(x) paste0(tumor_type[tt],'.',x)))
#colname_T  = as.vector(apply(matrix(colname_T,ncol=1 ),2,function(x) paste0(tumor_type[tt],'.',x)))
#colname_F  = as.vector(apply(matrix(colname_F,ncol=1 ),2,function(x) paste0(tumor_type[tt],'.',x)))
colnames(fp_mat_F)  = colname_F
fpredictor_list$'Frontal'         = fp_mat_F


fp_mat_F$id  = rownames(fp_mat_F )
tmp_dat_F  = dplyr::left_join(gbm_dat[,c('id','status','os','gender','age','ntumor','ks','roi_F' )], fp_mat_F , by=c('id'='id'))




#########   Penalized Cox   ########
glmnet.fcox.fit_F = glmnet::glmnet(x=as.matrix(tmp_dat_F[,-(1:3)]),
                                   y=Surv(time = tmp_dat_F$os, event = tmp_dat_F$status),
                                   penalty.factor = c(0,0,0,0,0,rep(1,ncol(tmp_dat_F)-8)),
                                   family = 'cox', alpha=1, standardize=FALSE,
                                   lambda = lambda_F)


coef_list$'Frontal'         = coef(glmnet.fcox.fit_F )

lambda_list$'Frontal'         = lambda_F



##   Functional coefficients
fcoef.F  = vector(mode = 'list',length=3)
names(fcoef.F) = c('dim0','dim1','dim2')

fcoef.F.df  = vector(mode = 'list',length=3)
names(fcoef.F.df) = c('dim0','dim1','dim2')

for( dd in 0:2){  # topological dim
  fcoef.F[[dd+1]]  = vector(mode = 'list',length=3)
  names(fcoef.F[[dd+1]]) = c('roi', 'non_roi', 'diff')
  
  fcoef.F.df[[dd+1]]  = vector(mode = 'list',length=3)
  names(fcoef.F.df[[dd+1]]) = c('roi', 'non_roi' , 'diff')
  
  fcoef.F.df[[dd+1]]$roi  = data.frame(x= gbm_ps_F$grid[[dd+1]][,1] , y = gbm_ps_F$grid[[dd+1]][,2], z = 0 )
  fcoef.F.df[[dd+1]]$non_roi  = data.frame(x= gbm_ps_F$grid[[dd+1]][,1] , y = gbm_ps_F$grid[[dd+1]][,2], z = 0 )
  fcoef.F.df[[dd+1]]$diff  = data.frame(x= gbm_ps_F$grid[[dd+1]][,1] , y = gbm_ps_F$grid[[dd+1]][,2], z = 0 )
  
  tmp_colname_F.w  = apply(cbind( rep(paste0('ps.w.dim',dd),hh.F[dd+1]) , 1:hh.F[dd+1]  ), 1, paste0,collapse = '.')
  tmp_colname_F.i  = apply(cbind( rep(paste0('ps.i.dim',dd),hh.F[dd+1]) , 1:hh.F[dd+1]  ), 1, paste0,collapse = '.')
  
  coef.w.F = as.vector(gbm_fp_F$Eigenfunc[[dd+1]][ ,1:hh.F[dd+1]]  %*% coef(glmnet.fcox.fit_F) [tmp_colname_F.w,,drop=F]  )
  coef.i.F = as.vector(gbm_fp_F$Eigenfunc[[dd+1]][ ,1:hh.F[dd+1]]  %*% coef(glmnet.fcox.fit_F) [tmp_colname_F.i,,drop=F]  )
  
  fcoef.F.df[[dd+1]]$non_roi$z[fcoef.F.df[[dd+1]]$non_roi$x<=fcoef.F.df[[dd+1]]$non_roi$y]    = coef.w.F
  
  fcoef.F.df[[dd+1]]$roi$z[fcoef.F.df[[dd+1]]$roi$x<=fcoef.F.df[[dd+1]]$roi$y]    = coef.w.F + coef.i.F
  
  fcoef.F.df[[dd+1]]$diff$z[fcoef.F.df[[dd+1]]$diff$x<=fcoef.F.df[[dd+1]]$diff$y]    = coef.i.F
}

lim.F  = max(sapply(fcoef.F.df,function(x){ max( c(max(abs(x$roi$z)),max(abs(x$non_roi$z)),max(abs(x$diff$z)))  )}  ) )

for( dd in 0:2){  # topological dim
  
  
  fcoef.F[[paste0('dim',dd)]]$roi<- ggplot(fcoef.F.df[[dd+1]]$roi,aes(x,y,fill=z))+
    geom_raster()+
    theme_tufte()+
    labs(x="birth",y="death",title = NULL) +
    geom_hline(yintercept=0)+
    geom_vline(xintercept=0)+
    geom_abline(slope=1,intercept=0)+
    scale_fill_gradient2(limits=c( -lim.F , lim.F ) ) +
    theme(legend.position = c(0.9, 0.25))+
    theme(plot.title = element_text(hjust = 0.4))
  
  
  
  fcoef.F[[paste0('dim',dd)]]$non_roi<- ggplot(fcoef.F.df[[dd+1]]$non_roi,aes(x,y,fill=z))+
    geom_raster()+
    theme_tufte()+
    labs(x="birth",y="death",title = NULL) +
    geom_hline(yintercept=0)+
    geom_vline(xintercept=0)+
    geom_abline(slope=1,intercept=0)+
    scale_fill_gradient2(limits=c( -lim.F , lim.F ) ) +
    theme(legend.position = c(0.9, 0.25))+
    theme(plot.title = element_text(hjust = 0.4))
  
  fcoef.F[[paste0('dim',dd)]]$diff<- ggplot(fcoef.F.df[[dd+1]]$diff,aes(x,y,fill=z))+
    geom_raster()+
    theme_tufte()+
    labs(x="birth",y="death",title = NULL) +
    geom_hline(yintercept=0)+
    geom_vline(xintercept=0)+
    geom_abline(slope=1,intercept=0)+
    scale_fill_gradient2(limits=c( -lim.F , lim.F ) ) +
    theme(legend.position = c(0.9, 0.25))+
    theme(plot.title = element_text(hjust = 0.4))
}

EVratio_tab_list$'Frontal'= t(sapply(gbm_fp_F$Eigenvalue, function(x) cumsum(x)[1:10]/sum(x)))


########## LOOCV for evaluating predictive risk  ##########
tmp.risk_FT = tmp.risk_T = tmp.risk_F = 0
for(ll in 1:nrow(gbm_dat)){
  
  #FPC order for each dim (threshold = 90%, upper bound = 10)
  ###    Entire PS
  gbm_ps_F.train  = gbm_ps_F
  
  gbm_ps_F.test = vector(mode='list',length = 3)
  for(dd in 1:3){
    gbm_ps_F.train$surface[[dd]]   = gbm_ps_F$surface[[dd]][-ll,]
    
    gbm_ps_F.test[[dd]]   = gbm_ps_F$surface[[dd]][ll,,drop=F]  - colMeans(gbm_ps_F$surface[[dd]][-ll,])
  }
  gbm_fp_F.train  = E_score(gbm_ps_F.train$surface, max_dim = 2)
  hh.F.cv  =  sapply(gbm_fp_F.train$propEV, function(x){ min(min(which(x>0.9)),fpc.max.order) })
  
  #functional feature (test)
  for(dd in 1:3){
    gbm_ps_F.test[[dd]]  = gbm_ps_F.test[[dd]]  %*% gbm_fp_F.train$Eigenfunc[[dd]][,1:hh.F.cv[dd]]
  }
  
  #functional features (train)
  fp_mat_train_F = data.frame(gbm_fp_F.train$Eigenscore[[1]][,1:hh.F.cv[1]], gbm_fp_F.train$Eigenscore[[2]][,1:hh.F.cv[2]], gbm_fp_F.train$Eigenscore[[3]][,1:hh.F.cv[3]])
  fp_mat_train_F  = data.frame(fp_mat_train_F, fp_mat_train_F*as.data.frame( matrix( rep(gbm_dat$'roi_F'[-ll] ,ncol(fp_mat_train_F)), ncol = ncol(fp_mat_train_F), nrow = nrow(fp_mat_train_F))) )
  
  fp_mat_train_F$id  = rownames(fp_mat_train_F)
  
  #  functional feature data frame
  train_dat_F  = dplyr::right_join(gbm_dat[,c('id','status','os','gender','age','ntumor','ks','roi_F')],  fp_mat_train_F , by=c('id'='id'))
  
  
  
  loocv.fit_F = glmnet::glmnet(x=as.matrix(train_dat_F[,-(1:3)]),
                               y=Surv(time = train_dat_F$os, event = train_dat_F$status),
                               penalty.factor = c(0,0,0,0,0,rep(1,ncol(train_dat_F)-8)),
                               #penalty.factor = c(0,0,0,0,rep(1,ncol(train_dat_F)-7)),
                               family = 'cox',
                               standardize = FALSE,
                               lambda = lambda_F,
                               alpha=1  )
  
  test_dat_F   = cbind( gbm_dat[ll, c('gender','age','ntumor','ks','roi_F' ),drop=F] ,matrix(c(unlist(lapply(gbm_ps_F.test, unlist)) ,unlist(lapply(gbm_ps_F.test , unlist))*gbm_dat$'roi_F'[ll] ), nrow=1))
  
  tmp.risk_F[ll]   = predict(loocv.fit_F ,  newx = as.matrix(test_dat_F) ,type = 'link') 
}

dat_km_F  = data.frame(tmp_dat_F[,c('os','status')] ,group = tmp.risk_F  < median(tmp.risk_F ))

loocv.pvalue_F  = survival::survdiff(Surv(os, status) ~ group, data = dat_km_F)$pvalue


km_plot_F  = survminer::ggsurvplot(survfit(Surv(os, status) ~ group, data = dat_km_F), 
                                   conf.int = TRUE,  pval = sprintf("p-value = %.0e", loocv.pvalue_F), pval.coord = c(1150, 0.9), 
                                   xlim = c(0, max(dat_km_F$os)),
                                   legend.labs = c("High","Low"),title=NULL)

pred_risk_list$'Frontal' = tmp.risk_F

survplot_list = list(0)
survplot_list$'Frontal' = km_plot_F
survplot_list$entire = km_plot_clinical

pvalue_list = c(loocv.pvalue_F)
names(pvalue_list) = c('Frontal')
group_tab_list = matrix(c(unname(table(dat_km_F$group))), ncol=2, byrow = T)
colnames(group_tab_list) = c('high_risk','low_risk')
rownames(group_tab_list) = c('Frontal')

#save.image(paste0('./GBM/GBM_PH-FCox_fitting.RData'))