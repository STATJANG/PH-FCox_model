rm(list=ls())
simul_size = 300
n_size = 70 # grp_A 70, grp_B 70, a total of 140

## Simulation under  constraint of effects model
rm(list=ls())
load('./GBM_R/gbm_simul/gbm_simul_simuldf.RData')
load('./GBM_R/gbm_simul/gbm_simul_fp.RData')

store_coef_cox = store_coef_cox2 = matrix(0,ncol=3,nrow=simul_size) 
store_coef_Fcox = store_coef_Fcox2 = vector(mode='list',length = simul_size)
f_ind = rep( c(rep(1, round(n_size*0.3) ),rep(0,n_size - round(n_size*0.3) )), 2) #frontal lobe

set.seed(700)

for(simul in 1:simul_size){
  simul_dat_f = data.frame(f = (f_ind+0)-(!f_ind+0) , fp_mat[[simul]] )
  
  simul.Fcox.cv  = glmnet::cv.glmnet( x = as.matrix(simul_dat_f), 
                                      y = Surv(time = simul_df[[simul]]$os, event = simul_df[[simul]]$status),
                                      penalty.factor = c(0, rep(1,ncol(fp_mat[[simul]])) ),
                                      family = 'cox', nfolds = 10, standardize = FALSE, alpha=1) 
  simul.Fcox2.cv = glmnet::cv.glmnet( x = as.matrix(simul_dat_f), 
                                      y = Surv(time = simul_df[[simul]]$os.null, event = simul_df[[simul]]$status ),
                                      penalty.factor = c(0, rep(1,ncol(fp_mat[[simul]])) ),
                                      family = 'cox', nfolds = 10, standardize = FALSE, alpha=1) 
  
  simul.Fcox  = glmnet::glmnet( x = as.matrix(simul_dat_f), 
                                y = Surv(time = simul_df[[simul]]$os, event = simul_df[[simul]]$status),
                                penalty.factor = c(0, rep(1,ncol(fp_mat[[simul]])) ),
                                family = 'cox', standardize = FALSE, alpha=1, lambda = simul.Fcox.cv$lambda.min) 
  simul.Fcox2 = glmnet::glmnet( x = as.matrix(simul_dat_f), 
                                y = Surv(time = simul_df[[simul]]$os.null, event = simul_df[[simul]]$status ),
                                penalty.factor = c(0, rep(1,ncol(fp_mat[[simul]])) ),
                                family = 'cox', standardize = FALSE, alpha=1, lambda = simul.Fcox2.cv$lambda.min) 
  
  simul.cox = coxph(Surv(time = os, event=status)~f.B+nf.B+f.A ,data= simul_df[[simul]])
  simul.cox2 = coxph(Surv(time = os.null, event=status)~f.B+nf.B+f.A, data= simul_df[[simul]])
  
  store_coef_Fcox[[simul]] = coef(simul.Fcox)
  store_coef_Fcox2[[simul]] = coef(simul.Fcox2)
  store_coef_cox[simul,] = coef(simul.cox)
  store_coef_cox2[simul,] = coef(simul.cox2)
  if (simul%%20==0) cat('\n','simul=',simul)
}
colMeans(store_coef_cox)
colMeans(store_coef_cox2)
save(simul_df, simul_dat_f, store_coef_Fcox, store_coef_Fcox2, store_coef_cox, store_coef_cox2,file = './GBM_R/gbm_simul/gbm_simul1_CornerStone.RData')
