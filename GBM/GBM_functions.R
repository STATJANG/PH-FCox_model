library(spatstat)
library(tidyverse)
library(tools)
library(survival)
library(survminer)
library(ggthemes)
library(glmnet)

# maximum distance weight (maximum of distance, absolute birth, and absolute death)
mdw=function(pd) {
  mdw=apply(cbind(pd$death-pd$birth,abs(pd$birth),abs(pd$death)),1,max)
  return(mdw)
}


##Here, I need to make a persistent surface value for each patient and for each dimension.
#calculate weights. # f(u) in the notation of Adams et. al
surface_weight = function(pd_input, w_type,C=1,d=1){ 
  #PD = matrix/DataFrame of n by 2.
  #C and d are only needed for pwgk
  if(w_type == "linear") { return( pd_input$death-pd_input$birth)}
  if(w_type == "pwgk") { return(atan(C*(pd_input$death-pd_input$birth)^d) )}
  if(w_type == "mdw") { return(apply(cbind(pd_input$death-pd_input$birth,abs(pd_input$birth),abs(pd_input$death)),1,max) )}
}


persistent_range= function(PD_list, max_dim = 2){
  #PD_list[[patient]] needs to contain (but not necessarily limited to) the three columns (dim, birth, death)
  minmax_table = matrix(NA, ncol=2, nrow = max_dim+1)
  for(dd in 0:max_dim){
    minmax = t(sapply(PD_list,function(pp){
      if(sum(pp$dim==dd)>0){  result = c(min(pp$birth[pp$dim==dd]),max(pp$death[pp$dim==dd])) 
      }else{
        result = c(NA,NA)
      }
      return(result)
    }))
    #minmax[is.na(minmax)] = NA
    minmax_table[dd+1,] = c( floor(min(minmax[,1],na.rm = TRUE)), ceiling(max(minmax[,2],na.rm = TRUE)) )
  }
  return(minmax_table)
}


persistent_surface = function(PD, minmax, par_sd, weight_type, pd_dim){
  #PD is a matrix or data frame that contains (but not necessrily limited to) three columns (dimension, birth, death)
  
  pixnum = minmax[2] - minmax[1]
  PD_grd_x = PD_grd_y = seq(minmax[1]+0.5, minmax[2]-0.5,length.out = pixnum)
  PD_grd = cbind( rep(PD_grd_x,each=pixnum), rep(PD_grd_y,pixnum) )
  
  PD_bd = PD[PD$dim==pd_dim, c("birth","death") ]
  #PD_bd = PD_bd[is.finite(PD_bd[,2]),]
  #PD_bd = PD_bd[is.finite(PD_bd[,1]),]
  weight_tmp = surface_weight(pd_input = PD_bd, w_type = weight_type)
  
  #print(c(dim(PD_bd)[1],length(weight_tmp)))
  surface_z = density(ppp(PD_bd[,1], PD_bd[,2], minmax,minmax),weights = weight_tmp, sigma = par_sd,
                      dimyx = c(pixnum,pixnum))$v #kernel = "gaussian" (default, to be omitted)
  surface_z = surface_z[PD_grd[,1]<=PD_grd[,2]]
  surface_z[surface_z<0] = 0
  return(list(surface = as.vector(surface_z), grid = PD_grd))
}


#surface[[topological dim+1]] & grid[[topological dim +1]]
persistent_surface_list = function(PD_list, minmax_table, sd_vector = c(1,1,1), dim_vector=c(0,1,2), weight_type){ #minmax_table is a matrix
  
  PD_grid = list(0)
  for(dd in 1:nrow(minmax_table)){
    minmax = minmax_table[dd,]
    pixnum = minmax[2] - minmax[1]
    PD_grd_x = PD_grd_y = seq(minmax[1]+0.5, minmax[2]-0.5,length.out = pixnum)
    PD_grid[[dd]] = cbind( rep(PD_grd_x,each=pixnum), rep(PD_grd_y,pixnum) )
  }
  
  surface_result = list(0)
  grid_result = list(0)
  for(pp in dim_vector){
    surface_tmp = 0
    for(pat in 1:length(PD_list)){
      #print(c(pp,pat))
      
      result = persistent_surface(PD_list[[pat]], minmax = minmax_table[pp+1,],
                                  pd_dim = pp, par_sd = sd_vector[pp+1], weight_type = weight_type) 
      surface_tmp = rbind(surface_tmp, result$surface)
    }
    grid_result[[pp+1]] = result$grid
    surface_tmp = surface_tmp[-1,]
    rownames(surface_tmp) = names(PD_list)
    surface_result[[pp+1]] = surface_tmp
  }
  return(list(surface = surface_result, grid = grid_result))
}

E_score = function(X, max_dim=2, trunc_order = NULL){ #X=PS, #trunc_order = c(3,3,3)
  center_ps = list(0)
  eigenfun_ps = list(0)
  eigenval_ps = list(0)
  propEV = list(0) #proportion of explained variance
  eigenscore_ps = list(0)
  for(dd in 0:max_dim){ 
    center_ps[[dd+1]] = scale(X[[dd+1]],center=TRUE,scale=FALSE) 
    svd.tmp = svd(center_ps[[dd+1]])
    eigenfun_ps[[dd+1]] = svd.tmp$v
    eigenval_ps[[dd+1]] = svd.tmp$d^2
    propEV[[dd+1]] = cumsum(eigenval_ps[[dd+1]])/sum(eigenval_ps[[dd+1]])
    eigenscore_ps[[dd+1]] =  center_ps[[dd+1]]%*%eigenfun_ps[[dd+1]]#[,1:trunc_order[dd+1], drop = F]
    if(!is.null(trunc_order)){ eigenscore_ps[[dd+1]] =  center_ps[[dd+1]]%*%eigenfun_ps[[dd+1]][,1:trunc_order[dd+1], drop = F]  }
  }
  return(list(Eigenscore = eigenscore_ps, Eigenvalue = eigenval_ps, Eigenfunc = eigenfun_ps, center_x = center_ps, propEV = propEV ))
}
