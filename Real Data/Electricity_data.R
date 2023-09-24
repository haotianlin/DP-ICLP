require(extraDistr)
require(pracma)
require(caret)  
require(RandomFieldsUtils)
require(gss)
require(diffpriv)
require(fda)


source("ICLP_functions.R")


args <-  as.numeric(commandArgs(trailingOnly=TRUE))
set.seed(2021*args)

Data <- mondaydemand$y
Data_s <- stand_curve(Data = Data, grid = NULL)
Data <- Data_s
N.size <- dim(Data)[2]
grid.size <- dim(Data)[1]



grid <- list(pt = seq(0.001,0.999,length.out = dim(Data)[1]),
             wt = rep(1/grid.size,grid.size))


eps_candi <- c(1/8, 1/4, 1/2, 1, 2, 4)


All_res = list()


for(kernel_type in c("M3/2", "M5/2")){
  ##### get the eigenpaires #####
  e_pairs = ICLP_kernel(grid = grid, kernel_type = kernel_type, rho = 1/10)
  e_val = e_pairs$e_val
  e_vec = e_pairs$e_vec
  
  if(kernel_type == "M3/2"){
    h = 2
  }else{
    h = 3
  }
  
  norms <- X_norm(X = Data, grid = grid, e_val = e_val, e_vec = e_vec)
  l1_tau <- norms$l1_tau
  l2_tau <- norms$l2_tau
  mean_model <- mean_est(t(Data))
  sample_mean <- mean_model(grid$pt)
  
  
  MSE_Eletricity <- matrix(NA, nrow = length(eps_candi), ncol = 4,
                           dimnames = list(c(eps_candi), c("iidLaplace","bernstein", "Soft","RKHS")))
  MU_tilde = list()
  
  for(i in 1:length(eps_candi)){
    eps <- eps_candi[i]
    ############### i.i.d. Laplace ###############
    h <- 2
    nocv_M_candi <- ceiling( N.size^(1/(2*h)) - 1 ):floor( (N.size)^{1/(3)} )
    MSE_IID <- rep(0,length(nocv_M_candi))
    for(nocv.M1 in nocv_M_candi){
      TMIID_res_plug = IID_Lap(X = Data, eps = eps, M = nocv.M1, grid = grid, e_val = e_val, e_vec = e_vec, l1.bound = F,
                               Rep = 100, ifMSE = T, mu0 = sample_mean)
      MSE_IID[which(nocv.M1 == nocv_M_candi)] <- TMIID_res_plug$MSE[2]
      # print(TMIID_res_plug$MSE)
    }
    nocv.M1 <- nocv_M_candi[which.min(MSE_IID)] 
    TMIID_res_plug = IID_Lap(X = Data, eps = eps, M = nocv.M1, grid = grid, e_val = e_val, e_vec = e_vec, l1.bound = F,
                             Rep = 100, ifMSE = T, mu0 = sample_mean)
    MSE_Eletricity[i,1] <- TMIID_res_plug$MSE[3]
    
    ########## bernstein ################
    bern_res = bern_DP(X = Data, eps = eps, N = N.size, K = 20, e_val = e_val, e_vec = e_vec, mu0 = sample_mean)
    MSE_Eletricity[i,2] <- bern_res$MSE[3]
    
    
    ############ ICLP_l1 ###############
    eta_ll <- 2*(1+1/h)
    nocv.psi <- l2_tau*(1/N.size)
    TMAUTO_res_plug <- ICLP_l1(X = Data, eps = eps, grid = grid, e_val = e_val, e_vec = e_vec, psi = nocv.psi, J_tau = NA,
                               eta = eta_ll, l1.bound = F, ifMSE = T, Rep = 100, mu0 = sample_mean)
    MSE_Eletricity[i,3] <- TMAUTO_res_plug$MSE[3]
    
    
    
    ############ RKHS ####################
    eta_rl <- 1 + 1/(2*h)
    nocv.psi = 2*(N.size*(1)/l2_tau)^{-eta_rl}
    RKHS_res_plug = ICLP_RKHS(X = Data, eps = eps, grid = grid, e_val = e_val, e_vec = e_vec,
                              psi = nocv.psi, eta = eta_rl, l1.bound = T, ifMSE = T, mu0 = sample_mean)
    MSE_Eletricity[i,4] <- RKHS_res_plug$MSE[3]
    
    
    MU_tilde[[i]] = cbind(TMIID_res_plug$mu_tilde_mean, bern_res$mu_tilde_mean,
                          TMAUTO_res_plug$mu_tilde_mean, RKHS_res_plug$mu_tilde_mean)
    
  }
  
  All_res[[kernel_type]] = list(MSE = MSE_Eletricity,
                                MU_tilde = MU_tilde)
  
}


save(All_res, file = paste("All_res",args,".RData",sep=""))

