require(extraDistr)
require(pracma)
require(caret)  
require(RandomFieldsUtils)
require(gss)
require(diffpriv)

source("ICLP_functions.R")


set.seed(2021)


grid.size = 101
grid = list(pt = seq(0,1,length = grid.size),
            wt = rep(1/grid.size,grid.size))
N.candi = c(seq(100,1000, by  = 300), seq(1500,5000, by  = 500), seq(6000, 10000, by = 1000))
Eps <- 1   # Privacy Budget

method.list <- c("IID_Lap","ICLP_l1","ICLP_RKHS", "Bernstein")
rep.time <- 100


MSE_all = list()

Kernels = c("M3/2", "M5/2")

for(kernel_type in Kernels){
  ##### get the eigenpaires #####
  e_pairs = ICLP_kernel(grid = grid, kernel_type = kernel_type, rho = 1/10)
  e_val = e_pairs$e_val
  e_vec = e_pairs$e_vec
  
  if(kernel_type == "M3/2"){
    h = 2
  }else{
    h = 3
  }
  
  #### setting different mu0 ##############
  mu1 = exp(-grid$pt)*10*grid$pt
  mu2 = 0.3*dnorm(grid$pt, mean = 0.3, sd = .05) + 0.7* dnorm(grid$pt, mean = 0.8, sd = 0.05)
  mu3 = 0.2*(dnorm(grid$pt, mean = 0, sd = 0.03) +  dnorm(grid$pt, mean = 0.2, sd = 0.05) + dnorm(grid$pt, mean = 0.5, sd = 0.05) - dnorm(grid$pt, mean = 0.75, sd = 0.03) + dnorm(grid$pt, mean = 1, sd = 0.03))
  mu4 = rowSums( e_vec[,1:25]%*%diag( runif(n = 25, min = -.5,max = .5) ))
  
  mu0_candi <- matrix(cbind(mu1,mu2,mu3,mu4), nrow = grid.size, ncol = 4)

  
  MSE.total <- list()
  
  for(k in 1:4){
    MSE_plug = Stat_plug = matrix(NA, nrow = length(N.candi),ncol = length(method.list), 
                       dimnames = list(c(N.candi), c(method.list)) )
    MSE_pcv = Stat_pcv = matrix(NA, nrow = length(N.candi),ncol = length(method.list), 
                      dimnames = list(c(N.candi), c(method.list)) )
    
    for(i in 1:length(N.candi)){
      print(i)
      N <- N.candi[i]
      
      mu0 <- mu0_candi[,k]
      X <- generate_X(N = N, grid = grid, e_val = e_val, e_vec = e_vec, tau = 1, mu0 = mu0)
      
      norms <- X_norm(X = X, grid = grid, e_val = e_val, e_vec = e_vec)
      l1_tau <- norms$l1_tau
      l2_tau <- norms$l2_tau
      
      ############ IID ####################
      nocv_M_candi <- sort(ceiling( N^(1/(2*h)) - 1 ):floor( (N)^{1/(3)} ))
      MSE_IID <- rep(0,length(nocv_M_candi))
      for(M in nocv_M_candi){
        iid_lap_res = IID_Lap(X = X, eps = Eps, M = M, grid = grid, e_val = e_val, e_vec = e_vec, l1.bound = F,
                              Rep = 100, ifMSE = T, mu0 = mu0)
        MSE_IID[which(M == nocv_M_candi)] <- iid_lap_res$MSE[3]
      }
      plug_M = nocv_M_candi[which.min(MSE_IID)]
      iid_lap_plug_res = IID_Lap(X = X, eps = Eps, M = plug_M, grid = grid, e_val = e_val, e_vec = e_vec, l1.bound = F,
                                 Rep = 100, ifMSE = T, mu0 = mu0)

      iid_lap_M = pcv_IID_Lap(X = X, eps = Eps, grid = grid, M = seq( max(min(nocv_M_candi)-2,1) , min(max(nocv_M_candi)+5,20), by = 1),
                      fold = 10, e_val = e_val, e_vec = e_vec, l1.bound = F)
      iid_lap_pcv_res = IID_Lap(X = X, eps = Eps, M = iid_lap_M$M, grid = grid, e_val = e_val, e_vec = e_vec, l1.bound = F,
                                 Rep = 100, ifMSE = T, mu0 = mu0)
      
      
      ############ ICLP_l1 ###############
      eta_ll <- 2*(1+1/h)
      lower_psi <- .5*l2_tau*(1/N); upper_psi <- 1*(1/N);
      iclp_l1_psi = pcv_ICLP_l1(X = X, eps = Eps, grid = grid, e_val = e_val, e_vec = e_vec, psi = exp(seq(log(lower_psi),log(upper_psi),length.out = 20)), 
                                fold = 10, J_tau = NA, eta = eta_ll, l1.bound = F)
      iclp_l1_pcv_res <- ICLP_l1(X = X, eps = Eps, grid = grid, e_val = e_val, e_vec = e_vec, psi = iclp_l1_psi$psi, J_tau = NA,
                                 eta = eta_ll, l1.bound = F, ifMSE = T, Rep = 100, mu0 = mu0)
      
      nocv.psi <- l2_tau*(1/N)
      iclp_l1_plug_res <- ICLP_l1(X = X, eps = Eps, grid = grid, e_val = e_val, e_vec = e_vec, psi = nocv.psi, J_tau = NA,
                                  eta = eta_ll, l1.bound = F, ifMSE = T, Rep = 100, mu0 = mu0)
      
      
      ############ RKHS ####################
      eta_rl <- 1 + 1/(2*h)
      lower_psi <- .5*(N*(Eps^2)/l2_tau)^{-eta_rl}; upper_psi <- 2*(1/N)
      iclp_rkhs_psi = pcv_ICLP_RKHS(X = X, eps = Eps, grid = grid, e_val = e_val, e_vec = e_vec, psi = exp(seq(log(lower_psi),log(upper_psi),length.out = 20)),
                                    fold = 10, eta = eta_rl, l1.bound = T)
      iclp_rkhs_pcv_res = ICLP_RKHS(X = X, eps = Eps, grid = grid, e_val = e_val, e_vec = e_vec,
                                    psi = iclp_rkhs_psi$psi, eta = eta_rl, l1.bound = T, ifMSE = T, mu0 = mu0)
      nocv.psi = (N*(Eps^2)/l2_tau)^{-eta_rl}
      iclp_rkhs_plug_res = ICLP_RKHS(X = X, eps = Eps, grid = grid, e_val = e_val, e_vec = e_vec,
                                     psi = nocv.psi, eta = eta_rl, l1.bound = T, ifMSE = T, mu0 = mu0)
      
      ########## bernstein ################
      if(k == 1){
        KKK = 10
      }else{
        KKK = 20
      }
      bern_res = bern_DP(X = X, eps = Eps, N = N, K = KKK, e_val = e_val, e_vec = e_vec, mu0 = mu0)
      
      
      Stat_plug[i,] <- c(iid_lap_plug_res$MSE[1], iclp_l1_plug_res$MSE[1], iclp_rkhs_plug_res$MSE[1], bern_res$MSE[1])
      MSE_plug[i,] <- c(iid_lap_plug_res$MSE[3], iclp_l1_plug_res$MSE[3], iclp_rkhs_plug_res$MSE[3], bern_res$MSE[3])
      Stat_pcv[i,] <- c(iid_lap_pcv_res$MSE[1], iclp_l1_pcv_res$MSE[1], iclp_rkhs_pcv_res$MSE[1], bern_res$MSE[1])
      MSE_pcv[i,] <- c(iid_lap_pcv_res$MSE[3], iclp_l1_pcv_res$MSE[3], iclp_rkhs_pcv_res$MSE[3], bern_res$MSE[3])
    }
    MSE.total[[k]] <- cbind(MSE_plug, MSE_pcv, Stat_plug, Stat_pcv)
  }
  names(MSE.total) <- paste("mu0",1:4,sep = "")
  MSE_all[[kernel_type]] = MSE.total
}


names(MSE_all) = Kernels


save(MSE_all, file = paste("MSE.res",args,".RData",sep=""))










