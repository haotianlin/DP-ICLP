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
N = 1000
Eps_candi = c(1/32, 1/16, 1/8, 1/4, 1/2, 1, 2, 3, 4)


method.list <- c("IID_Lap","ICLP_l1","ICLP_RKHS", "Bernstein")
rep.time <- 100

MSE_all = list()



# for(kernel_type in c("M3/2", "M5/2")){
for(kernel_type in c("M3/2")){
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
  
  MSE.total = MU_tilde = list()
  
  for(k in 1:4){
    error_plug <- array(data = 0, dim = c(length(Eps_candi), length(method.list), 3),
                        dimnames = list(c(Eps_candi), c(method.list), c("Statistical Error","Privacy Error","MSE"))  )
    for(i in 1:length(Eps_candi)){
      print(i)
      Eps <- Eps_candi[i]
      
      mu0 <- mu0_candi[,k]
      X <- generate_X(N = N, grid = grid, e_val = e_val, e_vec = e_vec, tau = 4, mu0 = mu0)
      
      norms <- X_norm(X = X, grid = grid, e_val = e_val, e_vec = e_vec)
      l1_tau <- norms$l1_tau
      l2_tau <- norms$l2_tau
      
      ############ IID ####################
      plug_M = floor( (N)^{1/(3)} )
      iid_lap_plug_res = IID_Lap(X = X, eps = Eps, M = plug_M, grid = grid, e_val = e_val, e_vec = e_vec, l1.bound = F,
                                 Rep = 100, ifMSE = T, mu0 = mu0)
      
      ############ ICLP_l1 ###############
      eta_ll <- 2*(1+1/h)
      nocv.psi <- l2_tau*(1/N)
      iclp_l1_plug_res <- ICLP_l1(X = X, eps = Eps, grid = grid, e_val = e_val, e_vec = e_vec, psi = nocv.psi, J_tau = NA,
                                  eta = eta_ll, l1.bound = F, ifMSE = T, Rep = 100, mu0 = mu0)
      
      
      ############ RKHS ####################
      eta_rl <- 1 + 1/(2*h)
      nocv.psi = (N*(1)/l2_tau)^{-eta_rl}
      iclp_rkhs_plug_res = ICLP_RKHS(X = X, eps = Eps, grid = grid, e_val = e_val, e_vec = e_vec,
                                     psi = nocv.psi, eta = eta_rl, l1.bound = T, ifMSE = T, mu0 = mu0)
      
      ########## bernstein ################
      bern_res = bern_DP(X = X, eps = Eps, N = N, K = 20, e_val = e_val, e_vec = e_vec, mu0 = mu0)
      
      error_plug[i,1,] = iid_lap_plug_res$MSE[1:3]
      error_plug[i,2,] = iclp_l1_plug_res$MSE[1:3]
      error_plug[i,3,] = iclp_rkhs_plug_res$MSE[1:3]
      error_plug[i,4,] = bern_res$MSE[1:3]
    }
    MU_tilde[[k]] = cbind(iid_lap_plug_res$mu_tilde, bern_res$mu_tilde,
                          iclp_l1_plug_res$mu_tilde, iclp_rkhs_plug_res$mu_tilde)
    MSE.total[[k]] = error_plug
    
  }
  names(MSE.total) <- paste("mu0",1:4,sep = "")
  MSE_all[[kernel_type]] = MSE.total
}


names(MSE_all) = c("M3/2","M5/2")


save(MSE_all, file = paste("MSE.res",args,".RData",sep=""))
