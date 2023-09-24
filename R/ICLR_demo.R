require(extraDistr)
require(pracma)
require(caret)  
require(RandomFieldsUtils)
require(gss)
require(diffpriv)

source("ICLP_functions.R")

#' @param N Sample size of the dataset
#' @param Eps Privacy Budget
#' @param h Decay rate of the eigenvalues, which is tied to the Matern Kernel 

N = 1000  
Eps = 1


###### general setting for functional data #######
grid.size = 101
grid = list(pt = seq(0,1,length = grid.size),
            wt = rep(1/grid.size,grid.size))

############## eigen-pairs from existing kernels ##############
e_pairs = ICLP_kernel(grid = grid, kernel_type = "M3/2", rho = 1/10)
e_val = e_pairs$e_val
e_vec = e_pairs$e_vec
h = 2


############# eigen-pairs self-defined ###############
# grid.size <- 100; mm <- 100
# grid <- gauss.quad(size = grid.size, interval = c(0,1))
# h = 2
# e_vec = matrix(0,nrow = grid.size, ncol = mm)
# for(KK in 1:mm){
#   e_vec[,KK] = self_e_vec(x = grid$pt, k = KK)
# }
# e_val = self_e_val(k = 1:mm, h = h)



####### choose the groud truth of the mean function #####
# mu0 = (grid$pt)^{-1} *sin(2*pi*grid$pt) + 3
mu0 = 5/(1 + 5*grid$pt)
# mu0 = exp(-grid$pt)*10*grid$pt
# mu0 = rowSums(   e_vec[,1:25]%*%diag( runif(n = 25, min = -.5,max = .5) ) )



X = generate_X(N = N, grid = grid, e_val = e_val, e_vec = e_vec, tau = 1, mu0 = mu0)

##### calculate the L^2 and L^1 norm for dataset #####
norms = X_norm(X = X, grid = grid, e_val = e_val, e_vec = e_vec)
l1_tau = norms$l1_tau
l2_tau = norms$l2_tau


### IID ##########
nocv_M_candi = ceiling( N^(1/(2*h)) - 1 ):floor( (N)^{1/(3)} )
MSE_IID = rep(0,length(nocv_M_candi))
for(nocv.M1 in nocv_M_candi){
  iid_lap_res = IID_Lap(X = X, eps = Eps, M = nocv.M1, grid = grid, e_val = e_val, e_vec = e_vec, l1.bound = F,
                        Rep = 100, ifMSE = T, mu0 = mu0)
  MSE_IID[which(nocv.M1 == nocv_M_candi)] = iid_lap_res$MSE[3]
}
nocv.M1 = nocv_M_candi[which.min(MSE_IID)]
iid_lap_res = IID_Lap(X = X, eps = Eps, M = nocv.M1, grid = grid, e_val = e_val, e_vec = e_vec, l1.bound = F,
                      Rep = 100, ifMSE = T, mu0 = mu0)

### ICLP l1 ######
eta_ll = 2*(1+1/h)
nocv.psi = l2_tau*(1/N)
# nocv.Jtau = ceiling((nocv.psi/l2_tau)^{-1/(eta_ll*h)})
nocv.Jtau = NA
iclp_l1_res = ICLP_l1(X = X, eps = Eps, grid = grid, e_val = e_val, e_vec = e_vec, psi = nocv.psi, J_tau = nocv.Jtau,
                      eta = eta_ll, l1.bound = F, ifMSE = T, Rep = 100, mu0 = mu0)


### ICLP_RKHS ########
eta_rl = 1 + 1/(2*h)
nocv.psi = (N*(Eps)/l2_tau)^{-eta_rl}
iclp_rkhs_res = ICLP_RKHS(X = X, eps = Eps, grid = grid, e_val = e_val, e_vec = e_vec,
                          psi = nocv.psi, eta = eta_rl, l1.bound = T, ifMSE = T, mu0 = mu0)



########## bernstein ################
bern_res = bern_DP(X = X, eps = Eps, N = N, K = 10, e_val = e_val, e_vec = e_vec, mu0 = mu0)



######### draw one santized mean function from each  mechanism ###############
plot(x=grid$pt, y=iid_lap_res$mu_tilde,type="l",col="blue",
     ylim= range(iid_lap_res$mu_tilde, iclp_l1_res$mu_tilde,
                 iclp_rkhs_res$mu_tilde, bern_res$mu_tilde, mu0),
     lwd=2, lty=3, ylab="", xlab="", main = "Sanitized Curve")
for(i in 1:ceiling(1*dim(X)[2])){
  points(y=X[,i],x=grid$pt,type="l",col="grey",lwd=1,lty=1)
}
points(y=iid_lap_res$mu_tilde, x=grid$pt,type="l",col="blue",lty=1,lwd=2)
points(y=iclp_l1_res$mu_tilde, x=grid$pt,type="l",col="red",bty="n",lty=2,lwd=2)
points(y=iclp_rkhs_res$mu_tilde, x=grid$pt,type="l",col="green",bty="n",lty=2,lwd=2)
points(y=bern_res$mu_tilde,x=grid$pt,type="l",col="purple",bty="n",lty=3,lwd=2)
points(y=mu0,x=grid$pt,type="l",col="black",bty="n",lty=3,lwd=2)

legend(x="topright",legend=c("true Mean","IID_Laplace", "ICLP_l1", "ICLP_RKHS", "Bernstein"),
       col=c("black","blue","red", "green", "purple"), lty=c(1,2,3,4,5), lwd=c(rep(2,5)), cex=1, bty="n",
       seg.len=1, text.font = 1, text.width = .2)



######## report total MSE, statistical error, privacy cost for each 
method.list = c("IIDLaplace","l1","RKHS","bernstein")
MSE = array(data = 0, dim = c( length(method.list), 3),
            dimnames = list( c(method.list), c("Statistical Error","Privacy Error","MSE"))  )
MSE[1,] = iid_lap_res$MSE[1:3]
MSE[2,] = iclp_l1_res$MSE[1:3]
MSE[3,] = iclp_rkhs_res$MSE[1:3]
MSE[4,] = bern_res$MSE[1:3]
MSE
MSE[,1]/MSE[,3]




mu_hat = matrix(0,grid.size,1)
mu_hat_coef = rep(0,50)
for(j in 1:25){
  mu_hat_coef[j] = sum(rowMeans(X)*e_vec[,j]*grid$wt)
  mu_hat = mu_hat + mu_hat_coef[j]*e_vec[,j]
}


plot(grid$pt, bern_res$mu_hat, type = "l")
points(grid$pt, iid_lap_res$mu_hat, type = "l", col = "blue")
points(grid$pt, iclp_l1_res$mu_hat, type = "l", col = "red")
points(grid$pt, iclp_rkhs_res$mu_hat, type = "l", col = "green")

