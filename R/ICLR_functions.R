##################### pick ICLP kernel and calculate eigen_pairs #####################
ICLP_kernel = function(grid, kernel_type = "Exp", rho = 1/4){
  if(kernel_type=="Exp"){
    C=function(t,s,rho){      # kernel of covariance operator of Exponential Process
      return(exp(-abs(t-s)/rho))
    }
  }
  if(kernel_type=="M3/2"){           # kernel of covariance operator of Matern Process nu=3/2
    C=function(t,s,rho){     
      return((1+sqrt(3)*abs(t-s)/rho)*exp(-sqrt(3)*abs(t-s)/rho))
    }
  }
  if(kernel_type=="M5/2"){           # kernel of covariance operator of Matern Process nu=5/2
    C=function(t,s,rho){    
      return((1+(sqrt(5)*abs(t-s)/rho)+(5*(abs(t-s)^2)/(3*rho^2)))*exp(-sqrt(5)*abs(t-s)/rho))
    }
  }
  if(kernel_type=="Gau"){          # kernel of covariance operator of Gaussian Process
    C=function(t,s,rho){    
      return(exp(-(abs(t-s)^2)/rho))
    }
  }
  
  Sig_RKHS = outer(grid$pt,grid$pt, C, rho = rho)
  
  n = length(grid$pt)
  gamma_RKHS=eigen(Sig_RKHS)$values
  e_val_RKHS=gamma_RKHS[1:n]/n
  # e_vec_RKHS=eigen(Sig_RKHS)$vectors[,1:n]
  e_vec_RKHS=eigen(Sig_RKHS)$vectors[,1:n]*sqrt(n)
  m = length(which(e_val_RKHS>0))
  
  return( list( kernel_type = kernel_type,
                e_val = e_val_RKHS,
                e_vec = e_vec_RKHS) )
}

##################### self-defined eigen-pairs #####################
self_e_vec <- function(x,k){
  if(k == 1){
    return(1)
  }else{
    return( sqrt(2)*cos((k-1)*pi*x) )
  }
}
self_e_vec <- Vectorize(self_e_vec)

self_e_val <- function(k,h = 1){
  if(k==1){
    return(1/(pi^{2*h}))
  }else{
    return( 1/( (k*pi)^{2*h} ) )
  }
}
self_e_val <- Vectorize(self_e_val)



##################### generate X ###########################
generate_X = function(N, grid, e_val, e_vec, tau = 1, mu0){
  n = length(grid$pt)
  m = length(e_val)
  Xm = matrix(NA,n,N)
  for(s in 1:N){
    X = matrix(0,n,1)
    U = matrix(NA,m,1)
    for(i in 1:m){
      U[i] = runif(1,min = -tau,max = tau)
      X = X+sqrt(e_val[i])*U[i]*e_vec[,i]
    }
    Xm[,s] = X + mu0
  } 
  return(Xm)
}

################ calculate X-l1 norm #################
X_norm = function(X, grid, e_val, e_vec){
  N = dim(X)[2]
  m = length(e_val)
  n = length(grid$pt)
  
  l1 = l2 = matrix(0,N,1)
  for(i in 1:N){
    for(j in 1:m){
      l1[i,1] = l1[i,1] + abs( sum(X[,i]*e_vec[,j]*grid$wt) )
      l2[i,1] = l2[i,1] + (sum(X[,i]*e_vec[,j]*grid$wt))^2
    }
    l2[i,1] = sqrt(l2[i,1])
  }
  l1_tau = max(l1)
  l2_tau = max(l2)
  
  out = list(l1_tau = l1_tau,
             l2_tau = l2_tau)
  return(out)
  
}


#################### IID + Laplace mechanism #################
IID_Lap = function(X, eps, grid, M, e_val, e_vec, l1.bound = F, tau = NA,
                   Rep = 100, ifMSE = TRUE, mu0 = NA){
  x = X
  n = dim(X)[1]
  N = dim(X)[2]
  
  X_bar = rowMeans(x)
  mu_hat = matrix(0,n,1)
  mu_hat_coef = rep(0,M)
  for(j in 1:M){
    # coef
    mu_hat_coef[j] = sum(X_bar*e_vec[,j]*grid$wt)
    # mu.hat
    mu_hat = mu_hat + mu_hat_coef[j]*e_vec[,j]
  }
  
  ## DP part
  # Calculate GS
  if(is.na(tau) == T){
    norms = X_norm(X = x, grid = grid, e_val = e_val, e_vec = e_vec)
    if(l1.bound == T){
      tau = norms$l1_tau
    }else{
      tau = norms$l2_tau
    }
  }
  
  Delta = abs( (2*tau*M/(N*eps) ) )
  # Adding Laplace noise to each coordinate
  mu_tilde_coef = rep(0,M)
  for(j in 1:M){
    mu_tilde_coef[j] = mu_hat_coef[j] + rlaplace(n = 1, mu = 0, sigma = Delta)
  }
  # Calculate mu_tilde
  mu_tilde = matrix(0,n,1)
  for(j in 1:M){
    mu_tilde = mu_tilde + mu_tilde_coef[j]*e_vec[,j]
  }
  
  
  # Calculate MSE
  MSE = NA
  if(ifMSE == T){
    MSE = rep(0,3)
    for(I in 1:Rep){
      # Adding Laplace noise to each coordinate
      mu_tilde_MSE = matrix(0,n,1)
      for(j in 1:M){
        mu_tilde_MSE = mu_tilde_MSE + rlaplace(n = 1, mu = 0, sigma = Delta)*e_vec[,j]
      }
      mu_tilde_MSE = mu_tilde_MSE + mu_hat
      
      MSE[1] = MSE[1] +  sum( (mu_hat - mu0)^2 *grid$wt )
      MSE[2] = MSE[2] +  sum( (mu_tilde_MSE - mu_hat)^2 *grid$wt )
      MSE[3] = MSE[3] +  sum( (mu_tilde_MSE - mu0)^2 *grid$wt )
    }
    MSE = MSE/Rep
    names(MSE) = c("statistical error","privacy error","MSE")
  }
  out = list(mu_hat = mu_hat,
             mu_tilde = mu_tilde,
             sigma = Delta,
             M = M,
             MSE = MSE,
             tau = tau,
             Delta = Delta)
  return(out)
}


pcv_IID_Lap = function(X, eps, grid, M, fold, e_val, e_vec,
                       l1.bound = F, tau = NA,
                       Rep = 100, ifMSE = TRUE){
  CV = matrix(NA, nrow = length(M),ncol = 1)
  count = 0
  Data = X
  n = length(grid$pt)
  flds = createFolds(c(1:(dim(Data)[2])), k = fold, list = TRUE, returnTrain = FALSE)
  
  
  ################ CV for hetero ########################
  for(I in 1:length(M)){
    MF = matrix(0,nrow = fold,1)
    
    for(i in 1:fold){
      X_train = Data[,-flds[[i]]]
      X_test = Data[,flds[[i]]]
      
      N_train = dim(X_train)[2]
      
      # Calculate mu.hat, mu.hat.coef
      mu_hat_train = matrix(0,n,1)
      mu_hat_coef_train = rep(0,M[I])
      
      X_bar_train = rowMeans(X_train)
      
      for(j in 1:M[I]){
        # coef
        mu_hat_coef_train[j] = sum(X_bar_train*e_vec[,j]*grid$wt)
        # mu.hat
        mu_hat_train = mu_hat_train + mu_hat_coef_train[j]*e_vec[,j]
      }
      
      ## DP part
      # Calculate GS
      if(is.na(tau) == T){
        norms = X_norm(X = X_train, grid = grid, e_val = e_val, e_vec = e_vec)
        if(l1.bound == T){
          tau = norms$l1_tau
        }else{
          tau = norms$l2_tau
        }
      }
      
      Delta = abs( (2*tau*M[I]/(N_train*eps) ) )
      
      mu_tilde = matrix(0,n,Rep)
      Q = 0
      for(l in 1:Rep){
        for(j in 1:M[I]){
          mu_tilde[,l] = mu_tilde[,l] + rlaplace(n = 1, mu = 0, sigma = Delta)*e_vec[,j]
        }
        mu_tilde[,l] = mu_tilde[,l] + mu_hat_train
        
        # test set
        for(k in 1:length(flds[[i]])){
          Q = Q + sum( grid$wt * (mu_tilde[,l]-X_train[,k])^2 )
        }
      }
      MF[i] = Q/Rep
    } # fold loop
    
    CV[I,1] = mean(MF)
    
  } # M loop
  arg.min.cv = which(CV == min(CV), arr.ind=TRUE)
  rownames(CV) = c(paste( "M=",M,sep=""))
  
  out = list(M = M[arg.min.cv[1]],
             CV = CV)
  return(out)
  
}




#################### l1 regularization ######################
ICLP_l1 = function(X, eps, grid, e_val, e_vec, 
                   psi, J_tau = NA, rho = NA, tau = NA, eta = NA, l1.bound = F,
                   ifMSE = T, Rep = 500, mu0 = NA){
  x = X
  n = dim(X)[1]
  N = dim(X)[2]
  m = length(e_val)
  
  # Calculate mu.hat, mu.hat.coef
  mu_hat = matrix(0,n,1)
  mu_hat_coef = rep(0,m)
  
  X_bar = rowMeans(x)
  b = sqrt(abs(e_val))/sqrt(2)
  b_eta = sqrt(abs(e_val)^{eta})/sqrt(2)
  
  for(j in 1:m){
    # coef
    coefs = sum(X_bar*e_vec[,j]*grid$wt)
    mu_hat_coef[j] = sign(coefs)*(abs(coefs)-(psi/b_eta[j]))*((abs(coefs)-(psi/b_eta[j]))>0)
    # mu.hat
    mu_hat = mu_hat + mu_hat_coef[j]*e_vec[,j]
  }
  
  A = c()
  for(j in 1:m){
    if(mu_hat_coef[j] !=0){
      A = c(A,j)
    }
  }
  
  ## DP part
  # Calculate GS
  if(is.na(tau) == T){
    norms = X_norm(X = x, grid = grid, e_val = e_val, e_vec = e_vec)
    if(l1.bound == T){
      tau = norms$l1_tau
    }else{
      tau = norms$l2_tau
    }
  }
  if(is.na(J_tau) == T){
    J_tau = min(which(tau <= psi/b_eta))
  }
  
  
  Delta =  rep(0,J_tau)
  for(j in 1:J_tau){
    Delta[j] = (2*tau*sum(1/b[1:J_tau])*b[j])/(N*eps)
  }
  
  
  # Adding Laplace noise to each coordinate
  mu_tilde_coef = rep(0,J_tau)
  for(j in 1:J_tau){
    mu_tilde_coef[j] = mu_hat_coef[j] + rlaplace(n = 1, mu = 0, sigma = Delta[j])
  }
  # Calculate mu.tilde
  mu_tilde = matrix(0,n,1)
  for(j in 1:J_tau){
    mu_tilde = mu_tilde + mu_tilde_coef[j]*e_vec[,j]
  }
  
  # Calculate MSE
  MSE = NA
  if(ifMSE == T){
    MSE = rep(0,3)
    for(I in 1:Rep){
      # Adding Laplace noise to each coordinate
      mu_tilde_MSE = matrix(0,n,1)
      for(j in 1:J_tau){
        mu_tilde_MSE = mu_tilde_MSE + rlaplace(n = 1, mu = 0, sigma = Delta[j])*e_vec[,j]
      }
      mu_tilde_MSE = mu_tilde_MSE + mu_hat
      
      MSE[1] = MSE[1] +  sum( (mu_hat - mu0)^2 *grid$wt )
      MSE[2] = MSE[2] +  sum( (mu_tilde_MSE - mu_hat)^2 *grid$wt )
      MSE[3] = MSE[3] +  sum( (mu_tilde_MSE - mu0)^2 *grid$wt )
    }
    MSE = MSE/Rep
    names(MSE) = c("statistical error","privacy error","MSE")
  }
  out = list(mu_hat = mu_hat,
             mu_tilde = mu_tilde,
             MSE = MSE,
             psi = psi,
             Delta = Delta,
             J_tau = J_tau)
  return(out)
}


pcv_ICLP_l1 = function(X, eps, grid, e_val, e_vec, 
                       psi, fold = 10, J_tau = NA, rho = NA, tau = NA, eta = NA, l1.bound = F,
                       ifMSE = T, Rep = 500, mu0 = NA){
  CV = matrix(NA, nrow = length(psi),ncol = 1)
  count = 0
  Data = X
  n = length(grid$pt)
  flds = createFolds(c(1:(dim(Data)[2])), k = fold, list = TRUE, returnTrain = FALSE)
  b = sqrt(e_val/2)
  b_eta = sqrt(e_val^{eta}/2)
  m = length(e_val)
  
  
  for(I in 1:length(psi)){
    ## Trainning set
    MF = matrix(0,nrow = fold,1)
    
    for(i in 1:fold){
      X_train = Data[,-flds[[i]]]
      X_test = Data[,flds[[i]]]
      N_train = dim(X_train)[2]
      
      # Calculate mu.hat, mu.hat.coef
      mu_hat_train = matrix(0,n,1)
      mu_hat_coef_train = rep(0,m)
      X_bar_train = rowMeans(X_train)
      
      for(j in 1:m){
        # coef
        coefs = sum(X_bar_train*e_vec[,j]*grid$wt)
        mu_hat_coef_train[j] = sign(coefs)*(abs(coefs)-(psi[I]/b_eta[j]))*((abs(coefs)-(psi[I]/b_eta[j]))>0)
        # mu.hat
        mu_hat_train = mu_hat_train + mu_hat_coef_train[j]*e_vec[,j]
      }
      
      A_train = c()
      for(j in 1:m){
        if(mu_hat_coef_train[j] !=0){
          A_train = c(A_train,j)
        }
      }
      
      ## DP part
      # Calculate GS
      if(is.na(tau) == T){
        norms = X_norm(X = X_train, grid = grid, e_val = e_val, e_vec = e_vec)
        if(l1.bound == T){
          tau = norms$l1_tau
        }else{
          tau = norms$l2_tau
        }
      }
      if(is.na(J_tau) == T){
        J_tau = min(which(tau <= psi[I]/b_eta))
      }
      
      
      Delta =  rep(0,J_tau)
      for(j in 1:J_tau){
        Delta[j] = (2*tau*sum(1/b[1:J_tau])*b[j])/(N_train*eps)
      }
      
      
      f_hat_train = rowMeans(Data[,flds[[i]]])
      mu_tilde = matrix(0,n,Rep)
      Q = 0
      for(l in 1:Rep){
        mu_tilde_coef = rep(0,J_tau)
        for(j in 1:J_tau){
          mu_tilde_coef[j] = mu_hat_coef_train[j] + rlaplace(n = 1, mu = 0, sigma = Delta[j])
        }
        # Calculate mu.tilde for l th Rep
        for(j in 1:J_tau){
          mu_tilde[,l] = mu_tilde[,l] + mu_tilde_coef[j]*e_vec[,j]
        }
        for(k in 1:length(flds[[i]])){
          Q = Q + sum( grid$wt * (mu_tilde[,l]-X_test[,k])^2 )
        }
      }
      MF[i] = Q/Rep
      CV[I,1] = mean(MF)
    } # fold loop
  } # psi loop
  
  
  
  arg.min.cv = which(CV == min(CV), arr.ind=TRUE)
  rownames(CV) = c(paste( "psi=",psi,sep=""))
  
  out = list(psi = psi[arg.min.cv[1]],
             CV = CV)
  return(out)
}






##################### RKHS regularization #####################
ICLP_RKHS = function(X, eps = 1, grid = NULL, e_val, e_vec,
                     psi = 0 , rho = NA, tau = NA, eta = NA, l1.bound = T, 
                     Rep = 100,ifMSE = TRUE, mu0 = NA){
  x = X
  n = dim(X)[1]
  N = dim(X)[2]
  m = length(e_val)
  
  if(is.na(eta)==T) stop("eta is required.")
  # Calculate mu.hat
  mu_hat = matrix(0,n,1)
  mu_hat_coef = rep(0,m)
  
  X_bar = rowMeans(x)
  for(j in 1:m){
    # coef
    mu_hat_coef[j] = (e_val[j]^{eta}/(e_val[j]^{eta} + psi))*sum(X_bar*e_vec[,j]*grid$wt)
    # mu_hat
    mu_hat = mu_hat + mu_hat_coef[j]*e_vec[,j]
  }
  
  ## generate laplace type process
  Z = matrix(0,nrow = n,ncol = 1)
  for(j in 1:m){
    Z = Z+sqrt(e_val[j]/2)*rlaplace(n = 1, mu = 0, sigma = 1)*e_vec[,j]
  }
  
  # Calculate GS
  if(is.na(tau) == T){
    norms = X_norm(X = x, grid = grid, e_val = e_val, e_vec = e_vec)
    if(l1.bound == T){
      tau = norms$l1_tau
    }else{
      tau = norms$l2_tau
    }
  }
  
  
  if(l1.bound==T){
    factor = max( e_val^{(eta-0.5)}/(e_val^{eta} + psi) )
    delta = factor*(2*tau)/(N*eps)
  }else{
    factor = sum( (e_val^{(eta - 0.5)})/( e_val^{eta} + psi )  )
    delta = factor*(2*tau)/(N*eps)
  }
  
  mu_tilde = mu_hat + delta*Z
  
  
  
  # Calculate MSE
  MSE = NA
  if(ifMSE == T){
    MSE = rep(0,3)
    for(I in 1:Rep){
      ## generate laplace type process
      Z_mse = matrix(0,nrow = n,ncol = 1)
      for(j in 1:m){
        Z_mse = Z_mse + sqrt(e_val[j]/2)*rlaplace(n = 1, mu = 0, sigma = 1)*e_vec[,j]
      }
      
      # Calculate mu.tilde
      mu_tilde_MSE = mu_hat + delta*Z_mse
      
      MSE[1] = MSE[1] +  sum( (mu_hat - mu0)^2 *grid$wt )
      MSE[2] = MSE[2] +  sum( (mu_tilde_MSE - mu_hat)^2 *grid$wt )
      MSE[3] = MSE[3] +  sum( (mu_tilde_MSE - mu0)^2 *grid$wt )
    }
    MSE = MSE/Rep
    names(MSE) = c("statistical error", "privacy error", "MSE")
  }
  
  out = list(mu_hat = mu_hat,
             mu_tilde = mu_tilde,
             psi = psi,
             sigma = delta,
             trace = factor,
             MSE = MSE,
             tau = tau)
  return(out)
}

pcv_ICLP_RKHS = function(X, eps = 1, grid = NULL, e_val, e_vec,
                         psi = 0 , fold = 10, rho = NA, tau = NA, eta = NA, l1.bound = T, 
                         Rep = 100,ifMSE = TRUE, mu0 = NA){
  CV = matrix(NA, nrow = length(psi),ncol = 1)
  count = 0
  Data = X
  n = length(grid$pt)
  flds = createFolds(c(1:(dim(Data)[2])), k = fold, list = TRUE, returnTrain = FALSE)
  m = length(e_val)
  
  for(I in 1:length(psi)){
    ## Trainning set
    MF = matrix(0,nrow = fold,1)
    
    for(i in 1:fold){
      X_train = Data[,-flds[[i]]]
      X_test = Data[,flds[[i]]]
      N_train = dim(X_train)[2]
      
      # Calculate mu.hat, mu.hat.coef
      mu_hat_train = matrix(0,n,1)
      mu_hat_coef_train = rep(0,m)
      X_bar_train = rowMeans(X_train)
      
      for(j in 1:m){
        # coef
        mu_hat_coef_train[j] = (e_val[j]^{eta}/(e_val[j]^{eta} + psi[I]))*sum(X_bar_train*e_vec[,j]*grid$wt)
        # mu_hat
        mu_hat_train = mu_hat_train + mu_hat_coef_train[j]*e_vec[,j]
      }
      # Calculate GS
      if(is.na(tau) == T){
        norms = X_norm(X = X_train, grid = grid, e_val = e_val, e_vec = e_vec)
        if(l1.bound == T){
          tau = norms$l1_tau
        }else{
          tau = norms$l2_tau
        }
      }
      if(l1.bound==T){
        factor = max( e_val^{(eta-0.5)}/(e_val^{eta} + psi[I]) )
        delta = factor*(2*tau)/(N*eps)
      }else{
        factor = sum( (e_val^{(eta - 0.5)})/( e_val^{eta} + psi[I] )  )
        delta = factor*(2*tau)/(N*eps)
      }
      
      
      
      f_hat_train = rowMeans(Data[,flds[[i]]])
      mu_tilde = matrix(0,n,Rep)
      Q = 0
      for(l in 1:Rep){
        ## generate laplace type process
        Z = matrix(0,nrow = n,ncol = 1)
        for(j in 1:m){
          Z = Z+sqrt(e_val[j]/2)*rlaplace(n = 1, mu = 0, sigma = 1)*e_vec[,j]
        }
        mu_tilde[,l] = mu_hat_train + delta*Z
        # test set
        # Q = Q + sum( grid$wt * (mu_tilde[,l]-f_hat_train)^2 )
        for(k in 1:length(flds[[i]])){
          Q = Q + sum( grid$wt * (mu_tilde[,l]-X_test[,k])^2 )
        }
      }
      MF[i] = Q/Rep
      CV[I,1] = mean(MF)
    } # fold loop
  } # psi loop
  
  
  
  arg.min.cv = which(CV == min(CV), arr.ind=TRUE)
  rownames(CV) = c(paste( "psi=",psi,sep=""))
  
  out = list(psi = psi[arg.min.cv[1]],
             CV = CV)
  return(out)
}





############ bernstein ####################
mean_est = function(X){
  n = length(grid$pt)
  X_bar = colMeans(X)
  m = length(e_val)
  mu_hat = matrix(0,n,1)
  mu_hat_coef = rep(0,m)
  
  for(j in 1:m){
    # coef
    mu_hat_coef[j] = sum(X_bar*e_vec[,j]*grid$wt)
    # mu.hat
    mu_hat = mu_hat + mu_hat_coef[j]*e_vec[,j]
  }
  
  predictor = function(t) {
    mu_hat[which( abs(t - grid$pt) == min( abs(t - grid$pt) ) ),]
  }
  return(predictor) 
}

bern_DP = function(X, eps, N, K = 20, e_val, e_vec,
                   Rep = 100, ifMSE = TRUE, mu0 = NA){
  MM1 = matrix(NA,N,1)
  for(ii in 1:dim(X)[2]){
    MM1[ii,1] = max(abs(X[,ii]))
  }
  tau1 = max(MM1)
  Data = t(X)
  model = mean_est(Data)
  
  Bern_M = DPMechBernstein(target=mean_est, latticeK=K, dims=1, sensitivity = 2*tau1/N)
  MSE.bernsterin = rep(0,100)
  priv.bernstein = rep(0,100)
  for(ii in 1:100){
    R = releaseResponse(Bern_M, privacyParams=DPParamsEps(epsilon=eps), X=Data)
    pmodel = R$response
    MSE.bernsterin[ii] = sum( (pmodel(grid$pt) - mu0)^2 * grid$wt  )
    priv.bernstein[ii] = sum( (pmodel(grid$pt) - model(grid$pt))^2 * grid$wt  )
  }
  MSE = c( sum( (model(grid$pt) - mu0)^2*grid$wt ), mean(priv.bernstein), mean(MSE.bernsterin) )
  names(MSE) = c("statistical error", "privacy error", "MSE")
  
  mu_hat = sapply(grid$pt, model)
  mu_tilde_bern = pmodel(grid$pt)
  
  out = list(mu_hat = mu_hat,
             mu_tilde = mu_tilde_bern,
             K = K,
             MSE = MSE)
  return(out)
  
}







comp_plot_func <- function(mu.re,X, mu0 = NA, grid = NULL,
                           text.main="",legend.cex=1,legend.loc="topright",
                           seg.lin=0.75,text.font=1,text.width=0.1,
                           xlab="",ylab="",extra_range=c(0,0),iflegend = T){
  x <- X
  cols = c('green', 'yellow', 'red', 'blue')
  
  plot(y=mu.re[,1],x=grid,type="l",col="green",ylim=range(x,mu.re)+c(0,0),lwd=3,
       lty=3,ylab=ylab,xlab=xlab,main=text.main,cex.lab=2,cex.main=2.4,cex.axis=1.5)
  
  for(i in 1:ceiling(1*dim(x)[2])){
    points(y=x[,i],x=grid,type="l",col="grey",lwd=1,lty=2)
  }
  for(i in 1:ncol(mu.re)){
    points(y=mu.re[,i],x=grid,type="l",col=cols[i],ylim=range(mu.re),lty=1,lwd=3)
  }
  if(unique(is.na(mu0)) == F){
    points(y=mu0,x=grid,type="l",col="black",bty="n",lty=1,lwd=3)
  }
  if(iflegend == T){
    legend(x=legend.loc,legend=c('IID Laplace', 'Bernstein', 'l1', 'RKHS'),
           col=c('green', 'yellow', 'red', 'blue', 'black'),
           lty=c(1,1,1,1,1),lwd=c(3,3,3),cex=legend.cex,bty="n",
           seg.len=seg.lin,text.font = text.font,text.width = text.width)
  }
}



mean.est <- function(X){
  m = length(e_val)
  n <- length(grid$pt)
  X.bar <- colMeans(X)
  mu.hat <- matrix(0,n,1)
  mu.hat.coef <- rep(0,m)
  
  phi <- matrix(0,nrow = n, ncol = m)
  for(i in 1:m){
    phi[,i] <- e_vec[,i]
  }
  
  for(j in 1:m){
    # coef
    mu.hat.coef[j] <- sum(grid$wt*X.bar*phi[,j])
    # mu.hat
    mu.hat <- mu.hat + mu.hat.coef[j]*phi[,j]
  }
  
  predictor <- function(t) {
    mu.hat[which( abs(t - grid$pt) == min( abs(t - grid$pt) ) ),]
  }
  return(predictor) 
}


stand_curve <- function(Data, grid = NULL){
  N_size <- dim(Data)[2]
  if(is.null(grid)){
    grid <- seq(0,1, length.out = dim(Data)[1])
  }
  M <- SD <- NULL
  for(i in 1:N_size){
    M <- c(M, trapz(grid, Data[,i]))
  }
  for(i in 1:N_size){
    SD <- c(SD, sqrt(trapz( grid, (Data[,i] - M[i])^2  ) ) )
  }
  
  Data_s <- NULL
  for(i in 1:N_size){
    Data_s <- cbind(Data_s, (Data[,i] - M[i]) / SD[i]  )
  }
  
  return(Data_s)
}
