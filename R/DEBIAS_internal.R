
DEBIAS_internal <- function(Y,Tx,tp,Xp,lambda,alphak = NULL, max_iter = 100, tol = 1e-4){
  #
  # solves the optimization problem once for a fixed lambda and for one non-negative outcome weighting
  #
  # Inputs:
  #   Y: outcome items, where column names should contain numeric time points. The number of outcome items per time point must be the same.
  #   Tx: treatment assignments, where column names should again contain numeric time points
  #   tp: time point of interest (variable p in the paper)
  #   Xp: covariates. Must not contain any missing values
  #   lambda: non-negative weight for confounding penalty (b)
  #   alphak: matrix of non-negative weights computed in past iterations for orthogonality penalty
  #   max_iter: the maximum number of gradient steps (default: 100)
  #   tol: convergence tolerance (default: 1e-4)
  #
  #   A row index should correspond to the same patient across all matrices Y, Tx, and Xp
  #
  # Outputs: A list with
  #   alpha: vector of non-negative weights
  #   main_correlations: main correlation term (a) for each time point
  #   confounding_correlations: square root of confounding penalty term (b) (i.e., its absolute value of correlation instead of squared correlation) for each time point
  #   confounding_pvalues: p-value of confounding penalty term (b) for each time point
  #   Mahalanobis_cosine_similarity: Mahalanobis cosine similarity for orthogonality penalty (c)
  #   lambda: lambda value
  #
  # Written by Eric V. Strobl 06/2025
  #
  
  
  tY = as.numeric(colnames(Y))
  tT = as.numeric(colnames(Tx))
  
  colnames(Y) = c(); colnames(Tx) = c()
  
  Tm = Tx; Tn = Tx
  Yn = Y; Ym = Y
  
  u_time = sort(unique(c(tT,tY)),decreasing=F)
  m = length(u_time)
  ip = which(u_time == tp)
  
  # compute residuals
  Cs = list() # correlation matrices for Y
  
  Tm = list(); Tn = list()
  Ym = list(); Yn = list()
  for (i in (ip+1):m){ # outer summation
    
    # residuals for main correlation term (a)
    ic = complete.cases(cbind(Xp,Tx[,tT==u_time[1]],Tx[,tT==tp],Y[,tY==u_time[i]])) # all variables for main correlation
    
    Tm[[i]] = lm.fit(cbind(1,Xp[ic,],Tx[ic,tT==u_time[1]]),Tx[ic,tT==tp])$residuals
    Tm[[i]] = normalize_uncor(Tm[[i]])
    
    Ym[[i]] = lm.fit(cbind(1,Xp[ic,],Tx[ic,tT==u_time[1]]),Y[ic,tY==u_time[i]])$residuals
    
    # correlation matrices of Y in orthogonality penalty (c)
    io = complete.cases(Y[,tY==u_time[i]])
    Cs[[i]] = cor(Y[io,tY==u_time[i]])
    
    # residuals for confounding penalty (b)
    Yn[[i]] = list()
    Tn[[i]] = list()
    for (j in 1:(ip-1)){
      
      ic = complete.cases(cbind(Xp,Tx[,tT==tp],Y[,tY==u_time[i]]),Tx[,tT==u_time[j]]) # all variables for main correlation
      
      Tn[[i]][[j]] = lm.fit(cbind(1,Xp[ic,],Tx[ic,tT==tp]),Tx[ic,tT==u_time[j]])$residuals
      Tn[[i]][[j]] = normalize_uncor(Tn[[i]][[j]])
      
      Yn[[i]][[j]] = lm.fit(cbind(1,Xp[ic,],Tx[ic,tT==tp]),Y[ic,tY==u_time[i]])$residuals
    }
    
  }
  
  # initial guess of uniform weights
  q <- ncol(Ym[[m]])
  alpha <- rep(1/q, q)
  # alpha = rnorm(q)
  
  for (iter in 1:max_iter) {
    alpha_old = alpha

    ### compute gradients
    grad = 0; grad_conf = 0; grad_cosine = 0 
    cosine_mean = 0
    for (i in (ip+1):m){
      
      # gradient of main correlation term (a)
      den = norm(Ym[[i]] %*% alpha, "2") + 1E-10
      grad = grad + t(Ym[[i]]) %*% Tm[[i]] / den -
        c(t(Tm[[i]]) %*% Ym[[i]] %*% alpha / den^3) * t(Ym[[i]]) %*% Ym[[i]] %*% alpha
      
      # gradient of confounding penalty (b)
      for (j in 1:(ip-1)){
        den = norm(Yn[[i]][[j]] %*% alpha, "2") + 1E-10
        U = c(t(Yn[[i]][[j]] %*% alpha) %*% Tn[[i]][[j]])
        grad_conf = grad_conf + (U * t(Yn[[i]][[j]]) %*% Tn[[i]][[j]] / den^2 - 
                                         ((U^2)/ den^4) * t(Yn[[i]][[j]]) %*% Yn[[i]][[j]] %*% alpha)
      
      }
      
      # gradient of orthogonality penalty (c)
      if (!is.null(alphak)){
        
        for (k in 1:ncol(alphak)){
          grad_cosine = grad_cosine + (Cs[[i]] %*% alphak[,k] / c(sqrt(t(alphak[,k]) %*% Cs[[i]] %*% alphak[,k]) * sqrt(t(alpha) %*% Cs[[i]] %*% alpha)+1E-10) - 
                                         c(t(alphak[,k]) %*% Cs[[i]] %*% alpha)/c(sqrt(t(alphak[,k]) %*% Cs[[i]] %*% alphak[,k]) *
                                                                                       (t(alpha) %*% Cs[[i]] %*% alpha)^(3/2)+1E-10) * Cs[[i]] %*% alpha)
        }
      }
      
    }
    
    # combine all gradients
    if (!is.null(alphak)){
      grad = grad - (2*lambda/(ip-1))*grad_conf - grad_cosine / ncol(alphak)
    } else{
      grad = grad - (2*lambda/(ip-1))*grad_conf
    }
    
    # backtracking line search with the Armijo condition
    eta = backtracking(alpha,grad,Ym,Tm,Yn,Tn,lambda,ip,alphak,Cs)
    
    # take gradient step (maximization)
    alpha = alpha + eta*grad
    
    # projection onto feasible set
    alpha = pmax(alpha,0)
    alpha = alpha / (sum(alpha)+1e-10)
    
    # check convergence
    if ( sum(abs(alpha - alpha_old)) < tol) {
      break
    }

  }
  
  # project back onto feasible set, again
  if (sum(alpha) == 0) alpha = rep(1/length(alpha),length(alpha))
  alpha = alpha / (sum(alpha)+1e-10)
  
  # compute confounding p-value
  penalty_pval = c()
  cors = c()
  conf_cors = c()
  for (i in (ip+1):m){
    cors = c(cors, cor(Ym[[i]] %*% alpha, Tm[[i]]))
    for (j in 1:(ip-1)){
      cor_abs = abs(cor(Yn[[i]][[j]] %*% alpha, Tn[[i]][[j]]))
      conf_cors = c(conf_cors,cor_abs)
      nsamps = nrow(Yn[[i]][[j]])
      stat = cor_abs*sqrt( (nsamps-2) / (1-cor_abs^2))
      penalty_pval = c(penalty_pval, 2 * (1 - pt(stat, df = nsamps - 2)))
    }
  }
  
  return(list(alpha = alpha, main_correlations = cors, confounding_correlations = conf_cors, 
              confounding_pvalues = penalty_pval, Mahalanobis_cosine_similarity = cosine_mean, lambda = lambda))
}

backtracking <- function(alpha,grad,Ym,Tm,Yn,Tn,lambda,ip,alphak=NULL,Cs=NULL,eta=1){
  ## backtracking line search with the Armijo condition
  
  for (i in 1:100){ # Armijo condition
    if (my_obj(alpha + eta*grad,Ym,Tm,Yn,Tn,lambda,ip,alphak,Cs) < (my_obj(alpha,Ym,Tm,Yn,Tn,lambda,ip,alphak,Cs) + 1E-10*eta*norm(grad,"2")^2)){
      eta = 0.5*eta
    } else{
      break
    }
    
  }

  return(eta)
  
}

my_obj <- function(alpha,Ym,Tm,Yn,Tn,lambda,ip,alphak = NULL,Cs=NULL){
  ## compute objective function for Armijo condition
  
  # project onto feasible set
  alpha = pmax(alpha,0)
  alpha = alpha / (sum(alpha)+1e-10)
  m = length(u_time)
  cors = 0; cors_conf = 0; cors_cosine = 0
  for (i in (ip+1):m){
    # main correlation term (a)
    cors = cors+ t(Ym[[i]] %*% alpha) %*% Tm[[i]] / (norm(Ym[[i]] %*% alpha,"2")+1E-10)
    
    # confounding penalty (b)
    for (j in 1:(ip-1)){
      cors_conf = cors_conf + (t(Yn[[i]][[j]] %*% alpha) %*% Tn[[i]][[j]] / (norm(Yn[[i]][[j]] %*% alpha,"2")+1E-10))^2
    }
    
    # orthogonality penalty (c)
    if (!is.null(alphak)){
      for (k in 1:ncol(alphak)){
        cors_cosine = cors_cosine + c(t(alphak[,k]) %*% Cs[[i]] %*% alpha)/c(sqrt(t(alphak[,k]) %*% Cs[[i]] %*% alphak[,k]) *
                                                                               sqrt(t(alpha) %*% Cs[[i]] %*% alpha)+1E-10)   
      }
    }
  }
  
  # put (a), (b), and (c) together
  if (!is.null(alphak)){
    cors = cors - cors_conf*lambda/(ip-1) - cors_cosine/ncol(alphak)
  } else{
    cors = cors - cors_conf*lambda/(ip-1)
  }
  
  return(cors)
}


normalize_uncor <- function(x){
  # normalize so that x is centered and sqrt(centered dot product) = 1
  
  x = x-mean(x)
  sd = sqrt(sum(x^2))
  x = x/sd
  # print(sqrt(sum((x - mean(x))^2)))
  x
}

