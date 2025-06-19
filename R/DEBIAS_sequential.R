DEBIAS_sequential <- function(Y, Tx, tp, Xp, Yte=NULL, Txte=NULL, Xpte=NULL, lambda=1, max_iterk=3) {
  #
  # solves the optimization problem max_iterk times for a fixed lambda to obtain max_iterk non-negative outcome weighting
  #
  # Inputs:
  #   Y: outcome items, where column names should contain numeric time points. The number of outcome items per time point must be the same.
  #   Tx: treatment assignments, where column names should again contain numeric time points
  #   tp: time point of interest (varible p in the paper)
  #   Xp: covariates. Must not contain any missing values
  #   For outer cross-validation:
  #     Yte: outcome items for test/validation set evaluation (default: NULL)
  #     Txte: treatment assignment for test/validation set evaluation (default: NULL)
  #     Xpte: outcome items for test/validation set evaluation (default: NULL)
  #   lambda: non-negative weight for confounding penalty (b)
  #   max_iterk: maximum number of non-negative weights to extract (variable s in the paper)
  #
  # Outputs:
  #   mods: a list of models with output of DEBIAS_internal.R and:
  #     main_correlations_testset: main correlation terms (a) on test/validation set
  #     confounding_pvalues_testset: confounding penalty (b) p-values on test/validation set
  #
  # Written by Eric V. Strobl 06/2025
  #
  
  
  mods = list()
  alphak = NULL
  for (k in 1:max_iterk){ # run DEBIAS sequentially max_iterk times
    mods[[k]] = DEBIAS_internal(Y,Tx,tp,Xp,lambda=lambda,alphak)
    alphak = cbind(alphak, mods[[k]]$alpha)
  }
  
  tY = as.numeric(colnames(Y))
  tT = as.numeric(colnames(Tx))
  
  colnames(Y) = c(); colnames(Tx) = c()
  
  u_time = sort(unique(c(tT,tY)),decreasing=F)
  m = length(u_time)
  ip = which(u_time == tp)
  
  if (!is.null(Xpte)){
    for (k in seq_len(length(mods))){
      
      for (i in (ip+1):m){
        
        ic = complete.cases(cbind(Xp,Tx[,tT==u_time[1]],Tx[,tT==tp],Y[,tY==u_time[i]])) # all variables for main correlation
        icte = complete.cases(cbind(Xpte,Txte[,tT==u_time[1]],Txte[,tT==tp],Yte[,tY==u_time[i]])) # all variables for main correlation
        
        residY = lm.fit(cbind(Xpte[icte,],Txte[icte,tT==u_time[1]],1),Yte[icte,tY==u_time[i]])$residuals
        residT = lm.fit(cbind(Xpte[icte,],Txte[icte,tT==u_time[1]],1),Txte[icte,tT==tp])$residuals
        
        new_corN = cor(residY %*% mods[[k]]$alpha, residT)
        if (is.na(new_corN)){ new_corN = 0}
        mods[[k]]$main_correlations_testset = c(mods[[k]]$main_correlations_testset,new_corN)
        
        # corr_abs = 0
        corr_abs = c()
        for (j in 1:(ip-1)) {
          ic = complete.cases(cbind(Xp,Tx[,tT==tp],Y[,tY==u_time[i]]),Tx[,tT==u_time[j]]) # all variables for main correlation
          icte = complete.cases(cbind(Xpte,Txte[,tT==tp],Yte[,tY==u_time[i]]),Txte[,tT==u_time[j]]) # all variables for main correlation
          
          residY = lm.fit(cbind(Xpte[icte,],Txte[icte,tT==tp],1),Yte[icte,tY==u_time[i]])$residuals ## must compute on test set for proper partial correlation
          residT = lm.fit(cbind(Xpte[icte,],Txte[icte,tT==tp],1),Txte[icte,tT==u_time[j]])$residuals## must compute on test set for proper partial correlation
          
          mods[[k]]$alpha[abs(mods[[k]]$alpha) < 0.001]=0
          
          new_corr_abs = abs(cor(residY%*% mods[[k]]$alpha,residT))
          
          if (is.na(new_corr_abs)){ new_corr_abs = 0}
          
          corr_abs <- c(corr_abs, cor_to_pval(new_corr_abs,sum(icte)))
        }
        mods[[k]]$confounding_pvalues_testset = c(mods[[k]]$confounding_pvalues_testset,min(corr_abs))
        
      }
      
    }  
  }
  
  
  return(mods)
}
