DEBIAS <- function(Y,Tx,tp,Xp,Yte,Txte,Xpte,lambdas = c(0,1,2,3,4,5,6,7,8,9,10),nf=5, max_iterk = 3){
  #
  # runs the full DEBIAS algorithm, including cross-validation of lambda
  #
  # Inputs:
  #   Y: outcome items, where column names should contain numeric time points. The number of outcome items per time point must be the same.
  #   Tx: treatment assignments, where column names should again contain numeric time points
  #   tp: time point of interest (variable p in the paper)
  #   Xp: covariates. Must not contain any missing values
  #   For outer cross-validation:
  #     Yte: outcome items for test/validation set evaluation (default: NULL)
  #     Txte: treatment assignment for test/validation set evaluation (default: NULL)
  #     Xpte: outcome items for test/validation set evaluation (default: NULL)
  #   lambdas: grid of lambda values for cross-validation (default: 0,1,...,10)
  #   nf: number of cross-validation folds (default: 5)
  #   max_iterk: maximum number of non-negative weights to extract (variable s in the paper)
  #
  # Outputs:
  #   mods: a list of max_iterk models with the following per model:
  #     alpha: vector of non-negative weights
  #     main_correlations: main correlation term (a) for each time point
  #     confounding_correlations: square root of confounding penalty term (b) (i.e., its absolute value of correlation instead of squared correlation) for each time point
  #     confounding_pvalues: p-value of confounding penalty term (b) for each time point
  #     lambda: optimal lambda value
  #     main_correlations_testset: main correlation terms (a) on test/validation set
  #     confounding_pvalues_testset: confounding penalty (b) p-values on test/validation set
  #
  # Written by Eric V. Strobl 06/2025
  #
  
  folds = rep(1:nf,length.out=nrow(Y))
  nL = length(lambdas)
  
  cor_vals = matrix(0,nL)
  conf_vals = matrix(0,nL)
  for (f in 1:nf){ # cross-validation
    itr = (folds!=f) # training fold
    ite = (folds==f) # validation fold
    
    for (l in 1:nL){
      mods = DEBIAS_sequential(Y[itr,],Tx[itr,],tp,Xp[itr,],
                               Y[ite,],Tx[ite,],Xp[ite,],
                               lambda=lambdas[l], max_iterk = max_iterk) # run DEBIAS max_iterk times, with evaluation on validation folds
      conf_all = c()
      for (k in 1:max_iterk){
        cor_vals[l] = cor_vals[l] + mean(mods[[k]]$main_correlations_testset)/(nf*max_iterk)
        conf_all = c(conf_all, min(log(mods[[k]]$confounding_pvalues_testset+1E-10)))
      }
      conf_vals[l] = conf_vals[l] + mean(conf_all)/nf
    }
  }
  
  # indices where geometric mean of minimum p-values is above exp(-3) = 0.05
  ok = which(conf_vals > -3)
  
  if (length(ok) > 0) {
    # among these, pick the index of max cor_vals
    best_idx = ok[which.max(cor_vals[ok])]
  } else {
    # if none, pick the index of conf_vals closest to -3
    best_idx = which.min(abs(conf_vals - (-3)))
  }
  
  # run DEBIAS sequentially max_iterk times on optimal lambda
  mods = DEBIAS_sequential(Y,Tx,tp,Xp, Yte,Txte,Xpte, lambda=lambdas[best_idx], max_iterk = max_iterk)
  
  return(mods)
  
}
