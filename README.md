# DEBIAS (Durable Effects with Backdoor-Invariant Aggregated Symptoms)

Causal inference in longitudinal biomedical data remains a central challenge, especially in psychiatry, where symptom heterogeneity and latent confounding frequently undermine classical estimators. Most existing methods for treatment effect estimation presuppose a fixed outcome variable and address confounding through observed covariate adjustment. However, the assumption of unconfoundedness may not hold for a fixed outcome in practice. To address this foundational limitation, DEBIAS directly optimizes the outcome definition to maximize causal identifiability. The algorithm learns non-negative, clinically interpretable weights for outcome aggregation, maximizing durable treatment effects and empirically minimizing both observed and latent confounding by leveraging the time-limited direct effects of prior treatments in psychiatric longitudinal data.

# Installation

> library(devtools)

> install_github("ericstrobl/DEBIAS")

> library(DEBIAS)

# Inputs

`Y` = outcome items, where column names should contain numeric time points. The number of outcome items per time point must be the same.

`Tx` = treatment assignments, where column names should again contain numeric time points

`tp` = time point of interest (variable p in the paper)
 
`Xp` = covariates. Must not contain any missing values

Optional:

`Yte` = outcome items for test/validation set evaluation (default: NULL)
  
`Txte` = treatment assignment for test/validation set evaluation (default: NULL)
  
`Xpte` = outcome items for test/validation set evaluation (default: NULL)
  
`lambdas` = grid of lambda values for cross-validation (default: 0,1,...,10)
  
`nf` = number of cross-validation folds (default: 5)
  
`max_iterk` = maximum number of non-negative weights to extract (default: 3, variable s in the paper)

# Run the Algorithm

> data = generate_synthetic_longitudinal_data(nsamps=1000) # generate synthetic longitudinal data

> mods = DEBIAS(Y,X,time,k=5) # run DEBIAS

# Outputs

A list of `max_iterk` models with the following per model:

`alpha` = vector of non-negative weights

`main_correlations` = main correlation term (a) for each time point

`confounding_correlations` = square root of confounding penalty term (b) for each time point

`confounding_correlations` = square root of confounding penalty term (b) (i.e., absolute value of correlation instead of squared correlation) for each time point

`confounding_pvalues` = p-value of confounding penalty term (b) for each time point

`Mahalanobis_cosine_similarity` = Mahalanobis cosine similarity for orthogonality penalty (c)

`lambda` = optimal lambda value

`main_correlations_testset` = main correlation term (a) on test/validation set (if present)

`confounding_pvalues_testset` = confounding penalty (b) p-value on test/validation set (if present)



