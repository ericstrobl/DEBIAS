# DEBIAS (Durable Effects with Backdoor-Invariant Aggregated Symptom)

Causal inference in longitudinal biomedical data remains a central challenge, especially in psychiatry, where symptom heterogeneity and latent confounding frequently undermine classical estimators. Most existing methods for treatment effect estimation presuppose a fixed outcome variable and address confounding through observed covariate adjustment. However, the assumption of unconfoundedness may not hold for a fixed outcome in practice. To address this foundational limitation, DEBIAS directly optimize sthe outcome definition to maximize causal identifiability. The algorithm learns non-negative, clinically interpretable weights for outcome aggregation, maximizing durable treatment effects and empirically minimizing both observed and latent confounding by leveraging the time-limited direct effects of prior treatments in psychiatric longitudinal data.

# Installation

> library(devtools)

> install_github("ericstrobl/DEBIAS")

> library(DEBIAS)

# Inputs

`X` = $n$ by $p$ matrix of $p$ predictors

`Y` = $n$ by $q$ matrix of $q$ rating scale items

`time` = column vector of length $n$ with follow-up times

`k` (optional) = number of summary scores to extract, default is $k=p$

# Run the Algorithm

> Y = matrix(runif(1000),100,10) # generate synthetic values of the clinical rating scale items

> X = matrix(rnorm(1000),100,10) # generate synthetic values of the predictors

> time = rep(c(1,2,3,4),length.out=100) # generate synthetic time points

> mods = SCORE(Y,X,time,k=5) # run SCORE

# Outputs

`mods` = list of $k$ models, where:

`$alpha` = length $q$ vector of non-negative weights of `Y`

`$beta` = $p$ by $m$ matrix of the (negative, zero or positive) weights of `X` 


