## insect index on sticky traps evaluated in BACI experiment
##
## Response: day pm with max increment
##
## First edit: 2020225
## Last edit: 2020225
##
## Author: Julian Klein

model{
  
  ## Likelihood: ---------------------------------------------------------------
  
  ## Process model:
  for(i in 1:nobs) {
    ## Stochastic model:
    pm[i] ~ dgamma(shape[i], rate[i])
    sim[i] ~ dgamma(shape[i], rate[i])
    ## Moment matching:
    shape[i] <- max(0.00001, mu[i]^2/sigma^2)
    rate[i] <- max(0.00001, mu[i]/sigma^2)
    ## Deterministic model:
    mu[i] <- a_cov[treat[i],exp[i]] + #e_site[site[i]] +
             b_2018*year_2018[i] + b_2019*year_2019[i]
  }
  
  # for(j in 1:nsites){e_site[j] ~ dnorm(0, 1/sigma_site^2)}
  
  ## Priors: -------------------------------------------------------------------
  
  ## Process model:
  for(m in 1:max(treat)){
    for(n in 1:max(exp)){
      a_cov[m,n] ~ dunif(1, 100)
  }}
  sigma ~ dgamma(0.001, 0.001)
  b_2018 ~ dnorm(0, 0.001)
  b_2019 ~ dnorm(0, 0.001)
  # sigma_site ~ dgamma(0.1, 0.1)
  
  ## Model validation: ---------------------------------------------------------

  ## Bayesian p-value:
  mean_obs <- mean(pm[])
  mean_sim <- mean(sim[])
  p_mean <- step(mean_sim - mean_obs)

  ## Coefficient of variation:
  cv_obs <- sd(pm[])/mean_obs
  cv_sim <- sd(sim[])/mean_sim
  p_cv <- step(cv_sim - cv_obs)

  ## Model fit:
  for(i in 1:nobs){
    sq[i] <- (pm[i] - shape[i]/rate[i])^2
    sq_sim[i] <- (sim[i] - shape[i]/rate[i])^2
  }

  fit <- sum(sq[])
  fit_sim <- sum(sq_sim[])
  p_fit <- step(fit_sim - fit)

  ## Posteriors: ---------------------------------------------------------------

  ## BACI indicators:
  for(o in eval){
    CI_div[o] <- abs(a_cov[o,2]-a_cov[ref,2]) - abs(a_cov[o,1]-a_cov[ref,1])
    CI_ctr[o] <- abs(a_cov[o,2]-a_cov[o,1]) - abs(a_cov[ref,2]-a_cov[ref,1])
    BACI[o] <- (a_cov[o,2]-a_cov[o,1]) - (a_cov[ref,2]-a_cov[ref,1])
  }
  
}

