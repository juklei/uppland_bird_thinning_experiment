## bird bpo model with a binomial process for detection
##
## First edit: 20191201
## Last edit: 20200428
##
## Author: Julian Klein

model{
  
  ## 1. Likelihood: ------------------------------------------------------------
  
  ## Observational model:
  for(i in 1:nobs){
    observed[i] ~ dbin(occ_true[species[i],year[i],site[i]]*pdet[i], nvisits[i])
    # sim[i] ~ dbin(occ_true[species[i],year[i],site[i]]*pdet[i], nvisits[i])
    logit(pdet[i]) <- a_pdet[species[i]] + OLRE[i] +
                      b_pdet_2018[species[i]]*ifelse(year[i]==2,1,0) +
                      b_pdet_2019[species[i]]*ifelse(year[i]==3,1,0) +
                      b_observer[species[i]]*observer[i]
    OLRE[i] ~ dnorm(0, 1/sd_OLRE^2) ## Observation level random effect
  }
  
  ## Ecological process model:
  for(k in 1:max(species)){
    for(y in 1:max(year)){
      for(p in 1:max(site)){
        occ_true[k,y,p] ~ dbern(pocc[k,y,p])
        logit(pocc[k,y,p]) <- a_pocc[k,treat[p],exp[y,p]] + 
                              e_block[k,block[p]] +
                              b_pocc_2018[k]*ifelse(y==2,1,0) +
                              b_pocc_2019[k]*ifelse(y==3,1,0) 
  }}}

  ## Group effects:
  for(k in 1:max(species)){
    for(b in 1:max(block)){e_block[k,b] ~ dnorm(0, 1/sd_block[k]^2)}
  }
  
  ## 2. Priors: ----------------------------------------------------------------
  
  ## Observational model:
  for(k in 1:max(species)){
    a_pdet[k] ~ dnorm(mu_a_pdet, 1/sd_a_pdet^2)
    b_pdet_2018[k] ~ dnorm(mu_b_pdet_2018, 1/sd_b_pdet_2018^2)
    b_pdet_2019[k] ~ dnorm(mu_b_pdet_2019, 1/sd_b_pdet_2019^2)
    b_observer[k] ~ dnorm(mu_obs, 1/sd_obs^2)
  }
  sd_OLRE ~ dt(0, pow(2.5,-2), 1)T(0,)

  ## Ecological process model:
  for(k in 1:max(species)){
    for(m in 1:max(treat)){
      for(n in 1:max(exp)){
        a_pocc[k,m,n] ~ dnorm(mu_a_pocc[m,n], 1/sd_a_pocc[m,n]^2)
      }}
    b_pocc_2018[k] ~ dnorm(mu_b_pocc_2018, 1/sd_b_pocc_2018^2)
    b_pocc_2019[k] ~ dnorm(mu_b_pocc_2019, 1/sd_b_pocc_2019^2)
    sd_block[k] ~ dt(0, pow(u_sd_block,-2), 1)T(0,)
  }
  
  ## Hyperpriors:
  
  ## Observational model:
  mu_a_pdet ~ dnorm(0, 0.1)
  sd_a_pdet ~ dt(0, pow(2.5,-2), 1)T(0,)
  mu_b_pdet_2018 ~ dnorm(0, 0.1)
  sd_b_pdet_2018 ~ dt(0, pow(2.5,-2), 1)T(0,)
  mu_b_pdet_2019 ~ dnorm(0, 0.1)
  sd_b_pdet_2019 ~ dt(0, pow(2.5,-2), 1)T(0,)
  mu_obs ~ dnorm(0, 0.1)
  sd_obs ~ dt(0, pow(2.5,-2), 1)T(0,)

  ## Ecological process model:
  for(m in 1:max(treat)){
    for(n in 1:max(exp)){
      mu_a_pocc[m,n] ~ dnorm(0, 0.1)
      sd_a_pocc[m,n] ~ dt(0, pow(2.5,-2), 1)T(0,)
    }}
  mu_b_pocc_2018 ~ dnorm(0, 0.1)
  sd_b_pocc_2018 ~ dt(0, pow(2.5,-2), 1)T(0,)
  mu_b_pocc_2019 ~ dnorm(0, 0.1)
  sd_b_pocc_2019 ~ dt(0, pow(2.5,-2), 1)T(0,)
  u_sd_block ~ dunif(0, 5)

  # ## 3. Model validation: ------------------------------------------------------
  # 
  # ## Bayesian p-value:
  # mean_obs <- mean(observed[])
  # mean_sim <- mean(sim[])
  # p_mean <- step(mean_sim - mean_obs)
  # 
  # ## Coefficient of variation:
  # cv_obs <- sd(observed[])/mean_obs
  # cv_sim <- sd(sim[])/mean_sim
  # p_cv <- step(cv_sim - cv_obs)
  # 
  # ## Model fit:
  # for(i in 1:nobs){
  #   sq[i] <- (observed[i] -
  #               occ_true[species[i],year[i],site[i]]*pdet[i]*nvisits[i])^2
  #   sq_sim[i] <- (sim[i] -
  #                   occ_true[species[i],year[i],site[i]]*pdet[i]*nvisits[i])^2
  # }
  # 
  # fit <- sum(sq[])
  # fit_sim <- sum(sq_sim[])
  # p_fit <- step(fit_sim - fit)

}


