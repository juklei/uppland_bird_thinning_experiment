## bird bpo model with a binomial process for detection
##
## First edit: 20191218
## Last edit: 20191218
##
## Author: Julian Klein

model{
  
  ## Likelihood: ---------------------------------------------------------------
  
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
  for(k in 1:nspecies){
    for(y in 1:nyears){
      for(p in 1:nsites){
        occ_true[k,y,p] ~ dbern(pocc[k,y,p])
        logit(pocc[k,y,p]) <- a_pocc[k] + 
                              b_pocc_2018[k]*ifelse(y==2,1,0) +
                              b_pocc_2019[k]*ifelse(y==3,1,0) +
                              #e_block[k,block[p]] +
                              b_umbr[k]*umbr[p,1] +
                              b_BA[k]*BA[p,1] +
                              b_BA_spr[k]*BA_spruce[p,1] +
                              b_BA_dec[k]*BA_dec[p,1] +
                              b_BA_dw[k]*BA_dw[p,1] +
                              b_vis[k]*vis[p,1]
  }}}
  
  # ## Group effects:
  # for(k in 1:nspecies){
  #   for(b in 1:nblocks){e_block[k,b] ~ dnorm(0, 1/sd_block[k]^2)}
  # }
  
  ## Priors: -------------------------------------------------------------------
  
  ## Observational model:
  for(k in 1:nspecies){
    a_pdet[k] ~ dnorm(mu_a_pdet, 1/sd_a_pdet^2)
    b_pdet_2018[k] ~ dnorm(mu_b_pdet_2018, 1/sd_b_pdet_2018^2)
    b_pdet_2019[k] ~ dnorm(mu_b_pdet_2019, 1/sd_b_pdet_2019^2)
    b_observer[k] ~ dnorm(mu_obs, 1/sd_obs^2)
  }
  sd_OLRE ~ dt(0, pow(2.5,-2), 1)T(0,)


  ## Ecological process model:
  for(k in 1:nspecies){
    a_pocc[k] ~ dnorm(mu_a_pocc, 1/sd_a_pocc^2)
    b_pocc_2018[k] ~ dnorm(mu_b_pocc_2018, 1/sd_b_pocc_2018^2)
    b_pocc_2019[k] ~ dnorm(mu_b_pocc_2019, 1/sd_b_pocc_2019^2)
    # sd_block[k] ~ dt(0, pow(u_sd_block,-2), 1)T(0,)
    b_umbr[k] ~ dnorm(mu_b_umbr, 1/sd_b_umbr^2)
    b_BA[k] ~ dnorm(mu_b_BA, 1/sd_b_BA^2)
    b_BA_spr[k] ~ dnorm(mu_b_BA_spr, 1/sd_b_BA_spr^2)
    b_BA_dec[k] ~ dnorm(mu_b_BA_dec, 1/sd_b_BA_dec^2)
    b_BA_dw[k] ~ dnorm(mu_b_BA_dw, 1/sd_b_BA_dw^2)
    b_vis[k] ~ dnorm(mu_b_vis, 1/sd_b_vis^2)
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
  mu_a_pocc ~ dnorm(0, 0.1)
  sd_a_pocc ~ dt(0, pow(2.5,-2), 1)T(0,)
  mu_b_pocc_2018 ~ dnorm(0, 0.1)
  sd_b_pocc_2018 ~ dt(0, pow(2.5,-2), 1)T(0,)
  mu_b_pocc_2019 ~ dnorm(0, 0.1)
  sd_b_pocc_2019 ~ dt(0, pow(2.5,-2), 1)T(0,)
  # u_sd_block ~ dunif(0, 5)
  mu_b_umbr ~ dnorm(0, 0.1)
  sd_b_umbr ~ dt(0, pow(2.5,-2), 1)T(0,)
  mu_b_BA ~ dnorm(0, 0.1)
  sd_b_BA ~ dt(0, pow(2.5,-2), 1)T(0,)
  mu_b_BA_spr ~ dnorm(0, 0.1)
  sd_b_BA_spr ~ dt(0, pow(2.5,-2), 1)T(0,)
  mu_b_BA_dec ~ dnorm(0, 0.1)
  sd_b_BA_dec ~ dt(0, pow(2.5,-2), 1)T(0,)
  mu_b_BA_dw ~ dnorm(0, 0.1)
  sd_b_BA_dw ~ dt(0, pow(2.5,-2), 1)T(0,)
  mu_b_vis ~ dnorm(0, 0.1)
  sd_b_vis ~ dt(0, pow(2.5,-2), 1)T(0,)
  
  # ## Model validation: ---------------------------------------------------------
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

  ## Posteriors: ---------------------------------------------------------------


}


