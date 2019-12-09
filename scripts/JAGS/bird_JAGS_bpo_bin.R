## bird bpo model with a binomial process for detection
##
## First edit: 20191201
## Last edit: 20191209
##
## Author: Julian Klein

model{
  
  ## Likelihood: ---------------------------------------------------------------
  
  ## Observational model:
  for(i in 1:nobs){
    observed[i] ~ dbin(occ_true[species[i],year[i],site[i]]*pdet[i], nvisits[i])
    sim[i] ~ dbin(occ_true[species[i],year[i],site[i]]*pdet[i], nvisits[i])
    logit(pdet[i]) <- a_pdet[species[i]] + OLRE[i] +
                      # e_pdet_year[species[i],year[i]] +
                      b_pdet_2018[species[i]]*ifelse(year[i]==2,1,0) +
                      b_pdet_2019[species[i]]*ifelse(year[i]==3,1,0) +
                      b_observer[species[i]]*observer[i]
    OLRE[i] ~ dnorm(0, 1/sd_OLRE^2)
  }
  
  ## Ecological process model:
  for(k in 1:nspecies){
    for(y in 1:nyears){
      for(p in 1:nsites){
        occ_true[k,y,p] ~ dbern(pocc[k,y,p])
        logit(pocc[k,y,p]) <- a_pocc[k,treat[p],exp[y,p]] + ## Change according to what you want to do here!
                              b_pocc_2018[k]*ifelse(y==2,1,0) +
                              b_pocc_2019[k]*ifelse(y==3,1,0) #+
                              # e_year[k,y] + e_site[k,p]
  }}}
  
  # ## Group effects:
  # for(k in 1:nspecies){
  #   for(y in 1:nyears){e_pdet_year[k,y] ~ dnorm(0, 1/sd_pdet_year[k]^2)}
  #   for(y in 1:nyears){e_year[k,y] ~ dnorm(0, 1/sd_year[k]^2)}
  #   for(p in 1:nsites){e_site[k,p] ~ dnorm(0, 1/sd_site[k]^2)}
  # }
  
  ## Priors: -------------------------------------------------------------------
  
  ## Observational model:
  for(k in 1:nspecies){
    a_pdet[k] ~ dnorm(mu_a_pdet, 1/sd_a_pdet^2)
    b_pdet_2018[k] ~ dnorm(mu_b_pdet_2018, 1/sd_b_pdet_2018^2)
    b_pdet_2019[k] ~ dnorm(mu_b_pdet_2019, 1/sd_b_pdet_2019^2)
    # sd_pdet_year[k] ~ dt(0, pow(2.5,-2), 1)T(0,)
    b_observer[k] ~ dnorm(mu_obs, 1/sd_obs^2)
  }
  sd_OLRE ~ dt(0, pow(2.5,-2), 1)T(0,)


  ## Ecological process model:
  for(k in 1:nspecies){
    for(m in 1:max(treat)){
      for(n in 1:max(exp)){
        a_pocc[k,m,n] ~ dnorm(mu_a_pocc[m,n], 1/sd_a_pocc[m,n]^2)
      }}
    b_pocc_2018[k] ~ dnorm(mu_b_pocc_2018, 1/sd_b_pocc_2018^2)
    b_pocc_2019[k] ~ dnorm(mu_b_pocc_2019, 1/sd_b_pocc_2019^2)
    # sd_year[k] ~ dt(0, pow(2.5,-2), 1)T(0,)
    # sd_site[k] ~ dunif(0, u_sd_site)
  }
  
  ## Hyperpriors:
  
  ## Observational model:
  mu_a_pdet ~ dnorm(0, 0.1)
  sd_a_pdet ~ dt(0, pow(2.5,-2), 1)T(0,)
  mu_b_pdet_2018 ~ dnorm(0, 0.1)
  sd_b_pdet_2018 ~ dt(0, pow(2.5,-2), 1)T(0,)
  mu_b_pdet_2019 ~ dnorm(0, 0.1)
  sd_b_pdet_2019 ~ dt(0, pow(2.5,-2), 1)T(0,)
  # u_sd_pdet_year ~ dunif(0, 5)
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
  # u_sd_year ~ dunif(0, 5)
  # u_sd_site ~ dt(0,1,1)T(0,)

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

  ## BACI indicators for mean community response:

  for(m in 1:max(treat)){ ## Backtransform to logit
    for(n in 1:max(exp)){
      logit(BACI_cm[m,n]) <- mu_a_pocc[m,n]
  }}

  CI_div_C_cm <- abs(BACI_cm[1,2]-BACI_cm[3,2]) - abs(BACI_cm[1,1]-BACI_cm[3,1])
  CI_ctr_C_cm <- abs(BACI_cm[1,2]-BACI_cm[1,1]) - abs(BACI_cm[3,2]-BACI_cm[3,1])
  BACI_C_cm <- (BACI_cm[1,2]-BACI_cm[1,1]) - (BACI_cm[3,2]-BACI_cm[3,1])

  CI_div_T_cm <- abs(BACI_cm[2,2]-BACI_cm[3,2]) - abs(BACI_cm[2,1]-BACI_cm[3,1])
  CI_ctr_T_cm <- abs(BACI_cm[2,2]-BACI_cm[2,1]) - abs(BACI_cm[3,2]-BACI_cm[3,1])
  BACI_T_cm <- (BACI_cm[2,2]-BACI_cm[2,1]) - (BACI_cm[3,2]-BACI_cm[3,1])

  CI_div_URT_cm <- abs(BACI_cm[4,2]-BACI_cm[3,2]) - abs(BACI_cm[4,1]-BACI_cm[3,1])
  CI_ctr_URT_cm <- abs(BACI_cm[4,2]-BACI_cm[4,1]) - abs(BACI_cm[3,2]-BACI_cm[3,1])
  BACI_URT_cm <- (BACI_cm[4,2]-BACI_cm[4,1]) - (BACI_cm[3,2]-BACI_cm[3,1])

  ## BACI indicators for species specific response:

  for(k in 1:nspecies){ ## Backtransform to logit
    for(m in 1:max(treat)){
      for(n in 1:max(exp)){
        logit(BACI[k,m,n]) <- a_pocc[k,m,n]
    }}

    CI_div_C[k] <-  abs(BACI[k,1,2]-BACI[k,3,2]) - abs(BACI[k,1,1]-BACI[k,3,1])
    CI_ctr_C[k] <- abs(BACI[k,1,2]-BACI[k,1,1]) - abs(BACI[k,3,2]-BACI[k,3,1])
    BACI_C[k] <- (BACI[k,1,2]-BACI[k,1,1]) - (BACI[k,3,2]-BACI[k,3,1])

    CI_div_T[k] <- abs(BACI[k,2,2]-BACI[k,3,2]) - abs(BACI[k,2,1]-BACI[k,3,1])
    CI_ctr_T[k] <- abs(BACI[k,2,2]-BACI[k,2,1]) - abs(BACI[k,3,2]-BACI[k,3,1])
    BACI_T[k] <- (BACI[k,2,2]-BACI[k,2,1]) - (BACI[k,3,2]-BACI[k,3,1])

    CI_div_URT[k] <- abs(BACI[k,4,2]-BACI[k,3,2]) - abs(BACI[k,4,1]-BACI[k,3,1])
    CI_ctr_URT[k] <- abs(BACI[k,4,2]-BACI[k,4,1]) - abs(BACI[k,3,2]-BACI[k,3,1])
    BACI_URT[k] <- (BACI[k,4,2]-BACI[k,4,1]) - (BACI[k,3,2]-BACI[k,3,1])

  }

  ## BACI indicators for species richness:

  ## Simulate species occupancy per species*treatment*experiment combination:
  for(m in 1:max(treat)){
    for(n in 1:max(exp)){
      for(k in 1:nspecies){
        occ_BACI[k,m,n]  <- ifelse(BACI[k,m,n] > 0.5, 1, 0)
        }
      BACI_r[m,n] <- sum(occ_BACI[,m,n])
  }}

  CI_div_C_r <- abs(BACI_r[1,2]-BACI_r[3,2]) - abs(BACI_r[1,1]-BACI_r[3,1])
  CI_ctr_C_r <- abs(BACI_r[1,2]-BACI_r[1,1]) - abs(BACI_r[3,2]-BACI_r[3,1])
  BACI_C_r <- (BACI_r[1,2]-BACI_r[1,1]) - (BACI_r[3,2]-BACI_r[3,1])

  CI_div_T_r <- abs(BACI_r[2,2]-BACI_r[3,2]) - abs(BACI_r[2,1]-BACI_r[3,1])
  CI_ctr_T_r <- abs(BACI_r[2,2]-BACI_r[2,1]) - abs(BACI_r[3,2]-BACI_r[3,1])
  BACI_T_r <- (BACI_r[2,2]-BACI_r[2,1]) - (BACI_r[3,2]-BACI_r[3,1])

  CI_div_URT_r <- abs(BACI_r[4,2]-BACI_r[3,2]) - abs(BACI_r[4,1]-BACI_r[3,1])
  CI_ctr_URT_r <- abs(BACI_r[4,2]-BACI_r[4,1]) - abs(BACI_r[3,2]-BACI_r[3,1])
  BACI_URT_r <- (BACI_r[4,2]-BACI_r[4,1]) - (BACI_r[3,2]-BACI_r[3,1])

  ## BACI indicators for beta diversity according to Whittaker:

  for(n in 1:max(exp)){
    for(k in 1:nspecies){
      BACI_gamma[k,n] <- ifelse(sum(occ_BACI[k,,n]) > 0, 1, 0)
    }
    BACI_gd[n] <- sum(BACI_gamma[,n])
    for(m in 1:max(treat)){
      BACI_bd[m,n] <- BACI_gd[n]/ifelse(BACI_r[m,n] < 1, 1, BACI_r[m,n])
  }}

  CI_div_C_bd <- abs(BACI_bd[1,2]-BACI_bd[3,2]) - abs(BACI_bd[1,1]-BACI_bd[3,1])
  CI_ctr_C_bd <- abs(BACI_bd[1,2]-BACI_bd[1,1]) - abs(BACI_bd[3,2]-BACI_bd[3,1])
  BACI_C_bd <- (BACI_bd[1,2]-BACI_bd[1,1]) - (BACI_bd[3,2]-BACI_bd[3,1])

  CI_div_T_bd <- abs(BACI_bd[2,2]-BACI_bd[3,2]) - abs(BACI_bd[2,1]-BACI_bd[3,1])
  CI_ctr_T_bd <- abs(BACI_bd[2,2]-BACI_bd[2,1]) - abs(BACI_bd[3,2]-BACI_bd[3,1])
  BACI_T_bd <- (BACI_bd[2,2]-BACI_bd[2,1]) - (BACI_bd[3,2]-BACI_bd[3,1])

  CI_div_URT_bd <- abs(BACI_bd[4,2]-BACI_bd[3,2]) - abs(BACI_bd[4,1]-BACI_bd[3,1])
  CI_ctr_URT_bd <- abs(BACI_bd[4,2]-BACI_bd[4,1]) - abs(BACI_bd[3,2]-BACI_bd[3,1])
  BACI_URT_bd <- (BACI_bd[4,2]-BACI_bd[4,1]) - (BACI_bd[3,2]-BACI_bd[3,1])
  
}


