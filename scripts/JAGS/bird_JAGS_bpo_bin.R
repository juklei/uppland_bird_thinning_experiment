## bird bpo model with a binomial process for detection
##
## First edit: 20191201
## Last edit: 20200325
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
  for(k in 1:max(species)){
    for(y in 1:max(year)){
      for(p in 1:max(site)){
        occ_true[k,y,p] ~ dbern(pocc[k,y,p])
        logit(pocc[k,y,p]) <- a_pocc[k,treat[p],exp[y,p]] + 
                              b_pocc_2018[k]*ifelse(y==2,1,0) +
                              b_pocc_2019[k]*ifelse(y==3,1,0) 
  }}}
  
  ## Priors: -------------------------------------------------------------------
  
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
      logit(mu_pocc_real[m,n]) <- mu_a_pocc[m,n]
  }}
  for(o in eval){
    CI_div_cm[o] <- abs(mu_pocc_real[o,2]-mu_pocc_real[ref,2]) -
                    abs(mu_pocc_real[o,1]-mu_pocc_real[ref,1])
    CI_ctr_cm[o] <- abs(mu_pocc_real[o,2]-mu_pocc_real[o,1]) -
                    abs(mu_pocc_real[ref,2]-mu_pocc_real[ref,1])
    BACI_cm[o] <- (mu_pocc_real[o,2]-mu_pocc_real[o,1]) -
                  (mu_pocc_real[ref,2]-mu_pocc_real[ref,1])
  }

  ## BACI indicators for species specific response:
  for(k in 1:max(species)){ ## Backtransform to logit
    for(m in 1:max(treat)){
      for(n in 1:max(exp)){
        logit(pocc_real[k,m,n]) <- a_pocc[k,m,n]
    }}
    for(o in eval){
      CI_div_sl[o,k] <- abs(pocc_real[k,o,2]-pocc_real[k,ref,2]) -
                        abs(pocc_real[k,o,1]-pocc_real[k,ref,1])
      CI_ctr_sl[o,k] <- abs(pocc_real[k,o,2]-pocc_real[k,o,1]) -
                        abs(pocc_real[k,ref,2]-pocc_real[k,ref,1])
      BACI_sl[o,k] <- (pocc_real[k,o,2]-pocc_real[k,o,1]) -
                      (pocc_real[k,ref,2]-pocc_real[k,ref,1])
  }}
  
  ## Compare the pre-experiment values of true control and mean treatments:
  for(k in 1:max(species)){
    compare_before[k] <- pocc_real[k,ref,1] - mean(pocc_real[k,eval,1])
  }
  
  ## BACI indicators for food type group response:
  for(o in eval){
    CI_div_insect[o] <- mean(CI_div_sl[o,insect])
    CI_ctr_insect[o] <- mean(CI_ctr_sl[o,insect])
    BACI_insect[o] <- mean(BACI_sl[o,insect])
    CI_div_omni[o] <- mean(CI_div_sl[o,omni])
    CI_ctr_omni[o] <- mean(CI_ctr_sl[o,omni])
    BACI_omni[o] <- mean(BACI_sl[o,omni])
  }

  ## BACI indicators for foraging group response:
  for(o in eval){
    CI_div_bark[o] <- mean(CI_div_sl[o,bark])
    CI_ctr_bark[o] <- mean(CI_ctr_sl[o,bark])
    BACI_bark[o] <- mean(BACI_sl[o,bark])
    CI_div_f_cpy[o] <- mean(CI_div_sl[o,f_cpy])
    CI_ctr_f_cpy[o] <- mean(CI_ctr_sl[o,f_cpy])
    BACI_f_cpy[o] <- mean(BACI_sl[o,f_cpy])
    CI_div_f_grd[o] <- mean(CI_div_sl[o,f_grd])
    CI_ctr_f_grd[o] <- mean(CI_ctr_sl[o,f_grd])
    BACI_f_grd[o] <- mean(BACI_sl[o,f_grd])
    CI_div_grd_cpy[o] <- mean(CI_div_sl[o,grd_cpy])
    CI_ctr_grd_cpy[o] <- mean(CI_ctr_sl[o,grd_cpy])
    BACI_grd_cpy[o] <- mean(BACI_sl[o,grd_cpy])
  }

  ## BACI indicators for nesting group response:
  for(o in eval){
    CI_div_n_grd[o] <- mean(CI_div_sl[o,n_grd])
    CI_ctr_n_grd[o] <- mean(CI_ctr_sl[o,n_grd])
    BACI_n_grd[o] <- mean(BACI_sl[o,n_grd])
    CI_div_n_cpy[o] <- mean(CI_div_sl[o,n_cpy])
    CI_ctr_n_cpy[o] <- mean(CI_ctr_sl[o,n_cpy])
    BACI_n_cpy[o] <- mean(BACI_sl[o,n_cpy])
    CI_div_hole[o] <- mean(CI_div_sl[o,hole])
    CI_ctr_hole[o] <- mean(CI_ctr_sl[o,hole])
    BACI_hole[o] <- mean(BACI_sl[o,hole])
  }

  ## BACI indicators for species richness:
  ## Predict richness per treatment*experiment combination:
  for(m in 1:max(treat)){
    for(n in 1:max(exp)){
      rich[m,n] <- sum(pocc_real[,m,n])
  }}
  for(o in eval){
    CI_div_r[o] <- abs(rich[o,2]-rich[ref,2]) - abs(rich[o,1]-rich[ref,1])
    CI_ctr_r[o] <- abs(rich[o,2]-rich[o,1]) - abs(rich[ref,2]-rich[ref,1])
    BACI_r[o] <- (rich[o,2]-rich[o,1]) - (rich[ref,2]-rich[ref,1])
  }

  ## BACI indicators for beta diversity (Jaccard index):
  ## Simulate 2 plots for each m*n:
  for(q in 1:2){
    for(m in 1:max(treat)){
      for(n in 1:max(exp)){
        for(k in 1:max(species)){
          sim_occ[k,m,n,q] ~ dbern(pocc_real[k,m,n])
  }}}}
  ## Calculate for every m*n the jaccard index:
  for(m in 1:max(treat)){
    for(n in 1:max(exp)){
      JI[m,n] <- sum(sim_occ[,m,n,1]*sim_occ[,m,n,2])/
                 (ifelse(sum(sim_occ[,m,n,1]) == 0, 1, sum(sim_occ[,m,n,1])) +
                 sum(sim_occ[,m,n,2]) -
                 sum(sim_occ[,m,n,1]*sim_occ[,m,n,2]))
      negJI[m,n] <- 1 - JI[m,n]
  }}
  for(o in eval){
    CI_div_bd[o] <- abs(negJI[o,2]-negJI[ref,2]) - abs(negJI[o,1]-negJI[ref,1])
    CI_ctr_bd[o] <- abs(negJI[o,2]-negJI[o,1]) - abs(negJI[ref,2]-negJI[ref,1])
    BACI_bd[o] <- (negJI[o,2]-negJI[o,1]) - (negJI[ref,2]-negJI[ref,1])
  }
  
}


