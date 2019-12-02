## bird bpo model with a binomial process for detection
##
## First edit: 20191201
## Last edit: 20191201
##
## Author: Julian Klein

model{
  
  ## Likelihood:
  
  ## Observational model:
  for(i in 1:nobs){
    observed[i] ~ dbin(occ_true[species[i],year[i],site[i]]*pdet[i], nvisits[i])
    logit(pdet[i]) <- a_pdet[species[i]] + e_pdet_year[species[i],year[i]] +
                      b_observer*observer[i]
  }
  
  ## Ecological process model:
  for(k in 1:nspecies){
    for(y in 1:nyears){
      for(p in 1:nsites){
        occ_true[k,y,p] ~ dbern(pocc[k,y,p])
        logit(pocc[k,y,p]) <- a_pocc[k,treat[p],exp[y,p]] + 
                              e_site[k,p] + e_year[k,y]
  }}}
  
  ## Group effects:
  for(k in 1:nspecies){
    for(y in 1:nyears){e_pdet_year[k,y] ~ dnorm(0, 1/sd_pdet_year[k]^2)}
    for(y in 1:nyears){e_year[k,y] ~ dnorm(0, 1/sd_year[k]^2)}
    for(p in 1:nsites){e_site[k,p] ~ dnorm(0, 1/sd_site[k]^2)}
  }
  
  ## Priors:
  
  ## Observational model:
  for(k in 1:nspecies){
    a_pdet[k] ~ dnorm(mu_a_pdet, 1/sd_a_pdet^2)
    sd_pdet_year[k] ~ dunif(0, u_sd_pdet_year)
  }
  b_observer ~ dnorm(0, 0.01)

  ## Ecological process model:
  for(k in 1:nspecies){
    for(m in 1:max(treat)){
      for(n in 1:max(exp)){
        a_pocc[k,m,n] ~ dnorm(mu_a_pocc[m,n], 1/sd_a_pocc[m,n]^2)
    }}
    sd_year[k] ~ dunif(0, u_sd_year)
    sd_site[k] ~ dunif(0, u_sd_site)
  }
  
  ## Hyperpriors:
  
  ## Observational model:
  mu_a_pdet ~ dnorm(0, 0.01)
  sd_a_pdet ~ dunif(0, 5)
  u_sd_pdet_year ~ dt(0,1,1)T(0,10)
  
  ## Ecological process model:
  for(m in 1:max(treat)){
    for(n in 1:max(exp)){
      mu_a_pocc[m,n] ~ dnorm(0, 0.01)
      sd_a_pocc[m,n]~ dunif(0, 5)
  }}
  u_sd_year ~ dt(0,1,1)T(0,10)
  u_sd_site ~ dt(0,1,1)T(0,10)

  ## Model validation:
  
  ## Posteriors:
  
  ## BACI indicators:
  
  ## Posterior distributions of the three measurements of impact:
  
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
  
  ## For mean community response: 
  
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

  ## For species richness: 
  
  ##...
  
  ## For beta diversirty:
  
  ##...
  
}


