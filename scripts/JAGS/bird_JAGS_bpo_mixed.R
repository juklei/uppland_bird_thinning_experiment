## bird bpo model with a bernoulli process for occurrence
##
## First edit: 20191030
## Last edit: 20191106
##
## Author: Julian Klein

model{
  
  ## Likelihood:
  
  ## Observational model:
  for(i in 1:nobs){
    observed[i] ~ dbern(occ_true[species[i],year[i],site[i]]*pdet[i])
    logit(pdet[i]) <- a_pdet[species[i]] + b_observer*observer[i] +
                      b_dpm[ldm[i]]*dpm[i,1] + b2_dpm[ldm[i]]*dpm[i,1]^2 +
                      b_mps*mps[i,1] + b2_mps*mps[i,1]^2
  }
  
  ## Ecological process model:
  for(k in 1:nspecies){
    for(y in 1:nyears){
      for(p in 1:nsites){
        occ_true[k,y,p] ~ dbern(pocc[k,y,p])
        logit(pocc[k,y,p]) <- a_pocc[k] + e_site[k,p] + #e_year[k,y] +
                              b_2018[k]*ifelse(y==2,1,0) + 
                              b_2019[k]*ifelse(y==3,1,0) + 
                              b_thinned[k]*thinned[y,p] +
                              b_control[k]*control[y,p] +
                              b_sdbh[k]*sdbh[y,p] +
                              b_sdbh_t[k]*thinned[y,p]*sdbh[y,p] +
                              b_sdbh_c[k]*control[y,p]*sdbh[y,p]
  }}}
  
  ## Group effects:
  for(k in 1:nspecies){
    # for(y in 1:nyears){e_year[k,y] ~ dnorm(0, 1/sd_year[k]^2)} 
    for(p in 1:nsites){e_site[k,p] ~ dnorm(0, 1/sd_site[k]^2)}
  }
  
  ## Priors:
  
  ## Observational model:
  for(k in 1:nspecies){a_pdet[k] ~ dnorm(mu_a_pdet, 1/sd_a_pdet^2)}
  b_observer ~ dnorm(0, 0.01)
  for(l in 1:max(ldm)){    
    b_dpm[l] ~ dnorm(0, 0.01)
    b2_dpm[l] ~ dnorm(0, 0.01)
  }
  b_mps ~ dnorm(0, 0.01)
  b2_mps ~ dnorm(0, 0.01)
  
  ## Ecological process model:
  for(k in 1:nspecies){
    a_pocc[k] ~ dnorm(mu_a_pocc, 1/sd_a_pocc^2)
    # sd_year[k] ~ dunif(0, u_sd_year)
    sd_site[k] ~ dunif(0, u_sd_site)
    b_2018[k] ~ dnorm(mu_b_2018, 1/sd_b_2018^2)
    b_2019[k] ~ dnorm(mu_b_2019, 1/sd_b_2019^2)
    b_thinned[k] ~ dnorm(mu_b_thinned, 1/sd_b_thinned^2)
    b_control[k] ~ dnorm(mu_b_control, 1/sd_b_control^2)
    b_sdbh[k] ~ dnorm(mu_b_sdbh, 1/sd_b_sdbh^2)
    b_sdbh_t[k] ~ dnorm(mu_b_sdbh_t, 1/sd_b_sdbh_t^2)
    b_sdbh_c[k]~ dnorm(mu_b_sdbh_c, 1/sd_b_sdbh_c^2)
  }
  
  ## Hyperpriors:
  
  ## Observational model:
  mu_a_pdet ~ dnorm(0, 0.01)
  sd_a_pdet ~ dunif(0, 5)
  # mu_b_dpm ~ dnorm(0, 0.01)
  # sd_b_dpm ~ dunif(0, 10)
  # mu_b2_dpm ~ dnorm(0, 0.01)
  # sd_b2_dpm ~ dunif(0, 10)
  # mu_b_mps ~ dnorm(0, 0.01)
  # sd_b_mps ~ dunif(0, 10)
  # mu_b2_mps ~ dnorm(0, 0.01)
  # sd_b2_mps ~ dunif(0, 10)
  
  ## Ecological process model:
  mu_a_pocc ~ dnorm(0, 0.01)
  sd_a_pocc ~ dunif(0, 5)
  # u_sd_year ~ dunif(0, 10)
  u_sd_site ~ dunif(0, 10)
  mu_b_2018 ~ dnorm(0, 0.1)
  sd_b_2018 ~ dunif(0, 5)
  mu_b_2019 ~ dnorm(0, 0.1)
  sd_b_2019 ~ dunif(0, 5)
  mu_b_thinned ~ dnorm(0, 0.1)
  sd_b_thinned ~ dunif(0, 5)
  mu_b_control ~ dnorm(0, 0.1)
  sd_b_control ~ dunif(0, 5)
  mu_b_sdbh ~ dnorm(0, 0.1)
  sd_b_sdbh ~ dunif(0, 5)
  mu_b_sdbh_t ~ dnorm(0, 0.1)
  sd_b_sdbh_t ~ dunif(0, 5)
  mu_b_sdbh_c ~ dnorm(0, 0.1) 
  sd_b_sdbh_c ~ dunif(0, 5)
  
  ## Model validation:
  
  ## Predictions:

  
}


