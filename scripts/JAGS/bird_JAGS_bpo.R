## bird bpo model with a bernoulli process for occurrence
##
## First edit: 20191030
## Last edit: 20191106
##
## Author: Julian Klein

model{
  
  ## Likelihood:
  
  ## Observational model:
  for(k in 1:nspecies){
    for(y in 1:nyears){
      for(i in 1:nsites){
        for(n in 1:nvisits){
          seen[k,y,i,n] ~ dbern(occ_true[k,y,i]*pdet[k,y,i,n])
          logit(pdet[k,y,i,n]) <- a_pdet[k] + 
                                  b_dpm[ldm[k]]*dpm[k,y,i,n] + 
                                  b2_dpm[ldm[k]]*dpm[k,y,i,n]^2 +
                                  b_mps[ldm[k]]*mps[k,y,i,n] +
                                  b2_mps[ldm[k]]*mps[k,y,i,n]^2
  }}}}
  
  ## Ecological process model:
  for(k in 1:nspecies){
    for(y in 1:nyears){
      for(i in 1:nsites){
        occ_true[k,y,i] ~ dbern(pocc[k,y,i])
        logit(pocc[k,y,i]) <- a_pocc[k] + e_year[k,y] + e_site[k,i] +
                              b_thinned[k]*thinned[y,i] +
                              b_control[k]*control[y,i] +
                              b_sdbh[k]*sdbh[y,i] +
                              b_sdbh_t[k]*thinned[y,i]*sdbh[y,i] +
                              b_sdbh_c[k]*control[y,i]*sdbh[y,i]
  }}}
  
  ## Group effects:
  for(k in 1:nspecies){
    for(y in 1:nyears){e_year[k,y] ~ dunif(0, sd_year[k])} 
    for(y in 1:nsites){e_site[k,y] ~ dunif(0, sd_site[k])}
  }
  
  ## Priors:
  
  ## Observational model:
  for(k in 1:nspecies){a_pdet[k] ~ dnorm(mu_a_pdet, 1/sd_a_pdet^2)}
  for(l in 1:max(ldm)){    
    b_dpm[l] ~ dnorm(0, 0.01) #dnorm(mu_b_dpm, 1/sd_b_dpm^2)
    b2_dpm[l] ~ dnorm(0, 0.01) #dnorm(mu_b2_dpm, 1/sd_b2_dpm^2)
    b_mps[l] ~ dnorm(0, 0.01) #dnorm(mu_b_mps, 1/sd_b_mps^2)
    b2_mps[l] ~ dnorm(0, 0.01) #dnorm(mu_b2_mps, 1/sd_b2_mps^2)
  }
  
  ## Ecological process model:
  for(k in 1:nspecies){
    a_pocc[k] ~ dnorm(mu_a_pocc, 1/sd_a_pocc^2)
    sd_year[k] ~ dunif(0, u_sd_year)
    sd_site[k] ~ dunif(0, u_sd_site)
    b_thinned[k] ~ dnorm(mu_b_thinned, 1/sd_b_thinned^2)
    b_control[k] ~ dnorm(mu_b_control, 1/sd_b_control^2)
    b_sdbh[k] ~ dnorm(mu_b_sdbh, 1/sd_b_sdbh^2)
    b_sdbh_t[k] ~ dnorm(mu_b_sdbh_t, 1/sd_b_sdbh_t^2)
    b_sdbh_c[k]~ dnorm(mu_b_sdbh_c, 1/sd_b_sdbh_c^2)
  }
  
  ## Hyperpriors:
  
  ## Observational model:
  mu_a_pdet ~ dnorm(0, 0.01)
  sd_a_pdet ~ dunif(0, 10)
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
  sd_a_pocc ~ dunif(0, 10)
  u_sd_year ~ dunif(0, 10)
  u_sd_site ~ dunif(0, 10)
  mu_b_thinned ~ dnorm(0, 0.01)
  sd_b_thinned ~ dunif(0, 10)
  mu_b_control ~ dnorm(0, 0.01)
  sd_b_control ~ dunif(0, 10)
  mu_b_sdbh ~ dnorm(0, 0.01)
  sd_b_sdbh ~ dunif(0, 10)
  mu_b_sdbh_t ~ dnorm(0, 0.01)
  sd_b_sdbh_t ~ dunif(0, 10)
  mu_b_sdbh_c ~ dnorm(0, 0.01) 
  sd_b_sdbh_c ~ dunif(0, 10)
  
  ## Model validation:
  
  ## Predictions:

  
}


