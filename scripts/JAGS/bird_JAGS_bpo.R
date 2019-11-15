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
    logit(pdet[i]) <- a_pdet[species[i]] + 
                      b_pdet_2018[species[i]]*ifelse(year[i]==2,1,0) + 
                      b_pdet_2019[species[i]]*ifelse(year[i]==3,1,0) +
                      b_observer*observer[i] +
                      b_dpm[ldm[i]]*dpm[i,1] + b2_dpm[ldm[i]]*dpm[i,1]^2 +
                      b_mps*mps[i,1] #+ b2_mps*mps[i,1]^2
  }
  
  ## Ecological process model:
  for(k in 1:nspecies){
    for(y in 1:nyears){
      for(p in 1:nsites){
        occ_true[k,y,p] ~ dbern(pocc[k,y,p])
        logit(pocc[k,y,p]) <- a_pocc[k,treat[p],exp[y,p]] #+ 
                              # e_site[k,p] + e_year[k,y] +
                              # b_pocc_2018[k]*ifelse(y==2,1,0) + 
                              # b_pocc_2019[k]*ifelse(y==3,1,0) + 
                              # b_sdbh[k]*sdbh[y,p] +
                              # b_dec[k]*dec[y,p] +
                              # b_umbr[k]*umbr[y,p] +
                              # b_lm[k]*lm[y,p] +
  }}}
  
  # ## Group effects:
  # for(k in 1:nspecies){
  #   # for(y in 1:nyears){e_year[k,y] ~ dnorm(0, 1/sd_year[k]^2)} 
  #   # for(p in 1:nsites){e_site[k,p] ~ dnorm(0, 1/sd_site[k]^2)}
  # }
  
  ## Priors:
  
  ## Observational model:
  for(k in 1:nspecies){
    a_pdet[k] ~ dnorm(mu_a_pdet, 1/sd_a_pdet^2)
    b_pdet_2018[k] ~ dnorm(mu_b_pdet_2018, 1/sd_b_pdet_2018^2)
    b_pdet_2019[k] ~ dnorm(mu_b_pdet_2019, 1/sd_b_pdet_2019^2)
  }
  b_observer ~ dnorm(0, 0.01)
  for(l in 1:max(ldm)){    
    b_dpm[l] ~ dnorm(0, 0.01)
    b2_dpm[l] ~ dnorm(0, 0.01)
  }
  b_mps ~ dnorm(0, 0.01)
  # b2_mps ~ dnorm(0, 0.01)
  
  ## Ecological process model:
  for(k in 1:nspecies){
    for(m in 1:max(treat)){
      for(n in 1:max(exp)){
        a_pocc[k,m,n] ~ dnorm(mu_a_pocc[m,n], 1/sd_a_pocc[m,n]^2)
    }}
    # sd_year[k] ~ dunif(0, u_sd_year)
    # sd_site[k] ~ dunif(0, u_sd_site)
    # b_pocc_2018[k] ~ dnorm(mu_b_pocc_2018, 1/sd_b_pocc_2018^2)
    # b_pocc_2019[k] ~ dnorm(mu_b_pocc_2019, 1/sd_b_pocc_2019^2)
    # b_sdbh[k] ~ dnorm(mu_b_sdbh, 1/sd_b_sdbh^2)
    # b_dec[k] ~ dnorm(mu_b_dec, 1/sd_b_dec^2)
    # b_umbr[k] ~ dnorm(mu_b_umbr, 1/sd_b_umbr^2)
    # b_lm[k] ~ dnorm(mu_b_lm, 1/sd_b_lm^2)
  }
  
  ## Hyperpriors:
  
  ## Observational model:
  mu_a_pdet ~ dnorm(0, 0.01)
  sd_a_pdet ~ dunif(0, 5)
  mu_b_pdet_2018 ~ dnorm(0, 0.1)
  sd_b_pdet_2018 ~ dunif(0, 5)
  mu_b_pdet_2019 ~ dnorm(0, 0.1)
  sd_b_pdet_2019 ~ dunif(0, 5)
  
  ## Ecological process model:
  for(m in 1:max(treat)){
    for(n in 1:max(exp)){
      mu_a_pocc[m,n] ~ dnorm(0, 0.01)
      sd_a_pocc[m,n]~ dunif(0, 5)
  }}
  # u_sd_year ~ dunif(0, 10)
  # u_sd_site ~ dunif(0, 10)
  # mu_b_pocc_2018 ~ dnorm(0, 0.1)
  # sd_b_pocc_2018 ~ dunif(0, 5)
  # mu_b_pocc_2019 ~ dnorm(0, 0.1)
  # sd_b_pocc_2019 ~ dunif(0, 5)
  # mu_b_sdbh ~ dnorm(0, 0.1)
  # sd_b_sdbh ~ dunif(0, 5)
  # mu_b_dec ~ dnorm(0, 0.1)
  # sd_b_dec ~ dunif(0, 5)
  # mu_b_umbr ~ dnorm(0, 0.1)
  # sd_b_umbr ~ dunif(0, 5)
  # mu_b_lm ~ dnorm(0, 0.1)
  # sd_b_lm ~ dunif(0, 5)

  ## Model validation:
  
  ## Posteriors:
  
  ## BACI indicators:
  
  # Posterior distribution of the three measures of impact
  CI.div=abs(aft.imp-aft.ctrl)-abs(bef.imp-bef.ctrl)
  CI.cont=abs(aft.imp-bef.imp)-abs(aft.ctrl-bef.ctrl)
  BACI=(aft.imp-bef.imp)-(aft.ctrl-bef.ctrl) 

  
}


