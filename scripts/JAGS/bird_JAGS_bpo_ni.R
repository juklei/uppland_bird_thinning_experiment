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
    observed[i] ~ dbern(occ_true[pys[i]]*pdet[i])
    logit(pdet[i]) <- a_pdet[sH1[i]] + 
                      b_dpm[ldm[i]]*dpm[i] + b2_dpm[ldm[i]]*dpm[i]^2 +
                      b_mps*mps[i] + b2_mps*mps[i]^2
  }
  
  ## Ecological process model:
  for(j in 1:npys){
    occ_true[j] ~ dbern(pocc[j])
    logit(pocc[j]) <- a_pocc[sH2[j]] + e_year[ys[j]] + e_site[ps[j]] + 
                      b_thinned[sH2[j]]*thinned[j] +
                      b_control[sH2[j]]*control[j] +
                      b_sdbh[sH2[j]]*sdbh[j] +
                      b_sdbh_t[sH2[j]]*thinned[j]*sdbh[j] +
                      b_sdbh_c[sH2[j]]*control[j]*sdbh[j]
  }

  ## Grouping effects:
  for(k in 1:nps){e_site[k] ~ dnorm(0, 1/sd_site[sH3p[k]]^2)}
  for(l in 1:nys){e_year[l] ~ dnorm(0, 1/sd_year[sH3y[l]]^2)}
  
  ## Priors:
  
  ## Observational model:
  for(m in 1:ns){a_pdet[m] ~ dnorm(mu_a_pdet, 1/sd_a_pdet^2)}
  for(n in 1:2){    
    b_dpm[n] ~ dnorm(0, 0.01) 
    b2_dpm[n] ~ dnorm(0, 0.01) 
  }
  b_mps ~ dnorm(0, 0.01) 
  b2_mps ~ dnorm(0, 0.01) 

  ## Ecological process model:
  for(m in 1:ns){
    a_pocc[m] ~ dnorm(mu_a_pocc, 1/sd_a_pocc^2)
    sd_year[m] ~ dunif(0, 5)
    sd_site[m] ~ dunif(0, 5)
    b_thinned[m] ~ dnorm(mu_b_thinned, 1/sd_b_thinned^2)
    b_control[m] ~ dnorm(mu_b_control, 1/sd_b_control^2)
    b_sdbh[m] ~ dnorm(mu_b_sdbh, 1/sd_b_sdbh^2)
    b_sdbh_t[m] ~ dnorm(mu_b_sdbh_t, 1/sd_b_sdbh_t^2)
    b_sdbh_c[m]~ dnorm(mu_b_sdbh_c, 1/sd_b_sdbh_c^2)
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


