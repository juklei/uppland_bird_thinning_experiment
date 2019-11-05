## bird bpo model with a bernoulli process for occurrence
##
## First edit: 20191030
## Last edit: 20191031
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
      pocc[k,y] ~ dunif(0.001, 0.999)
      for(i in 1:nsites){
        occ_true[k,y,i] ~ dbern(pocc[k,y])
        # logit(p_occ[k,y,i]) <- a_pocc[k,experiment[y,i]] + year_effect[k,y] #+ plot effect instead of intercept by expermient 
                               # b_sdbh[k]*sdbh[experiment[y,i],i]
  }}}
  
  # ## Year group effect:
  # for(k in 1:nspecies){
  #   for(y in 1:nyears){
  #     year_effect[k,y] ~ dnorm(0, 1/sd_year[k]^2)
  # }}
  
  ## Priors:
  
  ## Observational model:
  for(k in 1:nspecies){a_pdet[k] ~ dnorm(mu_a_pdet, 1/sd_a_pdet^2)}
  for(l in 1:max(ldm)){    
    b_dpm[l] ~ dnorm(0, 0.01) #dnorm(mu_b_dpm, 1/sd_b_dpm^2)
    b2_dpm[l] ~ dnorm(0, 0.01) #dnorm(mu_b2_dpm, 1/sd_b2_dpm^2)
    b_mps[l] ~ dnorm(0, 0.01) #dnorm(mu_b_mps, 1/sd_b_mps^2)
    b2_mps[l] ~ dnorm(0, 0.01) #dnorm(mu_b2_mps, 1/sd_b2_mps^2)
  }
  
  # ## Ecological process model:
  # for(k in 1:nspecies){
  #   for(y in 1:nyears) {
  #   a_pocc[k,experiment[y]] ~ dnorm(mu_a_pocc, 1/sd_a_pocc^2)
  #   }
  #   b_sdbh[k] ~ dnorm(mu_b_sdbh, 1/sd_b_sdbh^2)
  #   sd_year[k] ~ dgamma(a_sd_year, b_sd_year)
  # }
  
  ## Hyperpriors:
  mu_a_pdet ~ dnorm(0, 0.01)
  sd_a_pdet ~ dunif(0, 10)
  # mu_b_dpm ~ dnorm(0, 0.01)
  # sd_b_dpm ~ dgamma(0.001, 0.001)
  # mu_b2_dpm ~ dnorm(0, 0.01)
  # sd_b2_dpm ~ dgamma(0.001, 0.001)
  # mu_b_mps ~ dnorm(0, 0.01)
  # sd_b_mps ~ dgamma(0.001, 0.001)
  # mu_b2_mps ~ dnorm(0, 0.01)
  # sd_b2_mps ~ dgamma(0.001, 0.001)
  # mu_a_pocc ~ dnorm(0, 0.001)
  # sd_a_pocc ~ dgamma(0.001, 0.001)
  # mu_b_sdbh ~ dnorm(0, 0.001)
  # sd_b_sdbh ~ dgamma(0.001, 0.001)
  # a_sd_year ~ dunif(0, 1)
  # b_sd_year ~ dunif(0, 1)
  
  # Model validation:
  # 
  # ## Bayesian p-value:
  # mean_nseen <- mean(nseen[,,])
  # mean_nseen_sim <- mean(nseen_sim[,,])
  # p_mean <- step(mean_nseen_sim - mean_nseen)
  # 
  # ## Coefficient of variation:
  # cv_nseen <- sd(nseen[,,])/mean_nseen
  # cv_nseen_sim <- sd(nseen_sim[,,])/mean_nseen_sim
  # p_cv <- step(cv_nseen - cv_nseen_sim)
  # 
  # ## Model fit:
  # for(k in 1:nspecies){
  #   for(y in 1:nyears){
  #     for(i in 1:nsites){
  #       sq[i,k,y] <- (nseen[i,k,y] -
  #                     occ_true[i,k,y]*p_det[k,y]*nvisits[i,k,y])^2
  #       sq_sim[i,k,y] <- (nseen_sim[i,k,y] -
  #                         occ_true[i,k,y]*p_det[k,y]*nvisits[i,k,y])^2
  #     }
  #   }
  # }
  # 
  # fit <- sum(sq[,,])
  # fit_sim <- sum(sq_sim[,,])
  # p_fit <- step(fit_sim - fit)
  
  # ## Predictions:
  # 
  # ## Stand dbh:
  # for(m in 1:length(stand_dbh_pred)){
  #   for(k in 1:nspecies){
  #     occ_true_stand_dbh[m,k] ~ dbern(p_occ_stand_dbh[m,k])
  #     logit(p_occ_stand_dbh[m,k]) <- alpha[k] + 
  #                                    beta_stand_dbh[k]*stand_dbh_pred[m]
  #   }
  #   r_stand_dbh[m] <- sum(p_occ_stand_dbh[m,])
  # }
  # 
  # 
  # ## Plot level richness:
  # for(i in 1:nsites){
  #   for(y in 1:nyears){
  #     r_year[i,y] <- sum(occ_true[i,,y])
  #   }
  #   r_plot[i] <- sum(r_year[i,])/nyears
  # }
  
}


