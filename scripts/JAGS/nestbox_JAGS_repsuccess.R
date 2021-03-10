## Zero-Inflated Nestbox model for reproductive success in the great tit
##
## First edit: 20210208
## Last edit: 20210208
##
## Author: Julian Klein

model{
  
  ## 1. Likelihood: ------------------------------------------------------------
  
  for(i in 1:nobs){
    repsuc[i] ~ dpois((1 - z[i])*lambda[i])
    sim[i] ~ dpois((1 - z[i])*lambda[i])
    log(lambda[i]) <- exp_effect_l[treat[i],exp[i]] + 
                      y_effect_l[1]*ifelse(year[i] == 2018, 1, 0) +
                      y_effect_l[2]*ifelse(year[i] == 2019, 1, 0)
    z[i] ~ dbern(p[i])
    logit(p[i]) <- exp_effect_p[treat[i],exp[i]] + 
                   y_effect_p[1]*ifelse(year[i] == 2018, 1, 0) +
                   y_effect_p[2]*ifelse(year[i] == 2019, 1, 0)
  }

  ## 2. Priors: ----------------------------------------------------------------
  
  for(n in 1:max(treat)){
    for(o in 1:max(exp)){
      exp_effect_l[n,o] ~ dnorm(0, 0.001)
      exp_effect_p[n,o] ~ dnorm(0, 0.01)
  }}
  for(y in 1:2){
    y_effect_l[y] ~ dnorm(0, 0.001)
    y_effect_p[y] ~ dnorm(0, 0.01)
  }
  
  ## 3. Model validation: ------------------------------------------------------

  ## Bayesian p-value:
  mean_obs <- mean(repsuc[])
  mean_sim <- mean(sim[])
  ## Variation:
  sd_obs <- sd(repsuc[])
  sd_sim <- sd(sim[])
  ## Model fit:
  for(i in 1:nobs){
    sq[i] <- (repsuc[i] - (1 - z[i])*lambda[i])^2
    sq_sim[i] <- (sim[i] - (1 - z[i])*lambda[i])^2
  }
  fit <- sum(sq[])
  fit_sim <- sum(sq_sim[])
  ## Autocorrelation (Moran's I) for one year:
  mean_sq <- mean(sq[seq])
  for(q in 1:length(seq)){
    for(r in 1:length(seq)){
      ## To the nominator in Moran's I formula:
          nom[q,r] <- dm[q,r]*(sq[seq[q]] - mean_sq)*(sq[seq[r]] - mean_sq)
    }
    ## To the denominator in Moran's I formula:
    denom[q] <- (sq[seq[q]] - mean_sq)^2
  }
  ## Moran's I is calculated:
  I <- length(seq)*sum(nom[,])/(sum(dm[,])*sum(denom[]))

  ## 4. Posterior caclualtions: ------------------------------------------------

  for(n in 1:max(treat)){
    for(o in 1:max(exp)){
      log(lambda_post[n,o]) <- exp_effect_l[n,o] ## Prediction as if 2017
      logit(p_post[n,o]) <- exp_effect_p[n,o] ## Prediction as if 2017
  }}

}

  