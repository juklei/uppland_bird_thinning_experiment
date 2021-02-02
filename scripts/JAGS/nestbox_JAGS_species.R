## Nestbox model for nestbox occupancy by species
##
## First edit: 20200129
## Last edit: 20200201
##
## Author: Julian Klein

## Page 42 in: http://people.bu.edu/dietze/Bayes2018/Lesson21_GLM.pdf
## All: http://people.bu.edu/dietze/Bayes2018/GE509.htm

model{
  
  ## 1. Likelihood: ------------------------------------------------------------
  
  for(i in 1:nobs){
    occ[i,] ~ dmulti(p[i,], 1)
    sim[i,1:4] ~ dmulti(p[i,], 1) ## Define number of output categories!
    for(j in 1:3){
      ## Estimate the cummultative probabilities for j=1 to j=J-1:
      logit(cp[i,j]) <- alpha[j,year[i]] + exp_effect[j,treat[i],exp[i]]
    }
    ## Define the category-specific probabilities by substraction of 
    ## cummultative probabilities:
    p[i,1] <- cp[i,1]
    p[i,2] <- cp[i,2] - cp[i,1]
    p[i,3] <- cp[i,3] - cp[i,2]
    p[i,4] <- 1 - cp[i,3]
  }
  
  ## 2. Priors: ----------------------------------------------------------------
  
  for(m in 1:max(year)){
    alpha[1,m]~ dnorm(0, 0.001)T(, alpha[2,m]) ## Priors must not overlap!
    alpha[2,m]~ dnorm(0, 0.001)T(, alpha[3,m])
    alpha[3,m]~ dnorm(0, 0.001)
  }
  for(n in 1:max(treat)){
    for(o in 1:max(exp)){
      exp_effect[1,n,o] ~ dnorm(0, 0.001)T(, exp_effect[2,n,o])
      exp_effect[2,n,o] ~ dnorm(0, 0.001)T(, exp_effect[3,n,o])
      exp_effect[3,n,o] ~ dnorm(0, 0.001)
  }}
  
  ## 3. Model validation: ------------------------------------------------------
  
  ## Bayesian p-value:
  mean_obs <- mean(occ[,1:3])
  mean_sim <- mean(sim[,1:3])
  ## Coefficient of variation:
  cv_obs <- sd(occ[,1:3])/mean_obs
  cv_sim <- sd(sim[,1:3])/mean_sim
  ## Model fit:
  for(i in 1:nobs){
    sq[i,1:3] <- (occ[i,1:3] - p[i,1:3])^2
    sq_sim[i,1:3] <- (sim[i,1:3] - p[i,1:3])^2
  }
  fit <- sum(sq[,1:3])
  fit_sim <- sum(sq_sim[,1:3])
  ## Autocorrelation (Moran's I) for year 2017 for now:
  for(k in 1:4){
    for(q in 1:length(seq_2017)){
      sq_2017[q,k] <- (occ[seq_2017[q],k] - p[seq_2017[q],k])^2 ## Residuals
      for(r in 1:length(seq_2017)){
        ## To the nominator in Moran's I formula:
        nom[q,r,k] <- dm_2017[q,r]*(sq_2017[q,k] - mean(sq_2017[,k]))*(sq_2017[r,k] - mean(sq_2017[,k]))
      }
      ## To the denominator in Moran's I formula:
      denom[q,k] <- (sq_2017[q,k] - mean(sq_2017[,k]))^2
    }
    ## Moran's I is calculated:
    I[k] <- length(seq_2017)*sum(nom[,,k])/(sum(dm_2017[,])*sum(denom[,k]))
  }
  
  ## 4. Posterior caclualtions: ------------------------------------------------

  for(n in 1:max(treat)){
    for(o in 1:max(exp)){
      for(j in 1:3){
        logit(cp_post[j,n,o]) <- mean(alpha[j,]) + exp_effect[j,n,o]
      }
      p_post[1,n,o] <- cp_post[1,n,o]
      p_post[2,n,o] <- cp_post[2,n,o] - cp_post[1,n,o]
      p_post[3,n,o] <- cp_post[3,n,o] - cp_post[2,n,o]
      p_post[4,n,o] <- 1 - cp_post[3,n,o]
  }}

}

  