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
    for(j in 1:3){
      ## Estimate the cummultative probabilities for j=1 to j=J-1:
      logit(cp[i,j]) <- #alpha[j,block[i],year[i]] + 
                        exp_effect[j,treat[i],exp[i]]
    }
    ## Define the category-specific probabilities by substraction of 
    ## cummultative probabilities:
    p[i,1] <- cp[i,1]
    p[i,2] <- cp[i,2] - cp[i,1]
    p[i,3] <- cp[i,3] - cp[i,2]
    p[i,4] <- 1 - cp[i,3]
  }
  
  ## 2. Priors: ----------------------------------------------------------------
  
  # for(k in 1:max(block)){
  #   for(m in 1:max(year)){
  #     alpha[1,k,m]~ dnorm(0, 0.001)T(, alpha[2,k,m]) ## Priors must not overlap!
  #     alpha[2,k,m]~ dnorm(0, 0.001)T(, alpha[3,k,m])
  #     alpha[3,k,m]~ dnorm(0, 0.001)
  # }}
  for(n in 1:max(treat)){
    for(o in 1:max(exp)){
      exp_effect[1,n,o] ~ dnorm(0, 0.001)T(, exp_effect[2,n,o])
      exp_effect[2,n,o] ~ dnorm(0, 0.001)T(, exp_effect[3,n,o])
      exp_effect[3,n,o] ~ dnorm(0, 0.001)
  }}

  ## 3. Model validation: ------------------------------------------------------
  
  ## ...
  
  ## 4. Posterior caclualtions: ------------------------------------------------

  for(n in 1:max(treat)){
    for(o in 1:max(exp)){
      for(j in 1:3){
        logit(cp_post[j,n,o]) <- exp_effect[j,n,o]
      }
      p_post[1,n,o] <- cp_post[1,n,o]
      p_post[2,n,o] <- cp_post[2,n,o] - cp_post[1,n,o]
      p_post[3,n,o] <- cp_post[3,n,o] - cp_post[2,n,o]
      p_post[4,n,o] <- 1 - cp_post[3,n,o]
  }}

}

  