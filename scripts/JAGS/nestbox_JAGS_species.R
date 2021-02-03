## Nestbox model for nestbox occupancy by species
##
## First edit: 20200129
## Last edit: 20200203
##
## Author: Julian Klein

## Page 42 in: http://people.bu.edu/dietze/Bayes2018/Lesson21_GLM.pdf
## All: http://people.bu.edu/dietze/Bayes2018/GE509.htm

model{
  
  ## 1. Lijelihood: ------------------------------------------------------------
  
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
  for(j in 1:3){
    for(i in 1:nobs){
      sq[i,j] <- (occ[i,j] - p[i,j])^2
      sq_sim[i,j] <- (sim[i,j] - p[i,j])^2
    }
    fit[j] <- sum(sq[,j])
    fit_sim[j] <- sum(sq_sim[,j])
  }
  ## Autocorrelation (Moran's I) for one year:
  mean_sq <- mean(sq[seq,])
  for(q in 1:length(seq)){
    for(j in 1:3){
      for(r in 1:length(seq)){
        for(s in 1:3){
          ## To the nominator in Moran's I formula:
          nom[q,r,j,s] <- dm[q,r]*
                          (sq[seq[q],j] - mean_sq)*
                          (sq[seq[r],s] - mean_sq)
      }}
      ## To the denominator in Moran's I formula:
      denom[q,j] <- (sq[seq[q],j] - mean_sq)^2
  }}
  ## Moran's I is calculated:
  I <- length(seq)*sum(nom[,,,])/(sum(dm[,])*sum(denom[,]))
  
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

  