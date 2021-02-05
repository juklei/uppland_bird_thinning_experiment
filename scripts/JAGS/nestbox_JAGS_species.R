## Nestbox model for nestbox occupancy by species
##
## First edit: 20200129
## Last edit: 20200203
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
      logit(cp[i,j]) <- exp_effect[j,treat[i],exp[i]] + 
                        y_effect[j,1]*ifelse(year[i] == 2018, 1, 0) + 
                        y_effect[j,2]*ifelse(year[i] == 2019, 1, 0)
    }
    ## Define the category-specific probabilities by substraction of 
    ## cummultative probabilities:
    p[i,1] <- cp[i,1]
    p[i,2] <- cp[i,2] - cp[i,1]
    p[i,3] <- cp[i,3] - cp[i,2]
    p[i,4] <- 1 - cp[i,3]
  }
  
  ## 2. Priors: ----------------------------------------------------------------
  
  for(n in 1:max(treat)){
    for(o in 1:max(exp)){
      exp_effect[1,n,o] ~ dnorm(0, 0.001)T(, exp_effect[2,n,o])
      exp_effect[2,n,o] ~ dnorm(0, 0.001)T(, exp_effect[3,n,o])
      exp_effect[3,n,o] ~ dnorm(0, 0.001)
    }}
  for(y in 1:2){
    y_effect[1,y] ~ dnorm(0, 0.001)T(, y_effect[2,y])
    y_effect[2,y] ~ dnorm(0, 0.001)T(, y_effect[3,y])
    y_effect[3,y] ~ dnorm(0, 0.001)
  }
  
  ## 3. Model validation: ------------------------------------------------------

  ## Bayesian p-value:
  mean_obs <- mean(occ[,1:3])
  mean_sim <- mean(sim[,1:3])
  ## Coefficient of variation:
  sd_obs <- sd(occ[,1:3])
  sd_sim <- sd(sim[,1:3])
  ## Model fit:
  for(j in 1:3){
    for(i in 1:nobs){
      sq[i,j] <- (occ[i,j] - p[i,j])^2
      sq_sim[i,j] <- (sim[i,j] - p[i,j])^2
  }}
  fit <- sum(sq[,])
  fit_sim <- sum(sq_sim[,])
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
        logit(cp_post[j,n,o]) <- exp_effect[j,n,o] ## Prediction as if 2017 
      }
      p_post[1,n,o] <- cp_post[1,n,o]
      p_post[2,n,o] <- cp_post[2,n,o] - cp_post[1,n,o]
      p_post[3,n,o] <- cp_post[3,n,o] - cp_post[2,n,o]
      p_post[4,n,o] <- 1 - cp_post[3,n,o]
  }}

}

  