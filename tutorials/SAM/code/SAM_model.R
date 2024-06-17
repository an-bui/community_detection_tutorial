model{
  
  for(i in 1:n.data){
    
    #-------------------------------------## 
    # Likelihood ###
    #-------------------------------------##
    
    #d is dissimilarity and is proportional, so beta distribution works here
    d[i] ~ dbeta(alpha[i], beta[i])
    
    #var.process is scalar but could be made dependent on site/other variables
    #phi incorporates mu (mean estimate), var.estimate (which is from
    #"data" on standard deviation (squared) from the original detection 
    #correction model) 
    #and var.process is something we're trying to estimate,
    #basically, the rest of the variation not accounted for
    phi[i] <- (((1-mu[i])*mu[i])/(var.estimate[i] + var.process))-1
    
    #alpha and beta are based on mu and phi values
    #sometimes these values send alpha and beta outside
    #the domain, so we have extra code below to get them to
    #stay where they belong
    alphaX[i] <- mu[i] * phi[i]
    betaX[i] <- (1 - mu[i]) * phi[i]
    
    #here is where we get alpha and beta to stay in their
    #domain
    alpha[i] <- max(0.01, alphaX[i])
    beta[i] <- max(0.01, betaX[i])
    
    #to get a good estimate of a prior for var.process, we
    #track the difference between these two values for each
    #data point
    diff[i] <- (1-mu[i])*mu[i] - var.estimate[i]
    
    #Regression of mu, which is dependent on a hierrarchically-centered
    #site random effect and the weighted antecedent effects of two covariates,
    #plant biomass and temperature
    logit(mu[i]) <- b0.site[Site.ID[i]] +
      b[1]*AntPlant[i] +
      b[2]*AntTemp[i]
    
    #-------------------------------------## 
    # SAM summing ###
    #-------------------------------------##
    
    #summing the antecedent values
    AntPlant[i] <- sum(PlantTemp[i,]) #summing across the total number of antecedent years
    AntTemp[i] <- sum(TempTemp[i,]) #summing across the total num of antecedent months
    
    #Generating each year's weight to sum above
    for(t in 1:n.plantlag){ #number of time steps we're going back in the past
      PlantTemp[i,t] <- Plant[i,t]*wA[t] 
      
      #imputing missing data
      Plant[i,t] ~ dnorm(mu.plant, tau.plant)
    }
    
    #generating each month's weight to sum above
    for(t in 1:n.templag){ #number of time steps we're going back in the past
      TempTemp[i,t] <- Temp[i,t]*wB[t] 
      
      #missing data
      Temp[i,t] ~ dnorm(mu.temp, tau.temp)
    }
    
    #-------------------------------------## 
    # Goodness of fit parameters ###
    #-------------------------------------##
    # 
    #replicated data
    d.rep[i] ~ dbeta(alpha[i], beta[i])
    # 
    #residuals
    resid[i] <- d[i] - mu[i]
    
  }
  
  #-------------------------------------## 
  # Priors ###
  #-------------------------------------##
  
  # ANTECEDENT CLIMATE PRIORS
  #Sum of the weights for lag
  sumA <- sum(deltaA[]) #all the plant weights
  
  #Employing "delta trick" to give vector of weights dirichlet priors
  #this is doing the dirichlet in two steps 
  #see Ogle et al. 2015 SAM model paper in Ecology Letters
  for(t in 1:n.plantlag){ #for the total number of lags
    #the weights for kelp - getting the weights to sum to 1
    wA[t] <- deltaA[t]/sumA
    #and follow a relatively uninformative gamma prior
    deltaA[t] ~ dgamma(1,1)
    
    #to look at how weights accumulate through time
    cumm.plantwt[t] <- sum(wA[1:t])
  }
  
  #Sum of the weights for temp lag
  sumB <- sum(deltaB[]) #all the temp weights
  
  #Employing "delta trick" to give vector of weights dirichlet priors
  #this is doing the dirichlet in two steps 
  #see Ogle et al. 2015 SAM model paper in Ecology Letters
  for(t in 1:n.templag){ #for the total number of lags
    #the weights for kelp - getting the weights to sum to 1
    wB[t] <- deltaB[t]/sumB
    #and follow a relatively uninformative gamma prior
    deltaB[t] ~ dgamma(1,1)
    
    #to look at cummulative weigths through time
    cumm.tempwt[t] <- sum(wB[1:t])
  }
  
  #BETA PRIORS
  #HIERARCHICAL STRUCTURE PRIORS
  #hierarchical centering of sites on b0
  for(s in 1:n.sites){
    b0.site[s] ~ dnorm(b0, tau.site)
  }
  
  #overall intercept gets relatively uninformative prior
  b0 ~ dnorm(0, 1E-2)
  
  #for low # of levels, from Gelman paper - define sigma
  # as uniform and then precision in relation to this sigma
  sig.site ~ dunif(0, 10)
  tau.site <- 1/pow(sig.site,2)
  
  #covariate effects - again get relatively uninformative priors
  for(i in 1:2){
    b[i] ~ dnorm(0, 1E-2)
  }
  
  #PRior for process error
  var.process ~ dunif(0, min(diff[]))
  
  
  #MISSING DATA PRIORS
  mu.plant ~ dunif(-10, 10)
  sig.plant ~ dunif(0, 20)
  tau.plant <- pow(sig.plant, -2)
  mu.temp ~ dunif(-10, 10)
  sig.temp ~ dunif(0, 20)
  tau.temp <- pow(sig.temp, -2)
  
}