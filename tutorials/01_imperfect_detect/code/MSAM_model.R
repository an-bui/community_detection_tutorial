model{
  
  #Example MSAM model for tutorial with simulated data
  
  
  for(s in 1:n.species){ #species
    for(t in 1:n.transects){ #transects
      for(y in n.start[t]:n.end[t]){ #years: setting these start and end years to depend on transect allows them to vary by transect
        
        #BIOLOGICAL MODEL 
        
        #expected number of individuals of species s
        #in site t in time y is dependent on 
        #a rate parameter, lambda, for a poisson distribution
        N[s,t,y] ~ dpois(lambda[s,t,y])
        
        #this rate parameter is dependent on a regression
        #with a species intercept,
        #and species within site, and species within year 
        #random effects (post-sweeping code is below for this )
        log(lambda[s,t,y]) <- b0.species[s] + 
          eps.site[Site.ID[t], s] +
          eps.year[Year.ID[y], s]
        
        for(r in 1:n.rep[t,y]){ #for the number of surveys on each transect in each year
          
          # OBSERVATION MODEL
          
          #detection probability for species s in site t in
          #year y in replicate r
          logit(p[s,t,y,r]) <- a0[s] + #species-level intercept
            #a covariate that is related to the survey (e.g., weather, survey length)
            #dependent on site, year, and replicate, but not species
            a1*survey.covariate[t,y,r] #+
          #could also add species covariates (e.g., body size, color, call frequency)
            #a2*species.covariate[s]
          
          #abundance is binomial based on detection probability
          #and total true abundance at the site
          count[s,t,y,r] ~ dbin(p[s,t,y,r], N[s,t,y])
          
          #create replicate data based on model estimation to
          #look at goodness-of-fit (regress y.rep ~ y)
          count.rep[s,t,y,r] ~ dbin(p[s,t,y,r], N[s,t,y])
          
        }
      }
    }
    
    #SPECIES-LEVEL PRIORS:
    #Detection intercept and slopes for each species
    #are centered on community-level hyperpriors
    a0[s] ~ dnorm(mu.a0,tau.a0)
    
    #"baseline" detection at vis = 0 and size = 0 
    #on standardized scale
    p0[s] <- ilogit(a0[s])
    
    #species-level intercept - 
    #currently non-identifiable:
    b0.species[s] ~ dnorm(mu.b0species, tau.b0species)
    
    #compute the identifiable species-level intercept:
    #TRACK THIS ONE for convergence
    b0.star[s] <- b0.species[s] + ave.eps.site[s] + ave.eps.year[s]
    
  }
  
  #SITE W/IN SPECIES RANDOM EFFECTS
  #sites nested within species (sites sum to zero w/in each species)
  for(s in 1:n.species){
    for(t in 1:n.sites){
      #non-identifiable random effect
      eps.site[t,s] ~ dnorm(0, tau.eps.site)
      #identifiable site random effect (monitor this)
      eps.site.star[t,s] <- eps.site[t,s] - ave.eps.site[s]
    }
    #mean site level random effects within each species
    ave.eps.site[s] <- mean(eps.site[,s])
  }
  
  #YEARS W/IN SPECIES RANDOM EFFECTS
  #sites nested within species (sites sum to zero w/in each species)
  for(s in 1:n.species){
    for(y in 1:n.years){
      #non-identifiable random effect
      eps.year[y,s] ~ dnorm(0, tau.eps.year)
      #identifiable year random effect (monitor this)
      eps.year.star[y,s] <- eps.year[y,s] - ave.eps.year[s]
    }
    #mean year level random effects within each species
    ave.eps.year[s] <- mean(eps.year[,s])
  }
  
  #COMMUNITY HYPER PRIORS
  #initial occupancy
  #Detection intercept
  mu.a0 ~ dnorm(0, 0.001)
  tau.a0 <- pow(sig.a0, -2)
  sig.a0 ~ dunif(0, 50)
  
  #species-level abundance
  mu.b0species ~ dnorm(0, 0.001)
  tau.b0species <- pow(sig.b0species, -2)
  sig.b0species ~ dunif(0,50)
  
  #site and year variances
  sig.eps.site ~ dunif(0, 10)
  tau.eps.site <- pow(sig.eps.site, -2)
  sig.eps.year ~ dunif(0, 10)
  tau.eps.year <- pow(sig.eps.year, -2)
  
  #detection covariate effect priors
  a1 ~ dnorm(0, 0.001)
  #a2 ~ dnorm(0, 0.001)
  
  #IN CASE OF missing data 
  #If detection covariate data are missing, this will
  #impute based on mean and variance of that variable
  for(t in 1:n.transects){
    for(y in n.start[t]:n.end[t]){
      for(r in 1:n.rep[t,y]){
        
        survey.covariate[t,y,r] ~ dnorm(mu.surveycov, tau.surveycov)
      }
    }
  }
  
  #PRIORS FOR IMPUTING MISSING DATA
  #Priors for mean and tau of missing covariates in the model
  mu.surveycov ~ dunif(-10, 10)
  sig.surveycov ~ dunif(0, 20)
  tau.surveycov <- pow(sig.surveycov, -2)
  
  
  #DERIVED QUANTIIES
  
  
  #Bray-Curtis dissimilarity
  for(t in 1:n.transects){
    for(y in (n.start[t]+1):n.end[t]){
      for(s in 1:n.species){
        # num individuals in both time periods per species
        a[s,t,y] <- min(N[s,t,y-1], N[s,t,y])
        # num individuals only in first time point
        b[s,t,y] <- N[s,t,y-1] - a[s,t,y]
        # num individuals only in second time point
        c[s,t,y] <- N[s,t,y] - a[s,t,y]
      }
      #for all years 2 onward:
      #total number of shared individuals across time periods
      A[t,y] <- sum(a[,t,y])
      #total number of individuals in only first time period
      B[t,y] <- sum(b[,t,y])
      #total number of individuals in only second time period
      C[t,y] <- sum(c[,t,y])
      
      #total bray-curtis (B+C)/(2A+B+C)
      num[t,y] <- B[t,y] + C[t,y]
      denom1[t,y] <- 2*A[t,y]+B[t,y]+C[t,y]
      #if all values are zero - this just keeps the eqn. from
      #dividing by zero
      denom[t,y] <- ifelse(denom1[t,y]==0,1, denom1[t,y])
      
      #Calculate Bray-Curtis dissimiarlity
      bray[t,y] <- num[t,y]/denom[t,y]
      
    }
  }
}
  
  
  
  
  
  
  
  
  
  
  
  
  
  
