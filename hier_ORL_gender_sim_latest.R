hier_ORL_gender_sim <- function(payoff, nsubs, ntrials,
                                mu_a_rew, alpha_a_rew,
                                mu_a_pun, alpha_a_pun,
                                mu_K,
                                mu_omega_f, alpha_omega_f,
                                mu_omega_p,
                                sigma_a_rew, sigma_a_pun, sigma_K,
                                sigma_omega_f, sigma_omega_p) {
  
  gender <- c(rep(0, 98), rep(1, 41))  # 98 m, 41 f
  
  x <- array(NA, c(nsubs, ntrials[1]))    # Choices
  X <- array(NA, c(nsubs, ntrials[1]))    # Outcomes
  Ev <- array(NA, c(nsubs, ntrials[1], 4)) # Expected values
  
  for (s in 1:nsubs) {
    
    # parameters, including gender effects for hypothesized parameters
    a_rew <- rtruncnorm(1, 0, 1, mu_a_rew + alpha_a_rew*(gender[s]-0.5), sigma_a_rew)
    a_pun <- rtruncnorm(1, 0, 1, mu_a_pun + alpha_a_pun*(gender[s]-0.5), sigma_a_pun)
    K <- rtruncnorm(1, 0, 5, mu_K, sigma_K)
    omega_f <- rnorm(1, mu_omega_f + alpha_omega_f*(gender[s]-0.5), sigma_omega_f)
    omega_p <- rnorm(1, mu_omega_p, sigma_omega_p)
    theta <- 1  # Fixed
    
    # Initialize subject-specific arrays for the current subject
    Ef <- matrix(NA, ntrials[1], 4)
    PS <- matrix(NA, ntrials[1], 4)
    
    # Initialize values for the first trial (t=1) of the current subject
    Ev[s,1,] <- rep(0,4)
    Ef[1,] <- rep(0,4)
    PS[1,] <- rep(1,4) # can be set to 1 or 0, must match JAGS
    
    x[s,1] <- sample(1:4, 1, prob=rep(0.25,4))
    X[s,1] <- payoff[1, x[s,1]]
    
    for(t in 2:ntrials[s]) {
      
      signX <- ifelse(X[s,t-1] < 0, -1, 1)
      
      # Initialize for each trial
      exp_p <- numeric(4) 
      V <- numeric(4)
      
      # Deck Updating Loop
      for(d in 1:4) {
        # Value Update
        Ev_update <- if(X[s,t-1] >= 0) {
          Ev[s,t-1,d] + a_rew*(X[s,t-1] - Ev[s,t-1,d])
        } else {
          Ev[s,t-1,d] + a_pun*(X[s,t-1] - Ev[s,t-1,d])
        }
        Ev[s,t,d] <- ifelse(d == x[s,t-1], Ev_update, Ev[s,t-1,d])
        
        # Frequency Update
        Ef_cho <- if(X[s,t-1] >= 0) {
          Ef[t-1,d] + a_rew*(signX - Ef[t-1,d])
        } else {
          Ef[t-1,d] + a_pun*(signX - Ef[t-1,d])
        }
        Ef_not <- if(X[s,t-1] >= 0) {
          Ef[t-1,d] + a_pun*(-signX/3 - Ef[t-1,d])
        } else {
          Ef[t-1,d] + a_rew*(-signX/3 - Ef[t-1,d])
        }
        Ef[t,d] <- ifelse(d == x[s,t-1], Ef_cho, Ef_not)
        
        # Perseverance Update
        PS[t,d] <- ifelse(x[s,t-1] == d, 1/(1+K), PS[t-1,d]/(1+K))
        
        # Calculate Exponentiated Values
        V[d] <- Ev[s,t,d] + Ef[t,d]*omega_f + PS[t,d]*omega_p # V should be indexed by d
        exp_p[d] <- exp(theta * V[d]) # exp_p should be indexed by d
      }
      
      # Calculate Choice Probabilities
      p <- exp_p / sum(exp_p)
      stopifnot(length(p) == 4, all(is.finite(p)))  # Safety check
      
      x[s,t] <- sample(1:4, 1, prob = p)
      X[s,t] <- payoff[t, x[s,t]]
    }
  }
  
  # 4. Return Results
  return(list(x = x, X = X, Ev = Ev, gender = gender))
}