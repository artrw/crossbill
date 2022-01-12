I <- 500 # number of individuals
J <- 100 # number of sites
S <- 50  # number of species 
iterations <- 50 # number of periods in each simulation
split <- 20 # period in which final species becomes present at single site
turnover <- 50 # number of birders replaced each period
beta_params <- c(2, 8) # parameters for beta distribution of phi
run_count <- 50 # number of runs of the simulation for each phi distribution