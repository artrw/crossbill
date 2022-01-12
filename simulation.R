library(tidyverse)

source("simulation_parameters.R")

distance_matrix <- function(home_locations, site_locations){
  # Input: df of home locations with id column "indv" and coordinate columns
  # "x.i" and "y.i"; df of site locations with id column "site" and coordinate
  # columns "x.j" and "y.j"
  # output: distance matrix
  output <- cross_df(list(indv = home_locations$indv, 
                          site = site_locations$site)) %>%
    left_join(home_locations[,c("indv", "x.i", "y.i")], by = "indv") %>%
    left_join(site_locations[,c("site", "x.j", "y.j")], by = "site") %>%
    {.$distance <- apply(.[,c("x.i", "y.i", "x.j", "y.j")], 1,
                         function(x) dist(matrix(x, 
                                                 nrow = 2, 
                                                 byrow = TRUE))); .} %>%
    select(indv, site, distance) %>%
    pivot_wider(names_from = site,
                values_from = distance) %>%
    select(-indv) %>%
    rename_with(~ sub("site", "", .x))
  return(output)
}

simulation_olg <- function(iterations = 100, split = 20, seed = 94706, 
                           I = 500, turnover = 50, 
                           J = 100, S = 50, phi = "unif"){
  # I: number of birders
  # J: number of sites
  # S: number of species
  # iterations: number of periods in each simulation
  # split: period in which final species becomes present at single site
  # turnover: number of birders replaced each period
  set.seed(seed = seed)
  P <- matrix(rbinom(n = J*S, size = 1, prob = .5),
              nrow = J,
              ncol = S) # species presence/absence matrix 
  P[, S] <- rep.int(0, times = J) # cassia crossbill initially absent
  
  V <- matrix(runif(n = I*S, min = 0, max = 1),
              nrow = I,
              ncol = S) # individual values of seeing each species
  
  L <- matrix(1,
              nrow = I,
              ncol = S) # species seen, 1 means species not seen, 0 means seen
  
  if (phi == "unif"){phi <- runif(n = I, min = 0, max = 1)}# relative utility value of new species
  else if (phi == "beta"){phi <- rbeta(n = I, beta_params[1], beta_params[2])} 
  else {stop("phi distribution unrecognized")}
  
  home_locations <- data.frame(indv = paste0("indv", 1:I),
                               x.i = runif(I, 0, 1),
                               y.i = runif(I, 0, 1))
  
  site_locations <- data.frame(site = paste0("site", 1:J),
                               x.j = runif(J, 0, 1),
                               y.j = runif(J, 0, 1))
  
  visit_data <- data.frame(site = as.factor(1:J))
  for (i in 1:iterations){
    set.seed(seed + i)
    if (i == split){
      P[, S] <- c(rep.int(0, times = J-1), 1)
    }
    if (i > 1){
      V <- rbind(V[(turnover+1):I,], 
                 matrix(runif(n = turnover*S, min = 0, max = 1),
                        nrow = turnover, 
                        ncol = S))
      L <- rbind(L[(turnover+1):I,], 
                 matrix(1,
                        nrow = turnover, 
                        ncol = S))
      if (phi == "unif"){
        phi <- c(phi[(turnover+1):I], runif(n = turnover, min = 0, max = 1))
      }
      else if (phi == "beta"){
        phi <- c(phi[(turnover+1):I], rbeta(n = turnover, 
                                            beta_params[1], beta_params[2]))
      }
      home_locations <- rbind(home_locations[(turnover+1):I,],
                              data.frame(indv = paste0("indv", 
                                                       (1:turnover)+
                                                         (I+turnover*i)),
                                         x.i = runif(turnover, 0, 1),
                                         y.i = runif(turnover, 0, 1)))
      
    }
    D <- distance_matrix(home_locations, site_locations) # I x J matrix
    E <- matrix(rnorm(I * J, 0, .5), nrow = J, ncol = I)
    site_values <- log(V %*% t(P))+phi*(L %*% t(P))^(1/2)-D+E
    visit_matrix <- apply(site_values, 1, function(z){1 * (z == max(z))})
    visit_data <- cbind(visit_data, rowSums(visit_matrix))
    L <- L - t(visit_matrix) %*% P
    L[L < 0] <- 0
  }
  colnames(visit_data) <- c("site", 1:iterations)
  return(visit_data)
}

simdata_collected_unif <- data.frame()
for (run in 1:run_count){
  simdata <- simulation_olg(I = I, J = J, S = S, iterations = iterations,
                            split = split, turnover = turnover, 
                            seed = run*94706, phi = "unif") %>%
    pivot_longer(!site, names_to = "period", values_to = "visits") %>%
    mutate(period = as.integer(period),
           visits = as.integer(visits),
           treated = ifelse(as.numeric(site) == max(as.numeric(site)), 1, 0),
           post = ifelse(as.numeric(period) >= split, 1, 0),
           run = run)
  simdata_collected_unif <- rbind(simdata_collected_unif, simdata)
}

write_csv(simdata_collected_unif, "simdata_collected_unif.csv")

simdata_collected_beta <- data.frame()
for (run in 1:run_count){
  simdata <- simulation_olg(I = I, J = J, S = S, iterations = iterations,
                            split = split, turnover = turnover, 
                            seed = run*94706, phi = "beta") %>%
    pivot_longer(!site, names_to = "period", values_to = "visits") %>%
    mutate(period = as.integer(period),
           visits = as.integer(visits),
           treated = ifelse(as.numeric(site) == max(as.numeric(site)), 1, 0),
           post = ifelse(as.numeric(period) >= split, 1, 0),
           run = run)
  simdata_collected_beta <- rbind(simdata_collected_beta, simdata)
}
  
write_csv(simdata_collected_beta, "simdata_collected_beta.csv")
