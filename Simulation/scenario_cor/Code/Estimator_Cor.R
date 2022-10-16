#*******************************************************************************
#**********           Causal estimands estimation under               **********
#**********           the Incremental Propensiy Score policy          **********
#**********           (Scenario: Correctly specified)                 **********
#**********                                                           **********	
#**********           Written by:				                              **********
#**********           Chanhwa Lee - chanhwa@email.unc.edu             **********
#**********                                                           **********
#**********           Version: 1.0                                    **********
#**********           Jul 21, 2022                                    **********
#*******************************************************************************

############################################
# This file estimates the causal estimands under the IPS policy with correctly
# specified nuisance function scenario.
# The code will take a lot of time, so it is recommended to use parallel computing.
# To parallelize, submit jobs in for(simul.id in 1:D){...} separately
############################################

### Libraries ###
library(dplyr)
library(SuperLearner)

### Load estimator function ###
source("Simulation/scenario_cor/Code/Helpfunc_Cor.R")

time = proc.time()

D = 10                                                                           ## Number of simulations
m = 10                                                                           ## Number of clusters per simulation
K = 2                                                                           ## Number of sample splitting
r = 100                                                                           ## Number of binary vector sampling
S = 1                                                                           ## Number of sample splitting repetition
deltas <- c(0.5,1,2)

print("[Simulation setting]")
print(paste0("D: ", D))
print(paste0("m: ", m))
print(paste0("deltas: ", paste(signif(deltas, 4), collapse = ", ")))
print(paste0("K: ", K))
print(paste0("r: ", r))
print(paste0("S: ", S))

### Generate data set and Compute estimators D times ###
for(sim.id in 1:D){
  
  print(paste("Sim.id", sim.id))
  
  ## Simulation data generation
  data <- data.cor(m)
  
  ## Change format of data to input of imu_b function
  dat <- data %>% 
    select(id, Y, A) %>% 
    group_by(id) %>%
    mutate(g.A = (sum(A) - A) / (n()-1)) %>%
    ungroup()
  
  X.trt <- data %>% select(-c(id, Y, A))
  X.out <- data %>% select(-c(id, Y, A))
  
  print("Simulation data generated and format changed")
  
  ### Repeat computing estimator for robustness from sample splitting ###
  est.par.list <- list()
  se.par.list <- list()
  
  est.nonpar.list <- list()
  se.nonpar.list <- list()
  
  for(s in 1:S){
    
    result <- estimator.cor(dat, X.trt, X.out, deltas, K, r)
    
    est.par.list[[s]] <- result$est.par
    se.par.list[[s]] <- result$se.par
    
    est.nonpar.list[[s]] <- result$est.nonpar
    se.nonpar.list[[s]] <- result$se.nonpar
    
  }
  
  est.par.list <- bind_rows(est.par.list)
  se.par.list <- bind_rows(se.par.list)
  
  est.nonpar.list <- bind_rows(est.nonpar.list)
  se.nonpar.list <- bind_rows(se.nonpar.list)
  
  est.par <- se.par <- data.frame(delta = deltas,
                                  mu = 0, mu_1 = 0, mu_0 = 0,
                                  de = 0, se_1 = 0, se_0 = 0,
                                  oe = 0, te = 0)
  
  est.nonpar <- se.nonpar <- data.frame(delta = deltas,
                                        mu = 0, mu_1 = 0, mu_0 = 0,
                                        de = 0, se_1 = 0, se_0 = 0,
                                        oe = 0, te = 0)
  
  for(delta.idx in seq_len(length(deltas))){
    
    est.par[delta.idx, ] <- apply(est.par.list %>% filter(delta == deltas[delta.idx]), 2, median)
    
    se.par[delta.idx, ]  <- apply(se.par.list  %>% filter(delta == deltas[delta.idx]), 2, median)
    
    est.nonpar[delta.idx, ] <- apply(est.nonpar.list %>% filter(delta == deltas[delta.idx]), 2, median)
    
    se.nonpar[delta.idx, ]  <- apply(se.nonpar.list  %>% filter(delta == deltas[delta.idx]), 2, median)
  }
  
  
  ###-------- Save simulated estimator list as Rdata --------###
  
  ## Save output
  saveRDS(est.par, file = paste0("Simulation/scenario_cor/Data/estimate.par_id", sim.id,".rds"))
  saveRDS(se.par,  file = paste0("Simulation/scenario_cor/Data/se.par_id", sim.id,".rds"))
  
  saveRDS(est.nonpar, file = paste0("Simulation/scenario_cor/Data/estimate.nonpar_id", sim.id,".rds"))
  saveRDS(se.nonpar,  file = paste0("Simulation/scenario_cor/Data/se.nonpar_id", sim.id,".rds"))
  
  print("Estimates and SE estimates computed")
  print("")
}


### Stop timer and report total run time ###

script.time = proc.time() - time
print(paste0("Total run time was ", script.time[3], " seconds"))