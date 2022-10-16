#*******************************************************************************
#**********           Causal estimands computation under              **********
#**********           the Incremental Propensiy Score policy          **********
#**********           (Scenario: Mis-specified)                       **********
#**********                                                           **********	
#**********           Written by:				                              **********
#**********           Chanhwa Lee - chanhwa@email.unc.edu             **********
#**********                                                           **********
#**********           Version: 1.0                                    **********
#**********           Jul 21, 2022                                    **********
#*******************************************************************************

############################################
# This file computes the causal estimands under the IPS policy with mis-specified
# nuisance function scenario.
# The code will take a lot of time, so it is recommended to use parallel computing.
# To parallelize, submit jobs in for(simul.id in 1:D){...} separately
############################################

### Libraries ###
library(dplyr)

### Load estimand computation function ###
source("Simulation/scenario_mis/Code/Helpfunc_Mis.R")

time = proc.time()

D = 10                                                                         ## Number of simulations
m = 100                                                                      ## Number of clusters in one simulation
r = 100                                                                         ## Number of subsampling approximation
deltas <- c(0.5,1,2)                                                            ## delta values of the IPS policies

print("[Simulation setting]")
print(paste0("D: ", D))
print(paste0("m: ", m))
print(paste0("deltas: ", paste(signif(deltas, 4), collapse = ", ")))
print(paste0("r: ", r))

###----------- Estimand computation  ---------------###

for(simul.id in 1:D){
  
  ### Cluster sizes ###
  N = sample(x = 5:20, size = m, replace = T)
  
  ### Estimands computation ###
  estimands.list <- lapply(N, estimand.mis, deltas = deltas, r = r)
  
  estimands <- data.frame(delta = rep(0, length(deltas)),
                          mu = 0, mu_1 = 0, mu_0 = 0)
  
  for(i in 1:m){
    estimands = estimands + estimands.list[[i]]
  }
  
  estimands = estimands / m
  
  ### Causal effects computation using delta==1 as the standard ###
  standard = estimands %>% filter(delta == 1)
  
  estimands = estimands %>% mutate(de  = mu_1 - mu_0, 
                                   se_1 = mu_1 - standard$mu_1,
                                   se_0 = mu_0 - standard$mu_0,
                                   oe  = mu   - standard$mu,
                                   te  = mu_1 - standard$mu_0)
  
  ### Save output ###
  saveRDS(estimands, file = paste0("Simulation/scenario_mis/Data/estimand_id", simul.id,".rds"))
}

### Stop timer and report total run time ###
script.time = proc.time() - time
print(paste0("Total run time was ", script.time[3], " seconds"))