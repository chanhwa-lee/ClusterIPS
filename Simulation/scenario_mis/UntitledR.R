#*******************************************************************************
#**********           Causal estimands estimation under               **********
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
# This file requires simulated estimands and estimates under "Simulation/scenario_mis/Data_parallel/".
# to compute the simulation metrics and save the result tables in "Simulation/scenario_mis/".
# as csv files.
############################################

### Libraries ###
library(dplyr)

###------ Read estimand ------###
file.list <- list.files("Simulation/scenario_mis/Data_parallel", pattern = "estimand.*rds")
D <- length(file.list)

deltas <- c(0.5,1,2)
estimands <- data.frame(delta = rep(0, length(deltas)),
                        mu = 0, mu_1 = 0, mu_0 = 0,
                        de = 0, se_1 = 0, se_0 = 0,
                        oe = 0, te = 0)

for(file in file.list){
  estimands <- estimands + readRDS(paste0("Simulation/scenario_mis/Data_parallel/", file))
}

print(paste0(D, " estimand Rdata files were loaded"))

estimands = estimands / D

print(estimands)



###------ Read estimator ------###

est.par.file.list <- list.files("Simulation/scenario_mis/Data_parallel/", pattern = "estimate.par.*rds")
se.par.file.list <- list.files("Simulation/scenario_mis/Data_parallel/", pattern = "se.par.*rds")

est.nonpar.file.list <- list.files("Simulation/scenario_mis/Data_parallel/", pattern = "estimate.nonpar.*rds")
se.nonpar.file.list <- list.files("Simulation/scenario_mis/Data_parallel/", pattern = "se.nonpar.*rds")

D <- length(est.par.file.list)

for(i in 1:D){
  est.par.file = est.par.file.list[i]
  se.par.file  = se.par.file.list[i]
  est.nonpar.file = est.nonpar.file.list[i]
  se.nonpar.file  = se.nonpar.file.list[i]
  
  if(i == 1){
    est.par.list    = bind_rows(readRDS(paste0("Simulation/scenario_mis/Data_parallel/", est.par.file)))
    se.par.list     = bind_rows(readRDS(paste0("Simulation/scenario_mis/Data_parallel/", se.par.file)))
    est.nonpar.list = bind_rows(readRDS(paste0("Simulation/scenario_mis/Data_parallel/", est.nonpar.file)))
    se.nonpar.list  = bind_rows(readRDS(paste0("Simulation/scenario_mis/Data_parallel/", se.nonpar.file)))
  }else{
    est.par.list = rbind(est.par.list, bind_rows(readRDS(paste0("Simulation/scenario_mis/Data_parallel/", est.par.file))))
    se.par.list  = rbind(se.par.list,  bind_rows(readRDS(paste0("Simulation/scenario_mis/Data_parallel/", se.par.file))))
    est.nonpar.list = rbind(est.nonpar.list, bind_rows(readRDS(paste0("Simulation/scenario_mis/Data_parallel/", est.nonpar.file))))
    se.nonpar.list  = rbind(se.nonpar.list,  bind_rows(readRDS(paste0("Simulation/scenario_mis/Data_parallel/", se.nonpar.file))))
  }
}

print(paste0(D, " estimate Rdata files were loaded"))

### Simulation metrics tatble ###
for(delta.idx in seq_len(length(deltas))){
  
  delta = deltas[delta.idx]
  print(paste("delta =", delta, "/ Nonpara result / Para result"))
  
  est.par = est.par.list %>% filter(delta == deltas[delta.idx]) %>% select(-delta)
  se.par  = se.par.list  %>% filter(delta == deltas[delta.idx]) %>% select(-delta)
  
  est.nonpar = est.nonpar.list %>% filter(delta == deltas[delta.idx]) %>% select(-delta)
  se.nonpar  = se.nonpar.list  %>% filter(delta == deltas[delta.idx]) %>% select(-delta)
  
  estimand = as.numeric(estimands %>% filter(delta == deltas[delta.idx]) %>% select(-delta))
  
  result.table <- cbind(
    
    # True Estimand
    estimand,                                                                    
    
    # Nonpara result    
    apply(est.nonpar, 2, mean) - estimand,                                          # Bias
    
    sqrt(apply((est.nonpar - rep(1, nrow(est.nonpar)) %o% estimand)^2, 2, mean)),      # RMSE
    
    apply(se.nonpar, 2, mean)   ,                                                   # Avg. SE
    
    apply(est.nonpar, 2, sd)    ,                                                   # Emp. SE
    
    colMeans(
      (est.nonpar - 1.96 * se.nonpar < rep(1, nrow(est.nonpar)) %o% estimand) & 
        (rep(1, nrow(est.nonpar)) %o% estimand < est.nonpar + 1.96 * se.nonpar), 
      na.rm = T),                                                                  # 95% coverage
    
    # Para result
    apply(est.par, 2, mean) - estimand,                                       # Bias
    
    sqrt(apply((est.par - rep(1, nrow(est.par)) %o% estimand)^2, 2, mean)),      # RMSE
    
    apply(se.par, 2, mean)   ,                                                   # Avg. SE
    
    apply(est.par, 2, sd)    ,                                                   # Emp. SE
    
    colMeans(
      (est.par - 1.96 * se.par < rep(1, nrow(est.par)) %o% estimand) & 
        (rep(1, nrow(est.par)) %o% estimand < est.par + 1.96 * se.par), 
      na.rm = T),                                                                 # 95% coverage
    
    # RMSE ratio
    sqrt(apply((est.nonpar - rep(1, nrow(est.nonpar)) %o% estimand)^2, 2, mean)) /
      sqrt(apply((est.par - rep(1, nrow(est.par)) %o% estimand)^2, 2, mean))
    
  )
  
  colnames(result.table) <- c("Truth", 
                              "Bias", "RMSE", "ASE", "ESE", "Cov",               # Nonpar
                              "Bias", "RMSE", "ASE", "ESE", "Cov",               # Par
                              "RMSERatio")
  
  
  print(round(result.table, 3))
  
  print("")
  
  write.csv(x = result.table, file = paste0("Simulation/scenario_mis/Result_parallel/SimulResult_Mis_delta_",delta,".csv"))
  
}
