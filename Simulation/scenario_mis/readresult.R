suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(dplyr))

###------ Read estimand ------###
file.list <- list.files("./Rdata", pattern = "estimand.*rds")
M <- length(file.list)

deltas <- readRDS("deltas.rds")
estimands <- data.frame(delta = rep(0, length(deltas)),
                        mu = 0, mu_1 = 0, mu_0 = 0,
                        de = 0, se_1 = 0, se_0 = 0,
                        oe = 0, te = 0)

for(file in file.list){
  estimands <- estimands + readRDS(paste0("Rdata/", file))
}

print(paste0(M, " estimand Rdata files were loaded"))

estimands = estimands / M

print(estimands)



###------ Read estimator ------###

est.par.file.list <- list.files("./Rdata", pattern = "estimate.par.*rds")
se.par.file.list <- list.files("./Rdata", pattern = "se.par.*rds")

est.nonpar.file.list <- list.files("./Rdata", pattern = "estimate.nonpar.*rds")
se.nonpar.file.list <- list.files("./Rdata", pattern = "se.nonpar.*rds")

M <- length(est.par.file.list)

for(i in 1:M){
  est.par.file = est.par.file.list[i]
  se.par.file  = se.par.file.list[i]
  est.nonpar.file = est.nonpar.file.list[i]
  se.nonpar.file  = se.nonpar.file.list[i]
  
  if(i == 1){
    est.par.list = readRDS(paste0("Rdata/", est.par.file))
    se.par.list  = readRDS(paste0("Rdata/", se.par.file))
    est.nonpar.list = readRDS(paste0("Rdata/", est.nonpar.file))
    se.nonpar.list  = readRDS(paste0("Rdata/", se.nonpar.file))
  }else{
    est.par.list = c(est.par.list, readRDS(paste0("Rdata/", est.par.file)))
    se.par.list  = c(se.par.list,  readRDS(paste0("Rdata/", se.par.file)))
    est.nonpar.list = c(est.nonpar.list, readRDS(paste0("Rdata/", est.nonpar.file)))
    se.nonpar.list  = c(se.nonpar.list,  readRDS(paste0("Rdata/", se.nonpar.file)))
  }
}
print(paste0(M, " estimate Rdata files were loaded"))

est.par.list = bind_rows(est.par.list)
se.par.list  = bind_rows(se.par.list)
est.nonpar.list = bind_rows(est.nonpar.list)
se.nonpar.list  = bind_rows(se.nonpar.list)

for(delta.idx in seq_len(length(deltas))){
  
  delta = deltas[delta.idx]
  print(paste("delta = ", delta))
  
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
  
}
