#*******************************************************************************
#**********           < WASH effect on diarrhea incidence             **********
#**********           among children in Senegal >                     **********
#**********           WASH effect estimation                          **********
#**********                                                           **********	
#**********           Written by:				                              **********
#**********           Chanhwa Lee - chanhwa@email.unc.edu             **********
#**********                                                           **********
#**********           Version: 1.0                                    **********
#**********           Sep 19, 2022                                    **********
#*******************************************************************************

############################################
# This file requires "Code/Estimator.R" and "Data/DHS/HHData.Rdata".
# Senegal DHS data is analyzed to estimate the causal estimands under the IPS policy.
# The code will take a lot of time, so it is recommended to use parallel computing.
# To parallelize, submit jobs in for(s in 1:S){...} separately
# Estimates and SE estimates are saved at "Data/DHS/result.Rdata".
############################################

### Load estimator function ###
source("Code/Estimator.R")

### Super Learner libraries ###
library(SuperLearner)
library(caret)
library(nnet)
library(glmnet)
library(earth)
library(gam)
library(gbm)
library(xgboost)     
library(kernlab)
library(polspline)
library(ranger)

SL.library = c("SL.glm", "SL.glmnet", "SL.earth", 
               "SL.gam", "SL.xgboost",
               "SL.ranger", "SL.nnet")
  
### Load preprocessed dataset ###
load("Data/DHS/HHData.Rdata")

###------------------- Estimator computation ----------------------###

### Parameter setting ###
K = 2
r = 5000
S = 3
delta_l = 1/2
delta_u = 2
num.delta = 30

time = proc.time()
inter.time = proc.time()

### dat: A data.frame with columns cluster id 'id', outcome 'Y', 
### treatment 'A', leave-one-out within-cluster treatment proportion 'g.A'
dat <- HH.Data %>% 
  group_by(cid) %>%
  mutate(g.A = ifelse(n() == 1, 0, (sum(A) - A) / (n()-1))) %>%
  ungroup() %>%
  select(cid, Y, A, g.A) %>% 
  rename(id = cid)

### X.trt: A data.frame with covariates for propensity score model ###
X.trt <- HH.Data %>%
  select(-c(hhid, cid, Y, A))

### X.out: A data.frame with covariates for outcome regression model ###
X.out <- HH.Data %>%
  select(-c(hhid, cid, Y, A))

### delta values of interest ###
deltas <- c(exp(seq(from = log(delta_l), to = 0, length.out = num.delta))[-num.delta],
            1,
            exp(seq(from = 0, to = log(delta_u), length.out = num.delta))[-1])

### Repeat sample splitting and constructing estimators for sample-splitting robust estimator ###
est.list <- list()
se.list <- list()

for(s in 1:S){
  
  print(paste("s:",s))
  
  result <- estimator(dat, X.trt, X.out, deltas, K, r, SL.library)
  
  est.list[[s]] <- result$est
  se.list[[s]] <- result$se
  
  script.time = proc.time() - inter.time
  inter.time = proc.time()
  print(paste0("Elapsed time: ", script.time[3], " seconds"))
  
}

est.list.comb = bind_rows(est.list)
se.list.comb  = bind_rows(se.list)

est.med <- se.med <- data.frame(delta = deltas,
                                mu = 0, mu_1 = 0, mu_0 = 0,
                                de = 0, se_1 = 0, se_0 = 0,
                                oe = 0, te = 0)

for(delta.idx in seq_len(length(deltas))){
  
  est.med[delta.idx, ] <- apply(est.list.comb %>% filter(delta == deltas[delta.idx]), 2, median)
  se.med[delta.idx, ]  <- apply(se.list.comb  %>% filter(delta == deltas[delta.idx]), 2, median)
  
}

result = list(est.list = est.list, se.list = se.list,
              est.med = est.med, se.med = se.med)

### Save output ###
save(result, file = "Data/DHS/result.Rdata")

### Elapsed time ###
script.time = proc.time() - time
print(paste0("Total run time was ", script.time[3], " seconds"))