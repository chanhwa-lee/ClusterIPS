time = proc.time()

###------------------- Load libraries ----------------------###
suppressPackageStartupMessages(require(dplyr))
suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(SuperLearner))

###----------- Read arguments (simulation parameters) ---------------###
option_list = list(
  make_option(c("-n", "--n"), action = "store", default = NA, type = "integer",
              help = paste0("Number of clusters per each simulation")),
  make_option(c("-K", "--nsplits"), action = "store", default = NA, type = "integer",
              help = paste0("Number of splits for sample splitting estimator")),
  make_option(c("-r", "--r"), action = "store", default = NA, type = "integer",
              help = paste0("Number of binary vector sampling for outcome reg computation")),
  make_option(c("-s", "--s"), action = "store", default = NA, type = "integer",
              help = paste0("Number of sample splitting repetition"))
)

opt = parse_args(OptionParser(option_list = option_list))

n = opt$n              # Number of clusters per simulation
nsplits = opt$nsplits  # Number of sample splitting
r = opt$r              # Number of binary vector sampling
s = opt$s              # Number of sample splitting repetition
deltas <- readRDS("deltas.rds")

## Wrap with as.numeric to coerce from character variable
task_ID <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

print("[Simulation setting]")
print(paste0("n: ", n))
print(paste0("deltas: ", paste(signif(deltas, 4), collapse = ", ")))
print(paste0("K: ", nsplits))
print(paste0("r: ", r))
print(paste0("s: ", s))


###---------------------- Help Functions ----------------------###

### Simulation Data set generation ###
### Input:
###   n = number of clusters
###   N = number of individuals per cluster
### Output:
###   list of cluster-level dataframe of (id, Y, A, X1, X2)
data.sim <- function(n){
  
  data.list <- list()
  
  N = sample(x = 5:20, size = n, replace = T)
  
  for(i in 1:n){
    
    # Step1. Data generation
    X1 = rnorm(N[i], 0, 1)
    X2 = rbinom(N[i], 1, 0.5)
    C  = rnorm(1, 0, 1)      # Cluster-level covariate
    
    # Step2. Treatment model
    pi <- plogis(0.1 + 0.2*abs(X1) + 0.2*abs(X1)*X2 + 0.1*I(C>0))
    A = rbinom(N[i], 1, pi)
    g.A = (sum(A)-A)/(N[i]-1)
    
    # Step3. Outcome model
    p.Y = plogis(3 - 2*A - g.A - 1.5*abs(X1) + 2*X2 - 3*abs(X1)*X2 - 2*I(C>0))
    Y = rbinom(N[i], 1, p.Y)
    
    data.list[[i]] <- data.frame(id = i,
                                 Y = Y,
                                 A = A,
                                 X1 = X1,
                                 X2 = X2,
                                 C = C)
  }
  
  return(dplyr::bind_rows(data.list))
  
}


### Estimator Main Function ###
### this function requires the following inputs:
### dat: dataframe with cols cluster id 'i', outc 'Y', trt 'A', trt proportion 'g.A'
### X.trt: dataframe of covariates for treatment regression included in dat
### X.out: dataframe of covariates for outcome regression inculded in dat (excluding 'A')
### deltas: must include 1 (delta = 1 is the baseline value of delta, compared to other delta values)
### nsplits: number of sample splits
### r: number of random binary vector from A(N_i) w. uniform sampling, instead of summing all vectors in A(N)
### type: "par" (parametric) or "nonpar" (nonparametric)
### NOTE: dat, X.trt, X.out should all have the same number of rows


estimator <- function(dat, X.trt, X.out, deltas, nsplits, r){
  
  ## If deltas does not include 1, then add it to the list
  if(!(1 %in% deltas)) deltas <- c(deltas, 1)
  
  ## Split dat into lists to get N for each cluster
  dat.list <- split(dat, f = dat$id)
  
  ## setup storage
  n <- length(dat.list)        # Number of clusters
  N <- sapply(dat.list, nrow)  # Number of individuals in each cluster
  N.limit <- cumsum(c(0,N))
  
  est.par <- se.par <-  data.frame(delta = deltas,
                                   mu = 0, mu_1 = 0, mu_0 = 0,
                                   de = 0, se_1 = 0, se_0 = 0,
                                   oe = 0, te = 0)
  
  est.nonpar <- se.nonpar <-  data.frame(delta = deltas,
                                         mu = 0, mu_1 = 0, mu_0 = 0,
                                         de = 0, se_1 = 0, se_0 = 0,
                                         oe = 0, te = 0)
  
  pi.par <- numeric(sum(N))
  mu.par <- numeric(sum(N))
  
  pi.nonpar <- numeric(sum(N))
  mu.nonpar <- numeric(sum(N))
  
  ifvals.mu.par   <- matrix(0, nrow = length(deltas), ncol = n)   # IF for mu per cluster
  ifvals.mu_1.par <- matrix(0, nrow = length(deltas), ncol = n)   # IF for mu_1 per cluster
  ifvals.mu_0.par <- matrix(0, nrow = length(deltas), ncol = n)   # IF for mu_0 per cluster
  
  ifvals.mu.nonpar   <- matrix(0, nrow = length(deltas), ncol = n)   # IF for mu per cluster
  ifvals.mu_1.nonpar <- matrix(0, nrow = length(deltas), ncol = n)   # IF for mu_1 per cluster
  ifvals.mu_0.nonpar <- matrix(0, nrow = length(deltas), ncol = n)   # IF for mu_0 per cluster  
  
  fold <- sample(1:n, replace = F) %% nsplits + 1
  foldlong <- rep(fold, times = N)
  
  ## Fit treatment and outcome model at each split and evaluate IF values
  for (split in 1:nsplits) {
    
    print(paste("   Split:", split))
    
    train.idx = which(foldlong != split)
    eval.idx  = which(foldlong == split)
    
    ### 1. Main-effect only GLM (Logistic regression)
    pi.fit.par <- glm(A ~ ., data = cbind(A = dat$A, X.trt)[train.idx,], family = "binomial")
    
    pi.par[eval.idx] <- predict(pi.fit.par, newdata = X.trt[eval.idx,], type = "response")
    
    mu.fit.par <- glm(Y ~ ., data = cbind(Y = dat$Y, A = dat$A, g.A = dat$g.A, X.out)[train.idx,], family = "binomial")
    
    mu.par[eval.idx] <- predict(mu.fit.par, newdata = cbind(A = dat$A, g.A = dat$g.A, X.out)[eval.idx,], type = "response")
    
    ### 2. SuperLearner (RF & GLM)
    pi.fit.nonpar = SuperLearner(Y = dat$A[train.idx], X = X.trt[train.idx,],
                                 family = binomial(), SL.library = c("SL.glm", "SL.ranger", "SL.gam", "SL.nnet"))
    
    pi.nonpar[eval.idx] <- predict(pi.fit.nonpar, X.trt[eval.idx,], onlySL = TRUE)$pred
    
    mu.fit.nonpar <- SuperLearner(Y = dat$Y[train.idx], X = cbind(A = dat$A, g.A = dat$g.A, X.out)[train.idx,],
                                  family = binomial(), SL.library = c("SL.glm", "SL.ranger", "SL.gam", "SL.nnet"))
    
    mu.nonpar[eval.idx] <- predict(mu.fit.nonpar, cbind(A = dat$A, g.A = dat$g.A, X.out)[eval.idx,], onlySL = TRUE)$pred
    
    ## evaluate IF values for each cluster in eval fold
    for(i in which(fold == split)){
      
      # print(paste("   Cluster:", i))
      
      dat.sub.idx <- (N.limit[i]+1):N.limit[i+1]
      dat.sub <- dat[dat.sub.idx,]
      
      pi.sub.par  <- pi.par[dat.sub.idx]
      mu.sub.par  <- mu.par[dat.sub.idx]
      
      pi.sub.nonpar  <- pi.nonpar[dat.sub.idx]
      mu.sub.nonpar  <- mu.nonpar[dat.sub.idx]
      
      # random binary vector from A(N_i) w. uniform sampling, instead of summing all vectors in A(N)
      
      # Help function
      agg.fun <- function(a, fun){
        aggregate(a, list(rep(1:r, each = N[i])), fun)[,-1]
      }
      
      a <- sample(c(0,1), r*N[i], replace = T)                                   # rN[i] vector
      g.a = (rep(agg.fun(a, sum), each = N[i]) - a) / (N[i]-1)
      
      
      ### 1. Parametric
      
      ## Outcome regression
      mu.a.par <- predict(mu.fit.par,
                          newdata = cbind(data.frame(A = a, g.A = g.a),
                                          bind_rows(replicate(r, X.out[dat.sub.idx,], simplify = FALSE))),
                          type = "response")
      
      pi.delta.sub.par <- pi.sub.par %o% deltas / (pi.sub.par %o% deltas + 1 - pi.sub.par)       # N[i] x length(deltas) matrix
      
      cluster.prob.ratio.par <- deltas^(sum(dat.sub$A)) / 
        apply(pi.sub.par %o% deltas + 1 - pi.sub.par, 2, prod)                           # length(deltas) vector
      
      pi.delta.sub.rep.par <- rep(1, r) %x% pi.delta.sub.par                             # rN[i] x length(deltas) matrix (replicate pi.delta.sub r times)
      
      cluster.prob.par <- agg.fun(a*pi.delta.sub.rep.par + (1-a)*(1-pi.delta.sub.rep.par), prod)   # r x length(deltas) matrix
      
      add.factor.par <-
        (2*a-1) %o% deltas * rep(dat.sub$A - pi.sub.par, r) / 
        (rep(1, r) %x% (pi.sub.par %o% deltas + 1 - pi.sub.par)) /
        (a*rep(pi.sub.par, r) %o% deltas + (1-a)*(1-rep(pi.sub.par, r)))                 # rN[i] x length(deltas) matrix
      
      out.reg.mu.par <- apply(agg.fun(mu.a.par, mean) * (1 + agg.fun(add.factor.par, sum)) * 
                                2^(N[i]) * cluster.prob.par, 2, mean)                    # length(deltas) vector
      
      out.reg.mu_1.par <-
        apply(agg.fun( (mu.a.par * a) %o% rep(1, length(deltas)) / pi.delta.sub.rep.par * 
                         (1 + as.matrix(agg.fun(add.factor.par, sum)) %x% rep(1,N[i]) - add.factor.par), mean) *
                2^(N[i]) * cluster.prob.par, 2, mean)                                 # length(deltas) vector
      
      out.reg.mu_0.par <-
        apply(agg.fun( (mu.a.par * (1-a)) %o% rep(1, length(deltas)) / (1-pi.delta.sub.rep.par) * 
                         (1 + as.matrix(agg.fun(add.factor.par, sum)) %x% rep(1,N[i]) - add.factor.par), mean) *
                2^(N[i]) * cluster.prob.par, 2, mean)                                # length(deltas) vector
      
      ## Bias correction part for mu
      bias.cor.mu.par <- mean(dat.sub$Y - mu.sub.par, na.rm = T) * cluster.prob.ratio.par    # length(deltas) vector
      
      ## Bias correction part for mu_1
      bias.cor.mu_1.par <-
        apply( (dat.sub$Y - mu.sub.par) * I(dat.sub$A==1) / pi.delta.sub.par, 2, mean) * 
        cluster.prob.ratio.par                                                       # length(deltas) vector
      
      ## Bias correction part for mu_0
      bias.cor.mu_0.par <-
        apply( (dat.sub$Y - mu.sub.par) * I(dat.sub$A==0) / (1-pi.delta.sub.par), 2, mean) * 
        cluster.prob.ratio.par                                                       # length(deltas) vector
      
      ## IF for mu at cluster i
      ifvals.mu.par[,i]   <- out.reg.mu.par + bias.cor.mu.par
      
      ## IF for mu_1 at cluster i
      ifvals.mu_1.par[,i] <- out.reg.mu_1.par + bias.cor.mu_1.par
      
      ## IF for mu_0 at cluster i
      ifvals.mu_0.par[,i] <- out.reg.mu_0.par + bias.cor.mu_0.par
      
      
      ### 2. Nonparametric
      
      ## Outcome regression
      mu.a.nonpar <- predict(mu.fit.nonpar, 
                             cbind(data.frame(A = a, g.A = g.a), bind_rows(replicate(r, X.out[dat.sub.idx,], simplify = FALSE))),
                             onlySL = TRUE)$pred
      
      
      pi.delta.sub.nonpar <- pi.sub.nonpar %o% deltas / (pi.sub.nonpar %o% deltas + 1 - pi.sub.nonpar)       # N[i] x length(deltas) matrix
      
      cluster.prob.ratio.nonpar <- deltas^(sum(dat.sub$A)) / 
        apply(pi.sub.nonpar %o% deltas + 1 - pi.sub.nonpar, 2, prod)                           # length(deltas) vector
      
      pi.delta.sub.rep.nonpar <- rep(1, r) %x% pi.delta.sub.nonpar                             # rN[i] x length(deltas) matrix (replicate pi.delta.sub r times)
      
      cluster.prob.nonpar <- agg.fun(a*pi.delta.sub.rep.nonpar + (1-a)*(1-pi.delta.sub.rep.nonpar), prod)   # r x length(deltas) matrix
      
      add.factor.nonpar <-
        (2*a-1) %o% deltas * rep(dat.sub$A - pi.sub.nonpar, r) / 
        (rep(1, r) %x% (pi.sub.nonpar %o% deltas + 1 - pi.sub.nonpar)) /
        (a*rep(pi.sub.nonpar, r) %o% deltas + (1-a)*(1-rep(pi.sub.nonpar, r)))                 # rN[i] x length(deltas) matrix
      
      out.reg.mu.nonpar <- apply(agg.fun(mu.a.nonpar, mean) * (1 + agg.fun(add.factor.nonpar, sum)) * 
                                   2^(N[i]) * cluster.prob.nonpar, 2, mean)                    # length(deltas) vector
      
      out.reg.mu_1.nonpar <-
        apply(agg.fun( (mu.a.nonpar * a)[,rep(1, length(deltas))] / pi.delta.sub.rep.nonpar * 
                         (1 + as.matrix(agg.fun(add.factor.nonpar, sum)) %x% rep(1,N[i]) - add.factor.nonpar), mean) *
                2^(N[i]) * cluster.prob.nonpar, 2, mean)                                 # length(deltas) vector
      
      out.reg.mu_0.nonpar <-
        apply(agg.fun( (mu.a.nonpar * (1-a))[,rep(1, length(deltas))] / (1-pi.delta.sub.rep.nonpar) * 
                         (1 + as.matrix(agg.fun(add.factor.nonpar, sum)) %x% rep(1,N[i]) - add.factor.nonpar), mean) *
                2^(N[i]) * cluster.prob.nonpar, 2, mean)                                # length(deltas) vector
      
      ## Bias correction part for mu
      bias.cor.mu.nonpar <- mean(dat.sub$Y - mu.sub.nonpar, na.rm = T) * cluster.prob.ratio.nonpar    # length(deltas) vector
      
      ## Bias correction part for mu_1
      bias.cor.mu_1.nonpar <-
        apply( (dat.sub$Y - mu.sub.nonpar) * I(dat.sub$A==1) / pi.delta.sub.nonpar, 2, mean) * 
        cluster.prob.ratio.nonpar                                                       # length(deltas) vector
      
      ## Bias correction part for mu_0
      bias.cor.mu_0.nonpar <-
        apply( (dat.sub$Y - mu.sub.nonpar) * I(dat.sub$A==0) / (1-pi.delta.sub.nonpar), 2, mean) * 
        cluster.prob.ratio.nonpar                                                       # length(deltas) vector
      
      ## IF for mu at cluster i
      ifvals.mu.nonpar[,i]   <- out.reg.mu.nonpar + bias.cor.mu.nonpar
      
      ## IF for mu_1 at cluster i
      ifvals.mu_1.nonpar[,i] <- out.reg.mu_1.nonpar + bias.cor.mu_1.nonpar
      
      ## IF for mu_0 at cluster i
      ifvals.mu_0.nonpar[,i] <- out.reg.mu_0.nonpar + bias.cor.mu_0.nonpar
      
      
    }
    
  }
  
  ## helper function
  IF.to.est <- function(IF){
    est = mean(aggregate(IF, by = list(fold), mean)[,-1])
    se  = sqrt(1/n * (mean(aggregate(IF, by = list(fold), function(phi) mean(phi^2))[,-1]) - est^2))
    return(c(est, se))
  }
  
  ## Standard (delta = 1)
  delta1.idx = which(deltas == 1)
  
  for(delta.idx in seq_len(length(deltas))){
    
    
    ### Parametric
    # delta specific estimates and se estimates
    est.par[delta.idx, "mu"]   <- IF.to.est(ifvals.mu.par  [delta.idx,])[1]
    est.par[delta.idx, "mu_1"] <- IF.to.est(ifvals.mu_1.par[delta.idx,])[1]
    est.par[delta.idx, "mu_0"] <- IF.to.est(ifvals.mu_0.par[delta.idx,])[1]
    est.par[delta.idx, "de"]   <- IF.to.est(ifvals.mu_1.par[delta.idx,] - ifvals.mu_0.par[delta.idx,])[1]
    
    se.par[delta.idx, "mu"]   <- IF.to.est(ifvals.mu.par  [delta.idx,])[2]
    se.par[delta.idx, "mu_1"] <- IF.to.est(ifvals.mu_1.par[delta.idx,])[2]
    se.par[delta.idx, "mu_0"] <- IF.to.est(ifvals.mu_0.par[delta.idx,])[2]
    se.par[delta.idx, "de"]   <- IF.to.est(ifvals.mu_1.par[delta.idx,] - ifvals.mu_0.par[delta.idx,])[2]
    
    # delta versus delta==1 comparison estimates and se estimates
    est.par[delta.idx, "se_1"] <- IF.to.est(ifvals.mu_1.par[delta.idx,] - ifvals.mu_1.par[delta1.idx,])[1]
    est.par[delta.idx, "se_0"] <- IF.to.est(ifvals.mu_0.par[delta.idx,] - ifvals.mu_0.par[delta1.idx,])[1]
    est.par[delta.idx, "oe"]   <- IF.to.est(ifvals.mu.par  [delta.idx,] - ifvals.mu.par  [delta1.idx,])[1]
    est.par[delta.idx, "te"]   <- IF.to.est(ifvals.mu_1.par[delta.idx,] - ifvals.mu_0.par[delta1.idx,])[1]
    
    se.par[delta.idx, "se_1"] <- IF.to.est(ifvals.mu_1.par[delta.idx,] - ifvals.mu_1.par[delta1.idx,])[2]
    se.par[delta.idx, "se_0"] <- IF.to.est(ifvals.mu_0.par[delta.idx,] - ifvals.mu_0.par[delta1.idx,])[2]
    se.par[delta.idx, "oe"]   <- IF.to.est(ifvals.mu.par  [delta.idx,] - ifvals.mu.par  [delta1.idx,])[2]
    se.par[delta.idx, "te"]   <- IF.to.est(ifvals.mu_1.par[delta.idx,] - ifvals.mu_0.par[delta1.idx,])[2]
    
    ### Nonparametric
    # delta specific estimates and se estimates
    est.nonpar[delta.idx, "mu"]   <- IF.to.est(ifvals.mu.nonpar  [delta.idx,])[1]
    est.nonpar[delta.idx, "mu_1"] <- IF.to.est(ifvals.mu_1.nonpar[delta.idx,])[1]
    est.nonpar[delta.idx, "mu_0"] <- IF.to.est(ifvals.mu_0.nonpar[delta.idx,])[1]
    est.nonpar[delta.idx, "de"]   <- IF.to.est(ifvals.mu_1.nonpar[delta.idx,] - ifvals.mu_0.nonpar[delta.idx,])[1]
    
    se.nonpar[delta.idx, "mu"]   <- IF.to.est(ifvals.mu.nonpar  [delta.idx,])[2]
    se.nonpar[delta.idx, "mu_1"] <- IF.to.est(ifvals.mu_1.nonpar[delta.idx,])[2]
    se.nonpar[delta.idx, "mu_0"] <- IF.to.est(ifvals.mu_0.nonpar[delta.idx,])[2]
    se.nonpar[delta.idx, "de"]   <- IF.to.est(ifvals.mu_1.nonpar[delta.idx,] - ifvals.mu_0.nonpar[delta.idx,])[2]
    
    # delta versus delta==1 comparison estimates and se estimates
    est.nonpar[delta.idx, "se_1"] <- IF.to.est(ifvals.mu_1.nonpar[delta.idx,] - ifvals.mu_1.nonpar[delta1.idx,])[1]
    est.nonpar[delta.idx, "se_0"] <- IF.to.est(ifvals.mu_0.nonpar[delta.idx,] - ifvals.mu_0.nonpar[delta1.idx,])[1]
    est.nonpar[delta.idx, "oe"]   <- IF.to.est(ifvals.mu.nonpar  [delta.idx,] - ifvals.mu.nonpar  [delta1.idx,])[1]
    est.nonpar[delta.idx, "te"]   <- IF.to.est(ifvals.mu_1.nonpar[delta.idx,] - ifvals.mu_0.nonpar[delta1.idx,])[1]
    
    se.nonpar[delta.idx, "se_1"] <- IF.to.est(ifvals.mu_1.nonpar[delta.idx,] - ifvals.mu_1.nonpar[delta1.idx,])[2]
    se.nonpar[delta.idx, "se_0"] <- IF.to.est(ifvals.mu_0.nonpar[delta.idx,] - ifvals.mu_0.nonpar[delta1.idx,])[2]
    se.nonpar[delta.idx, "oe"]   <- IF.to.est(ifvals.mu.nonpar  [delta.idx,] - ifvals.mu.nonpar  [delta1.idx,])[2]
    se.nonpar[delta.idx, "te"]   <- IF.to.est(ifvals.mu_1.nonpar[delta.idx,] - ifvals.mu_0.nonpar[delta1.idx,])[2]
    
  }
  
  return(list(est.par = est.par, se.par = se.par, est.nonpar = est.nonpar, se.nonpar = se.nonpar))
  
}



### Simulation ####

m = 20

# ### Test ###
# m = 3
# n = 10
# deltas = c(0.5, 1, 2)
# nsplits = 2
# r = 100
# s = 1
# ############


estimates.par <- list()
ses.par <- list()

estimates.nonpar <- list()
ses.nonpar <- list()

for(sim.id in 1:m){

  print(paste("Sim.id", sim.id))
  
  ## Simulation data generation
  data <- data.sim(n)
  
  ## Change format of data to input of imu_b function
  dat <- data %>% 
    select(id, Y, A) %>% 
    group_by(id) %>%
    mutate(g.A = (sum(A) - A) / (n()-1)) %>%
    ungroup()
  
  X.trt <- data %>% select(-c(id, Y, A))
  X.out <- data %>% select(-c(id, Y, A))

  print("Data generated and format changed")
  
  ## Repeat computing estimator for robustness from sample splitting
  est.par.list <- list()
  se.par.list <- list()
  
  est.nonpar.list <- list()
  se.nonpar.list <- list()
  
  for(k in 1:s){
    
    result <- estimator(dat, X.trt, X.out, deltas, nsplits, r)
    
    est.par.list[[k]] <- result$est.par
    se.par.list[[k]] <- result$se.par
    
    est.nonpar.list[[k]] <- result$est.nonpar
    se.nonpar.list[[k]] <- result$se.nonpar
    
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
    
    se.par[delta.idx, ] <- sqrt(apply((se.par.list %>% filter(delta == deltas[delta.idx]))^2 +
                              (est.par.list %>% filter(delta == deltas[delta.idx]) - est.par[rep(delta.idx,s), ])^2, 
                              2, median))
    
    est.nonpar[delta.idx, ] <- apply(est.nonpar.list %>% filter(delta == deltas[delta.idx]), 2, median)
    
    se.nonpar[delta.idx, ] <- sqrt(apply((se.nonpar.list %>% filter(delta == deltas[delta.idx]))^2 +
                              (est.nonpar.list %>% filter(delta == deltas[delta.idx]) - est.nonpar[rep(delta.idx,s), ])^2, 
                              2, median))
  }
  
  estimates.par[[sim.id]] <- est.par
  ses.par[[sim.id]] <- se.par
  
  estimates.nonpar[[sim.id]] <- est.nonpar
  ses.nonpar[[sim.id]] <- se.nonpar
  
  print("Estimates and SE estimates computed")
  print("")
}

print(estimates.par)
print(ses.par)

print(estimates.nonpar)
print(ses.nonpar)

###-------- Save simulated estimator list as Rdata --------###

## Save output
saveRDS(estimates.par, file = paste0("Rdata/estimate.par_id", task_ID,".rds"))
saveRDS(ses.par, file = paste0("Rdata/se.par_id", task_ID,".rds"))

saveRDS(estimates.nonpar, file = paste0("Rdata/estimate.nonpar_id", task_ID,".rds"))
saveRDS(ses.nonpar, file = paste0("Rdata/se.nonpar_id", task_ID,".rds"))


###------------------------------------------------------------###

## Stop timer and report total run time

script.time = proc.time() - time
print(paste0("Total run time was ", script.time[3], " seconds"))