#*******************************************************************************
#**********           Causal estimands computation under              **********
#**********           the Incremental Propensiy Score policy          **********
#**********           (Scenario: Correctly specified)                 **********
#**********                                                           **********	
#**********           Written by:				                              **********
#**********           Chanhwa Lee - chanhwa@email.unc.edu             **********
#**********                                                           **********
#**********           Version: 1.0                                    **********
#**********           Jul 21, 2022                                    **********
#*******************************************************************************

### Libraries ###
library(dplyr)

###----------- Estimand computation function ---------------###

#' Causal estimands computation under the Incremental Propensity Score policy
#'
#' @param N An integer. Cluster size
#' @param deltas A numeric vector. Deltas of IPS policies. Should include 1 
#' (delta = 1 is the baseline value), otherwise 1 will be included automatically
#' @param r An integer. Number of subsampling approximation 
#' (number of random binary vector sampling)
#' @return Approximated causal estimands under the IPS policy over the delta values
#' @examples 
#' estimand.cor(N = 10, deltas = c(0.5,1,2), r = 100)

estimand.cor <- function(N, deltas, r){
  
  ### Add `delta` = 1 if not in `deltas` ###
  if(!(1 %in% deltas)) deltas <- c(deltas, 1)
  
  X1 = rnorm(N, 0, 1)
  X2 = rbinom(N, 1, 0.5)
  C  = rnorm(1, 0, 1)                                                           ## Cluster-level covariate
  
  pi <- plogis(0.1 + 0.2*X1 + 0.2*X2 + 0.1*C)
  
  estimands <- data.frame(delta = deltas,
                          mu = 0, mu_1 = 0, mu_0 = 0)
  
  for(i in seq_len(length(deltas))){
    
    delta = deltas[i]
    
    pi.delta <- delta*pi / (delta*pi + 1 - pi)
    
    ### Subsampling approximation ###
    A <- rbinom(n = r*N, size = 1, prob = rep(pi.delta, r))
    
    g.A = (rep(aggregate(A, list(rep(1:r, each = N)), sum)[,-1], each = N) - A) / (N-1)
    
    p.Y = plogis(3 - 2*A - g.A - 1.5*X1 + 2*X2 - 2*C)
    
    mu   = mean(p.Y, na.rm = T)
    
    mu_1 = mean(p.Y * A / rep(pi.delta, r), na.rm = T)
    
    mu_0 = mean(p.Y * (1-A) / rep(1-pi.delta, r), na.rm = T)
    
    estimands[i, ] = c(delta, mu, mu_1, mu_0)
    
  }
  
  return(estimands)
}

###----------- Simulation data set generating function ---------------###

#' Simulation data set generating function
#'
#' @param m An integer. Number of clusters
#' @return A list of cluster-level data.frame of (id, Y, A, X1, X2)
#' @examples 
#' data.cor(n = 500)

data.cor <- function(m){
  
  data.list <- list()
  
  N = sample(x = 5:20, size = m, replace = T)
  
  for(i in 1:m){
    
    # Step1. Data generation
    X1 = rnorm(N[i], 0, 1)
    X2 = rbinom(N[i], 1, 0.5)
    C  = rnorm(1, 0, 1)      # Cluster-level covariate
    
    # Step2. Treatment model
    pi <- plogis(0.1 + 0.2*X1 + 0.2*X2 + 0.1*C)
    A = rbinom(N[i], 1, pi)
    g.A = (sum(A)-A)/(N[i]-1)
    
    # Step3. Outcome model
    p.Y = plogis(3 - 2*A - g.A - 1.5*X1 + 2*X2 - 2*C)
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

###-------------- Simulation data set Estimator function --------------------###

#' Nonparametric efficient sample splitting estimator 
#'
#' @param dat A data.frame. Columns with cluster id 'id', outcome 'Y', 
#' treatment 'A', leave-one-out within-cluster treatment proportion 'g.A'
#' @param X.trt A data.frame. Columns with covariates for propensity score model.
#' Rows should be in the same order with `dat`
#' @param X.out A data.frame. Columns with covariates for outcome regression model.
#' Rows should be in the same order with `dat`
#' @param deltas A numeric vector. Deltas of IPS policies. Should include 1 
#' (delta = 1 is the baseline value), otherwise 1 will be included automatically
#' @param K An integer. Number of sample splitting groups
#' @param r An integer. Number of subsampling approximation 
#' (number of random binary vector sampling)
#' @return Estimate and standard error estimate of the causal estimands over the
#' delta values
#' @examples 
#' estimator.cor(dat, X.trt, X.out, deltas = seq(0.5,2,0.1), K = 2, r = 100)

estimator.cor <- function(dat, X.trt, X.out, deltas, K, r){
  
  ### If deltas does not include 1, then add it to the list ###
  if(!(1 %in% deltas)) deltas <- c(deltas, 1)
  
  ### Split dat into lists to get N for each cluster ###
  dat.list <- split(dat, f = dat$id)
  
  ### Storage setup ###
  n <- length(dat.list)                                                         ## Number of clusters
  N <- sapply(dat.list, nrow)                                                   ## Number of individuals in each cluster
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
  
  ifvals.mu.par   <- matrix(0, nrow = length(deltas), ncol = n)                 ## IF for mu per cluster
  ifvals.mu_1.par <- matrix(0, nrow = length(deltas), ncol = n)                 ## IF for mu_1 per cluster
  ifvals.mu_0.par <- matrix(0, nrow = length(deltas), ncol = n)                 ## IF for mu_0 per cluster
  
  ifvals.mu.nonpar   <- matrix(0, nrow = length(deltas), ncol = n)              ## IF for mu per cluster
  ifvals.mu_1.nonpar <- matrix(0, nrow = length(deltas), ncol = n)              ## IF for mu_1 per cluster
  ifvals.mu_0.nonpar <- matrix(0, nrow = length(deltas), ncol = n)              ## IF for mu_0 per cluster  
  
  fold <- sample(1:n, replace = F) %% K + 1
  foldlong <- rep(fold, times = N)
  
  ## Fit treatment and outcome model at each k and evaluate IF values
  for (k in 1:K) {
    
    print(paste("   Split:", k))
    
    train.idx = which(foldlong != k)
    eval.idx  = which(foldlong == k)
    
    ### 1. Main-effect only GLM (Logistic regression)
    pi.fit.par <- glm(A ~ ., data = cbind(A = dat$A, X.trt)[train.idx,], family = "binomial")
    
    pi.par[eval.idx] <- predict(pi.fit.par, newdata = X.trt[eval.idx,], type = "response")
    
    mu.fit.par <- glm(Y ~ ., data = cbind(Y = dat$Y, A = dat$A, g.A = dat$g.A, X.out)[train.idx,], family = "binomial")
    
    mu.par[eval.idx] <- predict(mu.fit.par, newdata = cbind(A = dat$A, g.A = dat$g.A, X.out)[eval.idx,], type = "response")
    
    ### 2. SuperLearner (GLM, RF, GAM, NNET)
    pi.fit.nonpar = SuperLearner(Y = dat$A[train.idx], X = X.trt[train.idx,],
                                 family = binomial(), SL.library = c("SL.glm", "SL.ranger", "SL.gam", "SL.nnet"))
    
    pi.nonpar[eval.idx] <- predict(pi.fit.nonpar, X.trt[eval.idx,], onlySL = TRUE)$pred
    
    mu.fit.nonpar <- SuperLearner(Y = dat$Y[train.idx], X = cbind(A = dat$A, g.A = dat$g.A, X.out)[train.idx,],
                                  family = binomial(), SL.library = c("SL.glm", "SL.ranger", "SL.gam", "SL.nnet"))
    
    mu.nonpar[eval.idx] <- predict(mu.fit.nonpar, cbind(A = dat$A, g.A = dat$g.A, X.out)[eval.idx,], onlySL = TRUE)$pred
    
    ## evaluate IF values for each cluster in eval fold
    for(i in which(fold == k)){
      
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

