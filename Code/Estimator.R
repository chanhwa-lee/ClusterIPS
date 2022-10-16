#*******************************************************************************
#**********           < WASH effect on diarrhea incidence             **********
#**********           among children in Senegal >                     **********
#**********           Nonparametric efficient sample splitting        **********
#**********           estimator function                              **********	
#**********                                                           **********
#**********           Written by:				                              **********
#**********           Chanhwa Lee - chanhwa@email.unc.edu             **********
#**********                                                           **********
#**********           Version: 1.0                                    **********
#**********           Sep 19, 2022                                    **********
#*******************************************************************************

###--------------------------- Estimator function ---------------------------###

#' Nonparametric efficient sample splitting estimator under 
#' the Incremental Propensity Score policy
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
#' @param SL.library A vector of strings. Names of data-adaptive methods available
#' in `SuperLearner` R package to be used in the ensemble estimation of nuisance
#' functions
#' @return Estimate and standard error estimate of the causal estimands over the
#' delta values
#' @examples 
#' estimator(dat, X.trt, X.out, deltas = seq(0.5,2,0.1), K = 2, r = 100, SL.library = c("SL.ranger"))


estimator <- function(dat, X.trt, X.out, deltas, K, r, SL.library){
  
  ### Add `delta` = 1 if not in `deltas` ###
  if(!(1 %in% deltas)) deltas <- c(deltas, 1)
  
  ### Split dat into lists with respect to cluster id to get cluster sizes ###
  dat.list <- split(dat, f = dat$id)
  
  ### Storage setup ###
  n <- length(dat.list)                                                         ## Number of clusters
  N <- sapply(dat.list, nrow)                                                   ## Cluster sizes
  N.limit <- cumsum(c(0,N))
  
  est <- se <-  data.frame(delta = deltas,                                      ## Estimates and SE estimates
                           mu = 0, mu_1 = 0, mu_0 = 0,
                           de = 0, se_1 = 0, se_0 = 0,
                           oe = 0, te = 0)
  
  pi <- numeric(sum(N))                                                         ## Propensity score estimates for clusters
  mu <- numeric(sum(N))                                                         ## Outcome regression function estimates for clusters
  
  ifvals.mu   <- matrix(0, nrow = length(deltas), ncol = n)                     ## IF for mu for clusters
  ifvals.mu_1 <- matrix(0, nrow = length(deltas), ncol = n)                     ## IF for mu_1 for clusters
  ifvals.mu_0 <- matrix(0, nrow = length(deltas), ncol = n)                     ## IF for mu_0 for clusters
  
  fold <- sample(1:n, replace = F) %% K + 1                                     ## Sample split index for clusters
  foldlong <- rep(fold, times = N)                                              ## Sample split index for units
  
  ### Estimate nuisance functions at groups and evaluate IF at clusters ###
  for (k in 1:K) {
    
    print(paste("   Group:", k))
    
    train.idx = which(foldlong != k)
    eval.idx  = which(foldlong == k)
    
    ### Propensity and Outcome models estimation using other groups' data ###
    ### and evaluate at group k                                           ###
    
    pi.fit = SuperLearner(Y = dat$A[train.idx], 
                          X = X.trt[train.idx,],
                          family = binomial(), 
                          SL.library = SL.library)
    
    pi[eval.idx] <- predict(pi.fit, 
                            X.trt[eval.idx,], 
                            onlySL = TRUE)$pred
    
    mu.fit <- SuperLearner(Y = dat$Y[train.idx], 
                           X = cbind(A = dat$A, g.A = dat$g.A, X.out)[train.idx,],
                           family = binomial(), 
                           SL.library = SL.library)
    
    mu[eval.idx] <- predict(mu.fit, 
                            cbind(A = dat$A, g.A = dat$g.A, X.out)[eval.idx,], 
                            onlySL = TRUE)$pred
    
    ### Evaluate IF values for clusters in group k ###
    for(i in which(fold == k)){
      
      dat.sub.idx <- (N.limit[i]+1):N.limit[i+1]
      dat.sub <- dat[dat.sub.idx,]
      
      pi.sub  <- pi[dat.sub.idx]
      mu.sub  <- mu[dat.sub.idx]
      
      ### Outcome regression ###
      
      if(2^N[i] <= r){
        ###### Summing all vectors in A(N) ######
        a = as.vector(t(expand.grid(replicate(N[i], 0:1, simplify = F))))       ## 2^N[i]*N[i] - vector
        rr = 2^N[i]                                                             ## rr: placeholder for r
      }else{
        ###### Subsampling approximation ######
        a <- sample(c(0,1), r*N[i], replace = T)                                ## r*N[i] - vector
        rr = r                                                                  ## rr: placeholder for r
      }
      
      ###### Aggregating help function ######
      agg.fun <- function(a, fun){
        aggregate(a, list(rep(1:rr, each = N[i])), fun)[,-1]
      }
      
      g.a = (rep(agg.fun(a, sum), each = N[i]) - a) / (N[i]-1)                  ## Leave-one-out within-cluster treatment proportion
      
      new.X.out.a = cbind(data.frame(A = a, g.A = g.a),                         ## New outcome regression model covariates with a and g.a
                          bind_rows(replicate(rr, X.out[dat.sub.idx,], 
                                              simplify = FALSE)))
      
      mu.a <- predict(mu.fit, new.X.out.a, onlySL = TRUE)$pred                  ## Outcomre regression model evaluated at a and g.a
      
      pi.delta.sub <- pi.sub %o% deltas / (pi.sub %o% deltas + 1 - pi.sub)      ## N[i] x length(deltas) matrix
      
      pi.delta.sub.rep <- rep(1, rr) %x% pi.delta.sub                           ## rr*N[i] x length(deltas) matrix (replicate `pi.delta.sub` rr times)
      
      cluster.prob <- agg.fun(a*pi.delta.sub.rep + (1-a)*(1-pi.delta.sub.rep), 
                              prod)                                             ## rr x length(deltas) matrix
      
      add.factor <-
        (2*a-1) %o% deltas * rep(dat.sub$A - pi.sub, rr) / 
        (rep(1, rr) %x% (pi.sub %o% deltas + 1 - pi.sub)) /
        (a*rep(pi.sub, rr) %o% deltas + (1-a)*(1-rep(pi.sub, rr)))              ## rr*N[i] x length(deltas) matrix
      
      ###### Outcome regression part for mu ######
      out.reg.mu <- apply(agg.fun(mu.a, mean) * (1 + agg.fun(add.factor, sum)) * 
                            2^(N[i]) * cluster.prob, 2, mean)                   ## length(deltas) vector
      
      ###### Outcome regression part for mu_1 ######
      out.reg.mu_1 <-
        apply(agg.fun( (mu.a * a)[,rep(1, length(deltas))] / 
                         pi.delta.sub.rep * 
                         (1 + 
                          as.matrix(agg.fun(add.factor, sum)) %x% rep(1,N[i]) - 
                          add.factor), 
                       mean) *
                2^(N[i]) * cluster.prob, 2, mean)                               ## length(deltas) vector
      
      ###### Outcome regression part for mu_0 ######
      out.reg.mu_0 <-
        apply(agg.fun( (mu.a * (1-a))[,rep(1, length(deltas))] / 
                         (1-pi.delta.sub.rep) * 
                         (1 + 
                          as.matrix(agg.fun(add.factor, sum)) %x% rep(1,N[i]) - 
                          add.factor), 
                       mean) *
                2^(N[i]) * cluster.prob, 2, mean)                               ## length(deltas) vector
      
      
      ### Bias correction ###
      cluster.prob.ratio <- deltas^(sum(dat.sub$A)) / 
        apply(pi.sub %o% deltas + 1 - pi.sub, 2, prod)                          ## length(deltas) vector
      
      ###### Bias correction part for mu ######
      bias.cor.mu <- mean(dat.sub$Y - mu.sub, na.rm = T) * cluster.prob.ratio   ## length(deltas) vector
      
      ###### Bias correction part for mu_1 ######
      bias.cor.mu_1 <-
        apply( (dat.sub$Y - mu.sub) * I(dat.sub$A==1) / pi.delta.sub, 2, mean) * 
        cluster.prob.ratio                                                      ## length(deltas) vector
      
      ###### Bias correction part for mu_0 ######
      bias.cor.mu_0 <-
        apply( (dat.sub$Y - mu.sub) * I(dat.sub$A==0) / (1-pi.delta.sub), 2, mean) * 
        cluster.prob.ratio                                                      ## length(deltas) vector
      
      
      ### IF for mu at cluster i ###
      ifvals.mu[,i]   <- out.reg.mu + bias.cor.mu
      
      ## IF for mu_1 at cluster i
      ifvals.mu_1[,i] <- out.reg.mu_1 + bias.cor.mu_1
      
      ## IF for mu_0 at cluster i
      ifvals.mu_0[,i] <- out.reg.mu_0 + bias.cor.mu_0
      
    }
    
  }
  
  ### Estimates and SE estimate helper function ###
  IF.to.est <- function(IF){
    est = mean(aggregate(IF, by = list(fold), mean)[,-1])
    se  = sqrt(1/n * (mean(aggregate(IF, by = list(fold), 
                                     function(phi) mean(phi^2))[,-1]) - est^2))
    return(c(est, se))
  }
  
  delta1.idx = which(deltas == 1)                                               ## Baseline delta == 1
  
  for(delta.idx in seq_len(length(deltas))){
    
    ### Causal estimands estimates and SE estimates for delta
    est[delta.idx, "mu"]   <- IF.to.est(ifvals.mu  [delta.idx,])[1]
    est[delta.idx, "mu_1"] <- IF.to.est(ifvals.mu_1[delta.idx,])[1]
    est[delta.idx, "mu_0"] <- IF.to.est(ifvals.mu_0[delta.idx,])[1]
    est[delta.idx, "de"]   <- IF.to.est(ifvals.mu_1[delta.idx,] - ifvals.mu_0[delta.idx,])[1]
    
    se[delta.idx, "mu"]   <- IF.to.est(ifvals.mu  [delta.idx,])[2]
    se[delta.idx, "mu_1"] <- IF.to.est(ifvals.mu_1[delta.idx,])[2]
    se[delta.idx, "mu_0"] <- IF.to.est(ifvals.mu_0[delta.idx,])[2]
    se[delta.idx, "de"]   <- IF.to.est(ifvals.mu_1[delta.idx,] - ifvals.mu_0[delta.idx,])[2]
    
    ### Causal effects when comparing delta versus delta == 1
    est[delta.idx, "se_1"] <- IF.to.est(ifvals.mu_1[delta.idx,] - ifvals.mu_1[delta1.idx,])[1]
    est[delta.idx, "se_0"] <- IF.to.est(ifvals.mu_0[delta.idx,] - ifvals.mu_0[delta1.idx,])[1]
    est[delta.idx, "oe"]   <- IF.to.est(ifvals.mu  [delta.idx,] - ifvals.mu  [delta1.idx,])[1]
    est[delta.idx, "te"]   <- IF.to.est(ifvals.mu_1[delta.idx,] - ifvals.mu_0[delta1.idx,])[1]
    
    se[delta.idx, "se_1"] <- IF.to.est(ifvals.mu_1[delta.idx,] - ifvals.mu_1[delta1.idx,])[2]
    se[delta.idx, "se_0"] <- IF.to.est(ifvals.mu_0[delta.idx,] - ifvals.mu_0[delta1.idx,])[2]
    se[delta.idx, "oe"]   <- IF.to.est(ifvals.mu  [delta.idx,] - ifvals.mu  [delta1.idx,])[2]
    se[delta.idx, "te"]   <- IF.to.est(ifvals.mu_1[delta.idx,] - ifvals.mu_0[delta1.idx,])[2]
    
  }
  
  return(list(est = est, se = se))
  
}