time = proc.time()

###------------------- Load libraries ----------------------###
suppressPackageStartupMessages(require(dplyr))
suppressPackageStartupMessages(require(optparse))


###----------- Read arguments (simulation parameters) ---------------###
option_list = list(
  make_option(c("-m", "--m"), action = "store", default = NA, type = "integer",
              help = paste0("Number of clusters for estimand simulation")),
  make_option(c("-r", "--r"), action = "store", default = NA, type = "integer",
              help = paste0("Number of binary vector sampling for outcome reg computation"))
)

opt = parse_args(OptionParser(option_list = option_list))

m = opt$m              # Number of clusters for target estimand computation
r = opt$r              # Number of binary vector sampling
deltas <- readRDS("deltas.rds")

print("[Simulation setting]")
print(paste0("m: ", m))
print(paste0("deltas: ", paste(signif(deltas, 4), collapse = ", ")))
print(paste0("r: ", r))


###----------- Estimand computation Main function ---------------###

estimand.sim <- function(N, deltas, r){
  
  X1 = rnorm(N, 0, 1)
  X2 = rbinom(N, 1, 0.5)
  C  = rnorm(1, 0, 1)      # Cluster-level covariate

  pi <- plogis(0.1 + 0.2*X1 + 0.2*X2 + 0.1*C)
  
  estimands <- data.frame(delta = deltas,
                          mu = 0, mu_1 = 0, mu_0 = 0)
  
  for(i in seq_len(length(deltas))){
    
    delta = deltas[i]
    
    pi.delta <- delta*pi / (delta*pi + 1 - pi)
    
    # random binary vector from A(N) w.p. P_delta(a_i|X_i), instead of summing all vectors in A(N)
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

###----------- Simulated estimand  ---------------###
# m = 1000
# deltas = c(0.5, 1, 2)
# r = 100

N = sample(x = 5:20, size = m, replace = T)

estimands.list <- lapply(N, estimand.sim, deltas = deltas, r = r)

estimands <- data.frame(delta = rep(0, length(deltas)),
                        mu = 0, mu_1 = 0, mu_0 = 0)

for(i in 1:m){
  estimands = estimands + estimands.list[[i]]
}

estimands = estimands / m

# Standard (delta = 1)
standard = estimands %>% filter(delta == 1)

estimands = estimands %>% mutate(de  = mu_1 - mu_0, 
                     se_1 = mu_1 - standard$mu_1,
                     se_0 = mu_0 - standard$mu_0,
                     oe  = mu   - standard$mu,
                     te  = mu_1 - standard$mu_0)

###-------- Save simulated estimand list as Rdata --------###
## Wrap with as.numeric to coerce from character variable
task_ID <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

## Save output
saveRDS(estimands, file = paste0("Rdata/estimand_id", task_ID,".rds"))


###------------------------------------------------------------###

## Stop timer and report total run time

script.time = proc.time() - time
print(paste0("Total run time was ", script.time[3], " seconds"))