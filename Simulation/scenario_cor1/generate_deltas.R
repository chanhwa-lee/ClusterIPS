### Generate and save deltas of interest ###
# Better if deltas includes 1 

deltas <- c(0.25, 0.5, 1, 2, 4)

saveRDS(deltas, file = "deltas.rds")