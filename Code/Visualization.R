#*******************************************************************************
#**********           < WASH effect on diarrhea incidence             **********
#**********           among children in Senegal >                     **********
#**********           Estimation result visualization                 **********
#**********                                                           **********	
#**********           Written by:				                              **********
#**********           Chanhwa Lee - chanhwa@email.unc.edu             **********
#**********                                                           **********
#**********           Version: 1.0                                    **********
#**********           Sep 19, 2022                                    **********
#*******************************************************************************

############################################
# This file requires "Data/DHS/result.Rdata".
# Causal estimands under the Cluster IPS policy estimation results are visualized.
# Generated figures are saved under "Data/DHS/".
############################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)

### Load estimation result ###
load("Data/DHS/result.Rdata")

### Estimation result ###
est.med = result$est.med
se.med = result$se.med

### Drawing arguments ###
limits = c(1/2,2)
breaks = c(1/2,1,2)
mu.limits.y = c(0.675, 0.85)
de.limits.y = c(-0.05, 0.125)
se1.limits.y = c(-0.025, 0.03)
se0.limits.y = c(-0.06, 0.05)
oe.limits.y = c(-0.05, 0.03)
te.limits.y = c(0, 0.10)

### 95% CI plots ###
plot.est <- list()

for(eff in c("mu", "mu_1", "mu_0", "de", "se_1", "se_0", "oe", "te")){
  
  data.eff = cbind(est.med %>% select(delta, eff) %>% rename(est.med = eff), 
                   se.med %>% select(eff) %>% rename(se.med = eff)) %>%
    mutate(lower.med = est.med - 1.96 * se.med,
           upper.med = est.med + 1.96 * se.med)
  
  plot.est[[eff]] <- 
    ggplot(data = data.eff) +
    geom_errorbar(aes(x = delta, ymin = lower.med, ymax = upper.med), width = 0.01) +
    geom_point(aes(x = delta, y = est.med), col = "red", size = 0.5) +
    ggtitle(eff) +
    xlab(expression(delta)) +
    ylab(NULL)
  
  if(eff %in% c("de", "se_1", "se_0", "oe", "te")){
    plot.est[[eff]] <- plot.est[[eff]] + geom_hline(yintercept = 0, color = "blue", size = 0.5, linetype = "dashed")
  }
  
}

### Arranged plots (all 8 plots) ###
plot_grid(
  plot_grid(
    plot.est[["mu"]] + 
      lims(y = mu.limits.y) + 
      scale_x_continuous(trans = 'log', limits = limits, breaks = breaks) + 
      ggtitle(expression(mu(delta))),
    plot.est[["mu_1"]] + lims(y = mu.limits.y) + scale_x_continuous(trans = 'log', limits = limits, breaks = breaks) + 
      ggtitle(expression(mu[1](delta))),
    plot.est[["mu_0"]] + lims(y = mu.limits.y) + scale_x_continuous(trans = 'log', limits = limits, breaks = breaks) + 
      ggtitle(expression(mu[0](delta))),
    nrow = 1)
  ,
  plot_grid(
    plot.est[["oe"]] + lims(y = oe.limits.y) + scale_x_continuous(trans = 'log', limits = limits, breaks = breaks) + 
      ggtitle(expression(OE(delta, 1))),
    plot.est[["se_1"]] + lims(y = se1.limits.y) + scale_x_continuous(trans = 'log', limits = limits, breaks = breaks) + 
      ggtitle(expression(SE[1](delta, 1))),
    plot.est[["se_0"]] + lims(y = se0.limits.y) + scale_x_continuous(trans = 'log', limits = limits, breaks = breaks) + 
      ggtitle(expression(SE[0](delta, 1))),
    nrow = 1)
  ,
  plot_grid(
    plot.est[["de"]] + lims(y = de.limits.y) + scale_x_continuous(trans = 'log', limits = limits, breaks = breaks) + 
      ggtitle(expression(DE(delta))),
    plot.est[["te"]] + lims(y = te.limits.y) + scale_x_continuous(trans = 'log', limits = limits, breaks = breaks) + 
      ggtitle(expression(TE(delta, 1))),
    nrow = 1)
  ,
  nrow = 3)

ggsave("Data/DHS/95CIs.pdf", width = 10, height = 10)

### Arranged plots (mu, mu1,mu0, oe, se1, se0 plots) ###
plot_grid(
  plot_grid(
    plot.est[["mu"]] + 
      lims(y = mu.limits.y) + 
      scale_x_continuous(trans = 'log', limits = limits, breaks = breaks) + 
      ggtitle(expression(mu(delta))),
    plot.est[["mu_1"]] + lims(y = mu.limits.y) + scale_x_continuous(trans = 'log', limits = limits, breaks = breaks) + 
      ggtitle(expression(mu[1](delta))),
    plot.est[["mu_0"]] + lims(y = mu.limits.y) + scale_x_continuous(trans = 'log', limits = limits, breaks = breaks) + 
      ggtitle(expression(mu[0](delta))),
    nrow = 1)
  ,
  plot_grid(
    plot.est[["oe"]] + lims(y = oe.limits.y) + scale_x_continuous(trans = 'log', limits = limits, breaks = breaks) + 
      ggtitle(expression(OE(delta, 1))),
    plot.est[["se_1"]] + lims(y = se1.limits.y) + scale_x_continuous(trans = 'log', limits = limits, breaks = breaks) + 
      ggtitle(expression(SE[1](delta, 1))),
    plot.est[["se_0"]] + lims(y = se0.limits.y) + scale_x_continuous(trans = 'log', limits = limits, breaks = breaks) + 
      ggtitle(expression(SE[0](delta, 1))),
    nrow = 1)
  ,
  nrow = 2)

ggsave("Data/DHS/95CIs_mu_se.pdf", width = 10, height = 6.5)

### Arranged plots (de and te plots) ###
plot_grid(
  plot.est[["de"]] + lims(y = de.limits.y) + scale_x_continuous(trans = 'log', limits = limits, breaks = breaks) + 
    ggtitle(expression(DE(delta))),
  plot.est[["te"]] + lims(y = te.limits.y) + scale_x_continuous(trans = 'log', limits = limits, breaks = breaks) + 
    ggtitle(expression(TE(delta, 1))),
  nrow = 1)

ggsave("Data/DHS/95CIs_de_te.pdf", width = 6.5, height = 3.5)

### Plots without se ###
est.med.mu = gather(est.med %>% select(delta:mu_0), type, Estimate, mu:mu_0)
est.med.effects = gather(est.med %>% select(delta, de:te), type, Estimate, de:te)

plot_grid(
  ggplot(data = est.med.mu, aes(x = delta, y = Estimate, color = type)) + 
    geom_line() + 
    xlab(expression(delta)) +
    scale_x_continuous(trans = 'log', breaks = c(0.5,1,2))
  ,
  ggplot(data = est.med.effects, aes(x = delta, y = Estimate)) + 
    geom_line(aes(color = type)) + 
    xlab(expression(delta)) +
    scale_x_continuous(trans = 'log', breaks = c(0.5,1,2))
  ,
  nrow = 1)

ggsave("Data/DHS/Effectsplots.pdf", width = 10, height = 4)
