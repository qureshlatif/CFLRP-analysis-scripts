library(jagsUI)
library(stringr)
library(dplyr)
library(R.utils)
library(ggplot2)
library(cowplot)

#setwd("/home/RMBO.LOCAL/quresh.latif/CFLRP")
setwd("C:/Users/Quresh.Latif/files/projects/FS/CFLRP/")
load("Data_compiled_RESQ.RData")

mods <- c("mod_RESQ_treatment_global", "mod_RESQ_treatment_reduced")

Y <- Y.dist
nsites <- length(Y)

cols <- c("WAIC", "DIC")
out <- matrix(NA, nrow = length(mods), ncol = length(cols),
              dimnames = list(mods, cols))

for(m in 1:length(mods)) {
  mod <- loadObject(mods[m])
  nsim <- mod$mcmc.info$n.samples
  ifelse(nsim > 10000, nsamp <- sample(nsim, 10000), nsamp <- 1:nsim)
  out[mods[m], "DIC"] <- mod$DIC
  lambda <- mod$sims.list$lambda[nsamp, ]
  a <- mod$sims.list$a[nsamp, ]
  b <- mod$sims.list$b[nsamp]
  p <- pi <- pic <- array(NA, dim = c(dim(a), nG))
  for(k in 1:nG) {
    for(j in 1:nsites) {
      p[, j, k] <- 1 - exp(-(((breaks[k] + breaks[k + 1]) / 2) / a[,j])^(-1*b))
      pi[, j, k] <- p[, j, k] * area.prop[k]
    }
  }
  pcap <- apply(pi, c(1, 2), sum)
  lpd <- pwaic <- rep(NA, nsites)
  for(i in 1:nsites) {
    lpd[i] <- log(mean(dpois(x = Y[i], lambda = lambda[, i] * pcap[, i])))
    pwaic[i] <- var(log(dpois(x = Y[i], lambda = lambda[, i] * pcap[, i])))
  }
  out[mods[m], "WAIC"] <- -2 * (sum(lpd) - sum(pwaic))
}

out <- as.data.frame(out) %>%
  mutate(mods = mods) %>%
  mutate(dWAIC = WAIC - min(WAIC)) %>%
  mutate(dDIC = DIC - min(DIC)) %>%
  arrange(dWAIC) %>%
  select(mods, WAIC, dWAIC, DIC, dDIC)

write.csv(out, "Model_selection_RESQ.csv", row.names = F)
