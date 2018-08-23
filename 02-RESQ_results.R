library(jagsUI)
library(stringr)
library(dplyr)
library(R.utils)
library(ggplot2)
library(cowplot)

#setwd("/home/RMBO.LOCAL/quresh.latif/CFLRP")
setwd("C:/Users/Quresh.Latif/files/projects/FS/CFLRP/")
load("Data_compiled_RESQ.RData")

#___ Inputs ___#
mod <- loadObject("mod_RESQ_treatment_global")
pars <- c("beta0.mean", "beta0.sd", "bl.trt", # Parameters of interest
          "bl.YST", "bd.TWIP", "bd.RDens",
          "a0", "a.Time", "a.Time2", "a.DOY",
          "a.DOY2", "a.trt", "a.YST", "b")
tab.out <- "Param_summ_RESQ.csv"
#______________#

# Parameter summary table #
sum.table <- mod$summary %>% tbl_df %>%
  mutate(parameter = dimnames(mod$summary)[[1]]) %>%
  filter(parameter %in% pars) %>%
  mutate(estimate = str_c(round(`50%`, digits = 3), "(", round(`2.5%`, digits = 3), ",", round(`97.5%`, digits = 3), ")")) %>%
  select(parameter, estimate, Rhat, n.eff, overlap0, f) %>%
  mutate(n.eff = n.eff %>% as.integer)

write.csv(sum.table, tab.out, row.names = F)

# Plot outbreak relationships #
Cov <- eval(as.name(str_c("Cov.", stratum)))
x.pctdead <- seq(0, 1, length.out = 20)
#z.pctdead <- (x.pctdead - mean(Cov[, "DeadConif"], na.rm = T)) / sd(Cov[, "DeadConif"], na.rm = T)
x.YSO <- seq(0, max(Cov[, "YSO"], na.rm = T))
#z.YSO <- (x.YSO - mean(Cov[, "YSO"], na.rm = T)) / sd(Cov[, "YSO"], na.rm = T)
plt.tbl <- expand.grid(x.pctdead, x.YSO) %>%
  rename(pctdead = Var1, YSO = Var2) %>%
  mutate(pctdead.z = (pctdead - mean(Cov[, "DeadConif"], na.rm = T)) / sd(Cov[, "DeadConif"], na.rm = T)) %>%
  mutate(YSO.z = (YSO - mean(Cov[, "YSO"], na.rm = T)) / sd(Cov[, "YSO"], na.rm = T)) %>%
  mutate(Outbreak = 1)
plt.tbl <- plt.tbl %>% rbind(c(mean(Cov[-ind, "DeadConif"], na.rm = T),
                               -2,
                               ((mean(Cov[-ind, "DeadConif"], na.rm = T) - mean(Cov[, "DeadConif"], na.rm = T)) / sd(Cov[, "DeadConif"], na.rm = T)),
                               0, 0))

pred.N <- matrix(NA, nrow = nrow(plt.tbl), ncol = mod$mcmc.info$n.samples)
for(j in 1:nrow(pred.N))
  pred.N[j, ] <- exp(mod$sims.list$beta0.mean +
                       mod$sims.list$bl.pdead * plt.tbl$pctdead.z[j] +
                       mod$sims.list$bl.pdead2 * (plt.tbl$pctdead.z[j]^2) +
                       mod$sims.list$bl.outbrk * plt.tbl$Outbreak[j] +
                       mod$sims.list$bl.YSO * plt.tbl$YSO.z[j] +
                       mod$sims.list$bl.YSO2 * (plt.tbl$YSO[j]^2) +
                       mod$sims.list$bl.pdXYSO * plt.tbl$pctdead.z[j] * plt.tbl$YSO.z[j])
plt.tbl$N.md <- apply(pred.N, 1, median)
plt.tbl$N.95lo <- apply(pred.N, 1, function(x) quantile(x, prob = 0.025, type = 8))
plt.tbl$N.95hi <- apply(pred.N, 1, function(x) quantile(x, prob = 0.975, type = 8))
