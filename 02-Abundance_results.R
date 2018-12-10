library(jagsUI)
library(stringr)
library(dplyr)
library(R.utils)
library(ggplot2)
library(cowplot)

#setwd("/home/RMBO.LOCAL/quresh.latif/CFLRP")
setwd("C:/Users/Quresh.Latif/files/projects/FS/CFLRP/")
load("Data_compiled_abundance.RData")

#___ Inputs ___#
spp <- "YRWA"
mod <- loadObject("abund_models/mod_YRWA_abundance_treatment_global")
pars <- c("beta0.mean", "beta0.sd", "bd.PctTrt", "bd.heatload", # Parameters of interest
          "bd.TWI", "bl.trt", "bl.YST", "bd.RDens",
          "a0", "a.Time", "a.Time2", "a.DOY",
          "a.DOY2", "a.trt", "a.YST", "b")
tab.out <- "Param_summ_abundance.csv"
#______________#

samp.acres <- sum(area.band.list[[which(Spp == spp)]]) * 0.000247105

# Parameter summary table #
sum.table <- mod$summary %>% tbl_df %>%
  mutate(parameter = dimnames(mod$summary)[[1]]) %>%
  filter(parameter %in% pars) %>%
  mutate(estimate = str_c(round(`50%`, digits = 3), "(", round(`2.5%`, digits = 3), ",", round(`97.5%`, digits = 3), ")")) %>%
  select(parameter, estimate, Rhat, n.eff, overlap0, f) %>%
  mutate(n.eff = n.eff %>% as.integer)

write.csv(sum.table, tab.out, row.names = F)

#### Plot abundance relationships ####
## Treatment ##
dat.obs <- Cov %>%
  as.data.frame() %>%
  mutate(N.md = apply(mod$sims.list$N, 2, median),
         N.lo = apply(mod$sims.list$N, 2, function(x) quantile(x, prob = 0.025, type = 8)),
         N.hi = apply(mod$sims.list$N, 2, function(x) quantile(x, prob = 0.975, type = 8)),
         Y = Y.dist)

dat.pred <- data.frame(x.trt = c(0, rep(1, 10)), x.YST = c(NA, 1:10)) %>%
  mutate(z.YST = (x.YST - mean(dat.obs$Trt_time, na.rm = T)) / sd(dat.obs$Trt_time, na.rm = T)) %>%
  mutate(x.YST = replace(x.YST, which(x.trt == 0), 0)) %>%
  mutate(z.YST = replace(z.YST, which(x.trt == 0), 0)) %>%
  mutate(pred = NA, pred.lo = NA, pred.hi = NA)

B.trt <- mod$sims.list$bl.trt
B.YST <- mod$sims.list$bl.YST
for(i in 1:nrow(dat.pred)) {
  beta0 <- mod$sims.list$beta0.mean
  lambda <- (exp(beta0 + B.trt*dat.pred$x.trt[i] + B.YST*dat.pred$z.YST[i]) / samp.acres) * 100 # Rescale to number per 100 acres
  dat.pred$pred[i] <- median(lambda)
  dat.pred$pred.lo[i] <- quantile(lambda, prob = 0.025, type = 8)
  dat.pred$pred.hi[i] <- quantile(lambda, prob = 0.975, type = 8)
}

dat.obs <- dat.obs %>%
  mutate(x = replace(Trt_time, which(is.na(Trt_time)), 0)) %>%
  dplyr::group_by(x) %>%
  summarise(N = (mean(N.md) / samp.acres) * 100)
p.trt <- ggplot(dat.pred, aes(x = x.YST, y = pred)) +
  #geom_point(data = dat.obs, aes(x = x, y = N), size = 1) +
  geom_errorbar(data = dat.pred %>% filter(x.trt == 0), aes(x = x.YST, ymin = pred.lo, ymax = pred.hi), size = 1, width = 0.2) +
  geom_point(data = dat.pred %>% filter(x.trt == 0), aes(x = x.YST, y = pred), size = 3) +
  geom_ribbon(data = dat.pred %>% filter(x.trt == 1), aes(x = x.YST, ymin = pred.lo, ymax = pred.hi), alpha = 0.3) +
  geom_line(data = dat.pred %>% filter(x.trt == 1), aes(x = x.YST, y = pred), size = 2) +
  scale_x_continuous(lim = c(-0.3, 10.3), breaks = 0:10, labels = c("pre", 1:10)) +
  ylim(0, 23) +
  xlab("Years since treatment") + ylab(NULL)

## TWI ##
dat.pred <- data.frame(x = seq(min(Cov[, "TWI"]), max(Cov[, "TWI"]), length.out = 10)) %>%
  mutate(z = (x - mean(Cov[, "TWI"], na.rm = T)) / sd(Cov[, "TWI"], na.rm = T)) %>%
  mutate(pred = NA, pred.lo = NA, pred.hi = NA)

B <- mod$sims.list$bd.TWI
for(i in 1:nrow(dat.pred)) {
  beta0 <- mod$sims.list$beta0.mean + B*dat.pred$z[i]
  lambda <- (exp(beta0) / samp.acres) * 100 # Rescale to number per 100 acres
  dat.pred$pred[i] <- median(lambda)
  dat.pred$pred.lo[i] <- quantile(lambda, prob = 0.025, type = 8)
  dat.pred$pred.hi[i] <- quantile(lambda, prob = 0.975, type = 8)
}

p.TWI <- ggplot(dat.pred, aes(x = x, y = pred)) +
  geom_ribbon(data = dat.pred, aes(x = x, ymin = pred.lo, ymax = pred.hi), alpha = 0.3) +
  geom_line(data = dat.pred, aes(x = x, y = pred), size = 2) +
  ylim(0, 23) +
  xlab("Topographic wetness index") + ylab(NULL)

p <- ggdraw() +
  draw_plot(p.trt, x = 0.05, y = 0, width = 0.475, height = 1) +
  draw_plot(p.TWI, x = 0.525, y = 0, width = 0.475, height = 1) +
  draw_plot_label("RESQ density per 100 acres", x = 0.02, y = 0.15, angle = 90, hjust = 0)

save_plot("figure_RESQ_abundance.tiff", p, ncol = 2, nrow = 1, dpi = 200)


#### Plot detection relationships??? ####
## Treatment ##
