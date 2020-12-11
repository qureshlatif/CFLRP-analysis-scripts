library(jagsUI)
library(stringr)
library(dplyr)
library(R.utils)
library(ggplot2)
library(cowplot)
library(lme4)
library(QSLpersonal)
library(MASS)
theme_set(theme_cowplot())

setwd("C:/Users/Quresh.Latif/files/projects/FS/CFLRP")
load("Data_compiled.RData")

mod <- loadObject("mod_treatment_d0yr")

#### Plot species richness ####

## Grid level ##
gridID <- Cov[, "gridIndex"]
yearID <- Cov[, "YearInd"]
PctTrt.d <- matrix(NA, nrow = max(gridID), ncol = max(yearID))
PctTrt.d[landscape_data %>% filter(YearInd == 1) %>% pull(gridIndex), 1] <-
  landscape_data %>% filter(YearInd == 1) %>% pull(PctTrt)
PctTrt.d[landscape_data %>% filter(YearInd == 2) %>% pull(gridIndex), 2] <-
  landscape_data %>% filter(YearInd == 2) %>% pull(PctTrt)
PctTrt.d[landscape_data %>% filter(YearInd == 3) %>% pull(gridIndex), 3] <-
  landscape_data %>% filter(YearInd == 3) %>% pull(PctTrt)
dat.pred <- data.frame(PctTrt.x = seq(min(PctTrt.d, na.rm = T), max(PctTrt.d, na.rm = T), length.out = 20)) %>%
  mutate(PctTrt.z = (PctTrt.x - mean(PctTrt.d, na.rm = T)) / sd(PctTrt.d, na.rm = T))
dat.x <- data.frame(PctTrt.z = PctTrt.d %>%
                      (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)) %>%
                      apply(1, mean, na.rm = T))

SPR <- mod$sims.list$SR.grid
  # Derive posterior samples for plotting spp richness trend #
Y <- matrix(NA, nrow = dim(SPR)[1], ncol = nrow(dat.pred))
B0.grid <- B1.grid <- B2.grid <- r <- numeric(length = dim(SPR)[1])
for(i in 1:dim(SPR)[1]) {
  dat = dat.x %>% mutate(N = SPR[i,,] %>% apply(1, mean))
  m <- suppressWarnings(glm(N ~ PctTrt.z + I(PctTrt.z^2), data = dat, family = "poisson"))
  mn <- m$coefficients
  vc <- vcov(m)
  cfs <- mvrnorm(1, mn, vc)
  B0.grid[i] <- cfs["(Intercept)"]
  B1.grid[i] <- cfs["PctTrt.z"]
  B2.grid[i] <- cfs["I(PctTrt.z^2)"]
  X <- cbind(1, dat.pred$PctTrt.z, dat.pred$PctTrt.z^2) %>% as.matrix
  Y[i, ] <- exp(X %*% cfs)
  #Y[i, ] <- predict(m, dat.pred, type = "response")
}
dat.pred <- dat.pred %>%
  mutate(Y.md = apply(Y, 2, median),
         Y.lo = apply(Y, 2, function(x) quantile(x, prob = 0.025, type = 8)),
         Y.hi = apply(Y, 2, function(x) quantile(x, prob = 0.975, type = 8)))

B0.grid <- str_c(median(B0.grid) %>% round(digits = 2),
                 " (",
                 quantile(B0.grid, prob = 0.025, type = 8) %>% round(digits = 2),
                 ",",
                 quantile(B0.grid, prob = 0.975, type = 8) %>% round(digits = 2),
                 ")")
B1.grid <- str_c(median(B1.grid) %>% round(digits = 2),
                 " (",
                 quantile(B1.grid, prob = 0.025, type = 8) %>% round(digits = 2),
                 ",",
                 quantile(B1.grid, prob = 0.975, type = 8) %>% round(digits = 2),
                 ")")
B2.grid <- str_c(median(B2.grid) %>% round(digits = 2),
                 " (",
                 quantile(B2.grid, prob = 0.025, type = 8) %>% round(digits = 2),
                 ",",
                 quantile(B2.grid, prob = 0.975, type = 8) %>% round(digits = 2),
                 ")")

dat.SR <- data.frame(X = (PctTrt.d %>% apply(1, mean, na.rm = T))) %>%
  mutate(Y = SPR %>% apply(c(1, 2), mean) %>% apply(2, median) %>% as.numeric,
         Y.lo = SPR %>% apply(c(1, 2), mean) %>%
           apply(2, function(x) quantile(x,prob=0.025,type=8)) %>% as.numeric,
         Y.hi = SPR %>% apply(c(1, 2), mean) %>%
           apply(2, function(x) quantile(x,prob=0.975,type=8)) %>% as.numeric)

p.ptrt <- ggplot(data = dat.SR, aes(x = X, y = Y)) + 
  geom_point(alpha = 0.3) +
  geom_errorbar(aes(ymin = Y.lo, ymax = Y.hi), width = 1, alpha = 0.3) +
  geom_ribbon(data = dat.pred, aes(x = PctTrt.x, ymin = Y.lo, ymax = Y.hi), alpha = 0.2, inherit.aes = F) +
  geom_line(data = dat.pred, aes(x = PctTrt.x, y = Y.md), colour = "blue", size = 1.5) + 
  labs(x= "Percent treated", y = "Species richness (grid)")

## Point level ##

Trt.b <- Cov[, "Trt_stat"] # Point-level values
YST.b <- Cov[, "Trt_time"] # Point-level values
SPR <- mod$sims.list$SR.point
dat.SR <- data.frame(Trt = (Trt.b %>% as.numeric),
                     YST = YST.b) %>%
  mutate(#X = X + runif(length(Trt.b), -0.4, 0.4),
    Y = apply(SPR, 2, median) %>% as.numeric,
    Y.lo = apply(SPR, 2,function(x) quantile(x,prob=0.025,type=8)) %>%
      as.numeric,
    Y.hi = apply(SPR, 2, function(x) quantile(x,prob=0.975,type=8)) %>%
      as.numeric)
psi <- expit(apply(mod$sims.list$d0, c(1, 2), mean))
b0 <- mod$sims.list$b0
omega <- mod$sims.list$omega

dat.pred = data.frame(x = c(0, 1))
b1 <- mod$sims.list$bb.trt
SR.pred <- matrix(NA, nrow = dim(b0)[1], ncol = nrow(dat.pred))
for(i in 1:dim(SR.pred)[2]) {
  theta <- expit(b0 + b1*dat.pred$x[i])
  SR.pred[, i] <- apply(psi * theta * omega, 1, sum)
}
dat.pred <- dat.pred %>%
  mutate(pred.md = apply(SR.pred, 2, median),
         pred.lo = apply(SR.pred, 2, function(x) quantile(x, prob = 0.025, type = 8)),
         pred.hi = apply(SR.pred, 2, function(x) quantile(x, prob = 0.975, type = 8)))
jitter = runif(length(Trt.b), -0.4, 0.4)
p.trt <- ggplot(data = dat.SR, aes(x = Trt, y = Y)) + 
  geom_point(aes(x = Trt + jitter), alpha = 0.1) + 
  geom_errorbar(aes(x = Trt + jitter, ymin = Y.lo, ymax = Y.hi), width = 0, alpha = 0.1) +
  geom_errorbar(data = dat.pred, aes(x = x, ymin = pred.lo, ymax = pred.hi), width = 0.3, size = 1, color = "blue", inherit.aes = F) +
  geom_point(data = dat.pred, aes(x = x, y = pred.md), size = 3, color = "blue", inherit.aes = F) + 
  scale_x_continuous(breaks = c(0, 1), labels = c("Untreated", "Treated")) +
  labs(x = "Treatment status", y = "Species richness (point)")

dat.pred = data.frame(x = seq(min(YST.b, na.rm = T), max(YST.b, na.rm = T), by = 1)) %>%
  mutate(z = (x - mean(YST.b, na.rm = T)) / sd(YST.b, na.rm = T))
b1 <- mod$sims.list$bb.trt
b2 <- mod$sims.list$bb.YST
SR.pred <- matrix(NA, nrow = dim(b0)[1], ncol = nrow(dat.pred))
for(i in 1:dim(SR.pred)[2]) {
  theta <- expit(b0 + b1 + b2*dat.pred$z[i])
  SR.pred[, i] <- apply(psi * theta * omega, 1, sum)
}
dat.pred <- dat.pred %>%
  mutate(pred.md = apply(SR.pred, 2, median),
         pred.lo = apply(SR.pred, 2, function(x) quantile(x, prob = 0.025, type = 8)),
         pred.hi = apply(SR.pred, 2, function(x) quantile(x, prob = 0.975, type = 8)))
p.yst <- ggplot(data = dat.SR %>% filter(!is.na(YST)), aes(x = YST, y = Y)) + 
  geom_point(aes(x = YST), alpha = 0.1) + 
  geom_errorbar(aes(x = YST, ymin = Y.lo, ymax = Y.hi), width = 0, alpha = 0.1) +
  geom_ribbon(data = dat.pred, aes(x = x, ymin = pred.lo, ymax = pred.hi), alpha = 0.3, inherit.aes = F) +
  geom_line(data = dat.pred, aes(x = x, y = pred.md), size = 1, color = "blue", inherit.aes = F) +
  scale_x_continuous(breaks = 1:10, labels = 1:10) +
  labs(x = "Year since treatment", y = NULL)

p <- ggdraw() + 
  draw_plot(p.ptrt, x = 0.0, y = 0.5, width = 0.525, height = 0.5) +
  draw_plot(p.trt, x = 0, y = 0, width = 0.525, height = 0.5) +
  draw_plot(p.yst, x = 0.5, y = 0, width = 0.475, height = 0.5)

#save_plot("Plot_richness_treatment.tiff", p, ncol = 2, nrow = 2, dpi = 200)
save_plot("manuscript/Figure4.tiff", p, ncol = 1.5, nrow = 2.5, dpi = 600)

#c(median(mod$sims.list$rho.ab), quantile(mod$sims.list$rho.ab, prob = 0.025, type = 8), quantile(mod$sims.list$rho.ab, prob = 0.975, type = 8))
#c(median(mod$sims.list$rho.bd), quantile(mod$sims.list$rho.bd, prob = 0.025, type = 8), quantile(mod$sims.list$rho.bd, prob = 0.975, type = 8))
