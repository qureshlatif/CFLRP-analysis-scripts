library(jagsUI)
library(stringr)
library(dplyr)
library(R.utils)
library(ggplot2)
library(cowplot)
library(lme4)

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
  B0.grid[i] <- m$coefficients["(Intercept)"]
  B1.grid[i] <- m$coefficients["PctTrt.z"]
  B2.grid[i] <- m$coefficients["I(PctTrt.z^2)"]
  Y[i, ] <- predict(m, dat.pred, type = "response")
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
           apply(2, function(x) quantile(x,prob=0.975,type=8)) %>% as.numeric) %>%
  filter(!is.na(X))

p.grid <- ggplot(data = dat.SR, aes(x = X, y = Y)) + 
  geom_point(alpha = 0.3) +
  geom_errorbar(aes(ymin = Y.lo, ymax = Y.hi), width = 1, alpha = 0.3) +
  geom_ribbon(data = dat.pred, aes(x = PctTrt.x, ymin = Y.lo, ymax = Y.hi), alpha = 0.2, inherit.aes = F) +
  geom_line(data = dat.pred, aes(x = PctTrt.x, y = Y.md), colour = "blue", size = 1.5) + 
  labs(x= "Percent treated (grid)", y = NULL)

## Point level ##

Trt.b <- Cov[, "Trt_stat"] # Point-level values
YST.b <- Cov[, "Trt_time"] %>%
  (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T))
YST.b <- replace(YST.b, which(is.na(YST.b)), 0) # Point-level values
ID <- Cov[, "gridIndex"] %>% as.factor
SPR <- mod$sims.list$SR.point

# dat.pred = data.frame(Trt = c(0, 1),
#                       YST = c(0, 0))
# Y <- matrix(NA, nrow = dim(SPR)[1], ncol = nrow(dat.pred))
# B0_pnt <- B_Trt <- B_YST <- B0_pntSD <- rep(NA, length = dim(SPR)[1])
# options(warn=1) # Change to warn=2 to investigate convergence warnings.
# for(i in 1:dim(SPR)[1]) {
#   dat <- data.frame(y = SPR[i,],
#                     Trt = Trt.b,
#                     YST = YST.b,
#                     ID = ID)
#   m <- glmer(y ~ Trt + YST + (1|ID), data = dat, family = "poisson")
#   B0_pnt[i] <- summary(m)$coefficients["(Intercept)", "Estimate"]
#   B_Trt[i] <- summary(m)$coefficients["Trt", "Estimate"]
#   B_YST[i] <- summary(m)$coefficients["YST", "Estimate"]
#   B0_pntSD[i] <- as.data.frame(VarCorr(m))$sdcor
#   Y[i, ] <- predict(m, dat.pred, re.form = NA, type = "response")
# }
# B0_pnt <- str_c(median(B0_pnt, na.rm = T) %>% round(digits = 2),
#                 " (",
#                 quantile(B0_pnt, prob = 0.025, type = 8, na.rm = T) %>% round(digits = 2),
#                 ",",
#                 quantile(B0_pnt, prob = 0.975, type = 8, na.rm = T) %>% round(digits = 2),
#                 ")")
# B_Trt <- str_c(median(B_Trt, na.rm = T) %>% round(digits = 2),
#                 " (",
#                 quantile(B_Trt, prob = 0.025, type = 8, na.rm = T) %>% round(digits = 2),
#                 ",",
#                 quantile(B_Trt, prob = 0.975, type = 8, na.rm = T) %>% round(digits = 2),
#                 ")")
# B_YST <- str_c(median(B_YST, na.rm = T) %>% round(digits = 2),
#                " (",
#                quantile(B_YST, prob = 0.025, type = 8, na.rm = T) %>% round(digits = 2),
#                ",",
#                quantile(B_YST, prob = 0.975, type = 8, na.rm = T) %>% round(digits = 2),
#                ")")
# B0_pntSD <- str_c(median(B0_pntSD, na.rm = T) %>% round(digits = 2),
#                " (",
#                quantile(B0_pntSD, prob = 0.025, type = 8, na.rm = T) %>% round(digits = 2),
#                ",",
#                quantile(B0_pntSD, prob = 0.975, type = 8, na.rm = T) %>% round(digits = 2),
#                ")")
# 
# dat.pred <- dat.pred %>%
#   mutate(Y.md = apply(Y, 2, median),
#          Y.lo = apply(Y, 2, function(x) quantile(x, prob = 0.025, type = 8)),
#          Y.hi = apply(Y, 2, function(x) quantile(x, prob = 0.975, type = 8)))
# write.csv(dat.pred, "Spp_richness_dat_pred_cache.csv", row.names = F)
# saveObject(B0_pnt, "B0_pnt_pnt_N_trt_cache")
# saveObject(B_Trt, "B_Trt_pnt_N_trt_cache")
# saveObject(B_YST, "B_YST_pnt_N_trt_cache")
# saveObject(B0_pntSD, "B0_pntSD_pnt_N_trt_cache")
# rm(i, m, Y)

dat.pred <- read.csv("Spp_richness_dat_pred_cache.csv", header = T)
B0_pnt <- loadObject("B0_pnt_pnt_N_trt_cache")
B_Trt <- loadObject("B_Trt_pnt_N_trt_cache")
B_YST <- loadObject("B_YST_pnt_N_trt_cache")
B0_pntSD <- loadObject("B0_pntSD_pnt_N_trt_cache")

dat.SR <- data.frame(X = (Trt.b %>% as.numeric)) %>%
  mutate(#X = X + runif(length(Trt.b), -0.4, 0.4),
    Y = apply(SPR, 2, median) %>% as.numeric,
    Y.lo = apply(SPR, 2,function(x) quantile(x,prob=0.025,type=8)) %>%
      as.numeric,
    Y.hi = apply(SPR, 2, function(x) quantile(x,prob=0.975,type=8)) %>%
      as.numeric) %>%
  filter(!is.na(X))

jitter = runif(length(Trt.b), -0.4, 0.4)
p.point <- ggplot(data = dat.SR, aes(x = X, y = Y)) + 
  geom_point(aes(x = X + jitter), alpha = 0.3) + 
  geom_errorbar(aes(x = X + jitter, ymin = Y.lo, ymax = Y.hi), width = 0, alpha = 0.3) +
  geom_errorbar(data = dat.pred, aes(x = Trt, ymin = Y.lo, ymax = Y.hi), width = 0.3, size = 1, color = "blue", inherit.aes = F) +
  geom_point(data = dat.pred, aes(x = Trt, y = Y.md), size = 3, color = "blue", inherit.aes = F) + 
  scale_x_continuous(breaks = c(0, 1), labels = c("Untreated", "Treated")) +
  labs(x = "Treatment status (point)", y = NULL)

p <- ggdraw() + 
  draw_plot(p.grid, x = 0.03, y = 0, width = 0.47, height = 1) +
  draw_plot(p.point, x = 0.53, y = 0, width = 0.47, height = 1) +
  draw_plot_label(c("N[psi]", "N[theta]"),
                  x = c(0, 0.5), y = c(0.5, 0.5),
                  size = c(20, 20), angle = c(90, 90),
                  hjust = c(0, 0), parse = T)

save_plot("Plot_richness_treatment.tiff", p, ncol = 2, nrow = 1, dpi = 200)

#c(median(mod$sims.list$rho.ab), quantile(mod$sims.list$rho.ab, prob = 0.025, type = 8), quantile(mod$sims.list$rho.ab, prob = 0.975, type = 8))
#c(median(mod$sims.list$rho.bd), quantile(mod$sims.list$rho.bd, prob = 0.025, type = 8), quantile(mod$sims.list$rho.bd, prob = 0.975, type = 8))
