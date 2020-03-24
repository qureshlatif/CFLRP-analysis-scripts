library(jagsUI)
library(stringr)
library(dplyr)
library(R.utils)
library(QSLpersonal)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())

setwd("C:/Users/Quresh.Latif/files/projects/FS/CFLRP")
load("Data_compiled.RData")

mod <- loadObject("mod_habitat_d0yr_reduced")

#### Plot species richness ####

## Grid level ##
gridID <- Cov[, "gridIndex"]
yearID <- Cov[, "YearInd"]
SPR <- mod$sims.list$SR.grid

  # Extent of canopy gaps (PACCGap) #
PACC10_3km.d <- matrix(NA, nrow = max(gridID), ncol = max(yearID))
PACC10_3km.d[landscape_data %>% filter(YearInd == 1) %>% pull(gridIndex), 1] <-
  landscape_data %>% filter(YearInd == 1) %>% pull(PACC10_3km)
PACC10_3km.d[landscape_data %>% filter(YearInd == 2) %>% pull(gridIndex), 2] <-
  landscape_data %>% filter(YearInd == 2) %>% pull(PACC10_3km)
PACC10_3km.d[landscape_data %>% filter(YearInd == 3) %>% pull(gridIndex), 3] <-
  landscape_data %>% filter(YearInd == 3) %>% pull(PACC10_3km)

PACC40_3km.d <- matrix(NA, nrow = max(gridID), ncol = max(yearID))
PACC40_3km.d[landscape_data %>% filter(YearInd == 1) %>% pull(gridIndex), 1] <-
  landscape_data %>% filter(YearInd == 1) %>% pull(PACC40_3km)
PACC40_3km.d[landscape_data %>% filter(YearInd == 2) %>% pull(gridIndex), 2] <-
  landscape_data %>% filter(YearInd == 2) %>% pull(PACC40_3km)
PACC40_3km.d[landscape_data %>% filter(YearInd == 3) %>% pull(gridIndex), 3] <-
  landscape_data %>% filter(YearInd == 3) %>% pull(PACC40_3km)

mnPerArRatio_Opn3km.d <- matrix(NA, nrow = max(gridID), ncol = max(yearID))
mnPerArRatio_Opn3km.d[landscape_data %>% filter(YearInd == 1) %>% pull(gridIndex), 1] <-
  landscape_data %>% filter(YearInd == 1) %>% pull(mnPerArRatio_Opn3km)
mnPerArRatio_Opn3km.d[landscape_data %>% filter(YearInd == 2) %>% pull(gridIndex), 2] <-
  landscape_data %>% filter(YearInd == 2) %>% pull(mnPerArRatio_Opn3km)
mnPerArRatio_Opn3km.d[landscape_data %>% filter(YearInd == 3) %>% pull(gridIndex), 3] <-
  landscape_data %>% filter(YearInd == 3) %>% pull(mnPerArRatio_Opn3km)

dat.pred.PACC10 <- data.frame(x = seq(min(PACC10_3km.d, na.rm = T), max(PACC10_3km.d, na.rm = T), length.out = 20)) %>%
  mutate(z = (x - mean(PACC10_3km.d, na.rm = T)) / sd(PACC10_3km.d, na.rm = T))
dat.pred.PACC40 <- data.frame(x = seq(min(PACC40_3km.d, na.rm = T), max(PACC40_3km.d, na.rm = T), length.out = 20)) %>%
  mutate(z = (x - mean(PACC40_3km.d, na.rm = T)) / sd(PACC40_3km.d, na.rm = T))
dat.pred.mnPAR <- data.frame(x = seq(min(mnPerArRatio_Opn3km.d, na.rm = T), max(mnPerArRatio_Opn3km.d, na.rm = T), length.out = 20)) %>%
  mutate(z = (x - mean(mnPerArRatio_Opn3km.d, na.rm = T)) / sd(mnPerArRatio_Opn3km.d, na.rm = T))
d0 <- apply(mod$sims.list$d0, c(1, 2), mean)

d1 <- mod$sims.list$bd.PACC10_3km
Y <- matrix(NA, nrow = dim(d0)[1], ncol = nrow(dat.pred.PACC10))
for(i in 1:dim(Y)[2]) {
  psi <- expit(d0 + d1*dat.pred.PACC10$z[i])
  Y[, i] <- apply(psi, 1, sum)
}
dat.pred.PACC10 <- dat.pred.PACC10 %>%
  mutate(Y.md = apply(Y, 2, median),
         Y.lo = apply(Y, 2, function(x) quantile(x, prob = 0.025, type = 8)),
         Y.hi = apply(Y, 2, function(x) quantile(x, prob = 0.975, type = 8)))

d1 <- mod$sims.list$bd.PACC40_3km
Y <- matrix(NA, nrow = dim(d0)[1], ncol = nrow(dat.pred.PACC40))
for(i in 1:dim(Y)[2]) {
  psi <- expit(d0 + d1*dat.pred.PACC40$z[i])
  Y[, i] <- apply(psi, 1, sum)
}
dat.pred.PACC40 <- dat.pred.PACC40 %>%
  mutate(Y.md = apply(Y, 2, median),
         Y.lo = apply(Y, 2, function(x) quantile(x, prob = 0.025, type = 8)),
         Y.hi = apply(Y, 2, function(x) quantile(x, prob = 0.975, type = 8)))

d1 <- mod$sims.list$bd.mnPerArRatio_Opn3km
Y <- matrix(NA, nrow = dim(d0)[1], ncol = nrow(dat.pred.mnPAR))
for(i in 1:dim(Y)[2]) {
  psi <- expit(d0 + d1*dat.pred.mnPAR$z[i])
  Y[, i] <- apply(psi, 1, sum)
}
dat.pred.mnPAR <- dat.pred.mnPAR %>%
  mutate(Y.md = apply(Y, 2, median),
         Y.lo = apply(Y, 2, function(x) quantile(x, prob = 0.025, type = 8)),
         Y.hi = apply(Y, 2, function(x) quantile(x, prob = 0.975, type = 8)))

dat.SR <- data.frame(PACC10 = (PACC10_3km.d %>% apply(1, mean, na.rm = T)),
                     PACC40 = (PACC40_3km.d %>% apply(1, mean, na.rm = T)),
                     mnPAR = (mnPerArRatio_Opn3km.d %>% apply(1, mean, na.rm = T)),
                     Y = SPR %>% apply(c(1, 2), mean) %>% apply(2, median) %>% as.numeric,
                     Y.lo = SPR %>% apply(c(1, 2), mean) %>%
                       apply(2, function(x) quantile(x,prob=0.025,type=8)) %>% as.numeric,
                     Y.hi = SPR %>% apply(c(1, 2), mean) %>%
                       apply(2, function(x) quantile(x,prob=0.975,type=8)) %>% as.numeric)

p.PACCGap <- ggplot(data = dat.SR, aes(x = PACC10, y = Y)) +
  geom_point(alpha = 0.3) +
  geom_errorbar(aes(ymin = Y.lo, ymax = Y.hi), width = 0, alpha = 0.3) +
  geom_ribbon(data = dat.pred.PACC10, aes(x = x, ymin = Y.lo, ymax = Y.hi), alpha = 0.2, inherit.aes = F) +
  geom_line(data = dat.pred.PACC10, aes(x = x, y = Y.md), colour = "blue", size = 1.5) +
  labs(x = "Extent of canopy gaps (PACCGap)", y = NULL)

p.PACCOpn <- ggplot(data = dat.SR, aes(x = PACC40, y = Y)) + 
  geom_point(alpha = 0.3) +
  geom_errorbar(aes(ymin = Y.lo, ymax = Y.hi), width = 0, alpha = 0.3) +
  geom_ribbon(data = dat.pred.PACC40, aes(x = x, ymin = Y.lo, ymax = Y.hi), alpha = 0.2, inherit.aes = F) +
  geom_line(data = dat.pred.PACC40, aes(x = x, y = Y.md), colour = "blue", size = 1.5) + 
  labs(x = "Extent of open forest (PACCOpn)", y = NULL)

p.PAROpn <- ggplot(data = dat.SR, aes(x = mnPAR, y = Y)) +
  geom_point(alpha = 0.3) +
  geom_errorbar(aes(ymin = Y.lo, ymax = Y.hi), width = 0, alpha = 0.3) +
  geom_ribbon(data = dat.pred.mnPAR, aes(x = x, ymin = Y.lo, ymax = Y.hi), alpha = 0.2, inherit.aes = F) +
  geom_line(data = dat.pred.mnPAR, aes(x = x, y = Y.md), colour = "blue", size = 1.5) +
  labs(x = "Open forest perimeter-area ratio (PAROpn)", y = NULL)

# p <- ggdraw() + 
#   draw_plot(p.PACCGap, x = 0.05, y = 0, width = 0.485, height = 1) +
#   draw_plot(p.PACCOpn, x = 0.525, y = 0, width = 0.485, height = 1) +
#   # draw_plot(p.PAROpn, x = 0.6766667, y = 0, width = 0.3233333, height = 1) +
#   draw_plot_label(c("N[psi]", "',k'"), x = c(0, 0.027),
#                   y = c(0.5, 0.58), size = c(20, 12),
#                   angle = c(90, 90), hjust = c(0, 0), parse = T) #****Need to figure out how to get k in there.
# 
# save_plot("Plot_richness_landscape.tiff", p, ncol = 2, nrow = 1, dpi = 200)

## Point level ##

SPR <- mod$sims.list$SR.point

# Canopy cover #
CanCov.b <- Cov[, "CanCov"]
CanHt.b <- Cov[, "CanHt"]
NSnag.b <- Cov[, "NumSnags"]
PIPO.b <- Cov[, "RCOV_PP"]
PSME.b <- Cov[, "RCOV_DF"]
POTR5.b <- Cov[, "RCOV_AS"]
ShrbVol.b <- Cov[, "ShrubVol"]
LadFuel.b <- Cov[, "RSCV_Ladder"]
Herb.b <- Cov[, "HerbGrassVol"]
ID <- Cov[, "gridIndex"] %>% as.factor

dat.SR <- data.frame(CanCov = CanCov.b,
                     CanHt = CanHt.b,
                     NSnag = NSnag.b,
                     PIPO = PIPO.b,
                     PSME = PSME.b,
                     POTR5 = POTR5.b,
                     ShrbVol = ShrbVol.b,
                     LadFuel = LadFuel.b,
                     Herb = Herb.b) %>%
  mutate(Y = apply(SPR, 2, median),
         Y.lo = apply(SPR, 2,function(x) quantile(x,prob=0.025,type=8)),
         Y.hi = apply(SPR, 2, function(x) quantile(x,prob=0.975,type=8))) %>%
  filter(!(is.na(CanCov) | is.na(CanHt) | is.na(NSnag) | is.na(PIPO) | is.na(PSME) |
             is.na(POTR5) | is.na(ShrbVol) | is.na(LadFuel) | is.na(Herb)))

psi <- expit(apply(mod$sims.list$d0, c(1, 2), mean))
b0 <- mod$sims.list$b0

v <- CanCov.b
vtit <- "Canopy cover (CanCov)"
b1 <- mod$sims.list$bb.CanCov
dat.pred <- data.frame(x = seq(min(v, na.rm = T), max(v, na.rm = T), length.out = 20)) %>%
  mutate(z = (x - mean(v, na.rm = T)) / sd(v, na.rm = T))
Y <- matrix(NA, nrow = dim(b0)[1], ncol = nrow(dat.pred))
for(i in 1:dim(Y)[2]) {
  theta <- expit(b0 + b1*dat.pred$z[i])
  Y[, i] <- apply(psi*theta, 1, sum)
}
dat.pred <- dat.pred %>%
  mutate(Y.md = apply(Y, 2, median),
         Y.lo = apply(Y, 2, function(x) quantile(x, prob = 0.025, type = 8)),
         Y.hi = apply(Y, 2, function(x) quantile(x, prob = 0.975, type = 8)))
p.CanCov <- ggplot(data = dat.SR, aes(x = CanCov, y = Y)) + 
  geom_point(alpha = 0.1) +
  geom_errorbar(aes(ymin = Y.lo, ymax = Y.hi), width = 0, alpha = 0.1) +
  geom_ribbon(data = dat.pred, aes(x = x, ymin = Y.lo, ymax = Y.hi), alpha = 0.2, inherit.aes = F) +
  geom_line(data = dat.pred, aes(x = x, y = Y.md), colour = "blue", size = 1.5, inherit.aes = F) + 
  labs(x = vtit, y = NULL)

v <- CanHt.b
vtit <- "Canopy height (CanHt)"
b1 <- mod$sims.list$bb.CanHt
dat.pred <- data.frame(x = seq(min(v, na.rm = T), max(v, na.rm = T), length.out = 20)) %>%
  mutate(z = (x - mean(v, na.rm = T)) / sd(v, na.rm = T))
Y <- matrix(NA, nrow = dim(b0)[1], ncol = nrow(dat.pred))
for(i in 1:dim(Y)[2]) {
  theta <- expit(b0 + b1*dat.pred$z[i])
  Y[, i] <- apply(psi*theta, 1, sum)
}
dat.pred <- dat.pred %>%
  mutate(Y.md = apply(Y, 2, median),
         Y.lo = apply(Y, 2, function(x) quantile(x, prob = 0.025, type = 8)),
         Y.hi = apply(Y, 2, function(x) quantile(x, prob = 0.975, type = 8)))
p.CanHt <- ggplot(data = dat.SR, aes(x = CanHt, y = Y)) +
  geom_point(alpha = 0.1) +
  geom_errorbar(aes(ymin = Y.lo, ymax = Y.hi), width = 0, alpha = 0.1) +
  geom_ribbon(data = dat.pred, aes(x = x, ymin = Y.lo, ymax = Y.hi), alpha = 0.2, inherit.aes = F) +
  geom_line(data = dat.pred, aes(x = x, y = Y.md), colour = "blue", size = 1.5, inherit.aes = F) +
  labs(x = vtit, y = NULL)

v <- NSnag.b
vtit <- "Number of snags (NSnag)"
b1 <- mod$sims.list$bb.NumSnags
dat.pred <- data.frame(x = seq(min(v, na.rm = T), max(v, na.rm = T), length.out = 20)) %>%
  mutate(z = (x - mean(v, na.rm = T)) / sd(v, na.rm = T))
Y <- matrix(NA, nrow = dim(b0)[1], ncol = nrow(dat.pred))
for(i in 1:dim(Y)[2]) {
  theta <- expit(b0 + b1*dat.pred$z[i])
  Y[, i] <- apply(psi*theta, 1, sum)
}
dat.pred <- dat.pred %>%
  mutate(Y.md = apply(Y, 2, median),
         Y.lo = apply(Y, 2, function(x) quantile(x, prob = 0.025, type = 8)),
         Y.hi = apply(Y, 2, function(x) quantile(x, prob = 0.975, type = 8)))
p.NSnag <- ggplot(data = dat.SR, aes(x = NSnag, y = Y)) + 
  geom_point(alpha = 0.1) +
  geom_errorbar(aes(ymin = Y.lo, ymax = Y.hi), width = 0, alpha = 0.1) +
  geom_ribbon(data = dat.pred, aes(x = x, ymin = Y.lo, ymax = Y.hi), alpha = 0.2, inherit.aes = F) +
  geom_line(data = dat.pred, aes(x = x, y = Y.md), colour = "blue", size = 1.5, inherit.aes = F) + 
  labs(x = vtit, y = NULL)

v <- PIPO.b
vtit <- "Ponderosa pine canopy (PIPO)"
b1 <- mod$sims.list$bb.RCOV_PP
dat.pred <- data.frame(x = seq(min(v, na.rm = T), max(v, na.rm = T), length.out = 20)) %>%
  mutate(z = (x - mean(v, na.rm = T)) / sd(v, na.rm = T))
Y <- matrix(NA, nrow = dim(b0)[1], ncol = nrow(dat.pred))
for(i in 1:dim(Y)[2]) {
  theta <- expit(b0 + b1*dat.pred$z[i])
  Y[, i] <- apply(psi*theta, 1, sum)
}
dat.pred <- dat.pred %>%
  mutate(Y.md = apply(Y, 2, median),
         Y.lo = apply(Y, 2, function(x) quantile(x, prob = 0.025, type = 8)),
         Y.hi = apply(Y, 2, function(x) quantile(x, prob = 0.975, type = 8)))
p.PIPO <- ggplot(data = dat.SR, aes(x = PIPO, y = Y)) +
  geom_point(alpha = 0.1) +
  geom_errorbar(aes(ymin = Y.lo, ymax = Y.hi), width = 0, alpha = 0.1) +
  geom_ribbon(data = dat.pred, aes(x = x, ymin = Y.lo, ymax = Y.hi), alpha = 0.2, inherit.aes = F) +
  geom_line(data = dat.pred, aes(x = x, y = Y.md), colour = "blue", size = 1.5, inherit.aes = F) +
  labs(x = vtit, y = NULL)

v <- PSME.b
vtit <- "Douglas fir canopy (PSME)"
b1 <- mod$sims.list$bb.RCOV_DF
dat.pred <- data.frame(x = seq(min(v, na.rm = T), max(v, na.rm = T), length.out = 20)) %>%
  mutate(z = (x - mean(v, na.rm = T)) / sd(v, na.rm = T))
Y <- matrix(NA, nrow = dim(b0)[1], ncol = nrow(dat.pred))
for(i in 1:dim(Y)[2]) {
  theta <- expit(b0 + b1*dat.pred$z[i])
  Y[, i] <- apply(psi*theta, 1, sum)
}
dat.pred <- dat.pred %>%
  mutate(Y.md = apply(Y, 2, median),
         Y.lo = apply(Y, 2, function(x) quantile(x, prob = 0.025, type = 8)),
         Y.hi = apply(Y, 2, function(x) quantile(x, prob = 0.975, type = 8)))
p.PSME <- ggplot(data = dat.SR, aes(x = PSME, y = Y)) + 
  geom_point(alpha = 0.1) +
  geom_errorbar(aes(ymin = Y.lo, ymax = Y.hi), width = 0, alpha = 0.1) +
  geom_ribbon(data = dat.pred, aes(x = x, ymin = Y.lo, ymax = Y.hi), alpha = 0.2, inherit.aes = F) +
  geom_line(data = dat.pred, aes(x = x, y = Y.md), colour = "blue", size = 1.5, inherit.aes = F) + 
  labs(x = vtit, y = NULL)

v <- POTR5.b
vtit <- "Aspen canopy (POTR5)"
b1 <- mod$sims.list$bb.RCOV_AS
dat.pred <- data.frame(x = seq(min(v, na.rm = T), max(v, na.rm = T), length.out = 20)) %>%
  mutate(z = (x - mean(v, na.rm = T)) / sd(v, na.rm = T))
Y <- matrix(NA, nrow = dim(b0)[1], ncol = nrow(dat.pred))
for(i in 1:dim(Y)[2]) {
  theta <- expit(b0 + b1*dat.pred$z[i])
  Y[, i] <- apply(psi*theta, 1, sum)
}
dat.pred <- dat.pred %>%
  mutate(Y.md = apply(Y, 2, median),
         Y.lo = apply(Y, 2, function(x) quantile(x, prob = 0.025, type = 8)),
         Y.hi = apply(Y, 2, function(x) quantile(x, prob = 0.975, type = 8)))
p.POTR5 <- ggplot(data = dat.SR, aes(x = POTR5, y = Y)) + 
  geom_point(alpha = 0.1) +
  geom_errorbar(aes(ymin = Y.lo, ymax = Y.hi), width = 0, alpha = 0.1) +
  geom_ribbon(data = dat.pred, aes(x = x, ymin = Y.lo, ymax = Y.hi), alpha = 0.2, inherit.aes = F) +
  geom_line(data = dat.pred, aes(x = x, y = Y.md), colour = "blue", size = 1.5, inherit.aes = F) + 
  labs(x = vtit, y = NULL)

v <- ShrbVol.b
vtit <- "Shrub-sapling volume (ShrbVol)"
b1 <- mod$sims.list$bb.shvol
dat.pred <- data.frame(x = seq(min(v, na.rm = T), max(v, na.rm = T), length.out = 20)) %>%
  mutate(z = (x - mean(v, na.rm = T)) / sd(v, na.rm = T))
Y <- matrix(NA, nrow = dim(b0)[1], ncol = nrow(dat.pred))
for(i in 1:dim(Y)[2]) {
  theta <- expit(b0 + b1*dat.pred$z[i])
  Y[, i] <- apply(psi*theta, 1, sum)
}
dat.pred <- dat.pred %>%
  mutate(Y.md = apply(Y, 2, median),
         Y.lo = apply(Y, 2, function(x) quantile(x, prob = 0.025, type = 8)),
         Y.hi = apply(Y, 2, function(x) quantile(x, prob = 0.975, type = 8)))
p.ShrbVol <- ggplot(data = dat.SR, aes(x = ShrbVol, y = Y)) + 
  geom_point(alpha = 0.1) +
  geom_errorbar(aes(ymin = Y.lo, ymax = Y.hi), width = 0, alpha = 0.1) +
  geom_ribbon(data = dat.pred, aes(x = x, ymin = Y.lo, ymax = Y.hi), alpha = 0.2, inherit.aes = F) +
  geom_line(data = dat.pred, aes(x = x, y = Y.md), colour = "blue", size = 1.5, inherit.aes = F) + 
  labs(x = vtit, y = NULL)

v <- LadFuel.b
vtit <- "Ladder fuels (LadFuel)"
b1 <- mod$sims.list$bb.RSCV_Ladder
dat.pred <- data.frame(x = seq(min(v, na.rm = T), max(v, na.rm = T), length.out = 20)) %>%
  mutate(z = (x - mean(v, na.rm = T)) / sd(v, na.rm = T))
Y <- matrix(NA, nrow = dim(b0)[1], ncol = nrow(dat.pred))
for(i in 1:dim(Y)[2]) {
  theta <- expit(b0 + b1*dat.pred$z[i])
  Y[, i] <- apply(psi*theta, 1, sum)
}
dat.pred <- dat.pred %>%
  mutate(Y.md = apply(Y, 2, median),
         Y.lo = apply(Y, 2, function(x) quantile(x, prob = 0.025, type = 8)),
         Y.hi = apply(Y, 2, function(x) quantile(x, prob = 0.975, type = 8)))
p.LadFuel <- ggplot(data = dat.SR, aes(x = LadFuel, y = Y)) + 
  geom_point(alpha = 0.1) +
  geom_errorbar(aes(ymin = Y.lo, ymax = Y.hi), width = 0, alpha = 0.1) +
  geom_ribbon(data = dat.pred, aes(x = x, ymin = Y.lo, ymax = Y.hi), alpha = 0.2, inherit.aes = F) +
  geom_line(data = dat.pred, aes(x = x, y = Y.md), colour = "blue", size = 1.5, inherit.aes = F) + 
  labs(x = vtit, y = NULL)

v <- Herb.b
vtit <- "Herbaceous volume (Herb)"
b1 <- mod$sims.list$bb.HerbGrassVol
dat.pred <- data.frame(x = seq(min(v, na.rm = T), max(v, na.rm = T), length.out = 20)) %>%
  mutate(z = (x - mean(v, na.rm = T)) / sd(v, na.rm = T))
Y <- matrix(NA, nrow = dim(b0)[1], ncol = nrow(dat.pred))
for(i in 1:dim(Y)[2]) {
  theta <- expit(b0 + b1*dat.pred$z[i])
  Y[, i] <- apply(psi*theta, 1, sum)
}
dat.pred <- dat.pred %>%
  mutate(Y.md = apply(Y, 2, median),
         Y.lo = apply(Y, 2, function(x) quantile(x, prob = 0.025, type = 8)),
         Y.hi = apply(Y, 2, function(x) quantile(x, prob = 0.975, type = 8)))
p.Herb <- ggplot(data = dat.SR, aes(x = Herb, y = Y)) + 
  geom_point(alpha = 0.1) +
  geom_errorbar(aes(ymin = Y.lo, ymax = Y.hi), width = 0, alpha = 0.1) +
  geom_ribbon(data = dat.pred, aes(x = x, ymin = Y.lo, ymax = Y.hi), alpha = 0.2, inherit.aes = F) +
  geom_line(data = dat.pred, aes(x = x, y = Y.md), colour = "blue", size = 1.5, inherit.aes = F) + 
  labs(x = vtit, y = NULL)

# Put it all together #
p <- ggdraw() + 
  draw_plot(p.PACCGap, x = 0.03, y = 0.75, width = 0.3233333, height = 0.25) +
  draw_plot(p.PACCOpn, x = 0.3533333, y = 0.75, width = 0.3233333, height = 0.25) +
  draw_plot(p.PAROpn, x = 0.6766667, y = 0.75, width = 0.3233333, height = 0.25) +
  draw_plot(p.CanCov, x = 0.03, y = 0.5, width = 0.3233333, height = 0.25) +
  draw_plot(p.CanHt, x = 0.3533333, y = 0.5, width = 0.3233333, height = 0.25) +
  draw_plot(p.NSnag, x = 0.6766667, y = 0.5, width = 0.3233333, height = 0.25) +
  draw_plot(p.PIPO, x = 0.03, y = 0.25, width = 0.3233333, height = 0.25) +
  draw_plot(p.PSME, x = 0.3533333, y = 0.25, width = 0.3233333, height = 0.25) +
  draw_plot(p.POTR5, x = 0.6766667, y = 0.25, width = 0.3233333, height = 0.25) +
  draw_plot(p.ShrbVol, x = 0.03, y = 0, width = 0.3233333, height = 0.25) +
  draw_plot(p.LadFuel, x = 0.3533333, y = 0, width = 0.3233333, height = 0.25) +
  draw_plot(p.Herb, x = 0.6766667, y = 0, width = 0.3233333, height = 0.25) +
  draw_plot_label(c("Species richness (point)", "Species richness (grid)"),
                  x = c(0, 0), y = c(0.29, 0.79),
                  angle = c(90, 90), hjust = c(0, 0))

#save_plot("Plot_richness_hab.tiff", p, ncol = 2.25, nrow = 4, dpi = 50)
save_plot("manuscript/Figure5.tiff", p, ncol = 2.25, nrow = 4, dpi = 300)
