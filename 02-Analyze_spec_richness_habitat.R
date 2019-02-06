library(jagsUI)
library(stringr)
library(dplyr)
library(R.utils)
library(ggplot2)
library(cowplot)
library(lme4)

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

dat.x <- data.frame(PACC10.z = PACC10_3km.d %>%
                      (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)) %>%
                      apply(1, mean, na.rm = T),
                    PACC40.z = PACC40_3km.d %>%
                      (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)) %>%
                      apply(1, mean, na.rm = T),
                    mnPAR.z = mnPerArRatio_Opn3km.d %>%
                      (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)) %>%
                      apply(1, mean, na.rm = T))

dat.pred.PACC10 <- data.frame(PACC10.x = seq(min(PACC10_3km.d, na.rm = T), max(PACC10_3km.d, na.rm = T), length.out = 20)) %>%
  mutate(PACC10.z = (PACC10.x - mean(PACC10_3km.d, na.rm = T)) / sd(PACC10_3km.d, na.rm = T),
         PACC40.z = 0, mnPAR.z = 0)
dat.pred.PACC40 <- data.frame(PACC40.x = seq(min(PACC40_3km.d, na.rm = T), max(PACC40_3km.d, na.rm = T), length.out = 20)) %>%
  mutate(PACC40.z = (PACC40.x - mean(PACC40_3km.d, na.rm = T)) / sd(PACC40_3km.d, na.rm = T),
         PACC10.z = 0, mnPAR.z = 0)
dat.pred.mnPAR <- data.frame(mnPAR.x = seq(min(mnPerArRatio_Opn3km.d, na.rm = T), max(mnPerArRatio_Opn3km.d, na.rm = T), length.out = 20)) %>%
  mutate(mnPAR.z = (mnPAR.x - mean(mnPerArRatio_Opn3km.d, na.rm = T)) / sd(mnPerArRatio_Opn3km.d, na.rm = T),
         PACC10.z = 0, PACC40.z = 0)

Y.PACC10 <- Y.PACC40 <- Y.mnPAR <- matrix(NA, nrow = dim(SPR)[1], ncol = nrow(dat.pred.PACC10))
B0.grid <- B.PACC10 <- B.PACC40 <- B.mnPAR <- numeric(length = dim(SPR)[1])
for(i in 1:dim(SPR)[1]) {
  dat = dat.x %>% mutate(N = SPR[i,,] %>% apply(1, mean))
  m <- suppressWarnings(glm(N ~ PACC10.z + PACC40.z + mnPAR.z, data = dat, family = "poisson"))
  B0.grid[i] <- m$coefficients["(Intercept)"]
  B.PACC10[i] <- m$coefficients["PACC10.z"]
  B.PACC40[i] <- m$coefficients["PACC40.z"]
  B.mnPAR[i] <- m$coefficients["mnPAR.z"]
  Y.PACC10[i, ] <- predict(m, dat.pred.PACC10, type = "response")
  Y.PACC40[i, ] <- predict(m, dat.pred.PACC40, type = "response")
  Y.mnPAR[i, ] <- predict(m, dat.pred.mnPAR, type = "response")
}

dat.pred.PACC10 <- dat.pred.PACC10 %>%
  mutate(Y.md = apply(Y.PACC10, 2, median),
         Y.lo = apply(Y.PACC10, 2, function(x) quantile(x, prob = 0.025, type = 8)),
         Y.hi = apply(Y.PACC10, 2, function(x) quantile(x, prob = 0.975, type = 8)))
dat.pred.PACC40 <- dat.pred.PACC40 %>%
  mutate(Y.md = apply(Y.PACC40, 2, median),
         Y.lo = apply(Y.PACC40, 2, function(x) quantile(x, prob = 0.025, type = 8)),
         Y.hi = apply(Y.PACC40, 2, function(x) quantile(x, prob = 0.975, type = 8)))
dat.pred.mnPAR <- dat.pred.mnPAR %>%
  mutate(Y.md = apply(Y.mnPAR, 2, median),
         Y.lo = apply(Y.mnPAR, 2, function(x) quantile(x, prob = 0.025, type = 8)),
         Y.hi = apply(Y.mnPAR, 2, function(x) quantile(x, prob = 0.975, type = 8)))

B0.grid <- str_c(median(B0.grid) %>% round(digits = 2),
                 " (",
                 quantile(B0.grid, prob = 0.025, type = 8) %>% round(digits = 2),
                 ",",
                 quantile(B0.grid, prob = 0.975, type = 8) %>% round(digits = 2),
                 ")")
B.PACC10 <- str_c(median(B.PACC10) %>% round(digits = 2),
                 " (",
                 quantile(B.PACC10, prob = 0.025, type = 8) %>% round(digits = 2),
                 ",",
                 quantile(B.PACC10, prob = 0.975, type = 8) %>% round(digits = 2),
                 ")")
B.PACC40 <- str_c(median(B.PACC40) %>% round(digits = 2),
                 " (",
                 quantile(B.PACC40, prob = 0.025, type = 8) %>% round(digits = 2),
                 ",",
                 quantile(B.PACC40, prob = 0.975, type = 8) %>% round(digits = 2),
                 ")")
B.mnPAR <- str_c(median(B.mnPAR) %>% round(digits = 3),
                  " (",
                  quantile(B.mnPAR, prob = 0.025, type = 8) %>% round(digits = 3),
                  ",",
                  quantile(B.mnPAR, prob = 0.975, type = 8) %>% round(digits = 3),
                  ")")

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
  geom_ribbon(data = dat.pred.PACC10, aes(x = PACC10.x, ymin = Y.lo, ymax = Y.hi), alpha = 0.2, inherit.aes = F) +
  geom_line(data = dat.pred.PACC10, aes(x = PACC10.x, y = Y.md), colour = "blue", size = 1.5) + 
  labs(x = "Extent of canopy gaps (PACCGap)", y = NULL)

p.PACCOpn <- ggplot(data = dat.SR, aes(x = PACC40, y = Y)) + 
  geom_point(alpha = 0.3) +
  geom_errorbar(aes(ymin = Y.lo, ymax = Y.hi), width = 0, alpha = 0.3) +
  geom_ribbon(data = dat.pred.PACC40, aes(x = PACC40.x, ymin = Y.lo, ymax = Y.hi), alpha = 0.2, inherit.aes = F) +
  geom_line(data = dat.pred.PACC40, aes(x = PACC40.x, y = Y.md), colour = "blue", size = 1.5) + 
  labs(x = "Extent of open forest (PACCOpn)", y = NULL)

# p.PAROpn <- ggplot(data = dat.SR, aes(x = mnPAR, y = Y)) + 
#   geom_point(alpha = 0.3) +
#   geom_errorbar(aes(ymin = Y.lo, ymax = Y.hi), width = 0, alpha = 0.3) +
#   geom_ribbon(data = dat.pred.mnPAR, aes(x = mnPAR.x, ymin = Y.lo, ymax = Y.hi), alpha = 0.2, inherit.aes = F) +
#   geom_line(data = dat.pred.mnPAR, aes(x = mnPAR.x, y = Y.md), colour = "blue", size = 1.5) + 
#   labs(x = "Open forest perimeter-area ratio (PAROpn)", y = NULL)

p <- ggdraw() + 
  draw_plot(p.PACCGap, x = 0.05, y = 0, width = 0.485, height = 1) +
  draw_plot(p.PACCOpn, x = 0.525, y = 0, width = 0.485, height = 1) +
  # draw_plot(p.PAROpn, x = 0.6766667, y = 0, width = 0.3233333, height = 1) +
  draw_plot_label(c("N[psi]", "',k'"), x = c(0, 0.027),
                  y = c(0.5, 0.58), size = c(20, 12),
                  angle = c(90, 90), hjust = c(0, 0), parse = T) #****Need to figure out how to get k in there.

save_plot("Plot_richness_landscape.tiff", p, ncol = 2, nrow = 1, dpi = 200)

## Point level ##

SPR <- mod$sims.list$SR.point
digts <- 5 # Set number of digits to save for tabulation summaries

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

dat.x <- data.frame(ID = ID,
                    CanCov.z = CanCov.b %>%
                      (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)),
                    CanHt.z = CanHt.b %>%
                      (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)),
                    NSnag.z = NSnag.b %>%
                      (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)),
                    PIPO.z = PIPO.b %>%
                      (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)),
                    PSME.z = PSME.b %>%
                      (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)),
                    POTR5.z = POTR5.b %>%
                      (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)),
                    ShrbVol.z = ShrbVol.b %>%
                      (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)),
                    LadFuel.z = LadFuel.b %>%
                      (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)),
                    Herb.z = Herb.b %>%
                      (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)))


# #__Run this chunk once and then load from cache__#
var.names <- c("CanCov", "CanHt", "NSnag", "PIPO", "PSME", "POTR5", "ShrbVol", "LadFuel", "Herb")
# 
# dat.pred.CanCov <- data.frame(CanCov.x = seq(min(CanCov.b, na.rm = T), max(CanCov.b, na.rm = T), length.out = 20)) %>%
#   mutate(CanCov.z = (CanCov.x - mean(CanCov.b, na.rm = T)) / sd(CanCov.b, na.rm = T),
#          CanHt.z = 0, NSnag.z = 0, PIPO.z = 0, PSME.z = 0, POTR5.z = 0, ShrbVol.z = 0, LadFuel.z = 0, Herb.z = 0)
# dat.pred.CanHt <- data.frame(CanHt.x = seq(min(CanHt.b, na.rm = T), max(CanHt.b, na.rm = T), length.out = 20)) %>%
#   mutate(CanHt.z = (CanHt.x - mean(CanHt.b, na.rm = T)) / sd(CanHt.b, na.rm = T),
#          CanCov.z = 0, NSnag.z = 0, PIPO.z = 0, PSME.z = 0, POTR5.z = 0, ShrbVol.z = 0, LadFuel.z = 0, Herb.z = 0)
# dat.pred.NSnag <- data.frame(NSnag.x = seq(min(NSnag.b, na.rm = T), max(NSnag.b, na.rm = T), length.out = 20)) %>%
#   mutate(NSnag.z = (NSnag.x - mean(NSnag.b, na.rm = T)) / sd(NSnag.b, na.rm = T),
#          CanCov.z = 0, CanHt.z = 0, PIPO.z = 0, PSME.z = 0, POTR5.z = 0, ShrbVol.z = 0, LadFuel.z = 0, Herb.z = 0)
# dat.pred.PIPO <- data.frame(PIPO.x = seq(min(PIPO.b, na.rm = T), max(PIPO.b, na.rm = T), length.out = 20)) %>%
#   mutate(PIPO.z = (PIPO.x - mean(PIPO.b, na.rm = T)) / sd(PIPO.b, na.rm = T),
#          CanCov.z = 0, CanHt.z = 0, NSnag.z = 0, PSME.z = 0, POTR5.z = 0, ShrbVol.z = 0, LadFuel.z = 0, Herb.z = 0)
# dat.pred.PSME <- data.frame(PSME.x = seq(min(PSME.b, na.rm = T), max(PSME.b, na.rm = T), length.out = 20)) %>%
#   mutate(PSME.z = (PSME.x - mean(PSME.b, na.rm = T)) / sd(PSME.b, na.rm = T),
#          CanCov.z = 0, CanHt.z = 0, NSnag.z = 0, PIPO.z = 0, POTR5.z = 0, ShrbVol.z = 0, LadFuel.z = 0, Herb.z = 0)
# dat.pred.POTR5 <- data.frame(POTR5.x = seq(min(POTR5.b, na.rm = T), max(POTR5.b, na.rm = T), length.out = 20)) %>%
#   mutate(POTR5.z = (POTR5.x - mean(POTR5.b, na.rm = T)) / sd(POTR5.b, na.rm = T),
#          CanCov.z = 0, CanHt.z = 0, NSnag.z = 0, PIPO.z = 0, PSME.z = 0, ShrbVol.z = 0, LadFuel.z = 0, Herb.z = 0)
# dat.pred.ShrbVol <- data.frame(ShrbVol.x = seq(min(ShrbVol.b, na.rm = T), max(ShrbVol.b, na.rm = T), length.out = 20)) %>%
#   mutate(ShrbVol.z = (ShrbVol.x - mean(ShrbVol.b, na.rm = T)) / sd(ShrbVol.b, na.rm = T),
#          CanCov.z = 0, CanHt.z = 0, NSnag.z = 0, PIPO.z = 0, PSME.z = 0, POTR5.z = 0, LadFuel.z = 0, Herb.z = 0)
# dat.pred.LadFuel <- data.frame(LadFuel.x = seq(min(LadFuel.b, na.rm = T), max(LadFuel.b, na.rm = T), length.out = 20)) %>%
#   mutate(LadFuel.z = (LadFuel.x - mean(LadFuel.b, na.rm = T)) / sd(LadFuel.b, na.rm = T),
#          CanCov.z = 0, CanHt.z = 0, NSnag.z = 0, PIPO.z = 0, PSME.z = 0, POTR5.z = 0, ShrbVol.z = 0, Herb.z = 0)
# dat.pred.Herb <- data.frame(Herb.x = seq(min(Herb.b, na.rm = T), max(Herb.b, na.rm = T), length.out = 20)) %>%
#   mutate(Herb.z = (Herb.x - mean(Herb.b, na.rm = T)) / sd(Herb.b, na.rm = T),
#          CanCov.z = 0, CanHt.z = 0, NSnag.z = 0, PIPO.z = 0, PSME.z = 0, POTR5.z = 0, ShrbVol.z = 0, LadFuel.z = 0)
# 
# Y.CanCov <- Y.CanHt <- Y.NSnag <- Y.PIPO <- Y.PSME <- Y.POTR5 <- Y.ShrbVol <- Y.LadFuel <- Y.Herb <-
#   matrix(NA, nrow = dim(SPR)[1], ncol = nrow(dat.pred.CanCov))
# B0.point <- B0.pntSD <- B.CanCov <- B.CanHt <- B.NSnag <- B.PIPO <- B.PSME <- B.POTR5 <- B.ShrbVol <- B.LadFuel <- B.Herb <-
#   numeric(length = dim(SPR)[1])
# options(warn=1) # Change to warn=2 to investigate convergence warnings.
# #***Note: convergence warnings arose, but for one iteration, the warning went away after removing the random effect,
#   # and slope parameters were similar with vs. without the random effect. So, ignoring convergence warnings.
# for(i in 1:dim(SPR)[1]) {
#   dat = dat.x %>% mutate(N = SPR[i,])
#   m <- suppressWarnings(glmer(N ~ CanCov.z + CanHt.z + NSnag.z + PIPO.z + PSME.z + POTR5.z + ShrbVol.z +
#                LadFuel.z + Herb.z + (1|ID), data = dat, family = "poisson"))
#   B0.point[i] <- summary(m)$coefficients["(Intercept)", "Estimate"]
#   B0.pntSD[i] <- as.data.frame(VarCorr(m))$sdcor
#   for(v in var.names) {
#     Y <- eval(as.name(str_c("Y.", v)))
#     B <- eval(as.name(str_c("B.", v)))
#     dat.pred <- eval(as.name(str_c("dat.pred.", v)))
#     Y[i, ] <- predict(m, dat.pred, re.form = NA, type = "response")
#     B[i] <- summary(m)$coefficients[str_c(v, ".z"), "Estimate"]
#     assign(str_c("Y.", v), Y)
#     assign(str_c("B.", v), B)
#   }
# }
# 
# for(v in var.names) {
#   Y <- eval(as.name(str_c("Y.", v)))
#   dat.pred <- eval(as.name(str_c("dat.pred.", v)))
#   dat.pred <- dat.pred %>%
#     mutate(Y.md = apply(Y, 2, median),
#            Y.lo = apply(Y, 2, function(x) quantile(x, prob = 0.025, type = 8)),
#            Y.hi = apply(Y, 2, function(x) quantile(x, prob = 0.975, type = 8)))
#   assign(str_c("dat.pred.", v), dat.pred)
# }
# 
# for(v in var.names) {
#   B <- eval(as.name(str_c("B.", v)))
#   saveObject(B, str_c("B", v, "_pnt_N_hab_cache"))
#   dat.pred <- eval(as.name(str_c("dat.pred.", v)))
#   saveObject(dat.pred, str_c("Spp_richness_dat_pred_", v, "_cache"))
# }
# saveObject(B0.point, "B0_pnt_N_hab_cache")
# saveObject(B0.pntSD, "B0_pntSD_N_hab_cache")
# #________________________________________________#
# rm(i, m, Y, B, dat.pred)

for(v in var.names) {
  assign(str_c("B.", v), loadObject(str_c("B", v, "_pnt_N_hab_cache")))
  assign(str_c("dat.pred.", v), loadObject(str_c("Spp_richness_dat_pred_", v, "_cache")))
}
B0.point <- loadObject("B0_pnt_N_hab_cache")
B0.pntSD <- loadObject("B0_pntSD_N_hab_cache")

for(v in var.names) {
  B <- eval(as.name(str_c("B.", v)))
  B <- str_c(median(B, na.rm = T) %>% round(digits = digts), " (",
             quantile(B, prob = 0.025, type = 8, na.rm = T) %>% round(digits = digts), ",",
             quantile(B, prob = 0.975, type = 8, na.rm = T) %>% round(digits = digts), ")")
  assign(str_c("B.", v), B)
}
B0.point <- str_c(median(B0.point, na.rm = T) %>% round(digits = digts), " (",
                  quantile(B0.point, prob = 0.025, type = 8, na.rm = T) %>% round(digits = digts), ",",
                  quantile(B0.point, prob = 0.975, type = 8, na.rm = T) %>% round(digits = digts), ")")
B0.pntSD <- str_c(median(B0.pntSD, na.rm = T) %>% round(digits = digts), " (",
                  quantile(B0.pntSD, prob = 0.025, type = 8, na.rm = T) %>% round(digits = digts), ",",
                  quantile(B0.pntSD, prob = 0.975, type = 8, na.rm = T) %>% round(digits = digts), ")")
rm(v, B)

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

p.CanCov <- ggplot(data = dat.SR, aes(x = CanCov, y = Y)) + 
  geom_point(alpha = 0.3) +
  geom_errorbar(aes(ymin = Y.lo, ymax = Y.hi), width = 0, alpha = 0.3) +
  geom_ribbon(data = dat.pred.CanCov, aes(x = CanCov.x, ymin = Y.lo, ymax = Y.hi), alpha = 0.2, inherit.aes = F) +
  geom_line(data = dat.pred.CanCov, aes(x = CanCov.x, y = Y.md), colour = "blue", size = 1.5, inherit.aes = F) + 
  labs(x = "Canopy cover (CanCov)", y = NULL)

p.CanHt <- ggplot(data = dat.SR, aes(x = CanHt, y = Y)) + 
  geom_point(alpha = 0.3) +
  geom_errorbar(aes(ymin = Y.lo, ymax = Y.hi), width = 0, alpha = 0.3) +
  geom_ribbon(data = dat.pred.CanHt, aes(x = CanHt.x, ymin = Y.lo, ymax = Y.hi), alpha = 0.2, inherit.aes = F) +
  geom_line(data = dat.pred.CanHt, aes(x = CanHt.x, y = Y.md), colour = "blue", size = 1.5, inherit.aes = F) + 
  labs(x = "Canopy height (CanHt)", y = NULL)

p.NSnag <- ggplot(data = dat.SR, aes(x = NSnag, y = Y)) + 
  geom_point(alpha = 0.3) +
  geom_errorbar(aes(ymin = Y.lo, ymax = Y.hi), width = 0, alpha = 0.3) +
  geom_ribbon(data = dat.pred.NSnag, aes(x = NSnag.x, ymin = Y.lo, ymax = Y.hi), alpha = 0.2, inherit.aes = F) +
  geom_line(data = dat.pred.NSnag, aes(x = NSnag.x, y = Y.md), colour = "blue", size = 1.5, inherit.aes = F) + 
  labs(x = "Number of snags (NSnag)", y = NULL)

# p.PIPO <- ggplot(data = dat.SR, aes(x = PIPO, y = Y)) + 
#   geom_point(alpha = 0.3) +
#   geom_errorbar(aes(ymin = Y.lo, ymax = Y.hi), width = 0, alpha = 0.3) +
#   geom_ribbon(data = dat.pred.PIPO, aes(x = PIPO.x, ymin = Y.lo, ymax = Y.hi), alpha = 0.2, inherit.aes = F) +
#   geom_line(data = dat.pred.PIPO, aes(x = PIPO.x, y = Y.md), colour = "blue", size = 1.5, inherit.aes = F) + 
#   labs(x = "Ponderosa pine canopy (PIPO)", y = NULL)
# 
# p.PSME <- ggplot(data = dat.SR, aes(x = PSME, y = Y)) + 
#   geom_point(alpha = 0.3) +
#   geom_errorbar(aes(ymin = Y.lo, ymax = Y.hi), width = 0, alpha = 0.3) +
#   geom_ribbon(data = dat.pred.PSME, aes(x = PSME.x, ymin = Y.lo, ymax = Y.hi), alpha = 0.2, inherit.aes = F) +
#   geom_line(data = dat.pred.PSME, aes(x = PSME.x, y = Y.md), colour = "blue", size = 1.5, inherit.aes = F) + 
#   labs(x = "Douglas fir canopy (PSME)", y = NULL)

p.POTR5 <- ggplot(data = dat.SR, aes(x = POTR5, y = Y)) + 
  geom_point(alpha = 0.3) +
  geom_errorbar(aes(ymin = Y.lo, ymax = Y.hi), width = 0, alpha = 0.3) +
  geom_ribbon(data = dat.pred.POTR5, aes(x = POTR5.x, ymin = Y.lo, ymax = Y.hi), alpha = 0.2, inherit.aes = F) +
  geom_line(data = dat.pred.POTR5, aes(x = POTR5.x, y = Y.md), colour = "blue", size = 1.5, inherit.aes = F) + 
  labs(x = "Aspen canopy (POTR5)", y = NULL)

p.ShrbVol <- ggplot(data = dat.SR, aes(x = ShrbVol, y = Y)) + 
  geom_point(alpha = 0.3) +
  geom_errorbar(aes(ymin = Y.lo, ymax = Y.hi), width = 0, alpha = 0.3) +
  geom_ribbon(data = dat.pred.ShrbVol, aes(x = ShrbVol.x, ymin = Y.lo, ymax = Y.hi), alpha = 0.2, inherit.aes = F) +
  geom_line(data = dat.pred.ShrbVol, aes(x = ShrbVol.x, y = Y.md), colour = "blue", size = 1.5, inherit.aes = F) + 
  labs(x = "Shrub volume (ShrbVol)", y = NULL)

p.LadFuel <- ggplot(data = dat.SR, aes(x = LadFuel, y = Y)) + 
  geom_point(alpha = 0.3) +
  geom_errorbar(aes(ymin = Y.lo, ymax = Y.hi), width = 0, alpha = 0.3) +
  geom_ribbon(data = dat.pred.LadFuel, aes(x = LadFuel.x, ymin = Y.lo, ymax = Y.hi), alpha = 0.2, inherit.aes = F) +
  geom_line(data = dat.pred.LadFuel, aes(x = LadFuel.x, y = Y.md), colour = "blue", size = 1.5, inherit.aes = F) + 
  labs(x = "Ladder fuels (LadFuel)", y = NULL)

# p.Herb <- ggplot(data = dat.SR, aes(x = Herb, y = Y)) + 
#   geom_point(alpha = 0.3) +
#   geom_errorbar(aes(ymin = Y.lo, ymax = Y.hi), width = 0, alpha = 0.3) +
#   geom_ribbon(data = dat.pred.Herb, aes(x = Herb.x, ymin = Y.lo, ymax = Y.hi), alpha = 0.2, inherit.aes = F) +
#   geom_line(data = dat.pred.Herb, aes(x = Herb.x, y = Y.md), colour = "blue", size = 1.5, inherit.aes = F) + 
#   labs(x = "Herbaceous volume (Herb)", y = NULL)

p <- ggdraw() + 
  draw_plot(p.CanCov, x = 0.03, y = 0.5, width = 0.3233333, height = 0.5) +
  draw_plot(p.CanHt, x = 0.3533333, y = 0.5, width = 0.3233333, height = 0.5) +
  draw_plot(p.NSnag, x = 0.6766667, y = 0.5, width = 0.3233333, height = 0.5) +
  # draw_plot(p.PIPO, x = 0.03, y = 0.3333, width = 0.3233333, height = 0.3333) +
  # draw_plot(p.PSME, x = 0.3533333, y = 0.3333, width = 0.3233333, height = 0.3333) +
  draw_plot(p.POTR5, x = 0.03, y = 0, width = 0.3233333, height = 0.5) +
  draw_plot(p.ShrbVol, x = 0.3533333, y = 0, width = 0.3233333, height = 0.5) +
  draw_plot(p.LadFuel, x = 0.6766667, y = 0, width = 0.3233333, height = 0.5) +
  # draw_plot(p.Herb, x = 0.6766667, y = 0, width = 0.3233333, height = 0.3333) +
  draw_plot_label("N[theta]", x = 0, y = 0.5, size = 20, angle = 90, hjust = 0, parse = T)

save_plot("Plot_richness_veg.tiff", p, ncol = 3, nrow = 2, dpi = 200)


# # Shrub-overstory height ratio #
# X.b <- Cov[, "SOHtRatio"] # Point-level values
# 
# X <- seq(min(X.b, na.rm = T), max(X.b, na.rm = T), length.out = 20)
# Y <- matrix(NA, nrow = dim(SPR)[1], ncol = length(X))
# r <- numeric(length = dim(SPR)[1])
# for(i in 1:dim(SPR)[1]) {
#   y <- as.numeric(SPR[i,])
#   x <- as.numeric(X.b)
#   y <- y[-which(is.na(x))]
#   x <- x[-which(is.na(x))]
#   m <- lm(y ~ x)
#   Y[i, ] <- predict(m, data.frame(x = X))
#   r[i] <- cor(x, y)
# }
# dat.pred <- data.frame(X = X,
#                        Y = apply(Y, 2, median),
#                        Y.lo = apply(Y, 2, function(x) quantile(x, prob = 0.025, type = 8)),
#                        Y.hi = apply(Y, 2, function(x) quantile(x, prob = 0.975, type = 8)))
# rm(i, y, x, m, X, Y)
# 
# r.SOHRatio <- str_c(median(r) %>% round(digits = 2),
#                     " (",
#                     quantile(r, prob = 0.025, type = 8) %>% round(digits = 2),
#                     ",",
#                     quantile(r, prob = 0.975, type = 8) %>% round(digits = 2),
#                     ")")
# 
# dat.SR <- data.frame(X = (X.b %>% as.numeric)) %>%
#   mutate(#X = X + runif(length(Trt.b), -0.4, 0.4),
#     Y = apply(SPR, 2, median) %>% as.numeric,
#     Y.lo = apply(SPR, 2,function(x) quantile(x,prob=0.025,type=8)) %>%
#       as.numeric,
#     Y.hi = apply(SPR, 2, function(x) quantile(x,prob=0.975,type=8)) %>%
#       as.numeric) %>%
#   filter(!is.na(X))
# 
# p.SOHRatio <- ggplot(data = dat.SR, aes(x = X, y = Y)) + 
#   geom_point(alpha = 0.3) +
#   geom_errorbar(aes(ymin = Y.lo, ymax = Y.hi), width = 0, alpha = 0.3) +
#   geom_ribbon(data = dat.pred, aes(x = X, ymin = Y.lo, ymax = Y.hi), alpha = 0.2) +
#   geom_line(data = dat.pred, aes(x = X, y = Y), colour = "blue", size = 1.5) + 
#   labs(x = "Shrub-overstory height ratio (SOHtRatio)", y = NULL) +
#   annotate("text", x = 0, y = 25, label = str_c("r = ", r.SOHRatio), hjust = 0, size = 6) # +
# #theme(axis.title.x=element_text(size=40)) +
# #theme(axis.text.x=element_text(size=30))

# # Shrub diversity #
# X.b <- Cov[, "ShrubDiv"] # Point-level values
# 
# X <- seq(min(X.b, na.rm = T), max(X.b, na.rm = T), length.out = 20)
# Y <- matrix(NA, nrow = dim(SPR)[1], ncol = length(X))
# r <- numeric(length = dim(SPR)[1])
# for(i in 1:dim(SPR)[1]) {
#   y <- as.numeric(SPR[i,])
#   x <- as.numeric(X.b)
#   y <- y[-which(is.na(x))]
#   x <- x[-which(is.na(x))]
#   m <- lm(y ~ x)
#   Y[i, ] <- predict(m, data.frame(x = X))
#   r[i] <- cor(x, y)
# }
# dat.pred <- data.frame(X = X,
#                        Y = apply(Y, 2, median),
#                        Y.lo = apply(Y, 2, function(x) quantile(x, prob = 0.025, type = 8)),
#                        Y.hi = apply(Y, 2, function(x) quantile(x, prob = 0.975, type = 8)))
# rm(i, y, x, m, X, Y)
# 
# r.ShrubDiv <- str_c(median(r) %>% round(digits = 2),
#                     " (",
#                     quantile(r, prob = 0.025, type = 8) %>% round(digits = 2),
#                     ",",
#                     quantile(r, prob = 0.975, type = 8) %>% round(digits = 2),
#                     ")")
# 
# dat.SR <- data.frame(X = (X.b %>% as.numeric)) %>%
#   mutate(#X = X + runif(length(Trt.b), -0.4, 0.4),
#     Y = apply(SPR, 2, median) %>% as.numeric,
#     Y.lo = apply(SPR, 2,function(x) quantile(x,prob=0.025,type=8)) %>%
#       as.numeric,
#     Y.hi = apply(SPR, 2, function(x) quantile(x,prob=0.975,type=8)) %>%
#       as.numeric) %>%
#   filter(!is.na(X))
# 
# p.ShrubDiv <- ggplot(data = dat.SR, aes(x = X, y = Y)) + 
#   geom_point(alpha = 0.3) +
#   geom_errorbar(aes(ymin = Y.lo, ymax = Y.hi), width = 0, alpha = 0.3) +
#   geom_ribbon(data = dat.pred, aes(x = X, ymin = Y.lo, ymax = Y.hi), alpha = 0.2) +
#   geom_line(data = dat.pred, aes(x = X, y = Y), colour = "blue", size = 1.5) + 
#   labs(x = "Shrub diversity (ShrubDiv)", y = NULL) +
#   annotate("text", x = 0, y = 25, label = str_c("r = ", r.ShrubDiv), hjust = 0, size = 6) # +
# #theme(axis.title.x=element_text(size=40)) +
# #theme(axis.text.x=element_text(size=30))
