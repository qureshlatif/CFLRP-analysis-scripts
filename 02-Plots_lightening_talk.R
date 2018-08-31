library(jagsUI)
library(stringr)
library(dplyr)
library(R.utils)
library(ggplot2)
library(cowplot)

setwd("C:/Users/Quresh.Latif/files/projects/FS/CFLRP")
load("Data_compiled.RData")

mod <- loadObject("mod_treatment_d0yr")

# Tabulate parameter estimates
cols <- c("bd.ptrt", "bd.ptrt.lo", "bd.ptrt.hi", "bb.trt", "bb.trt.lo", "bb.trt.hi")
tbl_pars <- matrix(NA, nrow = length(spp.list), ncol = length(cols), dimnames = list(NULL, cols))

parm <- mod$sims.list[["bd.ptrt"]]
tbl_pars[, "bd.ptrt"] <- apply(parm, 2, median)
tbl_pars[, "bd.ptrt.lo"] <- apply(parm, 2, function(x) quantile(x, prob = 0.025, type = 8))
tbl_pars[, "bd.ptrt.hi"] <- apply(parm, 2, function(x) quantile(x, prob = 0.975, type = 8))

parm <- mod$sims.list[["bb.trt"]]
tbl_pars[, "bb.trt"] <- apply(parm, 2, median)
tbl_pars[, "bb.trt.lo"] <- apply(parm, 2, function(x) quantile(x, prob = 0.025, type = 8))
tbl_pars[, "bb.trt.hi"] <- apply(parm, 2, function(x) quantile(x, prob = 0.975, type = 8))

rm(parm)

# Plot species treatment effects #
tbl_pars <- tbl_pars %>% tbl_df() %>%
  mutate(Spp = spp.list) %>%
  arrange(bd.ptrt) %>%
  mutate(index = row_number())

min.lo <- min(tbl_pars$bd.ptrt.lo)
max.hi <- max(tbl_pars$bd.ptrt.hi)

pd <- ggplot(dat = tbl_pars, aes(x = index, y = bd.ptrt)) +
  geom_errorbar(aes(ymin = bd.ptrt.lo, ymax = bd.ptrt.hi), size=1, width=0) +
  geom_point(size = 2.5) + 
  geom_hline(yintercept = 0) +
  coord_flip() +
  scale_x_continuous(breaks = NULL, labels = NULL, expand=c(0, 1)) +
  scale_y_continuous(lim = c(min.lo, max.hi)) +
  ylab(expression(hat(beta)["percent treated"])) + xlab(NULL) +
  theme(axis.title.y=element_text(size=30)) +
  theme(axis.title.x=element_text(size=30)) +
  theme(axis.text.x=element_text(size=15)) +
  theme(axis.text.y=element_text(size=15))

tbl_pars <- tbl_pars %>%
  arrange(bb.trt) %>%
  mutate(index = row_number())
min.lo <- min(tbl_pars$bb.trt.lo)
max.hi <- max(tbl_pars$bb.trt.hi)
pb <- ggplot(dat = tbl_pars, aes(x = index, y = bb.trt)) +
  geom_errorbar(aes(ymin = bb.trt.lo, ymax = bb.trt.hi), size=1, width=0) +
  geom_point(size = 2.5) + 
  geom_hline(yintercept = 0) +
  coord_flip() +
  scale_x_continuous(breaks = NULL, labels = NULL, expand=c(0, 1)) +
  scale_y_continuous(lim = c(min.lo, max.hi)) +
  ylab(expression(hat(beta)["treated"])) + xlab(NULL) +
  theme(axis.title.y=element_text(size=30)) +
  theme(axis.title.x=element_text(size=30)) +
  theme(axis.text.x=element_text(size=15)) +
  theme(axis.text.y=element_text(size=15))

# Plot species richness #
gridID <- Cov[, "gridIndex"]
yearID <- Cov[, "YearInd"]
PctTrt.d <- matrix(NA, nrow = max(gridID), ncol = max(yearID))
PctTrt.d[landscape_data %>% filter(YearInd == 1) %>% pull(gridIndex), 1] <-
  landscape_data %>% filter(YearInd == 1) %>% pull(PctTrt)
PctTrt.d[landscape_data %>% filter(YearInd == 2) %>% pull(gridIndex), 2] <-
  landscape_data %>% filter(YearInd == 2) %>% pull(PctTrt)
PctTrt.d[landscape_data %>% filter(YearInd == 3) %>% pull(gridIndex), 3] <-
  landscape_data %>% filter(YearInd == 3) %>% pull(PctTrt)

SPR <- mod$sims.list$SR.grid

dat.SR <- data.frame(X = (PctTrt.d %>% as.numeric)) %>%
  mutate(Y = apply(SPR,c(2,3),median) %>% as.numeric,
         Y.lo = apply(SPR,c(2,3),function(x) quantile(x,prob=0.025,type=8)) %>%
           as.numeric,
         Y.hi = apply(SPR,c(2,3),function(x) quantile(x,prob=0.975,type=8)) %>%
           as.numeric) %>%
  filter(!is.na(X))

p <- ggplot(data = dat.SR, aes(x = X, y = Y)) + 
  geom_point(alpha = 0.3) + stat_smooth(method = lm, colour = "blue", size = 1.5) + 
  geom_errorbar(aes(ymin = Y.lo, ymax = Y.hi), width = 0.1, alpha = 0.3) +
  labs(x= "Percent treated", y = "Species Richness") +
  theme(axis.title.x=element_text(size=40)) +
  theme(axis.title.y=element_text(size=40)) +
  theme(axis.text.x=element_text(size=30)) +
  theme(axis.text.y=element_text(size=30))


Trt.b <- Cov[, "Trt_stat"] # Point-level values
SPR <- mod$sims.list$SR.point

dat.SR <- data.frame(X = (Trt.b %>% as.numeric)) %>%
  mutate(#X = X + runif(length(Trt.b), -0.4, 0.4),
         Y = apply(SPR, 2, median) %>% as.numeric,
         Y.lo = apply(SPR, 2,function(x) quantile(x,prob=0.025,type=8)) %>%
           as.numeric,
         Y.hi = apply(SPR, 2, function(x) quantile(x,prob=0.975,type=8)) %>%
           as.numeric) %>%
  filter(!is.na(X))

jitter = runif(length(Trt.b), -0.4, 0.4)
p <- ggplot(data = dat.SR, aes(x = X, y = Y)) + 
  geom_point(aes(x = X + jitter), alpha = 0.3) + stat_smooth(method = lm, colour = "blue", size = 1.5) + 
  geom_errorbar(aes(x = X + jitter, ymin = Y.lo, ymax = Y.hi), width = 0.1, alpha = 0.3) +
  scale_x_continuous(breaks = c(0, 1), labels = c("Untreated", "Treated")) +
  labs(x= NULL, y = "Species Richness") +
  theme(axis.title.y=element_text(size=40)) +
  theme(axis.text.x=element_text(size=30)) +
  theme(axis.text.y=element_text(size=30))
