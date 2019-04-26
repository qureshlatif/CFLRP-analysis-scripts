library(jagsUI)
library(stringr)
library(dplyr)
library(R.utils)
library(ggplot2)
library(cowplot)
library(QSLpersonal)

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

SPR <- mod$sims.list$SR.grid
  # Derive posterior samples for plotting spp richness trend #
X <- seq(0:100)
Y <- matrix(NA, nrow = dim(SPR)[1], ncol = length(X))
for(i in 1:dim(SPR)[1]) {
  y <- as.numeric(SPR[i,,])
  x <- as.numeric(PctTrt.d)
  m <- lm(y ~ x + I(x^2))
  Y[i, ] <- predict(m, data.frame(x = X))
}
dat.pred <- data.frame(X = X,
                       Y = apply(Y, 2, median),
                       Y.lo = apply(Y, 2, function(x) quantile(x, prob = 0.025, type = 8)),
                       Y.hi = apply(Y, 2, function(x) quantile(x, prob = 0.975, type = 8)))
rm(i, y, x, m, X, Y)

dat.SR <- data.frame(X = (PctTrt.d %>% as.numeric)) %>%
  mutate(Y = apply(SPR,c(2,3),median) %>% as.numeric,
         Y.lo = apply(SPR,c(2,3),function(x) quantile(x,prob=0.025,type=8)) %>%
           as.numeric,
         Y.hi = apply(SPR,c(2,3),function(x) quantile(x,prob=0.975,type=8)) %>%
           as.numeric) %>%
  filter(!is.na(X))

p <- ggplot(data = dat.SR, aes(x = X, y = Y)) + 
  geom_point(alpha = 0.3) +
  geom_errorbar(aes(ymin = Y.lo, ymax = Y.hi), width = 0.1, alpha = 0.3) +
  geom_ribbon(data = dat.pred, aes(x = X, ymin = Y.lo, ymax = Y.hi), alpha = 0.2) +
  geom_line(data = dat.pred, aes(x = X, y = Y), colour = "blue", size = 1.5) + 
  labs(x= "Percent treated", y = "Species Richness (grid)") +
  theme(axis.title.x=element_text(size=40)) +
  theme(axis.title.y=element_text(size=40)) +
  theme(axis.text.x=element_text(size=30)) +
  theme(axis.text.y=element_text(size=30))

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

dat.pred = data.frame(x = c(0, 1))
b1 <- mod$sims.list$bb.trt
SR.pred <- matrix(NA, nrow = dim(b0)[1], ncol = nrow(dat.pred))
for(i in 1:dim(SR.pred)[2]) {
  theta <- expit(b0 + b1*dat.pred$x[i])
  SR.pred[, i] <- apply(psi * theta, 1, sum)
}
dat.pred <- dat.pred %>%
  mutate(pred.md = apply(SR.pred, 2, median),
         pred.lo = apply(SR.pred, 2, function(x) quantile(x, prob = 0.025, type = 8)),
         pred.hi = apply(SR.pred, 2, function(x) quantile(x, prob = 0.975, type = 8)))
jitter = runif(length(Trt.b), -0.4, 0.4)
p <- ggplot(data = dat.SR, aes(x = Trt, y = Y)) + 
  geom_point(aes(x = Trt + jitter), alpha = 0.1) + 
  geom_errorbar(aes(x = Trt + jitter, ymin = Y.lo, ymax = Y.hi), width = 0, alpha = 0.1) +
  geom_errorbar(data = dat.pred, aes(x = x, ymin = pred.lo, ymax = pred.hi), width = 0.3, size = 1.5, color = "blue", inherit.aes = F) +
  geom_point(data = dat.pred, aes(x = x, y = pred.md), size = 5, color = "blue", inherit.aes = F) + 
  scale_x_continuous(breaks = c(0, 1), labels = c("Untreated", "Treated")) +
  labs(x = NULL, y = "Species richness (point)") +
  theme(axis.title.y=element_text(size=40)) +
  theme(axis.text.x=element_text(size=30)) +
  theme(axis.text.y=element_text(size=30))
