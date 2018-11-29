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
cols <- c("bd.ptrt", "bd.ptrt.lo", "bd.ptrt.hi", "bb.trt", "bb.trt.lo", "bb.trt.hi", "bb.YST", "bb.YST.lo", "bb.YST.hi")
tbl_pars <- matrix(NA, nrow = length(spp.list), ncol = length(cols), dimnames = list(NULL, cols))

parm <- mod$sims.list[["bd.ptrt"]]
tbl_pars[, "bd.ptrt"] <- apply(parm, 2, median)
tbl_pars[, "bd.ptrt.lo"] <- apply(parm, 2, function(x) quantile(x, prob = 0.025, type = 8))
tbl_pars[, "bd.ptrt.hi"] <- apply(parm, 2, function(x) quantile(x, prob = 0.975, type = 8))

parm <- mod$sims.list[["bb.trt"]]
tbl_pars[, "bb.trt"] <- apply(parm, 2, median)
tbl_pars[, "bb.trt.lo"] <- apply(parm, 2, function(x) quantile(x, prob = 0.025, type = 8))
tbl_pars[, "bb.trt.hi"] <- apply(parm, 2, function(x) quantile(x, prob = 0.975, type = 8))

parm <- mod$sims.list[["bb.YST"]]
tbl_pars[, "bb.YST"] <- apply(parm, 2, median)
tbl_pars[, "bb.YST.lo"] <- apply(parm, 2, function(x) quantile(x, prob = 0.025, type = 8))
tbl_pars[, "bb.YST.hi"] <- apply(parm, 2, function(x) quantile(x, prob = 0.975, type = 8))

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

tbl_pars <- tbl_pars %>%
  arrange(bb.YST) %>%
  mutate(index = row_number())
min.lo <- min(tbl_pars$bb.YST.lo)
max.hi <- max(tbl_pars$bb.YST.hi)
pt <- ggplot(dat = tbl_pars, aes(x = index, y = bb.YST)) +
  geom_errorbar(aes(ymin = bb.YST.lo, ymax = bb.YST.hi), size=1, width=0) +
  geom_point(size = 2.5) + 
  geom_hline(yintercept = 0) +
  coord_flip() +
  scale_x_continuous(breaks = NULL, labels = NULL, expand=c(0, 1)) +
  scale_y_continuous(lim = c(min.lo, max.hi)) +
  ylab(expression(hat(beta)["years since treatment"])) + xlab(NULL) +
  theme(axis.title.y=element_text(size=30)) +
  theme(axis.title.x=element_text(size=30)) +
  theme(axis.text.x=element_text(size=15)) +
  theme(axis.text.y=element_text(size=15))

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
  labs(x= "Percent treated", y = "Species Richness") +
  theme(axis.title.x=element_text(size=40)) +
  theme(axis.title.y=element_text(size=40)) +
  theme(axis.text.x=element_text(size=30)) +
  theme(axis.text.y=element_text(size=30))

## Point level ##

Trt.b <- Cov[, "Trt_stat"] # Point-level values
SPR <- mod$sims.list$SR.point

# Derive posterior samples for plotting spp richness trend #
X <- c(0, 1)
Y <- matrix(NA, nrow = dim(SPR)[1], ncol = length(X))
for(i in 1:dim(SPR)[1]) {
  y <- as.numeric(SPR[i,])
  x <- as.numeric(Trt.b)
  m <- lm(y ~ x)
  Y[i, ] <- predict(m, data.frame(x = X))
}
dat.pred <- data.frame(X = X,
                       Y = apply(Y, 2, median),
                       Y.lo = apply(Y, 2, function(x) quantile(x, prob = 0.025, type = 8)),
                       Y.hi = apply(Y, 2, function(x) quantile(x, prob = 0.975, type = 8)))
rm(i, y, x, m, X, Y)

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
  geom_point(aes(x = X + jitter), alpha = 0.3) + 
  geom_errorbar(aes(x = X + jitter, ymin = Y.lo, ymax = Y.hi), width = 0, alpha = 0.3) +
  geom_errorbar(data = dat.pred, aes(x = X, ymin = Y.lo, ymax = Y.hi), width = 0.3, size = 2, color = "blue") +
  geom_point(data = dat.pred, aes(x = X, y = Y), size = 6, color = "blue") + 
  scale_x_continuous(breaks = c(0, 1), labels = c("Untreated", "Treated")) +
  labs(x= NULL, y = "Species Richness") +
  theme(axis.title.y=element_text(size=40)) +
  theme(axis.text.x=element_text(size=30)) +
  theme(axis.text.y=element_text(size=30))
