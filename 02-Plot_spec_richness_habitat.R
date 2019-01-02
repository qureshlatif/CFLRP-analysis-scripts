library(jagsUI)
library(stringr)
library(dplyr)
library(R.utils)
library(ggplot2)
library(cowplot)

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
X.d <- PACC10_3km.d

X <- seq(min(X.d, na.rm = T), max(X.d, na.rm = T), length.out = 20)
Y <- matrix(NA, nrow = dim(SPR)[1], ncol = length(X))
r <- numeric(length = dim(SPR)[1])
for(i in 1:dim(SPR)[1]) {
  y <- as.numeric(SPR[i,,])
  x <- as.numeric(X.d)
  y <- y[-which(is.na(x))]
  x <- x[-which(is.na(x))]
  m <- lm(y ~ x)
  Y[i, ] <- predict(m, data.frame(x = X))
  r[i] <- cor(x, y)
}
dat.pred <- data.frame(X = X,
                       Y = apply(Y, 2, median),
                       Y.lo = apply(Y, 2, function(x) quantile(x, prob = 0.025, type = 8)),
                       Y.hi = apply(Y, 2, function(x) quantile(x, prob = 0.975, type = 8)))
r.PACCGap <- str_c(median(r) %>% round(digits = 2),
              " (",
              quantile(r, prob = 0.025, type = 8) %>% round(digits = 2),
              ",",
              quantile(r, prob = 0.975, type = 8) %>% round(digits = 2),
              ")")
rm(i, y, x, m, X, Y)

dat.SR <- data.frame(X = (X.d %>% as.numeric)) %>%
  mutate(Y = apply(SPR, c(2,3),median) %>% as.numeric,
         Y.lo = apply(SPR, c(2,3), function(x) quantile(x, prob = 0.025, type=8)) %>%
           as.numeric,
         Y.hi = apply(SPR, c(2,3), function(x) quantile(x, prob = 0.975, type=8)) %>%
           as.numeric) %>%
  filter(!is.na(X))

p.PACCGap <- ggplot(data = dat.SR, aes(x = X, y = Y)) + 
  geom_point(alpha = 0.3) +
  geom_errorbar(aes(ymin = Y.lo, ymax = Y.hi), width = 0, alpha = 0.3) +
  geom_ribbon(data = dat.pred, aes(x = X, ymin = Y.lo, ymax = Y.hi), alpha = 0.2) +
  geom_line(data = dat.pred, aes(x = X, y = Y), colour = "blue", size = 1.5) + 
  labs(x = "Extent of canopy gaps (PACCGap)", y = NULL) +
  annotate("text", x = 0, y = 50, label = str_c("r = ", r.PACCGap), hjust = 0, size = 6) # +
  #theme(axis.title.x=element_text(size=40)) +
  #theme(axis.text.x=element_text(size=30))

# Extent of open forest (PACCGap) #
PACC40_3km.d <- matrix(NA, nrow = max(gridID), ncol = max(yearID))
PACC40_3km.d[landscape_data %>% filter(YearInd == 1) %>% pull(gridIndex), 1] <-
  landscape_data %>% filter(YearInd == 1) %>% pull(PACC40_3km)
PACC40_3km.d[landscape_data %>% filter(YearInd == 2) %>% pull(gridIndex), 2] <-
  landscape_data %>% filter(YearInd == 2) %>% pull(PACC40_3km)
PACC40_3km.d[landscape_data %>% filter(YearInd == 3) %>% pull(gridIndex), 3] <-
  landscape_data %>% filter(YearInd == 3) %>% pull(PACC40_3km)
X.d <- PACC40_3km.d

X <- seq(min(X.d, na.rm = T), max(X.d, na.rm = T), length.out = 20)
Y <- matrix(NA, nrow = dim(SPR)[1], ncol = length(X))
r <- numeric(length = dim(SPR)[1])
for(i in 1:dim(SPR)[1]) {
  y <- as.numeric(SPR[i,,])
  x <- as.numeric(X.d)
  y <- y[-which(is.na(x))]
  x <- x[-which(is.na(x))]
  m <- lm(y ~ x)
  Y[i, ] <- predict(m, data.frame(x = X))
  r[i] <- cor(x, y)
}
dat.pred <- data.frame(X = X,
                       Y = apply(Y, 2, median),
                       Y.lo = apply(Y, 2, function(x) quantile(x, prob = 0.025, type = 8)),
                       Y.hi = apply(Y, 2, function(x) quantile(x, prob = 0.975, type = 8)))
r.PACCOpn <- str_c(median(r) %>% round(digits = 2),
                   " (",
                   quantile(r, prob = 0.025, type = 8) %>% round(digits = 2),
                   ",",
                   quantile(r, prob = 0.975, type = 8) %>% round(digits = 2),
                   ")")
rm(i, y, x, m, X, Y)

dat.SR <- data.frame(X = (X.d %>% as.numeric)) %>%
  mutate(Y = apply(SPR, c(2,3),median) %>% as.numeric,
         Y.lo = apply(SPR, c(2,3), function(x) quantile(x, prob = 0.025, type=8)) %>%
           as.numeric,
         Y.hi = apply(SPR, c(2,3), function(x) quantile(x, prob = 0.975, type=8)) %>%
           as.numeric) %>%
  filter(!is.na(X))

p.PACCOpn <- ggplot(data = dat.SR, aes(x = X, y = Y)) + 
  geom_point(alpha = 0.3) +
  geom_errorbar(aes(ymin = Y.lo, ymax = Y.hi), width = 0, alpha = 0.3) +
  geom_ribbon(data = dat.pred, aes(x = X, ymin = Y.lo, ymax = Y.hi), alpha = 0.2) +
  geom_line(data = dat.pred, aes(x = X, y = Y), colour = "blue", size = 1.5) + 
  labs(x = "Extent of open forest (PACCOpn)", y = NULL) +
  annotate("text", x = 0, y = 50, label = str_c("r = ", r.PACCOpn), hjust = 0, size = 6) # +
#theme(axis.title.x=element_text(size=40)) +
#theme(axis.text.x=element_text(size=30))

# Perimeter-area ratio for open forest (PAROpn) #
mnPerArRatio_Opn3km.d <- matrix(NA, nrow = max(gridID), ncol = max(yearID))
mnPerArRatio_Opn3km.d[landscape_data %>% filter(YearInd == 1) %>% pull(gridIndex), 1] <-
  landscape_data %>% filter(YearInd == 1) %>% pull(mnPerArRatio_Opn3km)
mnPerArRatio_Opn3km.d[landscape_data %>% filter(YearInd == 2) %>% pull(gridIndex), 2] <-
  landscape_data %>% filter(YearInd == 2) %>% pull(mnPerArRatio_Opn3km)
mnPerArRatio_Opn3km.d[landscape_data %>% filter(YearInd == 3) %>% pull(gridIndex), 3] <-
  landscape_data %>% filter(YearInd == 3) %>% pull(mnPerArRatio_Opn3km)
X.d <- mnPerArRatio_Opn3km.d

X <- seq(min(X.d, na.rm = T), max(X.d, na.rm = T), length.out = 20)
Y <- matrix(NA, nrow = dim(SPR)[1], ncol = length(X))
r <- numeric(length = dim(SPR)[1])
for(i in 1:dim(SPR)[1]) {
  y <- as.numeric(SPR[i,,])
  x <- as.numeric(X.d)
  y <- y[-which(is.na(x))]
  x <- x[-which(is.na(x))]
  m <- lm(y ~ x)
  Y[i, ] <- predict(m, data.frame(x = X))
  r[i] <- cor(x, y)
}
dat.pred <- data.frame(X = X,
                       Y = apply(Y, 2, median),
                       Y.lo = apply(Y, 2, function(x) quantile(x, prob = 0.025, type = 8)),
                       Y.hi = apply(Y, 2, function(x) quantile(x, prob = 0.975, type = 8)))
r.PAROpn <- str_c(median(r) %>% round(digits = 2),
                   " (",
                   quantile(r, prob = 0.025, type = 8) %>% round(digits = 2),
                   ",",
                   quantile(r, prob = 0.975, type = 8) %>% round(digits = 2),
                   ")")
rm(i, y, x, m, X, Y)

dat.SR <- data.frame(X = (X.d %>% as.numeric)) %>%
  mutate(Y = apply(SPR, c(2,3),median) %>% as.numeric,
         Y.lo = apply(SPR, c(2,3), function(x) quantile(x, prob = 0.025, type=8)) %>%
           as.numeric,
         Y.hi = apply(SPR, c(2,3), function(x) quantile(x, prob = 0.975, type=8)) %>%
           as.numeric) %>%
  filter(!is.na(X))

p.PAROpn <- ggplot(data = dat.SR, aes(x = X, y = Y)) + 
  geom_point(alpha = 0.3) +
  geom_errorbar(aes(ymin = Y.lo, ymax = Y.hi), width = 0, alpha = 0.3) +
  geom_ribbon(data = dat.pred, aes(x = X, ymin = Y.lo, ymax = Y.hi), alpha = 0.2) +
  geom_line(data = dat.pred, aes(x = X, y = Y), colour = "blue", size = 1.5) + 
  labs(x = "Open forest perimeter-area ratio (PAROpn)", y = NULL) +
  annotate("text", x = 0.06, y = 50, label = str_c("r = ", r.PAROpn), hjust = 0, size = 6) # +
#theme(axis.title.x=element_text(size=40)) +
#theme(axis.text.x=element_text(size=30))

p <- ggdraw() + 
  draw_plot(p.PACCGap, x = 0.03, y = 0, width = 0.3233333, height = 1) +
  draw_plot(p.PACCOpn, x = 0.3533333, y = 0, width = 0.3233333, height = 1) +
  draw_plot(p.PAROpn, x = 0.6766667, y = 0, width = 0.3233333, height = 1) +
  draw_plot_label("Species richness", x = 0, y = 0.3, size = 20, angle = 90, hjust = 0)

save_plot("Plot_richness_landscape.tiff", p, ncol = 3, nrow = 1, dpi = 200)

## Point level ##

SPR <- mod$sims.list$SR.point

# Canopy cover #
X.b <- Cov[, "CanCov"] # Point-level values

X <- seq(min(X.b, na.rm = T), max(X.b, na.rm = T), length.out = 20)
Y <- matrix(NA, nrow = dim(SPR)[1], ncol = length(X))
r <- numeric(length = dim(SPR)[1])
for(i in 1:dim(SPR)[1]) {
  y <- as.numeric(SPR[i,])
  x <- as.numeric(X.b)
  y <- y[-which(is.na(x))]
  x <- x[-which(is.na(x))]
  m <- lm(y ~ x)
  Y[i, ] <- predict(m, data.frame(x = X))
  r[i] <- cor(x, y)
}
dat.pred <- data.frame(X = X,
                       Y = apply(Y, 2, median),
                       Y.lo = apply(Y, 2, function(x) quantile(x, prob = 0.025, type = 8)),
                       Y.hi = apply(Y, 2, function(x) quantile(x, prob = 0.975, type = 8)))
rm(i, y, x, m, X, Y)

r.CanCov <- str_c(median(r) %>% round(digits = 2),
                  " (",
                  quantile(r, prob = 0.025, type = 8) %>% round(digits = 2),
                  ",",
                  quantile(r, prob = 0.975, type = 8) %>% round(digits = 2),
                  ")")

dat.SR <- data.frame(X = (X.b %>% as.numeric)) %>%
  mutate(#X = X + runif(length(Trt.b), -0.4, 0.4),
    Y = apply(SPR, 2, median) %>% as.numeric,
    Y.lo = apply(SPR, 2,function(x) quantile(x,prob=0.025,type=8)) %>%
      as.numeric,
    Y.hi = apply(SPR, 2, function(x) quantile(x,prob=0.975,type=8)) %>%
      as.numeric) %>%
  filter(!is.na(X))

p.CanCov <- ggplot(data = dat.SR, aes(x = X, y = Y)) + 
  geom_point(alpha = 0.3) +
  geom_errorbar(aes(ymin = Y.lo, ymax = Y.hi), width = 0, alpha = 0.3) +
  geom_ribbon(data = dat.pred, aes(x = X, ymin = Y.lo, ymax = Y.hi), alpha = 0.2) +
  geom_line(data = dat.pred, aes(x = X, y = Y), colour = "blue", size = 1.5) + 
  labs(x = "Canopy cover (CanCov)", y = NULL) +
  annotate("text", x = 0.06, y = 25, label = str_c("r = ", r.CanCov), hjust = 0, size = 6) # +
#theme(axis.title.x=element_text(size=40)) +
#theme(axis.text.x=element_text(size=30))

# Canopy height #
X.b <- Cov[, "CanHt"] # Point-level values

X <- seq(min(X.b, na.rm = T), max(X.b, na.rm = T), length.out = 20)
Y <- matrix(NA, nrow = dim(SPR)[1], ncol = length(X))
r <- numeric(length = dim(SPR)[1])
for(i in 1:dim(SPR)[1]) {
  y <- as.numeric(SPR[i,])
  x <- as.numeric(X.b)
  y <- y[-which(is.na(x))]
  x <- x[-which(is.na(x))]
  m <- lm(y ~ x)
  Y[i, ] <- predict(m, data.frame(x = X))
  r[i] <- cor(x, y)
}
dat.pred <- data.frame(X = X,
                       Y = apply(Y, 2, median),
                       Y.lo = apply(Y, 2, function(x) quantile(x, prob = 0.025, type = 8)),
                       Y.hi = apply(Y, 2, function(x) quantile(x, prob = 0.975, type = 8)))
rm(i, y, x, m, X, Y)

r.CanHt <- str_c(median(r) %>% round(digits = 2),
                  " (",
                  quantile(r, prob = 0.025, type = 8) %>% round(digits = 2),
                  ",",
                  quantile(r, prob = 0.975, type = 8) %>% round(digits = 2),
                  ")")

dat.SR <- data.frame(X = (X.b %>% as.numeric)) %>%
  mutate(#X = X + runif(length(Trt.b), -0.4, 0.4),
    Y = apply(SPR, 2, median) %>% as.numeric,
    Y.lo = apply(SPR, 2,function(x) quantile(x,prob=0.025,type=8)) %>%
      as.numeric,
    Y.hi = apply(SPR, 2, function(x) quantile(x,prob=0.975,type=8)) %>%
      as.numeric) %>%
  filter(!is.na(X))

p.CanHt <- ggplot(data = dat.SR, aes(x = X, y = Y)) + 
  geom_point(alpha = 0.3) +
  geom_errorbar(aes(ymin = Y.lo, ymax = Y.hi), width = 0, alpha = 0.3) +
  geom_ribbon(data = dat.pred, aes(x = X, ymin = Y.lo, ymax = Y.hi), alpha = 0.2) +
  geom_line(data = dat.pred, aes(x = X, y = Y), colour = "blue", size = 1.5) + 
  labs(x = "Canopy height (CanHt)", y = NULL) +
  annotate("text", x = 3, y = 25, label = str_c("r = ", r.CanHt), hjust = 0, size = 6) # +
#theme(axis.title.x=element_text(size=40)) +
#theme(axis.text.x=element_text(size=30))

# Number of snags #
X.b <- Cov[, "NumSnags"] # Point-level values

X <- seq(min(X.b, na.rm = T), max(X.b, na.rm = T), length.out = 20)
Y <- matrix(NA, nrow = dim(SPR)[1], ncol = length(X))
r <- numeric(length = dim(SPR)[1])
for(i in 1:dim(SPR)[1]) {
  y <- as.numeric(SPR[i,])
  x <- as.numeric(X.b)
  y <- y[-which(is.na(x))]
  x <- x[-which(is.na(x))]
  m <- lm(y ~ x)
  Y[i, ] <- predict(m, data.frame(x = X))
  r[i] <- cor(x, y)
}
dat.pred <- data.frame(X = X,
                       Y = apply(Y, 2, median),
                       Y.lo = apply(Y, 2, function(x) quantile(x, prob = 0.025, type = 8)),
                       Y.hi = apply(Y, 2, function(x) quantile(x, prob = 0.975, type = 8)))
rm(i, y, x, m, X, Y)

r.NSnag <- str_c(median(r) %>% round(digits = 2),
                 " (",
                 quantile(r, prob = 0.025, type = 8) %>% round(digits = 2),
                 ",",
                 quantile(r, prob = 0.975, type = 8) %>% round(digits = 2),
                 ")")

dat.SR <- data.frame(X = (X.b %>% as.numeric)) %>%
  mutate(#X = X + runif(length(Trt.b), -0.4, 0.4),
    Y = apply(SPR, 2, median) %>% as.numeric,
    Y.lo = apply(SPR, 2,function(x) quantile(x,prob=0.025,type=8)) %>%
      as.numeric,
    Y.hi = apply(SPR, 2, function(x) quantile(x,prob=0.975,type=8)) %>%
      as.numeric) %>%
  filter(!is.na(X))

p.NSnag <- ggplot(data = dat.SR, aes(x = X, y = Y)) + 
  geom_point(alpha = 0.3) +
  geom_errorbar(aes(ymin = Y.lo, ymax = Y.hi), width = 0, alpha = 0.3) +
  geom_ribbon(data = dat.pred, aes(x = X, ymin = Y.lo, ymax = Y.hi), alpha = 0.2) +
  geom_line(data = dat.pred, aes(x = X, y = Y), colour = "blue", size = 1.5) + 
  labs(x = "Number of snags (NSnag)", y = NULL) +
  annotate("text", x = 0, y = 25, label = str_c("r = ", r.NSnag), hjust = 0, size = 6) # +
#theme(axis.title.x=element_text(size=40)) +
#theme(axis.text.x=element_text(size=30))

# PIPO dominance #
X.b <- Cov[, "RCOV_PP"] # Point-level values

X <- seq(min(X.b, na.rm = T), max(X.b, na.rm = T), length.out = 20)
Y <- matrix(NA, nrow = dim(SPR)[1], ncol = length(X))
r <- numeric(length = dim(SPR)[1])
for(i in 1:dim(SPR)[1]) {
  y <- as.numeric(SPR[i,])
  x <- as.numeric(X.b)
  y <- y[-which(is.na(x))]
  x <- x[-which(is.na(x))]
  m <- lm(y ~ x)
  Y[i, ] <- predict(m, data.frame(x = X))
  r[i] <- cor(x, y)
}
dat.pred <- data.frame(X = X,
                       Y = apply(Y, 2, median),
                       Y.lo = apply(Y, 2, function(x) quantile(x, prob = 0.025, type = 8)),
                       Y.hi = apply(Y, 2, function(x) quantile(x, prob = 0.975, type = 8)))
rm(i, y, x, m, X, Y)

r.PIPO <- str_c(median(r) %>% round(digits = 2),
                 " (",
                 quantile(r, prob = 0.025, type = 8) %>% round(digits = 2),
                 ",",
                 quantile(r, prob = 0.975, type = 8) %>% round(digits = 2),
                 ")")

dat.SR <- data.frame(X = (X.b %>% as.numeric)) %>%
  mutate(#X = X + runif(length(Trt.b), -0.4, 0.4),
    Y = apply(SPR, 2, median) %>% as.numeric,
    Y.lo = apply(SPR, 2,function(x) quantile(x,prob=0.025,type=8)) %>%
      as.numeric,
    Y.hi = apply(SPR, 2, function(x) quantile(x,prob=0.975,type=8)) %>%
      as.numeric) %>%
  filter(!is.na(X))

p.PIPO <- ggplot(data = dat.SR, aes(x = X, y = Y)) + 
  geom_point(alpha = 0.3) +
  geom_errorbar(aes(ymin = Y.lo, ymax = Y.hi), width = 0, alpha = 0.3) +
  geom_ribbon(data = dat.pred, aes(x = X, ymin = Y.lo, ymax = Y.hi), alpha = 0.2) +
  geom_line(data = dat.pred, aes(x = X, y = Y), colour = "blue", size = 1.5) + 
  labs(x = "Ponderosa dominance (PIPO)", y = NULL) +
  annotate("text", x = 0, y = 25, label = str_c("r = ", r.PIPO), hjust = 0, size = 6) # +
#theme(axis.title.x=element_text(size=40)) +
#theme(axis.text.x=element_text(size=30))

# PSME dominance #
X.b <- Cov[, "RCOV_DF"] # Point-level values

X <- seq(min(X.b, na.rm = T), max(X.b, na.rm = T), length.out = 20)
Y <- matrix(NA, nrow = dim(SPR)[1], ncol = length(X))
r <- numeric(length = dim(SPR)[1])
for(i in 1:dim(SPR)[1]) {
  y <- as.numeric(SPR[i,])
  x <- as.numeric(X.b)
  y <- y[-which(is.na(x))]
  x <- x[-which(is.na(x))]
  m <- lm(y ~ x)
  Y[i, ] <- predict(m, data.frame(x = X))
  r[i] <- cor(x, y)
}
dat.pred <- data.frame(X = X,
                       Y = apply(Y, 2, median),
                       Y.lo = apply(Y, 2, function(x) quantile(x, prob = 0.025, type = 8)),
                       Y.hi = apply(Y, 2, function(x) quantile(x, prob = 0.975, type = 8)))
rm(i, y, x, m, X, Y)

r.PSME <- str_c(median(r) %>% round(digits = 2),
                " (",
                quantile(r, prob = 0.025, type = 8) %>% round(digits = 2),
                ",",
                quantile(r, prob = 0.975, type = 8) %>% round(digits = 2),
                ")")

dat.SR <- data.frame(X = (X.b %>% as.numeric)) %>%
  mutate(#X = X + runif(length(Trt.b), -0.4, 0.4),
    Y = apply(SPR, 2, median) %>% as.numeric,
    Y.lo = apply(SPR, 2,function(x) quantile(x,prob=0.025,type=8)) %>%
      as.numeric,
    Y.hi = apply(SPR, 2, function(x) quantile(x,prob=0.975,type=8)) %>%
      as.numeric) %>%
  filter(!is.na(X))

p.PSME <- ggplot(data = dat.SR, aes(x = X, y = Y)) + 
  geom_point(alpha = 0.3) +
  geom_errorbar(aes(ymin = Y.lo, ymax = Y.hi), width = 0, alpha = 0.3) +
  geom_ribbon(data = dat.pred, aes(x = X, ymin = Y.lo, ymax = Y.hi), alpha = 0.2) +
  geom_line(data = dat.pred, aes(x = X, y = Y), colour = "blue", size = 1.5) + 
  labs(x = "Douglas fir dominance (PSME)", y = NULL) +
  annotate("text", x = 0, y = 25, label = str_c("r = ", r.PSME), hjust = 0, size = 6) # +
#theme(axis.title.x=element_text(size=40)) +
#theme(axis.text.x=element_text(size=30))

# Aspen dominance #
X.b <- Cov[, "RCOV_AS"] # Point-level values

X <- seq(min(X.b, na.rm = T), max(X.b, na.rm = T), length.out = 20)
Y <- matrix(NA, nrow = dim(SPR)[1], ncol = length(X))
r <- numeric(length = dim(SPR)[1])
for(i in 1:dim(SPR)[1]) {
  y <- as.numeric(SPR[i,])
  x <- as.numeric(X.b)
  y <- y[-which(is.na(x))]
  x <- x[-which(is.na(x))]
  m <- lm(y ~ x)
  Y[i, ] <- predict(m, data.frame(x = X))
  r[i] <- cor(x, y)
}
dat.pred <- data.frame(X = X,
                       Y = apply(Y, 2, median),
                       Y.lo = apply(Y, 2, function(x) quantile(x, prob = 0.025, type = 8)),
                       Y.hi = apply(Y, 2, function(x) quantile(x, prob = 0.975, type = 8)))
rm(i, y, x, m, X, Y)

r.POTR5 <- str_c(median(r) %>% round(digits = 2),
                " (",
                quantile(r, prob = 0.025, type = 8) %>% round(digits = 2),
                ",",
                quantile(r, prob = 0.975, type = 8) %>% round(digits = 2),
                ")")

dat.SR <- data.frame(X = (X.b %>% as.numeric)) %>%
  mutate(#X = X + runif(length(Trt.b), -0.4, 0.4),
    Y = apply(SPR, 2, median) %>% as.numeric,
    Y.lo = apply(SPR, 2,function(x) quantile(x,prob=0.025,type=8)) %>%
      as.numeric,
    Y.hi = apply(SPR, 2, function(x) quantile(x,prob=0.975,type=8)) %>%
      as.numeric) %>%
  filter(!is.na(X))

p.POTR5 <- ggplot(data = dat.SR, aes(x = X, y = Y)) + 
  geom_point(alpha = 0.3) +
  geom_errorbar(aes(ymin = Y.lo, ymax = Y.hi), width = 0, alpha = 0.3) +
  geom_ribbon(data = dat.pred, aes(x = X, ymin = Y.lo, ymax = Y.hi), alpha = 0.2) +
  geom_line(data = dat.pred, aes(x = X, y = Y), colour = "blue", size = 1.5) + 
  labs(x = "Aspen dominance (POTR5)", y = NULL) +
  annotate("text", x = 0, y = 25, label = str_c("r = ", r.POTR5), hjust = 0, size = 6) # +
#theme(axis.title.x=element_text(size=40)) +
#theme(axis.text.x=element_text(size=30))

# Shrub volume #
X.b <- Cov[, "ShrubVol"] # Point-level values

X <- seq(min(X.b, na.rm = T), max(X.b, na.rm = T), length.out = 20)
Y <- matrix(NA, nrow = dim(SPR)[1], ncol = length(X))
r <- numeric(length = dim(SPR)[1])
for(i in 1:dim(SPR)[1]) {
  y <- as.numeric(SPR[i,])
  x <- as.numeric(X.b)
  y <- y[-which(is.na(x))]
  x <- x[-which(is.na(x))]
  m <- lm(y ~ x)
  Y[i, ] <- predict(m, data.frame(x = X))
  r[i] <- cor(x, y)
}
dat.pred <- data.frame(X = X,
                       Y = apply(Y, 2, median),
                       Y.lo = apply(Y, 2, function(x) quantile(x, prob = 0.025, type = 8)),
                       Y.hi = apply(Y, 2, function(x) quantile(x, prob = 0.975, type = 8)))
rm(i, y, x, m, X, Y)

r.ShrubVol <- str_c(median(r) %>% round(digits = 2),
                 " (",
                 quantile(r, prob = 0.025, type = 8) %>% round(digits = 2),
                 ",",
                 quantile(r, prob = 0.975, type = 8) %>% round(digits = 2),
                 ")")

dat.SR <- data.frame(X = (X.b %>% as.numeric)) %>%
  mutate(#X = X + runif(length(Trt.b), -0.4, 0.4),
    Y = apply(SPR, 2, median) %>% as.numeric,
    Y.lo = apply(SPR, 2,function(x) quantile(x,prob=0.025,type=8)) %>%
      as.numeric,
    Y.hi = apply(SPR, 2, function(x) quantile(x,prob=0.975,type=8)) %>%
      as.numeric) %>%
  filter(!is.na(X))

p.ShrubVol <- ggplot(data = dat.SR, aes(x = X, y = Y)) + 
  geom_point(alpha = 0.3) +
  geom_errorbar(aes(ymin = Y.lo, ymax = Y.hi), width = 0, alpha = 0.3) +
  geom_ribbon(data = dat.pred, aes(x = X, ymin = Y.lo, ymax = Y.hi), alpha = 0.2) +
  geom_line(data = dat.pred, aes(x = X, y = Y), colour = "blue", size = 1.5) + 
  labs(x = "Shrub volume (ShrubVol)", y = NULL) +
  annotate("text", x = 0, y = 25, label = str_c("r = ", r.ShrubVol), hjust = 0, size = 6) # +
#theme(axis.title.x=element_text(size=40)) +
#theme(axis.text.x=element_text(size=30))

# Ladder fuels dominance #
X.b <- Cov[, "RSCV_Ladder"] # Point-level values

X <- seq(min(X.b, na.rm = T), max(X.b, na.rm = T), length.out = 20)
Y <- matrix(NA, nrow = dim(SPR)[1], ncol = length(X))
r <- numeric(length = dim(SPR)[1])
for(i in 1:dim(SPR)[1]) {
  y <- as.numeric(SPR[i,])
  x <- as.numeric(X.b)
  y <- y[-which(is.na(x))]
  x <- x[-which(is.na(x))]
  m <- lm(y ~ x)
  Y[i, ] <- predict(m, data.frame(x = X))
  r[i] <- cor(x, y)
}
dat.pred <- data.frame(X = X,
                       Y = apply(Y, 2, median),
                       Y.lo = apply(Y, 2, function(x) quantile(x, prob = 0.025, type = 8)),
                       Y.hi = apply(Y, 2, function(x) quantile(x, prob = 0.975, type = 8)))
rm(i, y, x, m, X, Y)

r.LadFuel <- str_c(median(r) %>% round(digits = 2),
                    " (",
                    quantile(r, prob = 0.025, type = 8) %>% round(digits = 2),
                    ",",
                    quantile(r, prob = 0.975, type = 8) %>% round(digits = 2),
                    ")")

dat.SR <- data.frame(X = (X.b %>% as.numeric)) %>%
  mutate(#X = X + runif(length(Trt.b), -0.4, 0.4),
    Y = apply(SPR, 2, median) %>% as.numeric,
    Y.lo = apply(SPR, 2,function(x) quantile(x,prob=0.025,type=8)) %>%
      as.numeric,
    Y.hi = apply(SPR, 2, function(x) quantile(x,prob=0.975,type=8)) %>%
      as.numeric) %>%
  filter(!is.na(X))

p.LadFuel <- ggplot(data = dat.SR, aes(x = X, y = Y)) + 
  geom_point(alpha = 0.3) +
  geom_errorbar(aes(ymin = Y.lo, ymax = Y.hi), width = 0, alpha = 0.3) +
  geom_ribbon(data = dat.pred, aes(x = X, ymin = Y.lo, ymax = Y.hi), alpha = 0.2) +
  geom_line(data = dat.pred, aes(x = X, y = Y), colour = "blue", size = 1.5) + 
  labs(x = "Ladder fuel dominance (LadFuel)", y = NULL) +
  annotate("text", x = 0, y = 25, label = str_c("r = ", r.LadFuel), hjust = 0, size = 6) # +
#theme(axis.title.x=element_text(size=40)) +
#theme(axis.text.x=element_text(size=30))

# Herbaceous volume #
X.b <- Cov[, "HerbGrassVol"] # Point-level values

X <- seq(min(X.b, na.rm = T), max(X.b, na.rm = T), length.out = 20)
Y <- matrix(NA, nrow = dim(SPR)[1], ncol = length(X))
r <- numeric(length = dim(SPR)[1])
for(i in 1:dim(SPR)[1]) {
  y <- as.numeric(SPR[i,])
  x <- as.numeric(X.b)
  y <- y[-which(is.na(x))]
  x <- x[-which(is.na(x))]
  m <- lm(y ~ x)
  Y[i, ] <- predict(m, data.frame(x = X))
  r[i] <- cor(x, y)
}
dat.pred <- data.frame(X = X,
                       Y = apply(Y, 2, median),
                       Y.lo = apply(Y, 2, function(x) quantile(x, prob = 0.025, type = 8)),
                       Y.hi = apply(Y, 2, function(x) quantile(x, prob = 0.975, type = 8)))
rm(i, y, x, m, X, Y)

r.Herb <- str_c(median(r) %>% round(digits = 2),
                   " (",
                   quantile(r, prob = 0.025, type = 8) %>% round(digits = 2),
                   ",",
                   quantile(r, prob = 0.975, type = 8) %>% round(digits = 2),
                   ")")

dat.SR <- data.frame(X = (X.b %>% as.numeric)) %>%
  mutate(#X = X + runif(length(Trt.b), -0.4, 0.4),
    Y = apply(SPR, 2, median) %>% as.numeric,
    Y.lo = apply(SPR, 2,function(x) quantile(x,prob=0.025,type=8)) %>%
      as.numeric,
    Y.hi = apply(SPR, 2, function(x) quantile(x,prob=0.975,type=8)) %>%
      as.numeric) %>%
  filter(!is.na(X))

p.Herb <- ggplot(data = dat.SR, aes(x = X, y = Y)) + 
  geom_point(alpha = 0.3) +
  geom_errorbar(aes(ymin = Y.lo, ymax = Y.hi), width = 0, alpha = 0.3) +
  geom_ribbon(data = dat.pred, aes(x = X, ymin = Y.lo, ymax = Y.hi), alpha = 0.2) +
  geom_line(data = dat.pred, aes(x = X, y = Y), colour = "blue", size = 1.5) + 
  labs(x = "Herbaceous volume (Herb)", y = NULL) +
  annotate("text", x = 0, y = 25, label = str_c("r = ", r.Herb), hjust = 0, size = 6) # +
#theme(axis.title.x=element_text(size=40)) +
#theme(axis.text.x=element_text(size=30))

p <- ggdraw() + 
  draw_plot(p.CanCov, x = 0.03, y = 0.6667, width = 0.3233333, height = 0.3333) +
  draw_plot(p.CanHt, x = 0.3533333, y = 0.6667, width = 0.3233333, height = 0.3333) +
  draw_plot(p.NSnag, x = 0.6766667, y = 0.6667, width = 0.3233333, height = 0.3333) +
  draw_plot(p.PIPO, x = 0.03, y = 0.3333, width = 0.3233333, height = 0.3333) +
  draw_plot(p.PSME, x = 0.3533333, y = 0.3333, width = 0.3233333, height = 0.3333) +
  draw_plot(p.POTR5, x = 0.6766667, y = 0.3333, width = 0.3233333, height = 0.3333) +
  draw_plot(p.ShrubVol, x = 0.03, y = 0, width = 0.3233333, height = 0.3333) +
  draw_plot(p.LadFuel, x = 0.3533333, y = 0, width = 0.3233333, height = 0.3333) +
  draw_plot(p.Herb, x = 0.6766667, y = 0, width = 0.3233333, height = 0.3333) +
  draw_plot_label("Species richness", x = 0, y = 0.4, size = 20, angle = 90, hjust = 0)

save_plot("Plot_richness_veg.tiff", p, ncol = 3, nrow = 3, dpi = 200)


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
