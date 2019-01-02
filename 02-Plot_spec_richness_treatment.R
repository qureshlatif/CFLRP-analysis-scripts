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
PctTrt.d <- matrix(NA, nrow = max(gridID), ncol = max(yearID))
PctTrt.d[landscape_data %>% filter(YearInd == 1) %>% pull(gridIndex), 1] <-
  landscape_data %>% filter(YearInd == 1) %>% pull(PctTrt)
PctTrt.d[landscape_data %>% filter(YearInd == 2) %>% pull(gridIndex), 2] <-
  landscape_data %>% filter(YearInd == 2) %>% pull(PctTrt)
PctTrt.d[landscape_data %>% filter(YearInd == 3) %>% pull(gridIndex), 3] <-
  landscape_data %>% filter(YearInd == 3) %>% pull(PctTrt)

SPR <- mod$sims.list$SR.grid
  # Derive posterior samples for plotting spp richness trend #
X <- 0:100
Y <- matrix(NA, nrow = dim(SPR)[1], ncol = length(X))
m1_B1 <- m1_B2 <- r <- numeric(length = dim(SPR)[1])
# Y2 <- matrix(NA, nrow = dim(SPR)[1], ncol = length(X))
# m2_B1 <- m2_B2 <- m2_B3 <- r <- numeric(length = dim(SPR)[1])
for(i in 1:dim(SPR)[1]) {
  y <- as.numeric(SPR[i,,])
  x <- as.numeric(PctTrt.d)
  m <- lm(y ~ x + I(x^2))
  m1_B1[i] <- m$coefficients[2]
  m1_B2[i] <- m$coefficients[3]
  Y[i, ] <- predict(m, data.frame(x = X))
  x0 <- (x > 0)*1
  # m <- lm(y ~ x0 + x + I(x^2))
  # m2_B1[i] <- m$coefficients[2]
  # m2_B2[i] <- m$coefficients[3]
  # m2_B3[i] <- m$coefficients[4]
  # Y2[i, ] <- predict(m, data.frame(x0 = (X > 0)*1, x = X))
  r[i] <- cor(x, y, use = "complete")
}
dat.pred <- data.frame(X = X,
                       Y = apply(Y, 2, median),
                       Y.lo = apply(Y, 2, function(x) quantile(x, prob = 0.025, type = 8)),
                       Y.hi = apply(Y, 2, function(x) quantile(x, prob = 0.975, type = 8))) #,
                       # Y2 = apply(Y2, 2, median),
                       # Y2.lo = apply(Y2, 2, function(x) quantile(x, prob = 0.025, type = 8)),
                       # Y2.hi = apply(Y2, 2, function(x) quantile(x, prob = 0.975, type = 8)))
r <- str_c(median(r) %>% round(digits = 2),
                   " (",
                   quantile(r, prob = 0.025, type = 8) %>% round(digits = 2),
                   ",",
                   quantile(r, prob = 0.975, type = 8) %>% round(digits = 2),
                   ")")
m1_B1 <- str_c(median(m1_B1) %>% round(digits = 2),
           " (",
           quantile(m1_B1, prob = 0.025, type = 8) %>% round(digits = 2),
           ",",
           quantile(m1_B1, prob = 0.975, type = 8) %>% round(digits = 2),
           ")")
m1_B2 <- str_c(median(m1_B2) %>% round(digits = 4),
               " (",
               quantile(m1_B2, prob = 0.025, type = 8) %>% round(digits = 4),
               ",",
               quantile(m1_B2, prob = 0.975, type = 8) %>% round(digits = 4),
               ")")
# m2_B1 <- str_c(median(m2_B1) %>% round(digits = 2),
#                " (",
#                quantile(m2_B1, prob = 0.025, type = 8) %>% round(digits = 2),
#                ",",
#                quantile(m2_B1, prob = 0.975, type = 8) %>% round(digits = 2),
#                ")")
# m2_B2 <- str_c(median(m2_B2) %>% round(digits = 2),
#                " (",
#                quantile(m2_B2, prob = 0.025, type = 8) %>% round(digits = 2),
#                ",",
#                quantile(m2_B2, prob = 0.975, type = 8) %>% round(digits = 2),
#                ")")
# m2_B3 <- str_c(median(m2_B3) %>% round(digits = 4),
#                " (",
#                quantile(m2_B3, prob = 0.025, type = 8) %>% round(digits = 4),
#                ",",
#                quantile(m2_B3, prob = 0.975, type = 8) %>% round(digits = 4),
#                ")")
rm(i, y, x, m, X, Y)

dat.SR <- data.frame(X = (PctTrt.d %>% as.numeric)) %>%
  mutate(Y = apply(SPR,c(2,3),median) %>% as.numeric,
         Y.lo = apply(SPR,c(2,3),function(x) quantile(x,prob=0.025,type=8)) %>%
           as.numeric,
         Y.hi = apply(SPR,c(2,3),function(x) quantile(x,prob=0.975,type=8)) %>%
           as.numeric) %>%
  filter(!is.na(X))

p.grid <- ggplot(data = dat.SR, aes(x = X, y = Y)) + 
  geom_point(alpha = 0.3) +
  geom_errorbar(aes(ymin = Y.lo, ymax = Y.hi), width = 0.1, alpha = 0.3) +
  geom_ribbon(data = dat.pred, aes(x = X, ymin = Y.lo, ymax = Y.hi), alpha = 0.2) +
  geom_line(data = dat.pred, aes(x = X, y = Y), colour = "blue", size = 1.5) + 
  labs(x= "Percent treated (grid)", y = NULL) #+
#  theme(axis.title.x=element_text(size=40)) +
#  theme(axis.title.y=element_text(size=40)) +
#  theme(axis.text.x=element_text(size=30)) +
#  theme(axis.text.y=element_text(size=30))

# p2 <- ggplot(data = dat.SR, aes(x = X, y = Y)) + 
#   geom_point(alpha = 0.3) +
#   geom_errorbar(aes(ymin = Y.lo, ymax = Y.hi), width = 0.1, alpha = 0.3) +
#   geom_point(data = (dat.pred %>% filter(X == 0)), aes(x = X, y = Y2), color = "blue", size = 3) +
#   geom_errorbar(data = (dat.pred %>% filter(X == 0)), aes(ymin = Y2.lo, ymax = Y2.hi), width = 5, color = "blue", size = 2) +
#   geom_ribbon(data = (dat.pred %>% filter(X > 0)), aes(x = X, ymin = Y2.lo, ymax = Y2.hi), alpha = 0.2) +
#   geom_line(data = (dat.pred %>% filter(X > 0)), aes(x = X, y = Y2), colour = "blue", size = 1.5) + 
#   labs(x= "Percent treated", y = "Species Richness") +
#   theme(axis.title.x=element_text(size=40)) +
#   theme(axis.title.y=element_text(size=40)) +
#   theme(axis.text.x=element_text(size=30)) +
#   theme(axis.text.y=element_text(size=30))

## Point level ##

Trt.b <- Cov[, "Trt_stat"] # Point-level values
SPR <- mod$sims.list$SR.point

# Derive posterior samples for plotting spp richness trend #
X <- c(0, 1)
Y <- matrix(NA, nrow = dim(SPR)[1], ncol = length(X))
B_pnt <- r.pnt <- numeric(length = dim(SPR)[1])
for(i in 1:dim(SPR)[1]) {
  y <- as.numeric(SPR[i,])
  x <- as.numeric(Trt.b)
  m <- lm(y ~ x)
  Y[i, ] <- predict(m, data.frame(x = X))
  B_pnt[i] <- m$coefficients[2]
  r.pnt[i] <- cor(x, y, use = "complete")
}
dat.pred <- data.frame(X = X,
                       Y = apply(Y, 2, median),
                       Y.lo = apply(Y, 2, function(x) quantile(x, prob = 0.025, type = 8)),
                       Y.hi = apply(Y, 2, function(x) quantile(x, prob = 0.975, type = 8)))
r.pnt <- str_c(median(r.pnt) %>% round(digits = 2),
           " (",
           quantile(r.pnt, prob = 0.025, type = 8) %>% round(digits = 2),
           ",",
           quantile(r.pnt, prob = 0.975, type = 8) %>% round(digits = 2),
           ")")
B_pnt <- str_c(median(B_pnt) %>% round(digits = 2),
               " (",
               quantile(B_pnt, prob = 0.025, type = 8) %>% round(digits = 2),
               ",",
               quantile(B_pnt, prob = 0.975, type = 8) %>% round(digits = 2),
               ")")
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
p.point <- ggplot(data = dat.SR, aes(x = X, y = Y)) + 
  geom_point(aes(x = X + jitter), alpha = 0.3) + 
  geom_errorbar(aes(x = X + jitter, ymin = Y.lo, ymax = Y.hi), width = 0, alpha = 0.3) +
  geom_errorbar(data = dat.pred, aes(x = X, ymin = Y.lo, ymax = Y.hi), width = 0.3, size = 1, color = "blue") +
  geom_point(data = dat.pred, aes(x = X, y = Y), size = 3, color = "blue") + 
  scale_x_continuous(breaks = c(0, 1), labels = c("Untreated", "Treated")) +
  labs(x = "Treatment status (point)", y = NULL) #+
#  theme(axis.title.y=element_text(size=40)) +
#  theme(axis.text.x=element_text(size=30)) +
#  theme(axis.text.y=element_text(size=30))

p <- ggdraw() + 
  draw_plot(p.grid, x = 0.03, y = 0, width = 0.485, height = 1) +
  draw_plot(p.point, x = 0.515, y = 0, width = 0.485, height = 1) +
  draw_plot_label("Species richness", x = 0, y = 0.3, size = 20, angle = 90, hjust = 0)

save_plot("Plot_richness_treatment.tiff", p, ncol = 2, nrow = 1, dpi = 200)
