library(QSLpersonal)
library(jagsUI)
library(stringr)
library(dplyr)
library(R.utils)
library(ggplot2)
library(cowplot)

setwd("C:/Users/Quresh.Latif/files/projects/FS/CFLRP")
load("Data_compiled_abundance.RData")

#_______Script inputs_______#
spp <- "WAVI"
ind.spp <- which(Spp == spp)
mod <- loadObject("abund_models/mod_WAVI_abundance_treatment_global")
#samp.acres <- sum(area.band.list[[ind.spp]]) * 0.000247105
#___________________________#


#_______ Functions _________#
dat.pred.grid.fn <- function(ind.spp, PctTrt, mod) {
  x <- seq(min(PctTrt), max(PctTrt), length.out = 11)
  z <- (x - mean(PctTrt, na.rm = T)) / sd(PctTrt, na.rm = T)
  dat.pred <- data.frame(x = x, z = z, B0 = NA, B0.lo = NA, B0.hi = NA)
  rm(x, z)
  
  ind.spp <- which(Spp == spp)
  B0.mean <- mod$sims.list$beta0.mean
  bd.PctTrt <- mod$sims.list$bd.PctTrt
  for(i in 1:nrow(dat.pred)) {
    B0 <- exp(B0.mean + bd.PctTrt*dat.pred$z[i])
    dat.pred$B0[i] <- median(B0)
    dat.pred$B0.lo[i] <- quantile(B0, prob = 0.025, type = 8)
    dat.pred$B0.hi[i] <- quantile(B0, prob = 0.975, type = 8)
  }
  return(dat.pred)
}

grd.obs.fn <- function(ind.spp, landscape_data, Y.dist, scaling_factor) {
  grd.obs <- landscape_data %>%
    select(gridIndex, YearInd, PctTrt) %>%
    left_join(as.data.frame(Cov) %>%
                select(gridIndex, YearInd) %>%
                mutate(Y = Y.dist) %>%
                dplyr::group_by(gridIndex, YearInd) %>%
                summarise(Y = mean(Y)), by = c("gridIndex", "YearInd"))
  grd.obs.sum <- grd.obs %>%
    select(PctTrt, Y) %>%
    rename(x = PctTrt) %>%
    mutate(x_class = cut(x, breaks = c(-Inf, 0.0001, 33, 67, Inf), labels = 1:4)) %>%
    group_by(x_class) %>%
    summarise(x = mean(x), Y = mean(Y) / scaling_factor)
  return(grd.obs.sum)
}

dat.pred.pnt.fn <- function(ind.spp, trt, yst, mod) {
  x.trt <- c(0, rep(1, 10))
  z.trt <- (x.trt - mean(trt, na.rm = T)) / sd(trt, na.rm = T)
  x.yst <- c(0, min(yst, na.rm = T):max(yst, na.rm = T))
  z.yst <- (x.yst - mean(yst, na.rm = T)) / sd(yst, na.rm = T)
  dat.pred <- data.frame(x.trt = x.trt, z.trt = z.trt, x.yst = x.yst, z.yst = z.yst,
                         lambda = NA, lambda.lo = NA, lambda.hi = NA) %>%
    mutate(x = 1:11)
  rm(x.trt, z.trt, x.yst, z.yst)
  
  B0 <- mod$sims.list$beta0.mean
  bl.trt <- mod$sims.list$bl.trt
  bl.YST <- mod$sims.list$bl.YST
  for(i in 1:nrow(dat.pred)) {
    lambda <- exp(B0 + bl.trt*dat.pred$z.trt[i] + bl.YST*dat.pred$z.yst[i])
    dat.pred$lambda[i] <- median(lambda)
    dat.pred$lambda.lo[i] <- quantile(lambda, prob = 0.025, type = 8)
    dat.pred$lambda.hi[i] <- quantile(lambda, prob = 0.975, type = 8)
  }
  return(dat.pred)
}

pnt.obs.fn <- function(ind.spp, Cov, Y.dist, scaling_factor) {
  pnt.obs <- Cov %>% as.data.frame %>%
    select(gridIndex, YearInd, Trt_stat, Trt_time) %>%
    mutate(Y = Y.dist) %>%
    mutate(x = replace(Trt_time, which(is.na(Trt_time)), 0)) %>%
    mutate(x_class = cut(x, breaks = c(-Inf, 0.0001, 2, 4, 6, 8, Inf), labels = 1:6)) %>%
    group_by(x_class) %>%
    summarise(x = mean(x), Y = mean(Y) / scaling_factor) %>%
    ungroup
  return(pnt.obs)
}
#___________________________#

# Grid level #
dat.pred <- dat.pred.grid.fn(ind.spp, landscape_data$PctTrt, mod)
grd.obs <- grd.obs.fn(ind.spp, landscape_data, eval(as.name(str_c("Y.", spp, ".dist"))), 0.5)
p <- ggplot(dat.pred, aes(x = x, y = B0)) +
  geom_ribbon(aes(ymin = B0.lo, ymax = B0.hi), alpha = 0.4) +
  geom_line(size = 2) +
  xlab("Percent landscape treated") + ylab("Mean grid-level density (per 100 ac)")
p <- p + geom_point(data = grd.obs, aes(x = x, y = Y), size = 3, color = "red")
p

# Point level #
dat.pred <- dat.pred.pnt.fn(ind.spp, Cov[, "Trt_stat"], Cov[, "Trt_time"], mod)
pnt.obs <- pnt.obs.fn(ind.spp, Cov, eval(as.name(str_c("Y.", spp, ".dist"))), 0.3)
p <- ggplot(dat.pred, aes(x = x, y = lambda)) +
  geom_ribbon(data = dat.pred %>% filter(x.trt == 1),
              aes(ymin = lambda.lo, ymax = lambda.hi), alpha = 0.4) +
  geom_line(data = dat.pred %>% filter(x.trt == 1), size = 2) +
  geom_errorbar(data = dat.pred %>% filter(x.trt == 0),
                aes(x = x, ymin = lambda.lo, ymax = lambda.hi), width = 0.2) +
  geom_point(data = dat.pred %>% filter(x.trt == 0), size = 3, shape = 15) +
  geom_point(data = pnt.obs, aes(x = x + 1, y = Y), size = 3, color = "red") +
  scale_x_continuous(breaks = c(1, 3, 5, 7, 9, 11), labels = c("pre-trt", "trt yr 2", "trt yr 4", "trt yr 6", "trt yr 8", "trt yr 10")) +
  xlab("Treatment status") + ylab("Point-level density (per 100 ac)")
p
