library(QSLpersonal)
library(jagsUI)
library(stringr)
library(dplyr)
library(R.utils)
library(ggplot2)
library(cowplot)

setwd("C:/Users/Quresh.Latif/files/projects/FS/CFLRP")
load("Data_compiled.RData")

#_______Script inputs_______#
plot.obs <- T # Set to T if apparent occupancy is desired in plots.
mod <- loadObject("mod_treatment_d0yr")
#___________________________#


#_______ Functions _________#
dat.pred.grid.fn <- function(ind.spp, PctTrt, mod) {
  x <- seq(min(PctTrt), max(PctTrt), length.out = 11)
  z <- (x - mean(PctTrt, na.rm = T)) / sd(PctTrt, na.rm = T)
  dat.pred <- data.frame(x = x, z = z, psi = NA, psi.lo = NA, psi.hi = NA)
  rm(x, z)
  
  ind.spp <- which(spp.list == spp)
  B0 <- apply(mod$sims.list$d0[, ind.spp, ], 1, mean)
  B1 <- mod$sims.list$bd.ptrt[, ind.spp]
  for(i in 1:nrow(dat.pred)) {
    psi <- expit(B0 + B1*dat.pred$z[i]) # Rescale to number per 100 acres
    dat.pred$psi[i] <- median(psi)
    dat.pred$psi.lo[i] <- quantile(psi, prob = 0.025, type = 8)
    dat.pred$psi.hi[i] <- quantile(psi, prob = 0.975, type = 8)
  }
  return(dat.pred)
}

grd.obs.fn <- function(ind.spp, landscape_data, Y.mat) {
  grd.obs <- landscape_data %>%
    select(gridIndex, YearInd, PctTrt) %>%
    left_join(as.data.frame(Cov) %>%
                select(gridIndex, YearInd) %>%
                mutate(Y = Y.mat[, ind.spp]) %>%
                dplyr::group_by(gridIndex, YearInd) %>%
                summarise(Y = any(Y > 0)*1), by = c("gridIndex", "YearInd"))
  grd.obs.sum <- grd.obs %>%
    select(PctTrt, Y) %>%
    rename(x = PctTrt) %>%
    mutate(x_class = cut(x, breaks = c(-Inf, 0.0001, 33, 67, Inf), labels = 1:4)) %>%
    group_by(x_class) %>%
    summarise(x = mean(x), Y = mean(Y))
  out <- list(raw = grd.obs, sum = grd.obs.sum)
  return(out)
}

dat.pred.pnt.fn <- function(ind.spp, trt, yst, mod) {
  x.trt <- c(0, rep(1, 10))
  z.trt <- (x.trt - mean(trt, na.rm = T)) / sd(trt, na.rm = T)
  x.yst <- c(0, min(yst, na.rm = T):max(yst, na.rm = T))
  z.yst <- (x.yst - mean(yst, na.rm = T)) / sd(yst, na.rm = T)
  dat.pred <- data.frame(x.trt = x.trt, z.trt = z.trt, x.yst = x.yst, z.yst = z.yst,
                         psi = NA, psi.lo = NA, psi.hi = NA) %>%
    mutate(x = 1:11)
  rm(x.trt, z.trt, x.yst, z.yst)
  
  B0 <- mod$sims.list$b0[, ind.spp]
  B.trt <- mod$sims.list$bb.trt[, ind.spp]
  B.yst <- mod$sims.list$bb.YST[, ind.spp]
  for(i in 1:nrow(dat.pred)) {
    psi <- expit(B0 + B.trt*dat.pred$z.trt[i] + B.yst*dat.pred$z.yst[i]) # Rescale to number per 100 acres
    dat.pred$psi[i] <- median(psi)
    dat.pred$psi.lo[i] <- quantile(psi, prob = 0.025, type = 8)
    dat.pred$psi.hi[i] <- quantile(psi, prob = 0.975, type = 8)
  }
  return(dat.pred)
}

pnt.obs.fn <- function(ind.spp, Cov, Y.mat, grd.obs) {
  pnt.obs <- Cov %>% tbl_df %>%
    select(gridIndex, YearInd, Trt_stat, Trt_time) %>%
    mutate(Y = Y.mat[, ind.spp]) %>%
    left_join(grd.obs %>% rename(Y.grd = Y), by = c("gridIndex", "YearInd")) %>%
    filter(Y.grd == 1) %>%
    select(-Y.grd) %>%
    mutate(x = replace(Trt_time, which(is.na(Trt_time)), 0)) %>%
    mutate(x_class = cut(x, breaks = c(-Inf, 0.0001, 3.5, 6.5, Inf), labels = 1:4)) %>%
    group_by(x_class) %>%
    summarise(x = mean(x), Y = mean(Y))
  return(pnt.obs)
}
#___________________________#

##____ Grid level ____##
signif.ptrt <- (apply(mod$sims.list$bd.ptrt, 2, function(x) quantile(x, prob = 0.025, type = 8)) > 0 |
                       apply(mod$sims.list$bd.ptrt, 2, function(x) quantile(x, prob = 0.975, type = 8)) < 0)
spp.plot <- spp.list[which(signif.ptrt)]
spp.plot <- sort(spp.plot)

#____WISA____#
# Plot grid-level relationship #
#spp <- "CONI"
#ind.spp <- which(spp.list == spp)

for(i in 1:length(spp.plot)) {
  spp <- spp.plot[i]
  ind.spp <- which(spp.list == spp)
  
  dat.pred <- dat.pred.grid.fn(ind.spp, landscape_data$PctTrt, mod)
  grd.obs <- grd.obs.fn(ind.spp, landscape_data, Y.mat)
  p.grid <- ggplot(dat.pred, aes(x = x, y = psi)) +
    geom_ribbon(aes(ymin = psi.lo, ymax = psi.hi), alpha = 0.4) +
    geom_line(size = 2) +
    ylim(0, 1) +
    xlab(NULL) + ylab(NULL)
  if(plot.obs == T) p.grid <- p.grid + geom_point(data = grd.obs$sum, aes(x = x, y = Y), shape = 1)

  # Plot point-level relationships #
  dat.pred <- dat.pred.pnt.fn(ind.spp, Cov[, "Trt_stat"], Cov[, "Trt_time"], mod)
  pnt.obs <- pnt.obs.fn(ind.spp, Cov, Y.mat, grd.obs$raw)
  p.pnt <- ggplot(dat.pred, aes(x = x, y = psi)) +
    geom_ribbon(data = dat.pred %>% filter(x.trt == 1),
                aes(ymin = psi.lo, ymax = psi.hi), alpha = 0.4) +
    geom_line(data = dat.pred %>% filter(x.trt == 1), size = 2) +
    geom_errorbar(data = dat.pred %>% filter(x.trt == 0),
                  aes(x = x, ymin = psi.lo, ymax = psi.hi), width = 0.2) +
    geom_point(data = dat.pred %>% filter(x.trt == 0), size = 3, shape = 15) +
    scale_x_continuous(breaks = c(1, 3, 5, 7, 9, 11), labels = c("pre-trt", "trt yr 2", "trt yr 4", "trt yr 6", "trt yr 8", "trt yr 10")) +
    ylim(0, 1) +
    xlab(NULL) + ylab(NULL)
  if(plot.obs == T) p.pnt <- p.pnt + geom_point(data = pnt.obs, aes(x = x + 1, y = Y), shape = 1)

  p.spp <- ggdraw() +
    draw_plot(p.grid, x = 0, y = 0, width = 0.475, height = 0.95) +
    draw_plot(p.pnt, x = 0.525, y = 0, width = 0.475, height = 0.95) +
    draw_plot_label(c(spp, spp), x = c(0.23, 0.77), y = c(1, 1))
  assign(str_c("p.spp", i), p.spp)
}

p <- ggdraw() +
  draw_plot(p.spp1, x = 0.05, y = 0.05, width = 0.2833333, height = 0.08636364) +
  draw_plot(p.spp2, x = 0.05, y = 0.1363636, width = 0.2833333, height = 0.08636364) +
  draw_plot(p.spp3, x = 0.05, y = 0.2227273, width = 0.2833333, height = 0.08636364) +
  draw_plot(p.spp4, x = 0.05, y = 0.3090909, width = 0.2833333, height = 0.08636364) +
  draw_plot(p.spp5, x = 0.05, y = 0.3954546, width = 0.2833333, height = 0.08636364) +
  draw_plot(p.spp6, x = 0.05, y = 0.4818182, width = 0.2833333, height = 0.08636364) +
  draw_plot(p.spp7, x = 0.05, y = 0.5681818, width = 0.2833333, height = 0.08636364) +
  draw_plot(p.spp8, x = 0.05, y = 0.6545455, width = 0.2833333, height = 0.08636364) +
  draw_plot(p.spp9, x = 0.05, y = 0.7409091, width = 0.2833333, height = 0.08636364) +
  draw_plot(p.spp10, x = 0.05, y = 0.8272728, width = 0.2833333, height = 0.08636364) +
  draw_plot(p.spp11, x = 0.05, y = 0.9136364, width = 0.2833333, height = 0.08636364) +
  draw_plot(p.spp12, x = 0.3833333, y = 0.05, width = 0.2833333, height = 0.08636364) +
  draw_plot(p.spp13, x = 0.3833333, y = 0.1363636, width = 0.2833333, height = 0.08636364) +
  draw_plot(p.spp14, x = 0.3833333, y = 0.2227273, width = 0.2833333, height = 0.08636364) +
  draw_plot(p.spp15, x = 0.3833333, y = 0.3090909, width = 0.2833333, height = 0.08636364) +
  draw_plot(p.spp16, x = 0.3833333, y = 0.3954546, width = 0.2833333, height = 0.08636364) +
  draw_plot(p.spp17, x = 0.3833333, y = 0.4818182, width = 0.2833333, height = 0.08636364) +
  draw_plot(p.spp18, x = 0.3833333, y = 0.5681818, width = 0.2833333, height = 0.08636364) +
  draw_plot(p.spp19, x = 0.3833333, y = 0.6545455, width = 0.2833333, height = 0.08636364) +
  draw_plot(p.spp20, x = 0.3833333, y = 0.7409091, width = 0.2833333, height = 0.08636364) +
  draw_plot(p.spp21, x = 0.3833333, y = 0.8272728, width = 0.2833333, height = 0.08636364) +
  draw_plot(p.spp22, x = 0.7166666, y = 0.05, width = 0.2833333, height = 0.08636364) +
  draw_plot(p.spp23, x = 0.7166666, y = 0.1363636, width = 0.2833333, height = 0.08636364) +
  draw_plot(p.spp24, x = 0.7166666, y = 0.2227273, width = 0.2833333, height = 0.08636364) +
  draw_plot(p.spp25, x = 0.7166666, y = 0.3090909, width = 0.2833333, height = 0.08636364) +
  draw_plot(p.spp26, x = 0.7166666, y = 0.3954546, width = 0.2833333, height = 0.08636364) +
  draw_plot(p.spp27, x = 0.7166666, y = 0.4818182, width = 0.2833333, height = 0.08636364) +
  draw_plot(p.spp28, x = 0.7166666, y = 0.5681818, width = 0.2833333, height = 0.08636364) +
  draw_plot(p.spp29, x = 0.7166666, y = 0.6545455, width = 0.2833333, height = 0.08636364) +
  draw_plot(p.spp30, x = 0.7166666, y = 0.7409091, width = 0.2833333, height = 0.08636364) +
  draw_plot(p.spp31, x = 0.7166666, y = 0.8272728, width = 0.2833333, height = 0.08636364)# +
#  draw_plot_label(c(spp, spp), x = c(0.23, 0.77), y = c(1, 1))

save_plot("Plot_spp_trt_rels_inspect_GOF.tiff", p, ncol = 8, nrow = 11, dpi = 200)


### Point level ###
signif.ptrt <- (apply(mod$sims.list$bd.ptrt, 2, function(x) quantile(x, prob = 0.025, type = 8)) > 0 |
                  apply(mod$sims.list$bd.ptrt, 2, function(x) quantile(x, prob = 0.975, type = 8)) < 0)
signif.trt <- (apply(mod$sims.list$bb.trt, 2, function(x) quantile(x, prob = 0.025, type = 8)) > 0 |
                 apply(mod$sims.list$bb.trt, 2, function(x) quantile(x, prob = 0.975, type = 8)) < 0)
signif.yst <- (apply(mod$sims.list$bb.YST, 2, function(x) quantile(x, prob = 0.025, type = 8)) > 0 |
                 apply(mod$sims.list$bb.YST, 2, function(x) quantile(x, prob = 0.975, type = 8)) < 0)
spp.plot <- spp.list[which(signif.ptrt | signif.trt | signif.yst)]
spp.plot <- sort(spp.plot)

for(i in 1:length(spp.plot)) {
  spp <- spp.plot[i]
  ind.spp <- which(spp.list == spp)
  
  dat.pred <- dat.pred.grid.fn(ind.spp, landscape_data$PctTrt, mod)
  grd.obs <- grd.obs.fn(ind.spp, landscape_data, Y.mat)
  p.grid <- ggplot(dat.pred, aes(x = x, y = psi)) +
    geom_ribbon(aes(ymin = psi.lo, ymax = psi.hi), alpha = 0.4) +
    geom_line(size = 2) +
    ylim(0, 1) +
    xlab(NULL) + ylab(NULL)
  if(plot.obs == T) p.grid <- p.grid + geom_point(data = grd.obs$sum, aes(x = x, y = Y), shape = 1)
  
  # Plot point-level relationships #
  dat.pred <- dat.pred.pnt.fn(ind.spp, Cov[, "Trt_stat"], Cov[, "Trt_time"], mod)
  pnt.obs <- pnt.obs.fn(ind.spp, Cov, Y.mat, grd.obs$raw)
  p.pnt <- ggplot(dat.pred, aes(x = x, y = psi)) +
    geom_ribbon(data = dat.pred %>% filter(x.trt == 1),
                aes(ymin = psi.lo, ymax = psi.hi), alpha = 0.4) +
    geom_line(data = dat.pred %>% filter(x.trt == 1), size = 2) +
    geom_errorbar(data = dat.pred %>% filter(x.trt == 0),
                  aes(x = x, ymin = psi.lo, ymax = psi.hi), width = 0.2) +
    geom_point(data = dat.pred %>% filter(x.trt == 0), size = 3, shape = 15) +
    scale_x_continuous(breaks = c(1, 3, 5, 7, 9, 11), labels = c("pre-trt", "trt yr 2", "trt yr 4", "trt yr 6", "trt yr 8", "trt yr 10")) +
    ylim(0, 1) +
    xlab(NULL) + ylab(NULL)
  if(plot.obs == T) p.pnt <- p.pnt + geom_point(data = pnt.obs, aes(x = x + 1, y = Y), shape = 1)
  
  p.spp <- ggdraw() +
    draw_plot(p.grid, x = 0, y = 0, width = 0.475, height = 0.95) +
    draw_plot(p.pnt, x = 0.525, y = 0, width = 0.475, height = 0.95) +
    draw_plot_label(c(spp, spp), x = c(0.23, 0.77), y = c(1, 1))
  assign(str_c("p.spp", i), p.spp)
}

p <- ggdraw() +
  draw_plot(p.spp1, x = 0.05, y = 0.05, width = 0.2833333, height = 0.08636364) +
  draw_plot(p.spp2, x = 0.05, y = 0.1363636, width = 0.2833333, height = 0.08636364) +
  draw_plot(p.spp3, x = 0.05, y = 0.2227273, width = 0.2833333, height = 0.08636364) +
  draw_plot(p.spp4, x = 0.05, y = 0.3090909, width = 0.2833333, height = 0.08636364) +
  draw_plot(p.spp5, x = 0.05, y = 0.3954546, width = 0.2833333, height = 0.08636364) +
  draw_plot(p.spp6, x = 0.05, y = 0.4818182, width = 0.2833333, height = 0.08636364) +
  draw_plot(p.spp7, x = 0.05, y = 0.5681818, width = 0.2833333, height = 0.08636364) +
  draw_plot(p.spp8, x = 0.05, y = 0.6545455, width = 0.2833333, height = 0.08636364) +
  draw_plot(p.spp9, x = 0.05, y = 0.7409091, width = 0.2833333, height = 0.08636364) +
  draw_plot(p.spp10, x = 0.05, y = 0.8272728, width = 0.2833333, height = 0.08636364) +
  draw_plot(p.spp11, x = 0.05, y = 0.9136364, width = 0.2833333, height = 0.08636364) +
  draw_plot(p.spp12, x = 0.3833333, y = 0.05, width = 0.2833333, height = 0.08636364) +
  draw_plot(p.spp13, x = 0.3833333, y = 0.1363636, width = 0.2833333, height = 0.08636364) +
  draw_plot(p.spp14, x = 0.3833333, y = 0.2227273, width = 0.2833333, height = 0.08636364) +
  draw_plot(p.spp15, x = 0.3833333, y = 0.3090909, width = 0.2833333, height = 0.08636364) +
  draw_plot(p.spp16, x = 0.3833333, y = 0.3954546, width = 0.2833333, height = 0.08636364) +
  draw_plot(p.spp17, x = 0.3833333, y = 0.4818182, width = 0.2833333, height = 0.08636364) +
  draw_plot(p.spp18, x = 0.3833333, y = 0.5681818, width = 0.2833333, height = 0.08636364) +
  draw_plot(p.spp19, x = 0.3833333, y = 0.6545455, width = 0.2833333, height = 0.08636364) +
  draw_plot(p.spp20, x = 0.3833333, y = 0.7409091, width = 0.2833333, height = 0.08636364) +
  draw_plot(p.spp21, x = 0.3833333, y = 0.8272728, width = 0.2833333, height = 0.08636364) +
  draw_plot(p.spp22, x = 0.7166666, y = 0.05, width = 0.2833333, height = 0.08636364) +
  draw_plot(p.spp23, x = 0.7166666, y = 0.1363636, width = 0.2833333, height = 0.08636364) +
  draw_plot(p.spp24, x = 0.7166666, y = 0.2227273, width = 0.2833333, height = 0.08636364) +
  draw_plot(p.spp25, x = 0.7166666, y = 0.3090909, width = 0.2833333, height = 0.08636364) +
  draw_plot(p.spp26, x = 0.7166666, y = 0.3954546, width = 0.2833333, height = 0.08636364) +
  draw_plot(p.spp27, x = 0.7166666, y = 0.4818182, width = 0.2833333, height = 0.08636364) +
  draw_plot(p.spp28, x = 0.7166666, y = 0.5681818, width = 0.2833333, height = 0.08636364) +
  draw_plot(p.spp29, x = 0.7166666, y = 0.6545455, width = 0.2833333, height = 0.08636364) +
  draw_plot(p.spp30, x = 0.7166666, y = 0.7409091, width = 0.2833333, height = 0.08636364) +
  draw_plot(p.spp31, x = 0.7166666, y = 0.8272728, width = 0.2833333, height = 0.08636364)# +
#  draw_plot_label(c(spp, spp), x = c(0.23, 0.77), y = c(1, 1))

save_plot("Plot_spp_trt_rels_inspect_GOF.tiff", p, ncol = 8, nrow = 11, dpi = 200)



