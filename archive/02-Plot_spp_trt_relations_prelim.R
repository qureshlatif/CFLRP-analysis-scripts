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
#___________________________#

##____ Species with positive grid-level relationships ____##
signif.ptrt <- (apply(mod$sims.list$bd.ptrt, 2, function(x) quantile(x, prob = 0.025, type = 8)) > 0 |
                  apply(mod$sims.list$bd.ptrt, 2, function(x) quantile(x, prob = 0.975, type = 8)) < 0)
spp.plot <- spp.list[which(signif.ptrt)]

for(i in 1:length(spp.plot)) {
  spp <- spp.plot[i]
  ind.spp <- which(spp.list == spp)
  
  # Grid level #
  dat.pred <- dat.pred.grid.fn(ind.spp, landscape_data$PctTrt, mod)
  p <- ggplot(dat.pred, aes(x = x, y = psi)) +
    geom_ribbon(aes(ymin = psi.lo, ymax = psi.hi), alpha = 0.4) +
    geom_line(size = 2) +
    ylim(0, 1) +
    xlab(NULL) + ylab(NULL)

  p <- ggdraw() +
    draw_plot(p, x = 0, y = 0, width = 1, height = 0.95) +
    draw_plot_label(spp, x = 0.5, y = 1)
  assign(str_c("pg", i), p)

  # Point level
  dat.pred <- dat.pred.pnt.fn(ind.spp, Cov[, "Trt_stat"], Cov[, "Trt_time"], mod)
  p <- ggplot(dat.pred, aes(x = x, y = psi)) +
    geom_ribbon(data = dat.pred %>% filter(x.trt == 1),
                aes(ymin = psi.lo, ymax = psi.hi), alpha = 0.4) +
    geom_line(data = dat.pred %>% filter(x.trt == 1), size = 2) +
    geom_errorbar(data = dat.pred %>% filter(x.trt == 0),
                  aes(x = x, ymin = psi.lo, ymax = psi.hi), width = 0.2) +
    geom_point(data = dat.pred %>% filter(x.trt == 0), size = 3, shape = 15) +
    scale_x_continuous(breaks = c(1, 3, 5, 7, 9, 11), labels = c("pre-trt", "trt yr 2", "trt yr 4", "trt yr 6", "trt yr 8", "trt yr 10")) +
    ylim(0, 1) +
    xlab(NULL) + ylab(NULL)

  p <- ggdraw() +
    draw_plot(p, x = 0, y = 0, width = 1, height = 0.95) +
    draw_plot_label(spp, x = 0.5, y = 1)
  assign(str_c("pp", i), p)
}

p <- ggdraw() +
  draw_plot(pg1, x = 0.03, y = 0.8416667, width = 0.1583333, height = 0.1583333) +
  draw_plot(pp1, x = 0.1883333, y = 0.8416667, width = 0.1583333, height = 0.1583333) +
  draw_plot(pg2, x = 0.03, y = 0.6833333, width = 0.1583333, height = 0.1583333) +
  draw_plot(pp2, x = 0.1883333, y = 0.6833333, width = 0.1583333, height = 0.1583333) +
  draw_plot(pg3, x = 0.03, y = 0.5250000, width = 0.1583333, height = 0.1583333) +
  draw_plot(pp3, x = 0.1883333, y = 0.5250000, width = 0.1583333, height = 0.1583333) +
  draw_plot(pg4, x = 0.03, y = 0.3666667, width = 0.1583333, height = 0.1583333) +
  draw_plot(pp4, x = 0.1883333, y = 0.3666667, width = 0.1583333, height = 0.1583333) +
  draw_plot(pg5, x = 0.03, y = 0.2083333, width = 0.1583333, height = 0.1583333) +
  draw_plot(pp5, x = 0.1883333, y = 0.2083333, width = 0.1583333, height = 0.1583333) +
  draw_plot(pg6, x = 0.03, y = 0.05, width = 0.1583333, height = 0.1583333) +
  draw_plot(pp6, x = 0.1883333, y = 0.05, width = 0.1583333, height = 0.1583333) +
  
  draw_plot(pg7, x = 0.3466667, y = 0.8416667, width = 0.1583333, height = 0.1583333) +
  draw_plot(pp7, x = 0.5050000, y = 0.8416667, width = 0.1583333, height = 0.1583333) +
  draw_plot(pg8, x = 0.3466667, y = 0.6833333, width = 0.1583333, height = 0.1583333) +
  draw_plot(pp8, x = 0.5050000, y = 0.6833333, width = 0.1583333, height = 0.1583333) +
  draw_plot(pg9, x = 0.3466667, y = 0.5250000, width = 0.1583333, height = 0.1583333) +
  draw_plot(pp9, x = 0.5050000, y = 0.5250000, width = 0.1583333, height = 0.1583333) +
  draw_plot(pg10, x = 0.3466667, y = 0.3666667, width = 0.1583333, height = 0.1583333) +
  draw_plot(pp10, x = 0.5050000, y = 0.3666667, width = 0.1583333, height = 0.1583333) +
  draw_plot(pg11, x = 0.3466667, y = 0.2083333, width = 0.1583333, height = 0.1583333) +
  draw_plot(pp11, x = 0.5050000, y = 0.2083333, width = 0.1583333, height = 0.1583333) +
  draw_plot(pg12, x = 0.3466667, y = 0.05, width = 0.1583333, height = 0.1583333) +
  draw_plot(pp12, x = 0.5050000, y = 0.05, width = 0.1583333, height = 0.1583333) +
  
  draw_plot(pg13, x = 0.6633333, y = 0.8416667, width = 0.1583333, height = 0.1583333) +
  draw_plot(pp13, x = 0.8216667, y = 0.8416667, width = 0.1583333, height = 0.1583333) +
  draw_plot(pg14, x = 0.6633333, y = 0.6833333, width = 0.1583333, height = 0.1583333) +
  draw_plot(pp14, x = 0.8216667, y = 0.6833333, width = 0.1583333, height = 0.1583333) +
  draw_plot(pg15, x = 0.6633333, y = 0.5250000, width = 0.1583333, height = 0.1583333) +
  draw_plot(pp15, x = 0.8216667, y = 0.5250000, width = 0.1583333, height = 0.1583333) +
  draw_plot(pg16, x = 0.6633333, y = 0.3666667, width = 0.1583333, height = 0.1583333) +
  draw_plot(pp16, x = 0.8216667, y = 0.3666667, width = 0.1583333, height = 0.1583333) +
  draw_plot(pg17, x = 0.6633333, y = 0.2083333, width = 0.1583333, height = 0.1583333) +
  draw_plot(pp17, x = 0.8216667, y = 0.2083333, width = 0.1583333, height = 0.1583333) +
  draw_plot(pg18, x = 0.6633333, y = 0.05, width = 0.1583333, height = 0.1583333) +
  draw_plot(pp18, x = 0.8216667, y = 0.05, width = 0.1583333, height = 0.1583333) +
  draw_plot_label(c("Occupancy", "Percent treated (grid)", "Treatment status (point)", "Percent treated (grid)",
                    "Treatment status (point)", "Percent treated (grid)", "Treatment status (point)"),
                  x = c(0, 0.07, 0.22, 0.38, 0.54, 0.71, 0.86), y = c(0.5, rep(0.05, 6)),
                  angle = c(90, rep(0, 6)), size = c(35, rep(15, 6)), hjust = rep(0, 7))

save_plot("Plot_spp_w_grid_relations.tiff", p, ncol = 5, nrow = 5, dpi = 200)


##____ Species without significant grid relationships ____##
signif.trt <- (apply(mod$sims.list$bb.trt, 2, function(x) quantile(x, prob = 0.025, type = 8)) > 0 |
                 apply(mod$sims.list$bb.trt, 2, function(x) quantile(x, prob = 0.975, type = 8)) < 0)
signif.yst <- (apply(mod$sims.list$bb.YST, 2, function(x) quantile(x, prob = 0.025, type = 8)) > 0 |
                 apply(mod$sims.list$bb.YST, 2, function(x) quantile(x, prob = 0.975, type = 8)) < 0)
spp.plot <- spp.list[which((signif.trt | signif.yst) & !signif.ptrt)]
spp.plot <- sort(spp.plot)

for(i in 1:length(spp.plot)) {
  spp <- spp.plot[i]
  ind.spp <- which(spp.list == spp)
  
  # Grid level #
  dat.pred <- dat.pred.grid.fn(ind.spp, landscape_data$PctTrt, mod)
  p <- ggplot(dat.pred, aes(x = x, y = psi)) +
    geom_ribbon(aes(ymin = psi.lo, ymax = psi.hi), alpha = 0.4) +
    geom_line(size = 2) +
    ylim(0, 1) +
    xlab(NULL) + ylab(NULL)
  
  p <- ggdraw() +
    draw_plot(p, x = 0, y = 0, width = 1, height = 0.95) +
    draw_plot_label(spp, x = 0.5, y = 1)
  assign(str_c("pg", i), p)
  
  # Point level
  dat.pred <- dat.pred.pnt.fn(ind.spp, Cov[, "Trt_stat"], Cov[, "Trt_time"], mod)
  p <- ggplot(dat.pred, aes(x = x, y = psi)) +
    geom_ribbon(data = dat.pred %>% filter(x.trt == 1),
                aes(ymin = psi.lo, ymax = psi.hi), alpha = 0.4) +
    geom_line(data = dat.pred %>% filter(x.trt == 1), size = 2) +
    geom_errorbar(data = dat.pred %>% filter(x.trt == 0),
                  aes(x = x, ymin = psi.lo, ymax = psi.hi), width = 0.2) +
    geom_point(data = dat.pred %>% filter(x.trt == 0), size = 3, shape = 15) +
    scale_x_continuous(breaks = c(1, 3, 5, 7, 9, 11), labels = c("pre-trt", "trt yr 2", "trt yr 4", "trt yr 6", "trt yr 8", "trt yr 10")) +
    ylim(0, 1) +
    xlab(NULL) + ylab(NULL)
  
  p <- ggdraw() +
    draw_plot(p, x = 0, y = 0, width = 1, height = 0.95) +
    draw_plot_label(spp, x = 0.5, y = 1)
  assign(str_c("pp", i), p)
}

p <- ggdraw() +
  draw_plot(pg1, x = 0.05, y = 0.8416667, width = 0.2375, height = 0.1583333) +
  draw_plot(pp1, x = 0.2875, y = 0.8416667, width = 0.2375, height = 0.1583333) +
  draw_plot(pg2, x = 0.05, y = 0.6833333, width = 0.2375, height = 0.1583333) +
  draw_plot(pp2, x = 0.2875, y = 0.6833333, width = 0.2375, height = 0.1583333) +
  draw_plot(pg3, x = 0.05, y = 0.5250000, width = 0.2375, height = 0.1583333) +
  draw_plot(pp3, x = 0.2875, y = 0.5250000, width = 0.2375, height = 0.1583333) +
  draw_plot(pg4, x = 0.05, y = 0.3666667, width = 0.2375, height = 0.1583333) +
  draw_plot(pp4, x = 0.2875, y = 0.3666667, width = 0.2375, height = 0.1583333) +
  draw_plot(pg5, x = 0.05, y = 0.2083333, width = 0.2375, height = 0.1583333) +
  draw_plot(pp5, x = 0.2875, y = 0.2083333, width = 0.2375, height = 0.1583333) +
  draw_plot(pg6, x = 0.05, y = 0.05, width = 0.2375, height = 0.1583333) +
  draw_plot(pp6, x = 0.2875, y = 0.05, width = 0.2375, height = 0.1583333) +
  
  draw_plot(pg7, x = 0.5250, y = 0.8416667, width = 0.2375, height = 0.1583333) +
  draw_plot(pp7, x = 0.7625, y = 0.8416667, width = 0.2375, height = 0.1583333) +
  draw_plot(pg8, x = 0.5250, y = 0.6833333, width = 0.2375, height = 0.1583333) +
  draw_plot(pp8, x = 0.7625, y = 0.6833333, width = 0.2375, height = 0.1583333) +
  draw_plot(pg9, x = 0.5250, y = 0.5250000, width = 0.2375, height = 0.1583333) +
  draw_plot(pp9, x = 0.7625, y = 0.5250000, width = 0.2375, height = 0.1583333) +
  draw_plot(pg10, x = 0.5250, y = 0.3666667, width = 0.2375, height = 0.1583333) +
  draw_plot(pp10, x = 0.7625, y = 0.3666667, width = 0.2375, height = 0.1583333) +
  draw_plot(pg11, x = 0.5250, y = 0.2083333, width = 0.2375, height = 0.1583333) +
  draw_plot(pp11, x = 0.7625, y = 0.2083333, width = 0.2375, height = 0.1583333) +
  
  draw_plot_label(c("Occupancy", "Percent treated (grid)", "Treatment status (point)",
                    "Percent treated (grid)", "Treatment status (point)"),
                  x = c(0, 0.12, 0.35, 0.6, 0.83), y = c(0.5, rep(0.05, 4)),
                  angle = c(90, rep(0, 4)), size = c(35, rep(15, 4)), hjust = rep(0, 5))

save_plot("Plot_spp_w_nogrid_relations.tiff", p, ncol = 4, nrow = 5, dpi = 200)



