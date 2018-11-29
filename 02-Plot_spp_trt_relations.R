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
  x.trt <- c(0, rep(1, 2))
  z.trt <- (x.trt - mean(trt, na.rm = T)) / sd(trt, na.rm = T)
  x.yst <- c(0, 1, 10)
  z.yst <- (x.yst - mean(yst, na.rm = T)) / sd(yst, na.rm = T)
  dat.pred <- data.frame(x.trt = x.trt, z.trt = z.trt, x.yst = x.yst, z.yst = z.yst,
                         psi = NA, psi.lo = NA, psi.hi = NA) %>%
    mutate(x = 1:n())
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

##____ Group 1 - positive relations at both scales ____##
spp.plot <- c("OSFL", "WEWP", "STJA", "PYNU", "CAFI", "RECR")

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
    draw_plot_label(spp, x = 0.4, y = 1)
  assign(str_c("pg", i), p)

  # Point level
  dat.pred <- dat.pred.pnt.fn(ind.spp, Cov[, "Trt_stat"], Cov[, "Trt_time"], mod)
  p <- ggplot(dat.pred, aes(x = x, y = psi)) +
    geom_line(data = dat.pred %>% filter(x.trt == 1), size = 1, linetype = "dashed") +
    geom_point(data = dat.pred %>% filter(x.trt == 1), size = 3, shape = 16) +
    geom_errorbar(data = dat.pred %>% filter(x.trt == 1),
                aes(ymin = psi.lo, ymax = psi.hi), width = 0.2) +
    geom_errorbar(data = dat.pred %>% filter(x.trt == 0),
                  aes(x = x, ymin = psi.lo, ymax = psi.hi), width = 0.2) +
    geom_point(data = dat.pred %>% filter(x.trt == 0), size = 3, shape = 15) +
    scale_x_continuous(breaks = c(1, 2, 3), labels = c("pre-trt", "trt yr 1", "trt yr 10")) +
    ylim(0, 1) +
    xlab(NULL) + ylab(NULL)

  p <- ggdraw() +
    draw_plot(p, x = 0, y = 0, width = 1, height = 0.95) +
    draw_plot_label(spp, x = 0.4, y = 1)
  assign(str_c("pp", i), p)
}

p1 <- ggdraw() +
  draw_plot(pg1, x = 0, y = 0.6667, width = 0.5, height = 0.3333) +
  draw_plot(pp1, x = 0.5, y = 0.6667, width = 0.5, height = 0.3333) +
  draw_plot(pg2, x = 0, y = 0.3333, width = 0.5, height = 0.3333) +
  draw_plot(pp2, x = 0.5, y = 0.3333, width = 0.5, height = 0.3333) +
  draw_plot(pg3, x = 0, y = 0, width = 0.5, height = 0.3333) +
  draw_plot(pp3, x = 0.5, y = 0, width = 0.5, height = 0.3333)

p2 <- ggdraw() +
  draw_plot(pg4, x = 0, y = 0.6667, width = 0.5, height = 0.3333) +
  draw_plot(pp4, x = 0.5, y = 0.6667, width = 0.5, height = 0.3333) +
  draw_plot(pg5, x = 0, y = 0.3333, width = 0.5, height = 0.3333) +
  draw_plot(pp5, x = 0.5, y = 0.3333, width = 0.5, height = 0.3333) +
  draw_plot(pg6, x = 0, y = 0, width = 0.5, height = 0.3333) +
  draw_plot(pp6, x = 0.5, y = 0, width = 0.5, height = 0.3333)

p <- ggdraw() +
  draw_plot(p1, x = 0.05, y = 0.05, width = 0.475, height = 0.95) +
  draw_plot(p2, x = 0.525, y = 0.05, width = 0.475, height = 0.95) +
  draw_plot_label(c("Occupancy", "Percent landscape treated", "Treatment status of point",
                    "Percent landscape treated", "Treatment status of point"),
                  x = c(0, 0.08, 0.33, 0.56, 0.8), y = c(0.47, 0.05, 0.05, 0.05, 0.05),
                  angle = c(90, 0, 0, 0, 0), size = c(20, 15, 15, 15, 15), hjust = c(0, 0, 0, 0, 0))

save_plot("Plot_spp_positive_grid&pnt_relations.tiff", p, ncol = 3, nrow = 2.5, dpi = 200)


##____ Group 2 - positive relations at grid only ____##
spp.plot <- c("CONI", "WISA", "GRAJ", "CLNU",
              "WBNU", "AMRO", "PIGR", "BHCO")

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
    draw_plot_label(spp, x = 0.4, y = 1)
  assign(str_c("pg", i), p)
  
  # Point level
  dat.pred <- dat.pred.pnt.fn(ind.spp, Cov[, "Trt_stat"], Cov[, "Trt_time"], mod)
  p <- ggplot(dat.pred, aes(x = x, y = psi)) +
    geom_line(data = dat.pred %>% filter(x.trt == 1), size = 1, linetype = "dashed") +
    geom_point(data = dat.pred %>% filter(x.trt == 1), size = 3, shape = 16) +
    geom_errorbar(data = dat.pred %>% filter(x.trt == 1),
                  aes(ymin = psi.lo, ymax = psi.hi), width = 0.2) +
    geom_errorbar(data = dat.pred %>% filter(x.trt == 0),
                  aes(x = x, ymin = psi.lo, ymax = psi.hi), width = 0.2) +
    geom_point(data = dat.pred %>% filter(x.trt == 0), size = 3, shape = 15) +
    scale_x_continuous(breaks = c(1, 2, 3), labels = c("pre-trt", "trt yr 1", "trt yr 10")) +
    ylim(0, 1) +
    xlab(NULL) + ylab(NULL)
  
  p <- ggdraw() +
    draw_plot(p, x = 0, y = 0, width = 1, height = 0.95) +
    draw_plot_label(spp, x = 0.4, y = 1)
  assign(str_c("pp", i), p)
}

p1 <- ggdraw() +
  draw_plot(pg1, x = 0, y = 0.75, width = 0.5, height = 0.25) +
  draw_plot(pp1, x = 0.5, y = 0.75, width = 0.5, height = 0.25) +
  draw_plot(pg2, x = 0, y = 0.5, width = 0.5, height = 0.25) +
  draw_plot(pp2, x = 0.5, y = 0.5, width = 0.5, height = 0.25) +
  draw_plot(pg3, x = 0, y = 0.25, width = 0.5, height = 0.25) +
  draw_plot(pp3, x = 0.5, y = 0.25, width = 0.5, height = 0.25) +
  draw_plot(pg4, x = 0, y = 0, width = 0.5, height = 0.25) +
  draw_plot(pp4, x = 0.5, y = 0, width = 0.5, height = 0.25)

p2 <- ggdraw() +
  draw_plot(pg5, x = 0, y = 0.75, width = 0.5, height = 0.25) +
  draw_plot(pp5, x = 0.5, y = 0.75, width = 0.5, height = 0.25) +
  draw_plot(pg6, x = 0, y = 0.5, width = 0.5, height = 0.25) +
  draw_plot(pp6, x = 0.5, y = 0.5, width = 0.5, height = 0.25) +
  draw_plot(pg7, x = 0, y = 0.25, width = 0.5, height = 0.25) +
  draw_plot(pp7, x = 0.5, y = 0.25, width = 0.5, height = 0.25) +
  draw_plot(pg8, x = 0, y = 0, width = 0.5, height = 0.25) +
  draw_plot(pp8, x = 0.5, y = 0, width = 0.5, height = 0.25)

p <- ggdraw() +
  draw_plot(p1, x = 0.05, y = 0.05, width = 0.475, height = 0.95) +
  draw_plot(p2, x = 0.525, y = 0.05, width = 0.475, height = 0.95) +
  draw_plot_label(c("Occupancy", "Percent landscape treated", "Treatment status of point",
                    "Percent landscape treated", "Treatment status of point"),
                  x = c(0, 0.08, 0.33, 0.56, 0.8), y = c(0.47, 0.05, 0.05, 0.05, 0.05),
                  angle = c(90, 0, 0, 0, 0), size = c(20, 15, 15, 15, 15), hjust = c(0, 0, 0, 0, 0))

save_plot("Plot_spp_positive_grid_only.tiff", p, ncol = 3, nrow = 3, dpi = 200)

##____ Group 3 - point relations only ____##
spp.plot <- c("HAWO", "VGSW", "HETH",
              "SPTO", "SOSP", "VIWA",
              "MGWA", "YEWA", "YRWA")

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
    draw_plot_label(spp, x = 0.4, y = 1)
  assign(str_c("pg", i), p)
  
  # Point level
  dat.pred <- dat.pred.pnt.fn(ind.spp, Cov[, "Trt_stat"], Cov[, "Trt_time"], mod)
  p <- ggplot(dat.pred, aes(x = x, y = psi)) +
    geom_line(data = dat.pred %>% filter(x.trt == 1), size = 1, linetype = "dashed") +
    geom_point(data = dat.pred %>% filter(x.trt == 1), size = 3, shape = 16) +
    geom_errorbar(data = dat.pred %>% filter(x.trt == 1),
                  aes(ymin = psi.lo, ymax = psi.hi), width = 0.2) +
    geom_errorbar(data = dat.pred %>% filter(x.trt == 0),
                  aes(x = x, ymin = psi.lo, ymax = psi.hi), width = 0.2) +
    geom_point(data = dat.pred %>% filter(x.trt == 0), size = 3, shape = 15) +
    scale_x_continuous(breaks = c(1, 2, 3), labels = c("pre-trt", "trt yr 1", "trt yr 10")) +
    ylim(0, 1) +
    xlab(NULL) + ylab(NULL)
  
  p <- ggdraw() +
    draw_plot(p, x = 0, y = 0, width = 1, height = 0.95) +
    draw_plot_label(spp, x = 0.4, y = 1)
  assign(str_c("pp", i), p)
}

p1 <- ggdraw() +
  draw_plot(pg1, x = 0, y = 0.6667, width = 0.5, height = 0.3333) +
  draw_plot(pp1, x = 0.5, y = 0.6667, width = 0.5, height = 0.3333) +
  draw_plot(pg2, x = 0, y = 0.3333, width = 0.5, height = 0.3333) +
  draw_plot(pp2, x = 0.5, y = 0.3333, width = 0.5, height = 0.3333) +
  draw_plot(pg3, x = 0, y = 0, width = 0.5, height = 0.3333) +
  draw_plot(pp3, x = 0.5, y = 0, width = 0.5, height = 0.3333)

p2 <- ggdraw() +
  draw_plot(pg4, x = 0, y = 0.6667, width = 0.5, height = 0.3333) +
  draw_plot(pp4, x = 0.5, y = 0.6667, width = 0.5, height = 0.3333) +
  draw_plot(pg5, x = 0, y = 0.3333, width = 0.5, height = 0.3333) +
  draw_plot(pp5, x = 0.5, y = 0.3333, width = 0.5, height = 0.3333) +
  draw_plot(pg6, x = 0, y = 0, width = 0.5, height = 0.3333) +
  draw_plot(pp6, x = 0.5, y = 0, width = 0.5, height = 0.3333)

p3 <- ggdraw() +
  draw_plot(pg7, x = 0, y = 0.6667, width = 0.5, height = 0.3333) +
  draw_plot(pp7, x = 0.5, y = 0.6667, width = 0.5, height = 0.3333) +
  draw_plot(pg8, x = 0, y = 0.3333, width = 0.5, height = 0.3333) +
  draw_plot(pp8, x = 0.5, y = 0.3333, width = 0.5, height = 0.3333) +
  draw_plot(pg9, x = 0, y = 0, width = 0.5, height = 0.3333) +
  draw_plot(pp9, x = 0.5, y = 0, width = 0.5, height = 0.3333)

p <- ggdraw() +
  draw_plot(p1, x = 0.05, y = 0.05, width = 0.3166667, height = 0.95) +
  draw_plot(p2, x = 0.3666667, y = 0.05, width = 0.3166667, height = 0.95) +
  draw_plot(p3, x = 0.6833333, y = 0.05, width = 0.3166667, height = 0.95) +
  draw_plot_label(c("Occupancy", "Percent landscape treated", "Treatment status of point",
                    "Percent landscape treated", "Treatment status of point",
                    "Percent landscape treated", "Treatment status of point"),
                  x = c(0, 0.08, 0.24, 0.4, 0.55, 0.71, 0.87),
                  y = c(0.47, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05),
                  angle = c(90, rep(0, 6)), size = c(20, rep(15, 6)),
                  hjust = rep(0, 7))

save_plot("Plot_spp_positive_point_only.tiff", p, ncol = 5, nrow = 3, dpi = 200)

##____ Group 4 - contrasting relations at different scales ____##
spp.plot <- c("COFL", "BRCR", "RCKI",
              "EVGR", "LISP", "DEJU")

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
    draw_plot_label(spp, x = 0.4, y = 1)
  assign(str_c("pg", i), p)
  
  # Point level
  dat.pred <- dat.pred.pnt.fn(ind.spp, Cov[, "Trt_stat"], Cov[, "Trt_time"], mod)
  p <- ggplot(dat.pred, aes(x = x, y = psi)) +
    geom_line(data = dat.pred %>% filter(x.trt == 1), size = 1, linetype = "dashed") +
    geom_point(data = dat.pred %>% filter(x.trt == 1), size = 3, shape = 16) +
    geom_errorbar(data = dat.pred %>% filter(x.trt == 1),
                  aes(ymin = psi.lo, ymax = psi.hi), width = 0.2) +
    geom_errorbar(data = dat.pred %>% filter(x.trt == 0),
                  aes(x = x, ymin = psi.lo, ymax = psi.hi), width = 0.2) +
    geom_point(data = dat.pred %>% filter(x.trt == 0), size = 3, shape = 15) +
    scale_x_continuous(breaks = c(1, 2, 3), labels = c("pre-trt", "trt yr 1", "trt yr 10")) +
    ylim(0, 1) +
    xlab(NULL) + ylab(NULL)
  
  p <- ggdraw() +
    draw_plot(p, x = 0, y = 0, width = 1, height = 0.95) +
    draw_plot_label(spp, x = 0.4, y = 1)
  assign(str_c("pp", i), p)
}

p1 <- ggdraw() +
  draw_plot(pg1, x = 0, y = 0.6667, width = 0.5, height = 0.3333) +
  draw_plot(pp1, x = 0.5, y = 0.6667, width = 0.5, height = 0.3333) +
  draw_plot(pg2, x = 0, y = 0.3333, width = 0.5, height = 0.3333) +
  draw_plot(pp2, x = 0.5, y = 0.3333, width = 0.5, height = 0.3333) +
  draw_plot(pg3, x = 0, y = 0, width = 0.5, height = 0.3333) +
  draw_plot(pp3, x = 0.5, y = 0, width = 0.5, height = 0.3333)

p2 <- ggdraw() +
  draw_plot(pg4, x = 0, y = 0.6667, width = 0.5, height = 0.3333) +
  draw_plot(pp4, x = 0.5, y = 0.6667, width = 0.5, height = 0.3333) +
  draw_plot(pg5, x = 0, y = 0.3333, width = 0.5, height = 0.3333) +
  draw_plot(pp5, x = 0.5, y = 0.3333, width = 0.5, height = 0.3333) +
  draw_plot(pg6, x = 0, y = 0, width = 0.5, height = 0.3333) +
  draw_plot(pp6, x = 0.5, y = 0, width = 0.5, height = 0.3333)

p <- ggdraw() +
  draw_plot(p1, x = 0.05, y = 0.05, width = 0.475, height = 0.95) +
  draw_plot(p2, x = 0.525, y = 0.05, width = 0.475, height = 0.95) +
  draw_plot_label(c("Occupancy", "Percent landscape treated", "Treatment status of point",
                    "Percent landscape treated", "Treatment status of point"),
                  x = c(0, 0.08, 0.33, 0.56, 0.8), y = c(0.47, 0.05, 0.05, 0.05, 0.05),
                  angle = c(90, 0, 0, 0, 0), size = c(20, 15, 15, 15, 15), hjust = c(0, 0, 0, 0, 0))

save_plot("Plot_spp_contrasting_relations.tiff", p, ncol = 3, nrow = 2.5, dpi = 200)

##____ RBNU ____##
spp <- "RBNU"
ind.spp <- which(spp.list == spp)

# Grid level #
dat.pred <- dat.pred.grid.fn(ind.spp, landscape_data$PctTrt, mod)
p.grd <- ggplot(dat.pred, aes(x = x, y = psi)) +
  geom_ribbon(aes(ymin = psi.lo, ymax = psi.hi), alpha = 0.4) +
  geom_line(size = 2) +
  ylim(0, 1) +
  xlab(NULL) + ylab(NULL)

p.grd <- ggdraw() +
  draw_plot(p.grd, x = 0, y = 0, width = 1, height = 0.95) +
  draw_plot_label(spp, x = 0.4, y = 1)

dat.pred <- dat.pred.pnt.fn(ind.spp, Cov[, "Trt_stat"], Cov[, "Trt_time"], mod)
p.pnt <- ggplot(dat.pred, aes(x = x, y = psi)) +
  geom_line(data = dat.pred %>% filter(x.trt == 1), size = 1, linetype = "dashed") +
  geom_point(data = dat.pred %>% filter(x.trt == 1), size = 3, shape = 16) +
  geom_errorbar(data = dat.pred %>% filter(x.trt == 1),
                aes(ymin = psi.lo, ymax = psi.hi), width = 0.2) +
  geom_errorbar(data = dat.pred %>% filter(x.trt == 0),
                aes(x = x, ymin = psi.lo, ymax = psi.hi), width = 0.2) +
  geom_point(data = dat.pred %>% filter(x.trt == 0), size = 3, shape = 15) +
  scale_x_continuous(breaks = c(1, 2, 3), labels = c("pre-trt", "trt yr 1", "trt yr 10")) +
  ylim(0, 1) +
  xlab(NULL) + ylab(NULL)

p.pnt <- ggdraw() +
  draw_plot(p.pnt, x = 0, y = 0, width = 1, height = 0.95) +
  draw_plot_label(spp, x = 0.4, y = 1)

p <- ggdraw() +
  draw_plot(p.grd, x = 0, y = 0, width = 0.5, height = 1) +
  draw_plot(p.pnt, x = 0.5, y = 0, width = 0.5, height = 1)

p <- ggdraw() +
  draw_plot(p, x = 0.05, y = 0.05, width = 0.95, height = 0.95) +
  draw_plot_label(c("Occupancy", "Percent landscape treated", "Treatment status of point"),
                  x = c(0, 0.15, 0.65), y = c(0.38, 0.08, 0.08),
                  angle = c(90, 0, 0), hjust = c(0, 0, 0))

save_plot("Plot_spp_RBNU.tiff", p, ncol = 2, nrow = 1, dpi = 200)
