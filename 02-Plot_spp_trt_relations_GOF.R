library(QSLpersonal)
library(jagsUI)
library(stringr)
library(dplyr)
library(R.utils)
library(ggplot2)
library(cowplot)

setwd("C:/Users/Quresh.Latif/files/projects/FS/CFLRP")
load("GOF_workspace.RData")

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

grd.obs.fn <- function(Y, landscape_data) {
  grd.obs <- landscape_data %>%
    select(gridIndex, YearInd, PctTrt) %>%
    left_join(as.data.frame(Cov) %>%
                select(gridIndex, YearInd) %>%
                mutate(Y = Y) %>%
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

pnt.obs.fn <- function(Y, Cov, grd.obs) {
  pnt.obs <- Cov %>% tbl_df %>%
    select(gridIndex, YearInd, Trt_stat, Trt_time) %>%
    mutate(Y = Y) %>%
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

##____ Tabulate posterior predictive reference values ____##
signif.ptrt <- (apply(mod$sims.list$bd.ptrt, 2, function(x) quantile(x, prob = 0.025, type = 8)) > 0 |
                  apply(mod$sims.list$bd.ptrt, 2, function(x) quantile(x, prob = 0.975, type = 8)) < 0)
signif.trt <- (apply(mod$sims.list$bb.trt, 2, function(x) quantile(x, prob = 0.025, type = 8)) > 0 |
                 apply(mod$sims.list$bb.trt, 2, function(x) quantile(x, prob = 0.975, type = 8)) < 0)
signif.yst <- (apply(mod$sims.list$bb.YST, 2, function(x) quantile(x, prob = 0.025, type = 8)) > 0 |
                 apply(mod$sims.list$bb.YST, 2, function(x) quantile(x, prob = 0.975, type = 8)) < 0)
spp.subset <- spp.list[which(signif.ptrt | signif.trt | signif.yst)]

for(i in 1:length(spp.subset)) {
  spp <- spp.subset[i]
  ind.spp <- which(spp.list == spp)
  
  z.grid <- array(NA, dim = c(mod$mcmc.info$n.samples, n.grid, n.year))
  for(g in 1:n.grid) for(t in 1:n.year) {
    psi <- expit(mod$sims.list$d0[, ind.spp, t] +
                   mod$sims.list$bd.ptrt[, ind.spp] * PctTrt.d[g, t] +
                   mod$sims.list$bd.YST[, ind.spp] * YST.d[g, t] +
                   mod$sims.list$bd.TWIP[, ind.spp] * TWIP.d[g, t] +
                   mod$sims.list$bd.Rdens[, ind.spp] * Rdens.d[g, t])
    z.grid[, g, t] <- rbinom(mod$mcmc.info$n.samples, 1, psi)
  }
  
  ypred <- matrix(NA, nrow = mod$mcmc.info$n.samples, ncol = nrow(Y.mat))
  for(j in 1:ncol(ypred)) {
    theta <- expit(mod$sims.list$b0[, ind.spp] +
                     mod$sims.list$bb.trt[, ind.spp] * Trt.b[j] +
                     mod$sims.list$bb.YST[, ind.spp] * YST.b[j])
    z.point <- rbinom(mod$mcmc.info$n.samples, 1,
                             theta * z.grid[, gridID[j], yearID[j]])
    p <- expit(mod$sims.list$a0[, ind.spp] +
                 mod$sims.list$ba.Time[, ind.spp] * Time.b[j] +
                 mod$sims.list$ba.Time2[, ind.spp] * (Time.b[j]^2) +
                 mod$sims.list$ba.DOY[, ind.spp] * DOY.b[j] +
                 mod$sims.list$ba.trt[, ind.spp] * Trt.b[j] +
                 mod$sims.list$ba.YST[, ind.spp] * YST.b[j])
    ypred[, j] <- (rbinom(mod$mcmc.info$n.samples, 6, p * z.point) > 0) * 1
  }
  
  
}

##____ Grid level plots ____##
signif.ptrt <- (apply(mod$sims.list$bd.ptrt, 2, function(x) quantile(x, prob = 0.025, type = 8)) > 0 |
                       apply(mod$sims.list$bd.ptrt, 2, function(x) quantile(x, prob = 0.975, type = 8)) < 0)
spp.plot <- spp.list[which(signif.ptrt)]
spp.plot <- sort(spp.plot)

for(i in 1:length(spp.plot)) {
  spp <- spp.plot[i]
  ind.spp <- which(spp.list == spp)
  
  # Grid level #
  dat.pred <- dat.pred.grid.fn(ind.spp, landscape_data$PctTrt, mod)
  grd.obs <- grd.obs.fn(Y.mat[, ind.spp], landscape_data)
  p <- ggplot(dat.pred, aes(x = x, y = psi)) +
    geom_ribbon(aes(ymin = psi.lo, ymax = psi.hi), alpha = 0.4) +
    geom_line(size = 2) +
    ylim(0, 1) +
    xlab(NULL) + ylab(NULL)
  if(plot.obs == T) p <- p + geom_point(data = grd.obs$sum, aes(x = x, y = Y), size = 3, color = "red")

  p <- ggdraw() +
    draw_plot(p, x = 0, y = 0, width = 1, height = 0.95) +
    draw_plot_label(spp, x = 0.5, y = 1)
  assign(str_c("p", i), p)
}

p <- ggdraw() +
  draw_plot(p1, x = 0.05, y = 0.81, width = 0.2375, height = 0.19) +
  draw_plot(p2, x = 0.05, y = 0.62, width = 0.2375, height = 0.19) +
  draw_plot(p3, x = 0.05, y = 0.43, width = 0.2375, height = 0.19) +
  draw_plot(p4, x = 0.05, y = 0.24, width = 0.2375, height = 0.19) +
  draw_plot(p5, x = 0.05, y = 0.05, width = 0.2375, height = 0.19) +
  draw_plot(p6, x = 0.2875, y = 0.62, width = 0.2375, height = 0.19) +
  draw_plot(p7, x = 0.2875, y = 0.43, width = 0.2375, height = 0.19) +
  draw_plot(p8, x = 0.2875, y = 0.24, width = 0.2375, height = 0.19) +
  draw_plot(p9, x = 0.2875, y = 0.05, width = 0.2375, height = 0.19) +
  draw_plot(p10, x = 0.525, y = 0.62, width = 0.2375, height = 0.19) +
  draw_plot(p11, x = 0.525, y = 0.43, width = 0.2375, height = 0.19) +
  draw_plot(p12, x = 0.525, y = 0.24, width = 0.2375, height = 0.19) +
  draw_plot(p13, x = 0.525, y = 0.05, width = 0.2375, height = 0.19) +
  draw_plot(p14, x = 0.7625, y = 0.62, width = 0.2375, height = 0.19) +
  draw_plot(p15, x = 0.7625, y = 0.43, width = 0.2375, height = 0.19) +
  draw_plot(p16, x = 0.7625, y = 0.24, width = 0.2375, height = 0.19) +
  draw_plot(p17, x = 0.7625, y = 0.05, width = 0.2375, height = 0.19) +
  draw_plot_label(c("Grid occupancy", "Percent treated"),
                  x = c(0, 0.4), y = c(0.4, 0.05), angle = c(90, 0), size = c(35, 35))

save_plot("Plot_spp_trt_grid_occupancy_inspect_GOF.tiff", p, ncol = 4, nrow = 5, dpi = 200)


##____ Point level plots ____##
signif.trt <- (apply(mod$sims.list$bb.trt, 2, function(x) quantile(x, prob = 0.025, type = 8)) > 0 |
                 apply(mod$sims.list$bb.trt, 2, function(x) quantile(x, prob = 0.975, type = 8)) < 0)
signif.yst <- (apply(mod$sims.list$bb.YST, 2, function(x) quantile(x, prob = 0.025, type = 8)) > 0 |
                 apply(mod$sims.list$bb.YST, 2, function(x) quantile(x, prob = 0.975, type = 8)) < 0)
spp.plot <- spp.list[which(signif.trt | signif.yst)]
spp.plot <- sort(spp.plot)

for(i in 1:length(spp.plot)) {
  spp <- spp.plot[i]
  ind.spp <- which(spp.list == spp)
  
  dat.pred <- dat.pred.pnt.fn(ind.spp, Cov[, "Trt_stat"], Cov[, "Trt_time"], mod)
  pnt.obs <- pnt.obs.fn(Y.mat[, ind.spp], Cov, grd.obs$raw)
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
  if(plot.obs == T) p <- p + geom_point(data = pnt.obs, aes(x = x + 1, y = Y), size = 3, color = "red")
  
  p <- ggdraw() +
    draw_plot(p, x = 0, y = 0, width = 1, height = 0.95) +
    draw_plot_label(spp, x = 0.5, y = 1)
  assign(str_c("p", i), p)
}

p <- ggdraw() +
  draw_plot(p1, x = 0.05, y = 0.8642858, width = 0.3166667, height = 0.1357143) +
  draw_plot(p2, x = 0.05, y = 0.7285715, width = 0.3166667, height = 0.1357143) +
  draw_plot(p3, x = 0.05, y = 0.5928572, width = 0.3166667, height = 0.1357143) +
  draw_plot(p4, x = 0.05, y = 0.4571429, width = 0.3166667, height = 0.1357143) +
  draw_plot(p5, x = 0.05, y = 0.3214286, width = 0.3166667, height = 0.1357143) +
  draw_plot(p6, x = 0.05, y = 0.1857143, width = 0.3166667, height = 0.1357143) +
  draw_plot(p7, x = 0.05, y = 0.05, width = 0.3166667, height = 0.1357143) +
  draw_plot(p8, x = 0.3666667, y = 0.8642858, width = 0.3166667, height = 0.1357143) +
  draw_plot(p9, x = 0.3666667, y = 0.7285715, width = 0.3166667, height = 0.1357143) +
  draw_plot(p10, x = 0.3666667, y = 0.5928572, width = 0.3166667, height = 0.1357143) +
  draw_plot(p11, x = 0.3666667, y = 0.4571429, width = 0.3166667, height = 0.1357143) +
  draw_plot(p12, x = 0.3666667, y = 0.3214286, width = 0.3166667, height = 0.1357143) +
  draw_plot(p13, x = 0.3666667, y = 0.1857143, width = 0.3166667, height = 0.1357143) +
  draw_plot(p14, x = 0.3666667, y = 0.05, width = 0.3166667, height = 0.1357143) +
  draw_plot(p15, x = 0.6833334, y = 0.8642858, width = 0.3166667, height = 0.1357143) +
  draw_plot(p16, x = 0.6833334, y = 0.7285715, width = 0.3166667, height = 0.1357143) +
  draw_plot(p17, x = 0.6833334, y = 0.5928572, width = 0.3166667, height = 0.1357143) +
  draw_plot(p18, x = 0.6833334, y = 0.4571429, width = 0.3166667, height = 0.1357143) +
  draw_plot(p19, x = 0.6833334, y = 0.3214286, width = 0.3166667, height = 0.1357143) +
  draw_plot(p20, x = 0.6833334, y = 0.1857143, width = 0.3166667, height = 0.1357143) +
  draw_plot(p21, x = 0.6833334, y = 0.05, width = 0.3166667, height = 0.1357143) +
  draw_plot_label(c("Point occupancy", "Treatment status & year"),
                  x = c(0, 0.3), y = c(0.4, 0.05), angle = c(90, 0), size = c(35, 35))

save_plot("Plot_spp_trt_point_occupancy_inspect_GOF.tiff", p, ncol = 4, nrow = 6, dpi = 200)



