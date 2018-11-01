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

pnt.obs.fn <- function(Y, Cov) {
  pnt.obs <- Cov %>% tbl_df %>%
    select(gridIndex, YearInd, Trt_stat, Trt_time) %>%
    mutate(Y = Y) %>%
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
  
  grd.ppred <- grd.obs.fn(Y.mat[, ind.spp], landscape_data)$sum %>%
    rename(Y.obs = Y) %>%
    mutate(Y.md = NA, Y.lo = NA, Y.hi = NA)
  grd.pmat <- matrix(NA, nrow = dim(ypred)[1], ncol = nrow(grd.ppred))
  
  pnt.ppred <- Y.mat[, ind.spp] %>%
    pnt.obs.fn(Cov) %>%
    rename(Y.obs = Y) %>%
    mutate(Y.md = NA, Y.lo = NA, Y.hi = NA)
  pnt.pmat <- matrix(NA, nrow = dim(ypred)[1], ncol = nrow(pnt.ppred))  
  
  for(s in 1:dim(ypred)[1]) {
    grd.pmat[s, ] <- grd.obs.fn(ypred[s, ], landscape_data)$sum$Y
    pnt.pmat[s, ] <- pnt.obs.fn(ypred[s, ], Cov)$Y
  }
  
  grd.ppred$Y.md <- apply(grd.pmat, 2, median)
  grd.ppred$Y.lo <- apply(grd.pmat, 2, function(x) quantile(x, prob = 0.025, type = 8))
  grd.ppred$Y.hi <- apply(grd.pmat, 2, function(x) quantile(x, prob = 0.975, type = 8))
  
  pnt.ppred$Y.md <- apply(pnt.pmat, 2, median)
  pnt.ppred$Y.lo <- apply(pnt.pmat, 2, function(x) quantile(x, prob = 0.025, type = 8))
  pnt.ppred$Y.hi <- apply(pnt.pmat, 2, function(x) quantile(x, prob = 0.975, type = 8))
  
  grd.ppred <- grd.ppred %>%
    mutate(test = (Y.obs < Y.lo | Y.obs > Y.hi)*1)
  pnt.ppred <- pnt.ppred %>%
    mutate(test = (Y.obs < Y.lo | Y.obs > Y.hi)*1)
  
  assign(str_c("GOF.grd.", spp), grd.ppred)
  assign(str_c("GOF.pnt.", spp), pnt.ppred)
}

GOF.grd <- GOF.grd.AMRO %>%
  mutate(spp = "AMRO") %>%
  bind_rows(GOF.grd.BRCR %>%
              mutate(spp = "BRCR")) %>%
  bind_rows(GOF.grd.BTLH %>%
              mutate(spp = "BTLH")) %>%
  bind_rows(GOF.grd.CAFI %>%
              mutate(spp = "CAFI")) %>%
  bind_rows(GOF.grd.CLNU %>%
              mutate(spp = "CLNU")) %>%
  bind_rows(GOF.grd.COFL %>%
              mutate(spp = "COFL")) %>%
  bind_rows(GOF.grd.CONI %>%
              mutate(spp = "CONI")) %>%
  bind_rows(GOF.grd.DEJU %>%
              mutate(spp = "DEJU")) %>%
  bind_rows(GOF.grd.EVGR %>%
              mutate(spp = "EVGR")) %>%
  bind_rows(GOF.grd.GRAJ %>%
              mutate(spp = "GRAJ")) %>%
  bind_rows(GOF.grd.GTTO %>%
              mutate(spp = "GTTO")) %>%
  bind_rows(GOF.grd.HAWO %>%
              mutate(spp = "HAWO")) %>%
  bind_rows(GOF.grd.HETH %>%
              mutate(spp = "HETH")) %>%
  bind_rows(GOF.grd.LISP %>%
              mutate(spp = "LISP")) %>%
  bind_rows(GOF.grd.MGWA %>%
              mutate(spp = "MGWA")) %>%
  bind_rows(GOF.grd.OSFL %>%
              mutate(spp = "OSFL")) %>%
  bind_rows(GOF.grd.PYNU %>%
              mutate(spp = "PYNU")) %>%
  bind_rows(GOF.grd.RBNU %>%
              mutate(spp = "RBNU")) %>%
  bind_rows(GOF.grd.RCKI %>%
              mutate(spp = "RCKI")) %>%
  bind_rows(GOF.grd.RECR %>%
              mutate(spp = "RECR")) %>%
  bind_rows(GOF.grd.SOSP %>%
              mutate(spp = "SOSP")) %>%
  bind_rows(GOF.grd.SPTO %>%
              mutate(spp = "SPTO")) %>%
  bind_rows(GOF.grd.STJA %>%
              mutate(spp = "STJA")) %>%
  bind_rows(GOF.grd.TOSO %>%
              mutate(spp = "TOSO")) %>%
  bind_rows(GOF.grd.VGSW %>%
              mutate(spp = "VGSW")) %>%
  bind_rows(GOF.grd.VIWA %>%
              mutate(spp = "VIWA")) %>%
  bind_rows(GOF.grd.WBNU %>%
              mutate(spp = "WBNU")) %>%
  bind_rows(GOF.grd.WEWP %>%
              mutate(spp = "WEWP")) %>%
  bind_rows(GOF.grd.WISA %>%
              mutate(spp = "WISA")) %>%
  bind_rows(GOF.grd.YEWA %>%
              mutate(spp = "YEWA")) %>%
  bind_rows(GOF.grd.YRWA %>%
              mutate(spp = "YRWA")) %>%
  select(spp, x_class:test)

GOF.pnt <- GOF.pnt.AMRO %>%
  mutate(spp = "AMRO") %>%
  bind_rows(GOF.pnt.BRCR %>%
              mutate(spp = "BRCR")) %>%
  bind_rows(GOF.pnt.BTLH %>%
              mutate(spp = "BTLH")) %>%
  bind_rows(GOF.pnt.CAFI %>%
              mutate(spp = "CAFI")) %>%
  bind_rows(GOF.pnt.CLNU %>%
              mutate(spp = "CLNU")) %>%
  bind_rows(GOF.pnt.COFL %>%
              mutate(spp = "COFL")) %>%
  bind_rows(GOF.pnt.CONI %>%
              mutate(spp = "CONI")) %>%
  bind_rows(GOF.pnt.DEJU %>%
              mutate(spp = "DEJU")) %>%
  bind_rows(GOF.pnt.EVGR %>%
              mutate(spp = "EVGR")) %>%
  bind_rows(GOF.pnt.GRAJ %>%
              mutate(spp = "GRAJ")) %>%
  bind_rows(GOF.pnt.GTTO %>%
              mutate(spp = "GTTO")) %>%
  bind_rows(GOF.pnt.HAWO %>%
              mutate(spp = "HAWO")) %>%
  bind_rows(GOF.pnt.HETH %>%
              mutate(spp = "HETH")) %>%
  bind_rows(GOF.pnt.LISP %>%
              mutate(spp = "LISP")) %>%
  bind_rows(GOF.pnt.MGWA %>%
              mutate(spp = "MGWA")) %>%
  bind_rows(GOF.pnt.OSFL %>%
              mutate(spp = "OSFL")) %>%
  bind_rows(GOF.pnt.PYNU %>%
              mutate(spp = "PYNU")) %>%
  bind_rows(GOF.pnt.RBNU %>%
              mutate(spp = "RBNU")) %>%
  bind_rows(GOF.pnt.RCKI %>%
              mutate(spp = "RCKI")) %>%
  bind_rows(GOF.pnt.RECR %>%
              mutate(spp = "RECR")) %>%
  bind_rows(GOF.pnt.SOSP %>%
              mutate(spp = "SOSP")) %>%
  bind_rows(GOF.pnt.SPTO %>%
              mutate(spp = "SPTO")) %>%
  bind_rows(GOF.pnt.STJA %>%
              mutate(spp = "STJA")) %>%
  bind_rows(GOF.pnt.TOSO %>%
              mutate(spp = "TOSO")) %>%
  bind_rows(GOF.pnt.VGSW %>%
              mutate(spp = "VGSW")) %>%
  bind_rows(GOF.pnt.VIWA %>%
              mutate(spp = "VIWA")) %>%
  bind_rows(GOF.pnt.WBNU %>%
              mutate(spp = "WBNU")) %>%
  bind_rows(GOF.pnt.WEWP %>%
              mutate(spp = "WEWP")) %>%
  bind_rows(GOF.pnt.WISA %>%
              mutate(spp = "WISA")) %>%
  bind_rows(GOF.pnt.YEWA %>%
              mutate(spp = "YEWA")) %>%
  bind_rows(GOF.pnt.YRWA %>%
              mutate(spp = "YRWA")) %>%
  select(spp, x_class:test)

GOF.grd %>%
  mutate(Scale = "grid") %>%
  bind_rows(GOF.pnt %>%
              mutate(Scale = "point")) %>%
  select(Scale, spp:test) %>%
  write.csv("GOF.csv", row.names = F)
