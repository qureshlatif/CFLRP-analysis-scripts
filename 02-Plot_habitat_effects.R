library(jagsUI)
library(stringr)
library(dplyr)
library(R.utils)
library(ggplot2)
library(cowplot)
library(QSLpersonal)

setwd("C:/Users/Quresh.Latif/files/projects/FS/CFLRP")
load("Data_compiled.RData")

mod <- loadObject("mod_habitat_d0yr_reduced")

spp_trt_effects <- c("WISA", "RECR", "WEWP", "CAFI", "CONI",
                     "AMRO", "GRAJ", "PIGR", "EVGR", "CLNU",
                     "PYNU", "LISP", "BRCR", "DEJU", "WBNU",
                     "OSFL", "STJA", "BHCO", "COFL", "HAWO",
                     "MGWA", "RCKI", "YEWA", "VGSW", "YRWA",
                     "SOSP", "HETH", "VIWA")

# Tabulate parameter estimates
pars <- c("psi", "bd.TWI", "bd.heatload", "bd.ForAR", "bd.PACC10_3km", "bd.PACC40_3km", "bd.mnPerArRatio_Opn3km", 
          "theta", "bb.CanCov", "bb.CanHt", "bb.NumSnags", "bb.RCOV_PP", "bb.RCOV_DF", "bb.RCOV_AS", "bb.shvol",
          "bb.RSCV_Ladder", "bb.HerbGrassVol")
cols <- (c("", ".lo", ".hi") %>%
           expand.grid(pars, stringsAsFactors = F) %>%
           select(Var2, Var1) %>%
           mutate(Var3 = str_c(Var2, Var1, sep = "")))$Var3
tbl_pars <- matrix(NA, nrow = length(spp.list), ncol = length(cols), dimnames = list(NULL, cols))

for(par in pars[-which(pars %in% c("psi", "theta"))]) {
  parm <- mod$sims.list[[par]]
  tbl_pars[, par] <- apply(parm, 2, median)
  tbl_pars[, str_c(par, ".lo")] <- apply(parm, 2, function(x) quantile(x, prob = 0.025, type = 8))
  tbl_pars[, str_c(par, ".hi")] <- apply(parm, 2, function(x) quantile(x, prob = 0.975, type = 8))
}
rm(par)

parm <- expit(mod$sims.list[["d0"]])
tbl_pars[, "psi"] <- apply(parm, 2, median)
tbl_pars[, "psi.lo"] <- apply(parm, 2, function(x) quantile(x, prob = 0.025, type = 8))
tbl_pars[, "psi.hi"] <- apply(parm, 2, function(x) quantile(x, prob = 0.975, type = 8))

parm <- expit(mod$sims.list[["b0"]])
tbl_pars[, "theta"] <- apply(parm, 2, median)
tbl_pars[, "theta.lo"] <- apply(parm, 2, function(x) quantile(x, prob = 0.025, type = 8))
tbl_pars[, "theta.hi"] <- apply(parm, 2, function(x) quantile(x, prob = 0.975, type = 8))

rm(parm)
tbl_pars_all <- tbl_pars
tbl_pars <- tbl_pars[, c(1:3, 13:51)]

beta.cols <- dimnames(tbl_pars)[[2]][-which(dimnames(tbl_pars)[[2]] %in% c("psi", "psi.lo", "psi.hi",
                                                                           "theta", "theta.lo", "theta.hi"))]
ind.spp <- c(which(spp.list %in% spp_trt_effects),
             tbl_pars[, beta.cols[beta.cols %>% str_detect(".lo") %>% which]] %>%
               apply(1, function(x) any(x > 0)) %>%
               which,
             tbl_pars[, beta.cols[beta.cols %>% str_detect(".hi") %>% which]] %>%
               apply(1, function(x) any(x < 0)) %>%
               which) %>% unique %>% sort
spp.plt <- spp.list[ind.spp]

dat.plt <- tbl_pars %>% tbl_df() %>%
  mutate(Spp = spp.list) %>%
  filter(Spp %in% spp.plt) %>%
  mutate(index = row_number()) %>%
  mutate(index = (max(index) - index) + 1)
dat.plt$Spp[which(dat.plt$Spp %in% spp_trt_effects)] <-
  dat.plt$Spp[which(dat.plt$Spp %in% spp_trt_effects)] %>% str_c("*")

cols <- pars[-c(1:4, 8)] %>% str_c(".supp")
dat.supp <- matrix("none", nrow = nrow(dat.plt), ncol = length(cols),
                   dimnames = list(NULL, cols))
for(i in 1:length(cols)) {
  col.chck <- str_sub(cols[i], 1, -6)
  chck <- dat.plt[, which(str_detect(names(dat.plt), col.chck))]
  dat.supp[which(chck[, 2] > 0), cols[i]] <- "pos"
  dat.supp[which(chck[, 3] < 0), cols[i]] <- "neg"
}
rm(col.chck, chck)

dat.plt <- dat.plt %>%
  bind_cols(dat.supp %>% as.data.frame)
rm(dat.supp)

## Grid level relationships ##
p.psi <- ggplot(dat = dat.plt, aes(x = index, y = psi)) +
  geom_errorbar(aes(ymin = psi.lo, ymax = psi.hi), size=1, width=0) +
  geom_point(size = 2.5) + 
  coord_flip() +
  scale_x_continuous(breaks = 1:nrow(dat.plt), labels = dat.plt$Spp %>% rev, expand=c(0, 1)) +
  scale_y_continuous(lim = c(0, 1)) +
  ylab(expression(hat(psi)["mean"])) + xlab(NULL) +
  theme(axis.title.x=element_text(size=25))

p.PACCGap <- ggplot(dat = dat.plt, aes(x = index, y = bd.PACC10_3km, color = bd.PACC10_3km.supp)) +
  geom_errorbar(aes(ymin = bd.PACC10_3km.lo, ymax = bd.PACC10_3km.hi, color = bd.PACC10_3km.supp), size=1, width=0) +
  geom_point(size = 2.5) + 
  geom_hline(yintercept = 0) +
  coord_flip() +
  scale_x_continuous(breaks = 1:nrow(dat.plt), labels = dat.plt$Spp %>% rev, expand=c(0, 1)) +
  scale_y_continuous(lim = c(min(dat.plt$bd.PACC10_3km.lo), max(dat.plt$bd.PACC10_3km.hi))) +
  scale_color_manual(values = c("#0072B2", "dark gray", "#D55E00")) +
  ylab(expression(hat(beta)["PACCGap"])) + xlab(NULL) +
  theme(axis.title.x=element_text(size=25)) +
  guides(color = F)

p.PACCOpn <- ggplot(dat = dat.plt, aes(x = index, y = bd.PACC40_3km, color = bd.PACC40_3km.supp)) +
  geom_errorbar(aes(ymin = bd.PACC40_3km.lo, ymax = bd.PACC40_3km.hi, color = bd.PACC40_3km.supp), size=1, width=0) +
  geom_point(size = 2.5) + 
  geom_hline(yintercept = 0) +
  coord_flip() +
  scale_x_continuous(breaks = 1:nrow(dat.plt), labels = dat.plt$Spp %>% rev, expand=c(0, 1)) +
  scale_y_continuous(lim = c(min(dat.plt$bd.PACC10_3km.lo), max(dat.plt$bd.PACC10_3km.hi))) +
  scale_color_manual(values = c("#0072B2", "dark gray", "#D55E00")) +
  ylab(expression(hat(beta)["PACCOpn"])) + xlab(NULL) +
  theme(axis.title.x=element_text(size=25)) +
  guides(color = F)

p.PAROpn <- ggplot(dat = dat.plt, aes(x = index, y = bd.mnPerArRatio_Opn3km, color = bd.mnPerArRatio_Opn3km.supp)) +
  geom_errorbar(aes(ymin = bd.mnPerArRatio_Opn3km.lo, ymax = bd.mnPerArRatio_Opn3km.hi, color = bd.mnPerArRatio_Opn3km.supp), size=1, width=0) +
  geom_point(size = 2.5) + 
  geom_hline(yintercept = 0) +
  coord_flip() +
  scale_x_continuous(breaks = 1:nrow(dat.plt), labels = dat.plt$Spp %>% rev, expand=c(0, 1)) +
  scale_y_continuous(lim = c(min(dat.plt$bd.PACC10_3km.lo), max(dat.plt$bd.PACC10_3km.hi))) +
  scale_color_manual(values = c("#0072B2", "dark gray", "#D55E00")) +
  ylab(expression(hat(beta)["PAROpn"])) + xlab(NULL) +
  theme(axis.title.x=element_text(size=25)) +
  guides(color = F)

p <- ggdraw() + 
  draw_plot(p.psi, x = 0.05, y = 0, width = 0.2375, height = 1) +
  draw_plot(p.PACCGap, x = 0.2875, y = 0, width = 0.2375, height = 1) +
  draw_plot(p.PACCOpn, x = 0.5250, y = 0, width = 0.2375, height = 1) +
  draw_plot(p.PAROpn, x = 0.7625, y = 0, width = 0.2375, height = 1) +
  draw_plot_label("Species", x = 0, y = 0.5, size = 40, angle = 90, hjust = 0)

save_plot("Plot_landscape_effects.tiff", p, ncol = 4, nrow = 3, dpi = 200)

## Point level relationships 1 ##
p.theta <- ggplot(dat = dat.plt, aes(x = index, y = theta)) +
  geom_errorbar(aes(ymin = theta.lo, ymax = theta.hi), size=1, width=0) +
  geom_point(size = 2.5) + 
  coord_flip() +
  scale_x_continuous(breaks = 1:nrow(dat.plt), labels = dat.plt$Spp %>% rev, expand=c(0, 1)) +
  scale_y_continuous(lim = c(0, 1)) +
  ylab(expression(hat(theta)["mean"])) + xlab(NULL) +
  theme(axis.title.x=element_text(size=25))

p.CanCov <- ggplot(dat = dat.plt, aes(x = index, y = bb.CanCov, color = bb.CanCov.supp)) +
  geom_errorbar(aes(ymin = bb.CanCov.lo, ymax = bb.CanCov.hi, color = bb.CanCov.supp), size=1, width=0) +
  geom_point(size = 2.5) + 
  geom_hline(yintercept = 0) +
  coord_flip() +
  scale_x_continuous(breaks = 1:nrow(dat.plt), labels = dat.plt$Spp %>% rev, expand=c(0, 1)) +
  scale_y_continuous(lim = c(min(dat.plt$bb.RCOV_DF.lo), max(dat.plt$bb.HerbGrassVol.hi))) +
#  scale_y_continuous(lim = c(min(dat.plt$bb.CanCov.lo), max(dat.plt$bb.CanCov.hi))) +
  scale_color_manual(values = c("#0072B2", "dark gray", "#D55E00")) +
  ylab(expression(hat(beta)["CanCov"])) + xlab(NULL) +
  theme(axis.title.x=element_text(size=25)) +
  guides(color = F)
  
p.CanHt <- ggplot(dat = dat.plt, aes(x = index, y = bb.CanHt, color = bb.CanHt.supp)) +
  geom_errorbar(aes(ymin = bb.CanHt.lo, ymax = bb.CanHt.hi, color = bb.CanHt.supp), size=1, width=0) +
  geom_point(size = 2.5) + 
  geom_hline(yintercept = 0) +
  coord_flip() +
  scale_x_continuous(breaks = 1:nrow(dat.plt), labels = dat.plt$Spp %>% rev, expand=c(0, 1)) +
  scale_y_continuous(lim = c(min(dat.plt$bb.RCOV_DF.lo), max(dat.plt$bb.HerbGrassVol.hi))) +
#  scale_y_continuous(lim = c(min(dat.plt$bb.CanHt.lo), max(dat.plt$bb.CanHt.hi))) +
  scale_color_manual(values = c("#0072B2", "dark gray", "#D55E00")) +
  ylab(expression(hat(beta)["CanHt"])) + xlab(NULL) +
  theme(axis.title.x=element_text(size=25)) +
  guides(color = F)
  
p.NSnag <- ggplot(dat = dat.plt, aes(x = index, y = bb.NumSnags, color = bb.NumSnags.supp)) +
  geom_errorbar(aes(ymin = bb.NumSnags.lo, ymax = bb.NumSnags.hi, color = bb.NumSnags.supp), size=1, width=0) +
  geom_point(size = 2.5) + 
  geom_hline(yintercept = 0) +
  coord_flip() +
  scale_x_continuous(breaks = 1:nrow(dat.plt), labels = dat.plt$Spp %>% rev, expand=c(0, 1)) +
  scale_y_continuous(lim = c(min(dat.plt$bb.RCOV_DF.lo), max(dat.plt$bb.HerbGrassVol.hi))) +
#  scale_y_continuous(lim = c(min(dat.plt$bb.NumSnags.lo), max(dat.plt$bb.NumSnags.hi))) +
  scale_color_manual(values = c("#0072B2", "dark gray", "#D55E00")) +
  ylab(expression(hat(beta)["NSnag"])) + xlab(NULL) +
  theme(axis.title.x=element_text(size=25)) +
  guides(color = F)

p <- ggdraw() + 
  draw_plot(p.theta, x = 0.0500, y = 0, width = 0.2375, height = 1) +
  draw_plot(p.CanCov, x = 0.2875, y = 0, width = 0.2375, height = 1) +
  draw_plot(p.CanHt, x = 0.5250, y = 0, width = 0.2375, height = 1) +
  draw_plot(p.NSnag, x = 0.7625, y = 0, width = 0.2375, height = 1) +
  draw_plot_label("Species", x = 0, y = 0.5, size = 40, angle = 90, hjust = 0)

save_plot("Plot_veg_canStruct_effects.tiff", p, ncol = 3, nrow = 3, dpi = 200)

p.PIPO <- ggplot(dat = dat.plt, aes(x = index, y = bb.RCOV_PP, color = bb.RCOV_PP.supp)) +
  geom_errorbar(aes(ymin = bb.RCOV_PP.lo, ymax = bb.RCOV_PP.hi, color = bb.RCOV_PP.supp), size=1, width=0) +
  geom_point(size = 2.5) + 
  geom_hline(yintercept = 0) +
  coord_flip() +
  scale_x_continuous(breaks = 1:nrow(dat.plt), labels = dat.plt$Spp %>% rev, expand=c(0, 1)) +
  scale_y_continuous(lim = c(min(dat.plt$bb.RCOV_DF.lo), max(dat.plt$bb.HerbGrassVol.hi))) +
#  scale_y_continuous(lim = c(min(dat.plt$bb.RCOV_PP.lo), max(dat.plt$bb.RCOV_PP.hi))) +
  scale_color_manual(values = c("#0072B2", "dark gray", "#D55E00")) +
  ylab(expression(hat(beta)["PIPO"])) + xlab(NULL) +
  theme(axis.title.x=element_text(size=25)) +
  guides(color = F)
  
p.PSME <- ggplot(dat = dat.plt, aes(x = index, y = bb.RCOV_DF, color = bb.RCOV_DF.supp)) +
  geom_errorbar(aes(ymin = bb.RCOV_DF.lo, ymax = bb.RCOV_DF.hi, color = bb.RCOV_DF.supp), size=1, width=0) +
  geom_point(size = 2.5) + 
  geom_hline(yintercept = 0) +
  coord_flip() +
  scale_x_continuous(breaks = 1:nrow(dat.plt), labels = dat.plt$Spp %>% rev, expand=c(0, 1)) +
  scale_y_continuous(lim = c(min(dat.plt$bb.RCOV_DF.lo), max(dat.plt$bb.HerbGrassVol.hi))) +
#  scale_y_continuous(lim = c(min(dat.plt$bb.RCOV_DF.lo), max(dat.plt$bb.RCOV_DF.hi))) +
  scale_color_manual(values = c("#0072B2", "dark gray", "#D55E00")) +
  ylab(expression(hat(beta)["PSME"])) + xlab(NULL) +
  theme(axis.title.x=element_text(size=25)) +
  guides(color = F)
  
p.POTR5 <- ggplot(dat = dat.plt, aes(x = index, y = bb.RCOV_AS, color = bb.RCOV_AS.supp)) +
  geom_errorbar(aes(ymin = bb.RCOV_AS.lo, ymax = bb.RCOV_AS.hi, color = bb.RCOV_AS.supp), size=1, width=0) +
  geom_point(size = 2.5) + 
  geom_hline(yintercept = 0) +
  coord_flip() +
  scale_x_continuous(breaks = 1:nrow(dat.plt), labels = dat.plt$Spp %>% rev, expand=c(0, 1)) +
  scale_y_continuous(lim = c(min(dat.plt$bb.RCOV_DF.lo), max(dat.plt$bb.HerbGrassVol.hi))) +
#  scale_y_continuous(lim = c(min(dat.plt$bb.RCOV_AS.lo), max(dat.plt$bb.RCOV_AS.hi))) +
  scale_color_manual(values = c("#0072B2", "dark gray", "#D55E00")) +
  ylab(expression(hat(beta)["POTR5"])) + xlab(NULL) +
  theme(axis.title.x=element_text(size=25)) +
  guides(color = F)
  
p.ShrbVol <- ggplot(dat = dat.plt, aes(x = index, y = bb.shvol, color = bb.shvol.supp)) +
  geom_errorbar(aes(ymin = bb.shvol.lo, ymax = bb.shvol.hi, color = bb.shvol.supp), size=1, width=0) +
  geom_point(size = 2.5) + 
  geom_hline(yintercept = 0) +
  coord_flip() +
  scale_x_continuous(breaks = 1:nrow(dat.plt), labels = dat.plt$Spp %>% rev, expand=c(0, 1)) +
  scale_y_continuous(lim = c(min(dat.plt$bb.RCOV_DF.lo), max(dat.plt$bb.HerbGrassVol.hi))) +
#  scale_y_continuous(lim = c(min(dat.plt$bb.shvol.lo), max(dat.plt$bb.shvol.hi))) +
  scale_color_manual(values = c("#0072B2", "dark gray", "#D55E00")) +
  ylab(expression(hat(beta)["ShrbVol"])) + xlab(NULL) +
  theme(axis.title.x=element_text(size=25)) +
  guides(color = F)

p.LadFuel <- ggplot(dat = dat.plt, aes(x = index, y = bb.RSCV_Ladder, color = bb.RSCV_Ladder.supp)) +
  geom_errorbar(aes(ymin = bb.RSCV_Ladder.lo, ymax = bb.RSCV_Ladder.hi, color = bb.RSCV_Ladder.supp), size=1, width=0) +
  geom_point(size = 2.5) + 
  geom_hline(yintercept = 0) +
  coord_flip() +
  scale_x_continuous(breaks = 1:nrow(dat.plt), labels = dat.plt$Spp %>% rev, expand=c(0, 1)) +
  scale_y_continuous(lim = c(min(dat.plt$bb.RCOV_DF.lo), max(dat.plt$bb.HerbGrassVol.hi))) +
#  scale_y_continuous(lim = c(min(dat.plt$bb.RSCV_Ladder.lo), max(dat.plt$bb.RSCV_Ladder.hi))) +
  scale_color_manual(values = c("#0072B2", "dark gray", "#D55E00")) +
  ylab(expression(hat(beta)["LadFuel"])) + xlab(NULL) +
  theme(axis.title.x=element_text(size=25)) +
  guides(color = F)

p.Herb <- ggplot(dat = dat.plt, aes(x = index, y = bb.HerbGrassVol, color = bb.HerbGrassVol.supp)) +
  geom_errorbar(aes(ymin = bb.HerbGrassVol.lo, ymax = bb.HerbGrassVol.hi, color = bb.HerbGrassVol.supp), size=1, width=0) +
  geom_point(size = 2.5) + 
  geom_hline(yintercept = 0) +
  coord_flip() +
  scale_x_continuous(breaks = 1:nrow(dat.plt), labels = dat.plt$Spp %>% rev, expand=c(0, 1)) +
  scale_y_continuous(lim = c(min(dat.plt$bb.RCOV_DF.lo), max(dat.plt$bb.HerbGrassVol.hi))) +
#  scale_y_continuous(lim = c(min(dat.plt$bb.HerbGrassVol.lo), max(dat.plt$bb.HerbGrassVol.hi))) +
  scale_color_manual(values = c("#0072B2", "dark gray", "#D55E00")) +
  ylab(expression(hat(beta)["Herb"])) + xlab(NULL) +
  theme(axis.title.x=element_text(size=25)) +
  guides(color = F)

p <- ggdraw() + 
  draw_plot(p.PIPO, x = 0.05, y = 0, width = 0.1583333, height = 1) +
  draw_plot(p.PSME, x = 0.2083333, y = 0, width = 0.1583333, height = 1) +
  draw_plot(p.POTR5, x = 0.3666667, y = 0, width = 0.1583333, height = 1) +
  draw_plot(p.ShrbVol, x = 0.5250000, y = 0, width = 0.1583333, height = 1) +
  draw_plot(p.LadFuel, x = 0.6833333, y = 0, width = 0.1583333, height = 1) +
  draw_plot(p.Herb, x = 0.8416667, y = 0, width = 0.1583333, height = 1) +
  draw_plot_label("Species", x = 0, y = 0.5, size = 40, angle = 90, hjust = 0)

save_plot("Plot_veg_CanComp&Understory_effects.tiff", p, ncol = 4, nrow = 3, dpi = 200)

p <- ggdraw() + 
  draw_plot(p.theta, x = 0.03, y = 0, width = 0.097, height = 1) +
  draw_plot(p.CanCov, x = 0.127, y = 0, width = 0.097, height = 1) +
  draw_plot(p.CanHt, x = 0.224, y = 0, width = 0.097, height = 1) +
  draw_plot(p.NSnag, x = 0.321, y = 0, width = 0.097, height = 1) +
  draw_plot(p.PIPO, x = 0.418, y = 0, width = 0.097, height = 1) +
  draw_plot(p.PSME, x = 0.515, y = 0, width = 0.097, height = 1) +
  draw_plot(p.POTR5, x = 0.612, y = 0, width = 0.097, height = 1) +
  draw_plot(p.ShrbVol, x = 0.709, y = 0, width = 0.097, height = 1) +
  draw_plot(p.LadFuel, x = 0.806, y = 0, width = 0.097, height = 1) +
  draw_plot(p.Herb, x = 0.903, y = 0, width = 0.097, height = 1) +
  draw_plot_label("Species", x = 0, y = 0.5, size = 40, angle = 90, hjust = 0)
save_plot("Plot_veg_all_effects.tiff", p, ncol = 6, nrow = 3, dpi = 200)
