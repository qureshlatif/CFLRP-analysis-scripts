library(jagsUI)
library(stringr)
library(dplyr)
library(R.utils)
library(ggplot2)
library(cowplot)

setwd("C:/Users/Quresh.Latif/files/projects/FS/CFLRP")
load("Data_compiled.RData")

mod <- loadObject("mod_habitat_d0yr_reduced")
mod_shrub <- loadObject("mod_habitat_d0yr_reduced_unconverged_temp")

spp_trt_effects <- c("WISA", "RECR", "WEWP", "CAFI", "CONI",
                     "AMRO", "GRAJ", "PIGR", "EVGR", "CLNU",
                     "PYNU", "LISP", "BRCR", "DEJU", "WBNU",
                     "OSFL", "STJA", "BHCO", "COFL", "HAWO",
                     "MGWA", "RCKI", "YEWA", "VGSW", "YRWA",
                     "SOSP", "HETH", "VIWA")

# Tabulate parameter estimates
pars <- c("bd.TWI", "bd.heatload", "bd.ForAR", "bd.PACC10_3km", "bd.PACC40_3km", "bd.mnPerArRatio_Opn3km", 
          "bb.CanCov", "bb.CanHt", "bb.NumSnags", "bb.RCOV_PP", "bb.RCOV_DF", "bb.RCOV_AS", "bb.shvol",
          "bb.RSCV_Ladder", "bb.HerbGrassVol")
cols <- (c("", ".lo", ".hi") %>%
           expand.grid(pars, stringsAsFactors = F) %>%
           select(Var2, Var1) %>%
           mutate(Var3 = str_c(Var2, Var1, sep = "")))$Var3
tbl_pars <- matrix(NA, nrow = length(spp.list), ncol = length(cols), dimnames = list(NULL, cols))

#for(i in 1:length(pars)) { # tenp off
parst <- pars[-which(pars %in% c("bd.ForAR", "bb.shvol"))] # temp on
for(i in 1:length(parst)) { # tenp on
  parm <- mod$sims.list[[parst[i]]]
  tbl_pars[, parst[i]] <- apply(parm, 2, median)
  tbl_pars[, str_c(parst[i], ".lo")] <- apply(parm, 2, function(x) quantile(x, prob = 0.025, type = 8))
  tbl_pars[, str_c(parst[i], ".hi")] <- apply(parm, 2, function(x) quantile(x, prob = 0.975, type = 8))
}
rm(parst)
parm <- mod_shrub$sims.list[["bb.shvol"]]
tbl_pars[, "bb.shvol"] <- apply(parm, 2, median)
tbl_pars[, "bb.shvol.lo"] <- apply(parm, 2, function(x) quantile(x, prob = 0.025, type = 8))
tbl_pars[, "bb.shvol.hi"] <- apply(parm, 2, function(x) quantile(x, prob = 0.975, type = 8))

rm(parm)
tbl_pars_all <- tbl_pars
tbl_pars <- tbl_pars[, 10:45]

ind.spp <- c(which(spp.list %in% spp_trt_effects),
             tbl_pars[, dimnames(tbl_pars)[[2]] %>% str_detect(".lo")] %>%
               apply(1, function(x) any(x > 0)) %>%
               which,
             tbl_pars[, dimnames(tbl_pars)[[2]] %>% str_detect(".hi")] %>%
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

cols <- pars[-c(1:3)] %>% str_c(".supp")
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
  draw_plot(p.PACCGap, x = 0.05, y = 0, width = 0.3166667, height = 1) +
  draw_plot(p.PACCOpn, x = 0.3666667, y = 0, width = 0.3166667, height = 1) +
  draw_plot(p.PAROpn, x = 0.6833333, y = 0, width = 0.3166667, height = 1) +
  draw_plot_label("Species", x = 0, y = 0.5, size = 40, angle = 90, hjust = 0)

save_plot("Plot_landscape_effects.tiff", p, ncol = 3, nrow = 3, dpi = 200)

## Point level relationships ##
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
  draw_plot(p.CanCov, x = 0.03, y = 0, width = 0.1077778, height = 1) +
  draw_plot(p.CanHt, x = 0.1377778, y = 0, width = 0.1077778, height = 1) +
  draw_plot(p.NSnag, x = 0.2455556, y = 0, width = 0.1077778, height = 1) +
  draw_plot(p.PIPO, x = 0.3533333, y = 0, width = 0.1077778, height = 1) +
  draw_plot(p.PSME, x = 0.4611111, y = 0, width = 0.1077778, height = 1) +
  draw_plot(p.POTR5, x = 0.5688889, y = 0, width = 0.1077778, height = 1) +
  draw_plot(p.ShrbVol, x = 0.6766667, y = 0, width = 0.1077778, height = 1) +
  draw_plot(p.LadFuel, x = 0.7844444, y = 0, width = 0.1077778, height = 1) +
  draw_plot(p.Herb, x = 0.8922222, y = 0, width = 0.1077778, height = 1) +
  draw_plot_label("Species", x = 0, y = 0.5, size = 40, angle = 90, hjust = 0)

save_plot("Plot_veg_effects.tiff", p, ncol = 6, nrow = 3, dpi = 200)
