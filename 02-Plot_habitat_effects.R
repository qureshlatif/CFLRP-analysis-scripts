library(jagsUI)
library(stringr)
library(dplyr)
library(R.utils)
library(ggplot2)
library(cowplot)

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
pars <- c("bd.PACC10_3km", "bd.PACC40_3km", "bd.mnPerArRatio_Opn3km", "bd.heatload", "bd.TWI",
          "bb.CanHt", "bb.NumSnags", "bb.RCOV_PP", "bb.RCOV_DF", "bb.RCOV_AS", "bb.RSCV_Ladder",
          "bb.HerbGrassVol")
cols <- (c("", ".lo", ".hi") %>%
           expand.grid(pars, stringsAsFactors = F) %>%
           select(Var2, Var1) %>%
           mutate(Var3 = str_c(Var2, Var1, sep = "")))$Var3
tbl_pars <- matrix(NA, nrow = length(spp.list), ncol = length(cols), dimnames = list(NULL, cols))

for(i in 1:length(pars)) {
  parm <- mod$sims.list[[pars[i]]]
  tbl_pars[, pars[i]] <- apply(parm, 2, median)
  tbl_pars[, str_c(pars[i], ".lo")] <- apply(parm, 2, function(x) quantile(x, prob = 0.025, type = 8))
  tbl_pars[, str_c(pars[i], ".hi")] <- apply(parm, 2, function(x) quantile(x, prob = 0.975, type = 8))
}
rm(parm)

#### Plot all species treatment effects ####
dat.plt <- tbl_pars %>% tbl_df() %>%
  mutate(Spp = spp.list) %>%
  arrange(bd.ptrt) %>%
  mutate(index = row_number())

dat.supp <- dat.plt %>% filter(bd.ptrt.lo > 0 |
                                 bd.ptrt.hi < 0 |
                                 bb.trt.lo > 0 |
                                 bb.trt.hi < 0 |
                                 bb.YST.lo > 0 |
                                 bb.YST.hi < 0) %>%
  arrange(bd.ptrt) %>%
  mutate(index = row_number()) %>%
  mutate(bd.ptrt.supp = "none") %>%
  mutate(bd.ptrt.supp = replace(bd.ptrt.supp, which(bd.ptrt.lo > 0), "pos")) %>%
  mutate(bb.trt.supp = "none") %>%
  mutate(bb.trt.supp = replace(bb.trt.supp, which(bb.trt.lo > 0), "pos")) %>%
  mutate(bb.trt.supp = replace(bb.trt.supp, which(bb.trt.hi < 0), "neg")) %>%
  mutate(bb.YST.supp = "none") %>%
  mutate(bb.YST.supp = replace(bb.YST.supp, which(bb.YST.hi < 0), "neg"))
  
pd.ptrt <- ggplot(dat = dat.plt, aes(x = index, y = bd.ptrt)) +
  geom_errorbar(aes(ymin = bd.ptrt.lo, ymax = bd.ptrt.hi), size=1, width=0) +
  geom_point(size = 2.5) + 
  geom_hline(yintercept = 0) +
  coord_flip() +
  scale_x_continuous(breaks = NULL, labels = NULL, expand=c(0, 1)) +
  scale_y_continuous(lim = c(min(dat.plt$bd.ptrt.lo), max(dat.plt$bd.ptrt.hi))) +
  ylab(NULL) + xlab(NULL) +
  theme(axis.title.y=element_text(size=30)) +
  theme(axis.title.x=element_text(size=30)) +
  theme(axis.text.x=element_text(size=15)) +
  theme(axis.text.y=element_text(size=15))

pd.ptrt.supp <- ggplot(dat = dat.supp, aes(x = index, y = bd.ptrt)) +
  geom_errorbar(aes(ymin = bd.ptrt.lo, ymax = bd.ptrt.hi, color = bd.ptrt.supp), size=1, width=0) +
  geom_point(size = 2.5) + 
  geom_hline(yintercept = 0) +
  coord_flip() +
  scale_x_continuous(breaks = 1:nrow(dat.supp), labels = dat.supp$Spp, expand=c(0, 1)) +
  scale_y_continuous(lim = c(min(dat.plt$bd.ptrt.lo), max(dat.plt$bd.ptrt.hi))) +
  scale_color_manual(values = c("#000000", "#D55E00")) +
  ylab(expression(hat(beta)["percent treated"])) + xlab(NULL) +
  theme(axis.title.y=element_text(size=30)) +
  theme(axis.title.x=element_text(size=30)) +
  theme(axis.text.x=element_text(size=15)) +
  theme(axis.text.y=element_text(size=15)) +
  guides(color = F)

pb.trt <- ggplot(dat = dat.plt, aes(x = index, y = bb.trt)) +
  geom_errorbar(aes(ymin = bb.trt.lo, ymax = bb.trt.hi), size=1, width=0) +
  geom_point(size = 2.5) + 
  geom_hline(yintercept = 0) +
  coord_flip() +
  scale_x_continuous(breaks = NULL, labels = NULL, expand=c(0, 1)) +
  scale_y_continuous(lim = c(min(dat.plt$bb.trt.lo), max(dat.plt$bb.trt.hi))) +
  ylab(NULL) + xlab(NULL) +
  theme(axis.title.y=element_text(size=30)) +
  theme(axis.title.x=element_text(size=30)) +
  theme(axis.text.x=element_text(size=15)) +
  theme(axis.text.y=element_text(size=15))

pb.trt.supp <- ggplot(dat = dat.supp, aes(x = index, y = bb.trt)) +
  geom_errorbar(aes(ymin = bb.trt.lo, ymax = bb.trt.hi, color = bb.trt.supp), size=1, width=0) +
  geom_point(size = 2.5) + 
  geom_hline(yintercept = 0) +
  coord_flip() +
  scale_x_continuous(breaks = 1:nrow(dat.supp), labels = dat.supp$Spp, expand=c(0, 1)) +
  scale_y_continuous(lim = c(min(dat.plt$bb.trt.lo), max(dat.plt$bb.trt.hi))) +
  scale_color_manual(values = c("#0072B2", "#000000", "#D55E00")) +
  ylab(expression(hat(beta)["treatment status"])) + xlab(NULL) +
  theme(axis.title.y=element_text(size=30)) +
  theme(axis.title.x=element_text(size=30)) +
  theme(axis.text.x=element_text(size=15)) +
  theme(axis.text.y=element_text(size=15)) +
  guides(color = F)

pb.YST <- ggplot(dat = dat.plt, aes(x = index, y = bb.YST)) +
  geom_errorbar(aes(ymin = bb.YST.lo, ymax = bb.YST.hi), size=1, width=0) +
  geom_point(size = 2.5) + 
  geom_hline(yintercept = 0) +
  coord_flip() +
  scale_x_continuous(breaks = NULL, labels = NULL, expand=c(0, 1)) +
  scale_y_continuous(lim = c(min(dat.plt$bb.YST.lo), max(dat.plt$bb.YST.hi))) +
  ylab(NULL) + xlab(NULL) +
  theme(axis.title.y=element_text(size=30)) +
  theme(axis.title.x=element_text(size=30)) +
  theme(axis.text.x=element_text(size=15)) +
  theme(axis.text.y=element_text(size=15))

pb.YST.supp <- ggplot(dat = dat.supp, aes(x = index, y = bb.YST)) +
  geom_errorbar(aes(ymin = bb.YST.lo, ymax = bb.YST.hi, color = bb.YST.supp), size=1, width=0) +
  geom_point(size = 2.5) + 
  geom_hline(yintercept = 0) +
  coord_flip() +
  scale_x_continuous(breaks = 1:nrow(dat.supp), labels = dat.supp$Spp, expand=c(0, 1)) +
  scale_y_continuous(lim = c(min(dat.plt$bb.YST.lo), max(dat.plt$bb.YST.hi))) +
  scale_color_manual(values = c("#0072B2", "#000000")) +
  ylab(expression(hat(beta)["years since treatment"])) + xlab(NULL) +
  theme(axis.title.y=element_text(size=30)) +
  theme(axis.title.x=element_text(size=30)) +
  theme(axis.text.x=element_text(size=15)) +
  theme(axis.text.y=element_text(size=15)) +
  guides(color = F)

p <- ggdraw() + 
  draw_plot(pd.ptrt, x = 0.05, y = 0.6333333, width = 0.3166667, height = 0.3166667) +
  draw_plot(pd.ptrt.supp, x = 0.05, y = 0, width = 0.3166667, height = 0.6333333) +
  draw_plot(pb.trt, x = 0.3666667, y = 0.6333333, width = 0.3166667, height = 0.3166667) +
  draw_plot(pb.trt.supp, x = 0.3666667, y = 0, width = 0.3166667, height = 0.6333333) +
  draw_plot(pb.YST, x = 0.6833333, y = 0.6333333, width = 0.3166667, height = 0.3166667) +
  draw_plot(pb.YST.supp, x = 0.6833333, y = 0, width = 0.3166667, height = 0.6333333) +
  draw_plot_label(c("Species", "Grid scale", "Point scale"),
                  x = c(0, 0.15, 0.63),
                  y = c(0.5, 0.98, 0.98),
                  size = c(40, 30, 30),
                  angle = c(90, 0, 0),
                  hjust = c(0, 0, 0))

save_plot("Plot_trt_effects.tiff", p, ncol = 3, nrow = 4.5, dpi = 200)
