library(QSLpersonal)
library(jagsUI)
library(stringr)
library(dplyr)
library(R.utils)
library(ggplot2)
library(cowplot)

setwd("C:/Users/Quresh.Latif/files/projects/FS/CFLRP")
load("Data_compiled.RData")

mod <- loadObject("mod_treat2_d0yr")

#_______Script inputs_______#
spp <- "RNSA"
#___________________________#

# Plot grid-level relationship #
x <- seq(min(landscape_data$PctTrt), max(landscape_data$PctTrt), length.out = 11)
z <- (x - mean(landscape_data$PctTrt, na.rm = T)) / sd(landscape_data$PctTrt, na.rm = T)
dat.pred <- data.frame(x = x, z = z, psi = NA, psi.lo = NA, psi.hi = NA)

ind.spp <- which(spp.list == spp)
B0 <- apply(mod$sims.list$d0[, ind.spp, ], 1, mean)
B1 <- mod$sims.list$bd.ptrt[, ind.spp]
B2 <- mod$sims.list$bd.ptrt2[, ind.spp]
for(i in 1:nrow(dat.pred)) {
  psi <- expit(B0 + B1*dat.pred$z[i] + B2*(dat.pred$z[i]^2)) # Rescale to number per 100 acres
  dat.pred$psi[i] <- median(psi)
  dat.pred$psi.lo[i] <- quantile(psi, prob = 0.025, type = 8)
  dat.pred$psi.hi[i] <- quantile(psi, prob = 0.975, type = 8)
}

dat.obs <- landscape_data %>%
  select(gridIndex, YearInd, PctTrt) %>%
  left_join(as.data.frame(Cov) %>%
              select(gridIndex, YearInd) %>%
              mutate(Y = Y.mat[, ind.spp]) %>%
              dplyr::group_by(gridIndex, YearInd) %>%
              summarise(Y = any(Y > 0)*1), by = c("gridIndex", "YearInd")) %>%
  select(PctTrt, Y) %>%
  rename(x = PctTrt)




#### Plot all species treatment effects ####
dat.plt <- tbl_pars %>% tbl_df() %>%
  mutate(Spp = spp.list) %>%
  arrange(bd.ptrt) %>%
  mutate(index = row_number())

dat.supp <- dat.plt %>% filter(bd.ptrt.lo > 0 |
                                 bd.ptrt.hi < 0 |
                                 bd.YST.lo > 0 |
                                 bd.YST.hi < 0 |
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

pd.YST <- ggplot(dat = dat.plt, aes(x = index, y = bd.YST)) +
  geom_errorbar(aes(ymin = bd.YST.lo, ymax = bd.YST.hi), size=1, width=0) +
  geom_point(size = 2.5) + 
  geom_hline(yintercept = 0) +
  coord_flip() +
  scale_x_continuous(breaks = NULL, labels = NULL, expand=c(0, 1)) +
  scale_y_continuous(lim = c(min(dat.plt$bd.YST.lo), max(dat.plt$bd.YST.hi))) +
  ylab(NULL) + xlab(NULL) +
  theme(axis.title.y=element_text(size=30)) +
  theme(axis.title.x=element_text(size=30)) +
  theme(axis.text.x=element_text(size=15)) +
  theme(axis.text.y=element_text(size=15))

pd.YST.supp <- ggplot(dat = dat.supp, aes(x = index, y = bd.YST)) +
  geom_errorbar(aes(ymin = bd.YST.lo, ymax = bd.YST.hi), size=1, width=0) +
  geom_point(size = 2.5) + 
  geom_hline(yintercept = 0) +
  coord_flip() +
  scale_x_continuous(breaks = 1:nrow(dat.supp), labels = dat.supp$Spp, expand=c(0, 1)) +
  scale_y_continuous(lim = c(min(dat.plt$bd.YST.lo), max(dat.plt$bd.YST.hi))) +
  ylab(expression(hat(beta)["mean years since treatment"])) + xlab(NULL) +
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
  draw_plot(pd.ptrt, x = 0.05, y = 0.6333333, width = 0.2375, height = 0.3166667) +
  draw_plot(pd.ptrt.supp, x = 0.05, y = 0, width = 0.2375, height = 0.6333333) +
  draw_plot(pd.YST, x = 0.2875, y = 0.6333333, width = 0.2375, height = 0.3166667) +
  draw_plot(pd.YST.supp, x = 0.2875, y = 0, width = 0.2375, height = 0.6333333) +
  draw_plot(pb.trt, x = 0.525, y = 0.6333333, width = 0.2375, height = 0.3166667) +
  draw_plot(pb.trt.supp, x = 0.525, y = 0, width = 0.2375, height = 0.6333333) +
  draw_plot(pb.YST, x = 0.7625, y = 0.6333333, width = 0.2375, height = 0.3166667) +
  draw_plot(pb.YST.supp, x = 0.7625, y = 0, width = 0.2375, height = 0.6333333) +
  draw_plot_label(c("Species", "Grid scale", "Point scale"),
                  x = c(0, 0.25, 0.7),
                  y = c(0.5, 0.98, 0.98),
                  size = c(40, 30, 30),
                  angle = c(90, 0, 0),
                  hjust = c(0, 0, 0))

save_plot("Plot_trt_effects.tiff", p, ncol = 4, nrow = 4.5, dpi = 200)
