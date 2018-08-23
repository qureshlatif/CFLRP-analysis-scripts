library(jagsUI)
library(stringr)
library(dplyr)
library(R.utils)
library(ggplot2)
library(cowplot)

setwd("C:/Users/Quresh.Latif/files/projects/FS/CFLRP")
load("Data_compiled.RData")

mod <- loadObject("mod_treatment_d0yr")

# Tabulate parameter estimates
cols <- c("bd.ptrt", "bd.ptrt.lo", "bd.ptrt.hi", "bb.trt", "bb.trt.lo", "bb.trt.hi")
tbl_pars <- matrix(NA, nrow = length(spp.list), ncol = length(cols), dimnames = list(NULL, cols))

parm <- mod$sims.list[["bd.ptrt"]]
tbl_pars[, "bd.ptrt"] <- apply(parm, 2, median)
tbl_pars[, "bd.ptrt.lo"] <- apply(parm, 2, function(x) quantile(x, prob = 0.025, type = 8))
tbl_pars[, "bd.ptrt.hi"] <- apply(parm, 2, function(x) quantile(x, prob = 0.975, type = 8))

parm <- mod$sims.list[["bb.trt"]]
tbl_pars[, "bb.trt"] <- apply(parm, 2, median)
tbl_pars[, "bb.trt.lo"] <- apply(parm, 2, function(x) quantile(x, prob = 0.025, type = 8))
tbl_pars[, "bb.trt.hi"] <- apply(parm, 2, function(x) quantile(x, prob = 0.975, type = 8))

rm(parm)

# Plot species treatment effects #
tbl_pars <- tbl_pars %>% tbl_df() %>%
  mutate(Spp = spp.list) %>%
  arrange(bd.ptrt) %>%
  mutate(index = row_number())

min.lo <- min(tbl_pars$bd.ptrt.lo)
max.hi <- max(tbl_pars$bd.ptrt.hi)

pd <- ggplot(dat = tbl_pars, aes(x = index, y = bd.ptrt)) +
  geom_errorbar(aes(ymin = bd.ptrt.lo, ymax = bd.ptrt.hi), size=1, width=0) +
  geom_point(size = 2.5) + 
  geom_hline(yintercept = 0) +
  coord_flip() +
  scale_x_continuous(breaks = NULL, labels = NULL, expand=c(0, 1)) +
  scale_y_continuous(lim = c(min.lo, max.hi)) +
  ylab(expression(hat(beta)["percent treated"])) + xlab(NULL) +
  theme(axis.title.y=element_text(size=30)) +
  theme(axis.title.x=element_text(size=30)) +
  theme(axis.text.x=element_text(size=15)) +
  theme(axis.text.y=element_text(size=15))

tbl_pars <- tbl_pars %>%
  arrange(bb.trt) %>%
  mutate(index = row_number())
min.lo <- min(tbl_pars$bb.trt.lo)
max.hi <- max(tbl_pars$bb.trt.hi)
pb <- ggplot(dat = tbl_pars, aes(x = index, y = bb.trt)) +
  geom_errorbar(aes(ymin = bb.trt.lo, ymax = bb.trt.hi), size=1, width=0) +
  geom_point(size = 2.5) + 
  geom_hline(yintercept = 0) +
  coord_flip() +
  scale_x_continuous(breaks = NULL, labels = NULL, expand=c(0, 1)) +
  scale_y_continuous(lim = c(min.lo, max.hi)) +
  ylab(expression(hat(beta)["treated"])) + xlab(NULL) +
  theme(axis.title.y=element_text(size=30)) +
  theme(axis.title.x=element_text(size=30)) +
  theme(axis.text.x=element_text(size=15)) +
  theme(axis.text.y=element_text(size=15))

# Plot species richness #
gridID <- Cov[, "gridIndex"]
PctTrt.d <- tapply(Trt.b, gridID, mean, na.rm = T) # Grid-level values

SPR <- mod$sims.list$SR.grid
SPR.median <- apply(SPR,c(2,3),median)
SPR.95lo <- apply(SPR,c(2,3),function(x) quantile(x,prob=0.025,type=8))
SPR.95hi <- apply(SPR,c(2,3),function(x) quantile(x,prob=0.975,type=8))

X <- as.numeric(E[,1,])
SPR <- as.numeric(SPR.median)
SPR.lo <- as.numeric(SPR.95lo)
SPR.hi <- as.numeric(SPR.95hi)
dat <- data.frame(cbind(X,SPR,SPR.lo,SPR.hi))

x.labs <- round((c(-2,-1,0,1,2,3)*dNBR.sd)+dNBR.mn,digits=1)
p <- ggplot(data = dat,aes(x=X,y=SPR)) + 
  geom_point(alpha=0.3) + stat_smooth(method=lm,level=F,colour="black",size=1.5) + 
  geom_errorbar(aes(ymin=SPR.lo,ymax=SPR.hi),width=0.1,alpha=0.3) +
  scale_y_continuous(breaks=c(10,20,30)) +
  scale_x_continuous(limits=c(-2,3.1),breaks=c(-2,-1,0,1,2,3),labels=x.labs) +
  labs(x=expression(paste(Delta,"NBR",sep="")),y=NULL) +
  theme(axis.title.x=element_text(size=40)) +
  theme(axis.title.y=element_blank()) +
  theme(axis.text.x=element_text(size=40)) +
  theme(axis.text.y=element_text(size=40)) +
  theme(plot.margin = unit(c(0.1,1,1,0.5),"cm")) +
  geom_text(aes(x=-2,y=33),label="(b)",size=10)


Trt.b <- Cov[, "Trt_stat"] # Point-level values
