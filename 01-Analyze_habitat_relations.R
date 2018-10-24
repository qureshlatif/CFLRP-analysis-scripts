library(jagsUI)
library(stringr)
library(dplyr)

#setwd("/home/RMBO.LOCAL/quresh.latif/CFLRP")
setwd("C:/Users/Quresh.Latif/files/projects/FS/CFLRP")
load("Data_compiled.RData")

#### Script inputs ####
model.file <- "CFLRP-analysis-scripts/model_habitat_d0yr.jags"
RESQ.model <- "mod_RESQ_treatment_reduced"

# MCMC values
nc <- 3 # number of chains
nb <- 5 #5000 # burn in
ni <- 10 #30000 # number of iterations
nt <- 1 #10 # thinning

save.out <- "mod_habitat_d0yr"
##########################

# Data objects to send to JAGS
data <- list("Y", "TPeriod", "gridID", "yearID", "n.grid", "n.year", "n.point_year", "n.spp",
             "ccov.b", "ccov.means", "ccov.sd", "ccov.b.missing",
             "CanHt.b","CanHt.means", "CanHt.sd", "CanHt.b.missing",
             "NumSnags.b","NumSnags.means", "NumSnags.sd", "NumSnags.b.missing",
             "RCOV_PP.b","RCOV_PP.means", "RCOV_PP.sd", "RCOV_PP.b.missing",
             "RCOV_DF.b","RCOV_DF.means", "RCOV_DF.sd", "RCOV_DF.b.missing",
             "RCOV_AS.b","RCOV_AS.means", "RCOV_AS.sd", "RCOV_AS.b.missing",
             "shvol.b", "shvol.means", "shvol.sd", "shvol.b.missing",
             "RSCV_Ladder.b", "RSCV_Ladder.means", "RSCV_Ladder.sd", "RSCV_Ladder.b.missing",
             "RSCV_Ber.b", "RSCV_Ber.means", "RSCV_Ber.sd", "RSCV_Ber.b.missing",
             "HerbGrassVol.b", "HerbGrassVol.means", "HerbGrassVol.sd", "HerbGrassVol.b.missing",
             "SOHtRatio.b", "SOHtRatio.means", "SOHtRatio.sd", "SOHtRatio.b.missing",
             "PACC10_3km.d", "mnPtchAr_Gap3km.d",
             "mnPerArRatio_Gap3km.d",
             "NNdist_Gap3km.d", "mnPerArRatio_Gap3km.d.missing", "mnPerArRatio_Gap3km.d.mean", "mnPerArRatio_Gap3km.d.sd",
             "NNdist_Gap3km.d.missing", "NNdist_Gap3km.d.min", "NNdist_Gap3km.d.max",
             "PACC40_3km.d", "mnPtchAr_Opn3km.d",
             "mnPerArRatio_Opn3km.d",
             "NNdist_Opn3km.d", "NNdist_Opn3km.d.missing", "NNdist_Opn3km.d.min", "NNdist_Opn3km.d.max",
             "TWIP.d","heatload.d","TWI.d","Rdens.d",
             "RESQ.b", "RESQ.wts", "RESQ.bm",
             "DOY.b", "Time.b")

# Stuff to save from JAGS
parameters <- c("omega", "rho.ab", "rho.bd",
                "alpha0", "sigma.a0", "beta0", "sigma.b0", "delta0", "sigma.d0",
                
                "Betad.yr", "sigma.Betad.yr",
                "Betad.PACC10_3km", "sigma.Betad.PACC10_3km",
                "Betad.PACC10_3km2", "sigma.Betad.PACC10_3km2",
                "Betad.mnPtchAr_Gap3km", "sigma.Betad.mnPtchAr_Gap3km",
                "Betad.mnPerArRatio_Gap3km", "sigma.Betad.mnPerArRatio_Gap3km",
                "Betad.NNdist_Gap3km", "sigma.Betad.NNdist_Gap3km",
                "Betad.PACC40_3km", "sigma.Betad.PACC40_3km",
                "Betad.PACC40_3km2", "sigma.Betad.PACC40_3km2",
                "Betad.mnPtchAr_Opn3km", "sigma.Betad.mnPtchAr_Opn3km",
                "Betad.mnPerArRatio_Opn3km", "sigma.Betad.mnPerArRatio_Opn3km",
                "Betad.NNdist_Opn3km", "sigma.Betad.NNdist_Opn3km",
                "Betad.TWIP", "sigma.Betad.TWIP","Betad.TWI", "sigma.Betad.TWI",
                "Betad.heatload", "sigma.Betad.heatload","Betad.Rdens", "sigma.Betad.Rdens",
                
                "Betab.CanCov", "sigma.Betab.CanCov",
                "Betab.CanCov2", "sigma.Betab.CanCov2",
                "Betab.CanHt", "sigma.Betab.CanHt",
                "Betab.NumSnags", "sigma.Betab.NumSnags",
                "Betab.RCOV_PP", "sigma.Betab.RCOV_PP",
                "Betab.RCOV_DF", "sigma.Betab.RCOV_DF",
                "Betab.RCOV_AS", "sigma.Betab.RCOV_AS",
                "Betab.ShrubVol", "sigma.Betab.ShrubVol",
                "Betab.RSCV_Ladder", "sigma.Betab.RSCV_Ladder",
                "Betab.RSCV_Ber", "sigma.Betab.RSCV_Ber",
                "Betab.HerbGrassVol", "sigma.Betab.HerbGrassVol",
                "Betab.SOHtRatio", "sigma.Betab.SOHtRatio",
                "Betaa.Time", "sigma.Betaa.Time", "Betaa.Time2", "sigma.Betaa.Time2",
                "Betaa.DOY", "sigma.Betaa.DOY",
                "Betaa.CCov", "sigma.Betaa.CCov",
                "Betaa.SHVol", "sigma.Betaa.SHVol",
                "Betab.RESQ", "sigma.Betab.RESQ",
                
                "d0", "b0", "a0",
                "bd.yr", # For year effect
                "bd.pers", # For persistence effect
                "bd.PACC10_3km", "bd.PACC10_3km2", "bd.mnPtchAr_Gap3km", "bd.mnPerArRatio_Gap3km", "bd.NNdist_Gap3km",
                "bd.PACC40_3km", "bd.PACC40_3km2", "bd.mnPtchAr_Opn3km", "bd.mnPerArRatio_Opn3km", "bd.NNdist_Opn3km",
                "bd.TWIP", "bd.heatload", "bd.TWI", "bd.Rdens",
                "bb.CanCov", "bb.CanCov2", "bb.CanHt", "bb.NumSnags", "bb.RCOV_PP", "bb.RCOV_DF",
                "bb.RCOV_AS", "bb.shvol", "bb.RSCV_Ladder", "bb.RSCV_Ber", "bb.HerbGrassVol", "bb.SOHtRatio",
                "bb.RESQ", "ba.Time", "ba.Time2", "ba.DOY", "ba.ccov", "ba.shvol",
                
                "SR.grid", "SR.point")

# Function for setting initial values in JAGS
inits <- function()
  list(z=z.init, u=u.init, w=w.init, tvar.sigma.a0 = rnorm(1), tvar.sigma.b0 = rnorm(1), tvar.sigma.d0 = rnorm(n.year),
       tvar.Betad.yr = rnorm(n.year - 1),
       tvar.Betad.PACC10_3km = rnorm(1), tvar.Betad.PACC10_3km2 = rnorm(1), tvar.Betad.mnPtchAr_Gap3km = rnorm(1),
       tvar.Betad.mnPerArRatio_Gap3km = rnorm(1), tvar.Betad.NNdist_Gap3km = rnorm(1),
       tvar.Betad.PACC40_3km = rnorm(1), tvar.Betad.PACC40_3km2 = rnorm(1), tvar.Betad.mnPtchAr_Opn3km = rnorm(1),
       tvar.Betad.mnPerArRatio_Opn3km = rnorm(1), tvar.Betad.NNdist_Opn3km = rnorm(1),
       tvar.Betad.TWIP = rnorm(1), tvar.Betad.Rdens = rnorm(1),
       tvar.Betad.TWI = rnorm(1), tvar.Betad.heatload = rnorm(1), 
       tvar.Betab.CanCov = rnorm(1), tvar.Betab.CanCov2 = rnorm(1), tvar.Betab.CanHt = rnorm(1),
       tvar.Betab.NumSnags = rnorm(1), tvar.Betab.RCOV_PP = rnorm(1),
       tvar.Betab.RCOV_DF = rnorm(1), tvar.Betab.RCOV_AS = rnorm(1),
       tvar.Betab.ShrubVol = rnorm(1),
       tvar.Betab.RSCV_Ladder = rnorm(1), tvar.Betab.RSCV_Ber = rnorm(1),
       tvar.Betab.HerbGrassVol = rnorm(1), tvar.Betab.SOHtRatio = rnorm(1), tvar.Betab.RESQ = rnorm(1),
       tvar.Betaa.Time = rnorm(1), tvar.Betaa.Time2 = rnorm(1),
       tvar.Betaa.DOY = rnorm(1),
       tvar.Betaa.CCov = rnorm(1), tvar.Betaa.SHVol = rnorm(1))

# Detection data #
Y <- Y.mat
TPeriod <- TR.mat
Cov <- Cov
gridID <- Cov[, "gridIndex"]
yearID <- Cov[, "YearInd"]
n.grid <- max(gridID)
n.year <- max(yearID)
n.point_year <- dim(Y)[1]
n.spp <- dim(Y)[2]

# Covariates #
ccov.b <- Cov[, "CanCov"] %>% (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)) # Point-level values
ccov.means <- tapply(ccov.b, gridID, mean, na.rm = T) # Grid-level means for imputing missing values
ccov.sd <- tapply(ccov.b, gridID, sd, na.rm = T) # Grid-level SDs for imputing missing values
ccov.sd[which(ccov.sd == 0)] <- sd(ccov.b, na.rm = T) # Zeros won't work for SD!
ccov.b.missing <- is.na(ccov.b) %>% as.integer # Index missing values to be imputed
ccov.b[is.na(ccov.b)] <- 0

CanHt.b <- Cov[, "CanHt"] %>% (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)) # Point-level values
CanHt.means <- tapply(CanHt.b, gridID, mean, na.rm = T) # Grid-level means for imputing missing values
CanHt.sd <- tapply(CanHt.b, gridID, sd, na.rm = T) # Grid-level SDs for imputing missing values
CanHt.sd[which(CanHt.sd == 0)] <- sd(CanHt.b, na.rm = T) # Zeros won't work for SD!
CanHt.b.missing <- is.na(CanHt.b) %>% as.integer # Index missing values to be imputed
CanHt.b[is.na(CanHt.b)] <- 0

NumSnags.b <- Cov[, "NumSnags"] %>% (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)) # Point-level values
NumSnags.means <- tapply(NumSnags.b, gridID, mean, na.rm = T) # Grid-level means for imputing missing values
NumSnags.sd <- tapply(NumSnags.b, gridID, sd, na.rm = T) # Grid-level SDs for imputing missing values
NumSnags.sd[which(NumSnags.sd == 0)] <- sd(NumSnags.b, na.rm = T) # Zeros won't work for SD!
NumSnags.b.missing <- is.na(NumSnags.b) %>% as.integer # Index missing values to be imputed
NumSnags.b[is.na(NumSnags.b)] <- 0

RCOV_PP.b <- Cov[, "RCOV_PP"] %>% (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)) # Point-level values
RCOV_PP.means <- tapply(RCOV_PP.b, gridID, mean, na.rm = T) # Grid-level means for imputing missing values
RCOV_PP.sd <- tapply(RCOV_PP.b, gridID, sd, na.rm = T) # Grid-level SDs for imputing missing values
RCOV_PP.sd[which(RCOV_PP.sd == 0)] <- sd(RCOV_PP.b, na.rm = T) # Zeros won't work for SD!
RCOV_PP.b.missing <- is.na(RCOV_PP.b) %>% as.integer # Index missing values to be imputed
RCOV_PP.b[is.na(RCOV_PP.b)] <- 0

RCOV_DF.b <- Cov[, "RCOV_DF"] %>% (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)) # Point-level values
RCOV_DF.means <- tapply(RCOV_DF.b, gridID, mean, na.rm = T) # Grid-level means for imputing missing values
RCOV_DF.sd <- tapply(RCOV_DF.b, gridID, sd, na.rm = T) # Grid-level SDs for imputing missing values
RCOV_DF.sd[which(RCOV_DF.sd == 0)] <- sd(RCOV_DF.b, na.rm = T) # Zeros won't work for SD!
RCOV_DF.b.missing <- is.na(RCOV_DF.b) %>% as.integer # Index missing values to be imputed
RCOV_DF.b[is.na(RCOV_DF.b)] <- 0

RCOV_AS.b <- Cov[, "RCOV_AS"] %>% (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)) # Point-level values
RCOV_AS.means <- tapply(RCOV_AS.b, gridID, mean, na.rm = T) # Grid-level means for imputing missing values
RCOV_AS.sd <- tapply(RCOV_AS.b, gridID, sd, na.rm = T) # Grid-level SDs for imputing missing values
RCOV_AS.sd[which(RCOV_AS.sd == 0)] <- sd(RCOV_AS.b, na.rm = T) # Zeros won't work for SD!
RCOV_AS.b.missing <- is.na(RCOV_AS.b) %>% as.integer # Index missing values to be imputed
RCOV_AS.b[is.na(RCOV_AS.b)] <- 0

shvol.b <- Cov[, "ShrubVol"] %>% # area covered X volume
  (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)) # Point-level values
shvol.means <- tapply(shvol.b, gridID, mean, na.rm = T) # Grid-level means for imputing missing values
shvol.sd <- tapply(shvol.b, gridID, sd, na.rm = T) # Grid-level SDs for imputing missing values
shvol.b.missing <- is.na(shvol.b) %>% as.integer # Index missing values to be imputed
shvol.b[is.na(shvol.b)] <- 0

RSCV_Ladder.b <- Cov[, "RSCV_Ladder"] %>% (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)) # Point-level values
RSCV_Ladder.means <- tapply(RSCV_Ladder.b, gridID, mean, na.rm = T) # Grid-level means for imputing missing values
RSCV_Ladder.sd <- tapply(RSCV_Ladder.b, gridID, sd, na.rm = T) # Grid-level SDs for imputing missing values
RSCV_Ladder.sd[which(RSCV_Ladder.sd == 0)] <- sd(RSCV_Ladder.b, na.rm = T) # Zeros won't work for SD!
RSCV_Ladder.b.missing <- is.na(RSCV_Ladder.b) %>% as.integer # Index missing values to be imputed
RSCV_Ladder.b[is.na(RSCV_Ladder.b)] <- 0

RSCV_Ber.b <- Cov[, "RSCV_Ber"] %>% (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)) # Point-level values
RSCV_Ber.means <- tapply(RSCV_Ber.b, gridID, mean, na.rm = T) # Grid-level means for imputing missing values
RSCV_Ber.sd <- tapply(RSCV_Ber.b, gridID, sd, na.rm = T) # Grid-level SDs for imputing missing values
RSCV_Ber.sd[which(RSCV_Ber.sd == 0)] <- sd(RSCV_Ber.b, na.rm = T) # Zeros won't work for SD!
RSCV_Ber.b.missing <- is.na(RSCV_Ber.b) %>% as.integer # Index missing values to be imputed
RSCV_Ber.b[is.na(RSCV_Ber.b)] <- 0

HerbGrassVol.b <- Cov[, "HerbGrassVol"] %>% (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)) # Point-level values
HerbGrassVol.means <- tapply(HerbGrassVol.b, gridID, mean, na.rm = T) # Grid-level means for imputing missing values
HerbGrassVol.sd <- tapply(HerbGrassVol.b, gridID, sd, na.rm = T) # Grid-level SDs for imputing missing values
HerbGrassVol.sd[which(HerbGrassVol.sd == 0)] <- sd(HerbGrassVol.b, na.rm = T) # Zeros won't work for SD!
HerbGrassVol.b.missing <- is.na(HerbGrassVol.b) %>% as.integer # Index missing values to be imputed
HerbGrassVol.b[is.na(HerbGrassVol.b)] <- 0

SOHtRatio.b <- Cov[, "SOHtRatio"] %>% (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)) # Point-level values
SOHtRatio.means <- tapply(SOHtRatio.b, gridID, mean, na.rm = T) # Grid-level means for imputing missing values
SOHtRatio.sd <- tapply(SOHtRatio.b, gridID, sd, na.rm = T) # Grid-level SDs for imputing missing values
SOHtRatio.sd[which(SOHtRatio.sd == 0)] <- sd(SOHtRatio.b, na.rm = T) # Zeros won't work for SD!
SOHtRatio.b.missing <- is.na(SOHtRatio.b) %>% as.integer # Index missing values to be imputed
SOHtRatio.b[is.na(SOHtRatio.b)] <- 0

PACC10_3km.d <- matrix(NA, nrow = max(gridID), ncol = max(yearID))
PACC10_3km.d[landscape_data %>% filter(YearInd == 1) %>% pull(gridIndex), 1] <-
  landscape_data %>% filter(YearInd == 1) %>% pull(PACC10_3km)
PACC10_3km.d[landscape_data %>% filter(YearInd == 2) %>% pull(gridIndex), 2] <-
  landscape_data %>% filter(YearInd == 2) %>% pull(PACC10_3km)
PACC10_3km.d[landscape_data %>% filter(YearInd == 3) %>% pull(gridIndex), 3] <-
  landscape_data %>% filter(YearInd == 3) %>% pull(PACC10_3km)
PACC10_3km.d <- PACC10_3km.d %>%
  (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T))
PACC10_3km.d[which(is.na(PACC10_3km.d))] <- 0

mnPtchAr_Gap3km.d <- matrix(NA, nrow = max(gridID), ncol = max(yearID))
mnPtchAr_Gap3km.d[landscape_data %>% filter(YearInd == 1) %>% pull(gridIndex), 1] <-
  landscape_data %>% filter(YearInd == 1) %>% pull(mnPtchAr_Gap3km)
mnPtchAr_Gap3km.d[landscape_data %>% filter(YearInd == 2) %>% pull(gridIndex), 2] <-
  landscape_data %>% filter(YearInd == 2) %>% pull(mnPtchAr_Gap3km)
mnPtchAr_Gap3km.d[landscape_data %>% filter(YearInd == 3) %>% pull(gridIndex), 3] <-
  landscape_data %>% filter(YearInd == 3) %>% pull(mnPtchAr_Gap3km)
mnPtchAr_Gap3km.d <- mnPtchAr_Gap3km.d %>%
  (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T))
mnPtchAr_Gap3km.d[which(is.na(mnPtchAr_Gap3km.d))] <- 0

mnPerArRatio_Gap3km.d <- matrix(NA, nrow = max(gridID), ncol = max(yearID))
mnPerArRatio_Gap3km.d[landscape_data %>% filter(YearInd == 1) %>% pull(gridIndex), 1] <-
  landscape_data %>% filter(YearInd == 1) %>% pull(mnPerArRatio_Gap3km)
mnPerArRatio_Gap3km.d[landscape_data %>% filter(YearInd == 2) %>% pull(gridIndex), 2] <-
  landscape_data %>% filter(YearInd == 2) %>% pull(mnPerArRatio_Gap3km)
mnPerArRatio_Gap3km.d[landscape_data %>% filter(YearInd == 3) %>% pull(gridIndex), 3] <-
  landscape_data %>% filter(YearInd == 3) %>% pull(mnPerArRatio_Gap3km)
mnPerArRatio_Gap3km.d <- mnPerArRatio_Gap3km.d %>%
  (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T))
mnPerArRatio_Gap3km.d.missing <- is.na(mnPerArRatio_Gap3km.d)*1
mnPerArRatio_Gap3km.d.mean <- mean(mnPerArRatio_Gap3km.d, na.rm = T)
mnPerArRatio_Gap3km.d.sd <- sd(mnPerArRatio_Gap3km.d, na.rm = T)
mnPerArRatio_Gap3km.d[which(is.na(mnPerArRatio_Gap3km.d))] <- 0

NNdist_Gap3km.d <- matrix(NA, nrow = max(gridID), ncol = max(yearID))
NNdist_Gap3km.d[landscape_data %>% filter(YearInd == 1) %>% pull(gridIndex), 1] <-
  landscape_data %>% filter(YearInd == 1) %>% pull(NNdist_Gap3km)
NNdist_Gap3km.d[landscape_data %>% filter(YearInd == 2) %>% pull(gridIndex), 2] <-
  landscape_data %>% filter(YearInd == 2) %>% pull(NNdist_Gap3km)
NNdist_Gap3km.d[landscape_data %>% filter(YearInd == 3) %>% pull(gridIndex), 3] <-
  landscape_data %>% filter(YearInd == 3) %>% pull(NNdist_Gap3km)
NNdist_Gap3km.d <- NNdist_Gap3km.d %>%
  (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T))
NNdist_Gap3km.d.missing <- is.na(NNdist_Gap3km.d)*1
NNdist_Gap3km.d.min <- ((3^2 + 3^2)^0.5) / 2
NNdist_Gap3km.d.max <- ((3^2 + 3^2)^0.5)
NNdist_Gap3km.d[which(is.na(NNdist_Gap3km.d))] <- 0

PACC40_3km.d <- matrix(NA, nrow = max(gridID), ncol = max(yearID))
PACC40_3km.d[landscape_data %>% filter(YearInd == 1) %>% pull(gridIndex), 1] <-
  landscape_data %>% filter(YearInd == 1) %>% pull(PACC40_3km)
PACC40_3km.d[landscape_data %>% filter(YearInd == 2) %>% pull(gridIndex), 2] <-
  landscape_data %>% filter(YearInd == 2) %>% pull(PACC40_3km)
PACC40_3km.d[landscape_data %>% filter(YearInd == 3) %>% pull(gridIndex), 3] <-
  landscape_data %>% filter(YearInd == 3) %>% pull(PACC40_3km)
PACC40_3km.d <- PACC40_3km.d %>%
  (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T))
PACC40_3km.d[which(is.na(PACC40_3km.d))] <- 0

mnPtchAr_Opn3km.d <- matrix(NA, nrow = max(gridID), ncol = max(yearID))
mnPtchAr_Opn3km.d[landscape_data %>% filter(YearInd == 1) %>% pull(gridIndex), 1] <-
  landscape_data %>% filter(YearInd == 1) %>% pull(mnPtchAr_Opn3km)
mnPtchAr_Opn3km.d[landscape_data %>% filter(YearInd == 2) %>% pull(gridIndex), 2] <-
  landscape_data %>% filter(YearInd == 2) %>% pull(mnPtchAr_Opn3km)
mnPtchAr_Opn3km.d[landscape_data %>% filter(YearInd == 3) %>% pull(gridIndex), 3] <-
  landscape_data %>% filter(YearInd == 3) %>% pull(mnPtchAr_Opn3km)
mnPtchAr_Opn3km.d <- mnPtchAr_Opn3km.d %>%
  (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T))
mnPtchAr_Opn3km.d[which(is.na(mnPtchAr_Opn3km.d))] <- 0

mnPerArRatio_Opn3km.d <- matrix(NA, nrow = max(gridID), ncol = max(yearID))
mnPerArRatio_Opn3km.d[landscape_data %>% filter(YearInd == 1) %>% pull(gridIndex), 1] <-
  landscape_data %>% filter(YearInd == 1) %>% pull(mnPerArRatio_Opn3km)
mnPerArRatio_Opn3km.d[landscape_data %>% filter(YearInd == 2) %>% pull(gridIndex), 2] <-
  landscape_data %>% filter(YearInd == 2) %>% pull(mnPerArRatio_Opn3km)
mnPerArRatio_Opn3km.d[landscape_data %>% filter(YearInd == 3) %>% pull(gridIndex), 3] <-
  landscape_data %>% filter(YearInd == 3) %>% pull(mnPerArRatio_Opn3km)
mnPerArRatio_Opn3km.d <- mnPerArRatio_Opn3km.d %>%
  (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T))
mnPerArRatio_Opn3km.d.missing <- is.na(mnPerArRatio_Opn3km.d)*1
mnPerArRatio_Opn3km.d.mean <- mean(mnPerArRatio_Opn3km.d, na.rm = T)
mnPerArRatio_Opn3km.d.sd <- sd(mnPerArRatio_Opn3km.d, na.rm = T)
mnPerArRatio_Opn3km.d[which(is.na(mnPerArRatio_Opn3km.d))] <- 0

NNdist_Opn3km.d <- matrix(NA, nrow = max(gridID), ncol = max(yearID))
NNdist_Opn3km.d[landscape_data %>% filter(YearInd == 1) %>% pull(gridIndex), 1] <-
  landscape_data %>% filter(YearInd == 1) %>% pull(NNdist_Opn3km)
NNdist_Opn3km.d[landscape_data %>% filter(YearInd == 2) %>% pull(gridIndex), 2] <-
  landscape_data %>% filter(YearInd == 2) %>% pull(NNdist_Opn3km)
NNdist_Opn3km.d[landscape_data %>% filter(YearInd == 3) %>% pull(gridIndex), 3] <-
  landscape_data %>% filter(YearInd == 3) %>% pull(NNdist_Opn3km)
NNdist_Opn3km.d <- NNdist_Opn3km.d %>%
  (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T))
NNdist_Opn3km.d.missing <- is.na(NNdist_Opn3km.d)*1
NNdist_Opn3km.d.min <- ((3^2 + 3^2)^0.5) / 2
NNdist_Opn3km.d.max <- ((3^2 + 3^2)^0.5)
NNdist_Opn3km.d[which(is.na(NNdist_Opn3km.d))] <- 0

#TWIP.b <- Cov[, "TWIP"] %>% (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)) # Point-level values
TWIP.d <- tapply(Cov[, "TWIP"], gridID, mean, na.rm = T) %>% # Grid-level values
  (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T))

heatload.d <- tapply(Cov[, "heatload"], gridID, mean, na.rm = T) %>% # Grid-level values
  (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T))

TWI.d <- tapply(Cov[, "TWI"], gridID, mean, na.rm = T) %>% # Grid-level values
  (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T))

Rdens.d <- tapply(Cov[, "Rdens"], gridID, mean, na.rm = T) %>% # Grid-level values
  (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T))

library(R.utils)
mod <- loadObject(RESQ.model)
RESQ.b <- mod$sims.list$N
rm(mod) #N, cl.size.samps
# z-score
RESQ.b <- RESQ.b %>%
  (function(x) (x - mean(x)) / sd(x))
RESQ.wts <- rep(1, dim(RESQ.b)[1])
RESQ.bm <- apply(RESQ.b, 2, median)

DOY.b <- Cov[, "DayOfYear"] %>% (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)) # Point-level values

Time.b <- Cov[, "Time"] %>% (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)) # Point-level values

# Assemble the initial values for JAGS.
u.init <- Y
dimnames(u.init) <- NULL
z.init <- array(NA, dim = c(n.grid, n.spp, n.year))
for(sp in 1:n.spp) {
  for(t in 1:n.year) {
    vals <- tapply(u.init[, sp], gridID, max)
    z.init[, sp, t] <- vals
  }
}
w.init <- apply(Y, 2, max) %>% as.integer

# Fit model
st.time <- Sys.time()
out <- jagsUI(data, inits, parameters.to.save = parameters, model.file, n.thin=nt, n.chains=nc,
              n.burnin=nb, n.iter=ni, parallel=TRUE)
end.time <- Sys.time()
run.time <- end.time - st.time
run.time
rm(st.time,end.time)

# Check the basics
max(out$summary[,"Rhat"])
#sort(out$summary[,"Rhat"], decreasing = T)[1:100]

min(out$summary[,"n.eff"])
#sort(out$summary[,"n.eff"])[1:50]

# Save output
library(R.utils)
saveObject(out, save.out)
