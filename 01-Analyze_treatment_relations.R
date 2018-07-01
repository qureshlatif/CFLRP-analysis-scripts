library(jagsUI)
library(stringr)
library(dplyr)

#setwd("/home/RMBO.LOCAL/quresh.latif/CFLRP")
setwd("C:/Users/Quresh.Latif/files/projects/FS/CFLRP")
load("Data_compiled.RData")

#### Script inputs ####
model.file <- "model_treatment.jags"

# Data objects to send to JAGS
data <- list("Y", "TPeriod", "gridID", "yearID", "n.grid", "n.year", "n.point_year", "n.spp", "Trt.b",
             "PctTrt.d", "YST.b", "YST.d", "TWIP.b", "TWIP.d",
             "DOY.b", "Time.b",
             "ccov.b", "ccov.means", "ccov.sd", "ccov.b.missing",
             "shvol.b", "shvol.means", "shvol.sd", "shvol.b.missing")

# Stuff to save from JAGS
parameters <- c("omega", "rho.ab", "rho.bd",
                "alpha0", "sigma.a0", "beta0", "sigma.b0", "delta0", "sigma.d0",
                
                "Betad.pers", "sigma.Betad.pers",
                "Betad.PctTrt", "sigma.Betad.PctTrt", "Betad.YST", "sigma.Betad.YST", "Betad.TWIP", "sigma.Betad.TWIP",
                "Betab.Trt", "sigma.Betab.Trt", "Betab.YST", "sigma.Betab.YST", "Betab.TWIP", "sigma.Betab.TWIP",
                "Betaa.Time", "sigma.Betaa.Time", "Betaa.Time2", "sigma.Betaa.Time2",
                "Betaa.DOY", "sigma.Betaa.DOY", "Betaa.DOY2", "sigma.Betaa.DOY2",
                "Betaa.CCov", "sigma.Betaa.CCov",
                "Betaa.SHVol", "sigma.Betaa.SHVol",
                
                "d0", "b0", "a0",
                "bd.pers", "bd.PctTrt", "bd.YST", "bd.TWIP",
                "bb.Trt", "bb.YST", "bb.TWIP",
                "ba.Time", "ba.Time2", "ba.DOY", "ba.DOY2", "ba.ccov", "ba.shvol")

# Function for setting initial values in JAGS
inits <- function()
  list(z=z.init, u=u.init, w=w.init, tvar.sigma.a0 = rnorm(1), tvar.sigma.b0 = rnorm(1), tvar.sigma.d0 = rnorm(1),
       tvar.Betad.pers = rnorm(1),
       tvar.Betad.PctTrt = rnorm(1), tvar.Betad.YST = rnorm(1), tvar.Betad.TWIP = rnorm(1),
       tvar.Betab.Trt = rnorm(1), tvar.Betab.YST = rnorm(1), tvar.Betab.TWIP = rnorm(1),
       tvar.Betaa.Time = rnorm(1), tvar.Betaa.Time2 = rnorm(1),
       tvar.Betaa.DOY = rnorm(1), tvar.Betaa.DOY2 = rnorm(1),
       tvar.Betaa.CCov = rnorm(1), tvar.Betaa.SHVol = rnorm(1))

# MCMC values
nc <- 3 # number of chains
nb <- 1000 # burn in
ni <- 15000 # number of iterations
nt <- 10 # thinning

save.out <- "mod_treatment"
##########################

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
Trt.b <- Cov[, "Trt_stat"] # Point-level values
PctTrt.d <- tapply(Trt.b, gridID, mean, na.rm = T) # Grid-level values

YST.b <- Cov[, "Trt_time"] %>% (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)) # Point-level values
YST.d <- tapply(YST.b, gridID, mean, na.rm = T) # Grid-level values
YST.b[is.na(YST.b)] <- 0
YST.d[is.na(YST.d)] <- 0

TWIP.b <- Cov[, "TWIP"] %>% (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)) # Point-level values
TWIP.d <- tapply(TWIP.b, gridID, mean, na.rm = T) # Grid-level values

DOY.b <- Cov[, "DayOfYear"] %>% (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)) # Point-level values

Time.b <- Cov[, "Time"] %>% (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)) # Point-level values

ccov.b <- Cov[, "CanCov"] %>% (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)) # Point-level values
ccov.means <- tapply(ccov.b, gridID, mean, na.rm = T) # Grid-level means for imputing missing values
ccov.sd <- tapply(ccov.b, gridID, sd, na.rm = T) # Grid-level SDs for imputing missing values
ccov.sd[which(ccov.sd == 0)] <- sd(ccov.b, na.rm = T) # Zeros won't work for SD!
ccov.b.missing <- is.na(ccov.b) %>% as.integer # Index missing values to be imputed
ccov.b[is.na(ccov.b)] <- 0

shvol.b <- (Cov[, "shrub_cover"] / 100) * (pi * 50^2) * Cov[, "ShrubHt"] %>% # area covered X volume
  (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)) # Point-level values
shvol.means <- tapply(shvol.b, gridID, mean, na.rm = T) # Grid-level means for imputing missing values
shvol.sd <- tapply(shvol.b, gridID, sd, na.rm = T) # Grid-level SDs for imputing missing values
shvol.b.missing <- is.na(shvol.b) %>% as.integer # Index missing values to be imputed
shvol.b[is.na(shvol.b)] <- 0

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