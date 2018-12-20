library(jagsUI)
library(stringr)
library(dplyr)

#setwd("/home/RMBO.LOCAL/quresh.latif/CFLRP")
setwd("C:/Users/Quresh.Latif/files/projects/FS/CFLRP")
load("Data_compiled.RData")

#_____ Script inputs _____#
model.file <- "model_treatment_d0yr.jags"

# MCMC values
nc <- 3 # number of chains
nb <- 1000 # burn in
ni <- 35000 # number of iterations
nt <- 10 # thinning

save.out <- "mod_treatment_d0yr"
#_________________________#

# Data objects to send to JAGS
data <- list("Y", "TPeriod", "gridID", "yearID", "n.grid", "n.year", "n.point_year", "n.spp",
             "Trt.b","PctTrt.d", "YST.b", "TWIP.d","TWI.d", "heatload.d", "ForAR.d", "Rdens.d", #"YST.d", 
             "DOY.b", "Time.b")

# Stuff to save from JAGS
parameters <- c("omega", "rho.ab", "rho.bd",
                "alpha0", "sigma.a0", "beta0", "sigma.b0", "delta0", "sigma.d0",
                
                "Betad.yr", "sigma.Betad.yr",
                "Betad.PctTrt", "sigma.Betad.PctTrt",
                "Betad.PctTrt2", "sigma.Betad.PctTrt2",
                "Betad.YST", "sigma.Betad.YST",
                "Betad.TWIP", "sigma.Betad.TWIP","Betad.Rdens", "sigma.Betad.Rdens",
                "Betad.TWI", "sigma.Betad.TWI", "Betad.heatload", "sigma.Betad.heatload",
                "Betad.ForAR", "sigma.Betad.ForAR",
                "Betab.Trt", "sigma.Betab.Trt", "Betab.YST", "sigma.Betab.YST",
                "Betaa.Time", "sigma.Betaa.Time", "Betaa.Time2", "sigma.Betaa.Time2",
                "Betaa.DOY", "sigma.Betaa.DOY",
                "Betaa.DOY2", "sigma.Betaa.DOY2",
                "Betaa.Trt", "sigma.Betaa.Trt", "Betaa.YST", "sigma.Betaa.YST",
                
                "d0", "b0", "a0",
                "bd.yr", # For year effect
                "bd.pers", # For persistence effect
                "bd.ptrt", "bd.ptrt2", "bd.YST",
                "bd.TWIP", "bd.TWI", "bd.heatload", "bd.ForAR", "bd.Rdens",
                "bb.trt", "bb.YST", "ba.Time", "ba.Time2",
                "ba.DOY", "ba.trt", "ba.YST", "ba.DOY2",
                
                "SR.grid", "SR.point")

# Function for setting initial values in JAGS
inits <- function()
  list(z=z.init, u=u.init, w=w.init, tvar.sigma.a0 = rnorm(1), tvar.sigma.b0 = rnorm(1), tvar.sigma.d0 = rnorm(n.year),
       tvar.Betad.yr = rnorm(n.year - 1),
       tvar.Betad.PctTrt = rnorm(1),
       tvar.Betad.PctTrt2 = rnorm(1),
       tvar.Betad.YST = rnorm(1),
       tvar.Betad.TWIP = rnorm(1), tvar.Betad.Rdens = rnorm(1),
       tvar.Betad.TWI = rnorm(1), tvar.Betad.heatload = rnorm(1), tvar.Betad.ForAR = rnorm(1),
       tvar.Betab.Trt = rnorm(1), tvar.Betab.YST = rnorm(1),
       tvar.Betaa.Time = rnorm(1), tvar.Betaa.Time2 = rnorm(1),
       tvar.Betaa.DOY = rnorm(1), tvar.Betaa.DOY2 = rnorm(1),
       tvar.Betaa.Trt = rnorm(1), tvar.Betaa.YST = rnorm(1))

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
PctTrt.d <- matrix(NA, nrow = max(gridID), ncol = max(yearID))
PctTrt.d[landscape_data %>% filter(YearInd == 1) %>% pull(gridIndex), 1] <-
  landscape_data %>% filter(YearInd == 1) %>% pull(PctTrt_1kmNB)
PctTrt.d[landscape_data %>% filter(YearInd == 2) %>% pull(gridIndex), 2] <-
  landscape_data %>% filter(YearInd == 2) %>% pull(PctTrt_1kmNB)
PctTrt.d[landscape_data %>% filter(YearInd == 3) %>% pull(gridIndex), 3] <-
  landscape_data %>% filter(YearInd == 3) %>% pull(PctTrt_1kmNB)
PctTrt.d <- PctTrt.d %>%
  (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T))
PctTrt.d[which(is.na(PctTrt.d))] <- 0

YST.b <- Cov[, "Trt_time"] %>% (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)) # Point-level values
YST.b[is.na(YST.b)] <- 0
# YST.d <- matrix(NA, nrow = max(gridID), ncol = max(yearID))
# YST.d[landscape_data %>% filter(YearInd == 1) %>% pull(gridIndex), 1] <-
#   landscape_data %>% filter(YearInd == 1) %>% pull(Trt_time)
# YST.d[landscape_data %>% filter(YearInd == 2) %>% pull(gridIndex), 2] <-
#   landscape_data %>% filter(YearInd == 2) %>% pull(Trt_time)
# YST.d[landscape_data %>% filter(YearInd == 3) %>% pull(gridIndex), 3] <-
#   landscape_data %>% filter(YearInd == 3) %>% pull(Trt_time)
# YST.d <- YST.d %>%
#   (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T))
# YST.d[is.na(YST.d)] <- 0

TWIP.d <- tapply(Cov[, "TWIP"], gridID, mean, na.rm = T) %>% # Grid-level values
  (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T))

heatload.d <- tapply(Cov[, "heatload"], gridID, mean, na.rm = T) %>% # Grid-level values
  (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T))

TWI.d <- tapply(Cov[, "TWI"], gridID, mean, na.rm = T) %>% # Grid-level values
  (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T))

ForAR.d <- tapply(Cov[, "ForAR"], gridID, mean, na.rm = T) %>% # Grid-level values
  (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T))

Rdens.d <- tapply(Cov[, "Rdens"], gridID, mean, na.rm = T) %>% # Grid-level values
  (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T))

DOY.b <- Cov[, "DayOfYear"] %>% (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)) # Point-level values

Time.b <- Cov[, "Time"] %>% (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)) # Point-level values

# Save workspace for GOF #
#save.image("GOF_workspace.RData")

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
