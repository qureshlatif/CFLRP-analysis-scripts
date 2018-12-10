library(jagsUI)
library(stringr)
library(dplyr)

#setwd("/home/RMBO.LOCAL/quresh.latif/CFLRP")
setwd("C:/Users/Quresh.Latif/files/projects/FS/CFLRP")
load("Data_compiled_abundance.RData")

#___ Script inputs ___#
model.file <- "CFLRP-analysis-scripts/model_abundance_treatment_global.jags"
spp <- "CHSP"
spp.ind <- which(Spp == spp)
#detection.model <- "HZ" # Set to "HZ" for hazard rate of "HN" for half-normal

# MCMC values.  Adjust as needed.
nc <- 3
nb <- 5000
ni <- 30000
nt <- 10

save.out <- str_c("abund_models/mod_", spp, "_abundance_treatment_global")
#_____________________#

area.band <- area.band.list[[spp.ind]]
area.prop <- area.prop.list[[spp.ind]]
breaks <- breaks.list[[spp.ind]]
area.circle <- area.circle.list[spp.ind]
cutoff <- cutoff.list[spp.ind]

data <- list("Y",
             "dclass", "mean.cl", "sd.cl", # needed for distance sampling
             "gridID", "yearID", "nGrid", "n.year", "nPoint",
             "nInd", "nG", "area.band", "area.prop", "breaks", # needed for distance sampling
             "PctTrt.d", "TWI.d", "heatload.d", "trt.b", "YST.b",
             "DOY.b", "Time.b")

parameters <- c("beta0.mean", "beta0.sd", #"N.mean", "p.mean", # Assemble the parameters vector for JAGS (What we want to track).
                "beta0", "N", "cl.size",
                "bl.trt",
                "bl.YST",
                "bd.PctTrt", "bd.TWI", "bd.heatload",
                ##___ Hazard rate parameters ___##
                "a0", "a.Time", "a.Time2",
                "a.DOY", "a.DOY2",
                "a.trt", "a.YST",
                "b",
                ##______________________________##
                "lambda", "a", # Needed for WAIC
                "test") # GOF

Y <- eval(as.name(str_c("Y.", spp,".dist")))
dclass <- eval(as.name(str_c("dclass.", spp)))
mean.cl <- max(1.001, mean(dclass[, "CL_Count"]))
sd.cl <- max(0.001, sd(dclass[, "CL_Count"]))
dimnames(dclass) <- NULL
nInd <- nrow(dclass)
nPoint <- length(Y)
inits <- function() # Setting these based on posterior distribution from an initial successful run.
  list(N = Y, beta0.mean = rnorm(1, 1, 0.1), beta0.sd = rnorm(1, 0.66, 0.07),
       a0 = rnorm(1, 3.5, 0.5), b = rnorm(1, 3.16, 0.5)) # for hazard rate model
#bt.0 = rnorm(1, 3.4, 0.03)) # for half-normal model

# Detection data #
gridID <- Cov[, "gridIndex"]
yearID <- Cov[, "YearInd"]
nGrid <- max(gridID)
n.year <- max(yearID)

# Grid-level covariates #
heatload.d <- Cov[, "heatload"]  %>%
  tapply(gridID, mean) %>%
  (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T))

TWI.d <- Cov[, "TWI"]  %>%
  tapply(gridID, mean) %>%
  (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T))

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

# Point-level covariates #
trt.b <- Cov[, "Trt_stat"]

YST.b <- Cov[, "Trt_time"] %>%
  (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)) %>%
  replace(., which(is.na(.)), 0)

DOY.b <- Cov[, "DayOfYear"] %>%
  (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T))

Time.b <- Cov[, "Time"] %>%
  (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T))

st.time <- Sys.time()
out <- jagsUI(data, inits, parameters.to.save = parameters, model.file, n.thin=nt, n.chains=nc,
              n.burnin=nb, n.iter=ni, parallel=TRUE)
end.time <- Sys.time()
run.time <- end.time - st.time
run.time
rm(st.time,end.time)

max(out$summary[which(!is.na(out$summary[ ,"Rhat"])) ,"Rhat"])
sort(out$summary[,"Rhat"], decreasing = T)[1:100]

min(out$summary[,"n.eff"])
sort(out$summary[,"n.eff"])[1:100]

sum(out$sims.list$test) / out$mcmc.info$n.samples

library(R.utils)
saveObject(out, save.out)
