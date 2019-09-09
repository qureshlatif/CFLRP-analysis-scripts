library(jagsUI)
library(stringr)
library(dplyr)
library(R.utils)

setwd("C:/Users/Quresh.Latif/files/projects/FS/CFLRP")
load("Data_compiled.RData")

mod.type <- c("treatment", "habitat")
mods <- c("mod_treatment_d0yr", "mod_habitat_d0yr_reduced")
pars <- c("omega", "rho.ab", "rho.bd",
          "alpha0", "sigma.a0", "beta0", "sigma.b0", "delta0", "sigma.d0",
          
          "Betad.PctTrt", "sigma.Betad.PctTrt",
          "Betad.PACC10_3km", "sigma.Betad.PACC10_3km",
          "Betad.PACC40_3km", "sigma.Betad.PACC40_3km",
          "Betad.mnPerArRatio_Opn3km", "sigma.Betad.mnPerArRatio_Opn3km",
          "Betad.TWI", "sigma.Betad.TWI", "Betad.heatload", "sigma.Betad.heatload",
          "Betad.ForAR", "sigma.Betad.ForAR",
          
          "Betab.Trt", "sigma.Betab.Trt",
          "Betab.YST", "sigma.Betab.YST",
          "Betab.CanCov", "sigma.Betab.CanCov",
          "Betab.CanHt", "sigma.Betab.CanHt",
          "Betab.NumSnags", "sigma.Betab.NumSnags",
          "Betab.RCOV_PP", "sigma.Betab.RCOV_PP",
          "Betab.RCOV_DF", "sigma.Betab.RCOV_DF",
          "Betab.RCOV_AS", "sigma.Betab.RCOV_AS",
          "Betab.ShrubVol", "sigma.Betab.ShrubVol",
          "Betab.RSCV_Ladder", "sigma.Betab.RSCV_Ladder",
          "Betab.HerbGrassVol", "sigma.Betab.HerbGrassVol",
          
          "Betaa.Time", "sigma.Betaa.Time", "Betaa.Time2", "sigma.Betaa.Time2",
          "Betaa.DOY", "sigma.Betaa.DOY",
          "Betaa.Trt", "sigma.Betaa.Trt", "Betaa.YST", "sigma.Betaa.YST",
          "Betaa.CCov", "sigma.Betaa.CCov", "Betaa.SHVol", "sigma.Betaa.SHVol")

out <- matrix(NA, nrow = length(pars), ncol = length(mods),
              dimnames = list(pars, mods))

for(m in 1:length(mods)) {
  mod <- loadObject(mods[m])
  for(p in which(pars %in% names(mod$sims.list)))
    out[p, m] <- str_c(
      median(mod$sims.list[[pars[p]]]) %>% round(digits = 2),
      " (",
      quantile(mod$sims.list[[pars[p]]],
               prob = 0.025, type = 8) %>% round(digits = 2),
      ",",
      quantile(mod$sims.list[[pars[p]]],
               prob = 0.975, type = 8) %>% round(digits = 2),
      ")")
}

out %>% write.csv("Hyper_parameters.csv", row.names = T)
