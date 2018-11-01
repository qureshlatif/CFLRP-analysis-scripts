library(jagsUI)
library(stringr)
library(dplyr)
library(R.utils)

setwd("C:/Users/Quresh.Latif/files/projects/FS/CFLRP")
load("Data_compiled.RData")

#__________ Script inputs _____________#
mod <- loadObject("mod_treatment_d0yr")
params <- c("d0.1", "d0.2", "d0.3", "b0", "a0", "bd.ptrt", "bd.ptrt2", "bd.YST",
            "bd.PACC10_3km", "bd.mnPtchAr_Gap3km", "bd.mnPerArRatio_Gap3km",
            "bd.NNdist_Gap3km", "bd.PACC40_3km", "bd.mnPtchAr_Opn3km",
            "bd.mnPerArRatio_Opn3km", "bd.NNdist_Opn3km",
            "bd.TWIP", "bd.Rdens",
            "bb.trt", "bb.YST",
            "bb.CanCov", "bb.CanCov2", "bb.CanHt", "bb.NumSnags", "bb.RCOV_PP",
            "bb.RCOV_DF", "bb.RCOV_AS", "bb.RSCV_Ladder", "bb.RSCV_Ber",
            "bb.HerbGrassVol", "bb.SOHtRatio", "bb.RESQ",
            "ba.DOY", "ba.trt", "ba.YST","ba.Time",
            "ba.Time2","ba.ccov", "ba.shvol")
out.vals <- c("est", "f")
params <- c("d0.1", "d0.2", "d0.3", params[which(params %in% names(mod$sims.list))])
#______________________________________#

cols <- (expand.grid(out.vals, params, stringsAsFactors = F) %>%
  select(Var2, Var1) %>%
  mutate(Var3 = str_c(Var2, Var1, sep = ".")))$Var3
out <- matrix(NA, nrow = length(spp.list), ncol = length(cols),
              dimnames = list(spp.list, cols))

for(i in 4:length(params)) {
  parm <- mod$sims.list[[params[i]]]
  if(!is.null(parm)) {
    out[, (i*2 - 1)] <- str_c(
      apply(parm, 2, median) %>% round(digits = 2),
      "(",
      apply(parm, 2, function(x) quantile(x, prob = 0.025, type = 8)) %>% round(digits = 2),
      ",",
      apply(parm, 2, function(x) quantile(x, prob = 0.975, type = 8)) %>% round(digits = 2),
      ")")
    out[, (i*2)] <- apply(parm, 2, function(x) max(c(sum(x > 0), sum(x < 0))) / length(x)) %>%
      round(digits = 2)
  }
}

for(i in 1:3) {
  parm <- mod$sims.list[["d0"]][,,i]
  out[, (i*2 - 1)] <- str_c(
    apply(parm, 2, median) %>% round(digits = 2),
    "(",
    apply(parm, 2, function(x) quantile(x, prob = 0.025, type = 8)) %>% round(digits = 2),
    ",",
    apply(parm, 2, function(x) quantile(x, prob = 0.975, type = 8)) %>% round(digits = 2),
    ")")
  out[, (i*2)] <- apply(parm, 2, function(x) max(c(sum(x > 0), sum(x < 0))) / length(x)) %>%
    round(digits = 2)
}

write.csv(out, "Parameter_est.csv", row.names = T)
