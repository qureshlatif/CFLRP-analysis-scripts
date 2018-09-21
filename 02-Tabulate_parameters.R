library(jagsUI)
library(stringr)
library(dplyr)
library(R.utils)

setwd("C:/Users/Quresh.Latif/files/projects/FS/CFLRP")
load("Data_compiled.RData")

#__________ Script inputs _____________#
mod <- loadObject("mod_treat2_d0yr")
params <- c("d0.1", "d0.2", "d0.3",
            "bd.ptrt",
            "bd.ptrt2",
            "bd.YST",
            "bd.TWIP",
            "bd.Rdens",
            "bb.trt",
            "bb.YST",
            "ba.Time",
            "ba.Time2",
            "ba.DOY",
            "ba.trt",
            "ba.YST")
out.vals <- c("est", "f")
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
