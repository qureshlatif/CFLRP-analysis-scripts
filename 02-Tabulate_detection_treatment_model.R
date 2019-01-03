library(jagsUI)
library(stringr)
library(dplyr)
library(R.utils)
library(ggplot2)
library(cowplot)
library(QSLpersonal)

setwd("C:/Users/Quresh.Latif/files/projects/FS/CFLRP")
load("Data_compiled.RData")

mod <- loadObject("mod_treatment_d0yr")

# Tabulate parameter estimates
pars <- c("p_star", "p_star_md", "ba.Time", "ba.Time2", "ba.DOY", "ba.trt", "ba.YST")
tbl_pars <- matrix(NA, nrow = length(spp.list), ncol = length(pars),
                   dimnames = list(NULL, pars))

for(i in 3:length(pars)) {
  parm <- mod$sims.list[[pars[i]]]
  med <- apply(parm, 2, median)
  lo <- apply(parm, 2, function(x) quantile(x, prob = 0.025, type = 8))
  hi <- apply(parm, 2, function(x) quantile(x, prob = 0.975, type = 8))
  if(any(lo > 0 | hi < 0)) {
    ind.support <- which(lo > 0 | hi < 0)
    tbl_pars[ind.support, pars[i]] <- str_c(med[ind.support] %>% round(digits = 2), " (",
                                            lo[ind.support] %>% round(digits = 2), ",",
                                            hi[ind.support] %>% round(digits = 2), ")*")
    tbl_pars[-ind.support, pars[i]] <- str_c(med[-ind.support] %>% round(digits = 2), " (",
                                             lo[-ind.support] %>% round(digits = 2), ",",
                                             hi[-ind.support] %>% round(digits = 2), ")")
  } else {
    tbl_pars[, pars[i]] <- str_c(med %>% round(digits = 2), " (",
                                 lo %>% round(digits = 2), ",",
                                 hi %>% round(digits = 2), ")")
  }
}

p_star <- 1 - (1 - expit(mod$sims.list$a0))^6
tbl_pars[, "p_star_md"] <- apply(p_star, 2, median)
tbl_pars[, "p_star"] <- str_c(apply(p_star, 2, median) %>% round(digits = 2), " (",
                              apply(p_star, 2, function(x) quantile(x, prob = 0.025, type = 8)) %>% round(digits = 2), ",",
                              apply(p_star, 2, function(x) quantile(x, prob = 0.975, type = 8)) %>% round(digits = 2), ")")

rm(parm, p_star, i)

tbl_pars <- tbl_pars %>% tbl_df() %>%
  mutate(Spp = spp.list) %>%
  mutate(NDET = apply(Y.mat, 2, sum) == 0) %>%
  mutate(p_star_md = as.numeric(p_star_md)) %>%
  arrange(p_star_md) %>%
  select(Spp, NDET, p_star, ba.Time:ba.YST)

write.csv(tbl_pars, "Detection_estimates_treatment.csv", row.names = F)
