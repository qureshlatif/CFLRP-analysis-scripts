require(QSLpersonal)
require(ggplot2)
require(cowplot)
theme_set(theme_cowplot())

setwd("C:/Users/Quresh.Latif/files/projects/FS/CFLRP")
load("Data_compiled.RData")

dat.point <- Cov %>% tbl_df %>%
  mutate(Point = row.names(Cov)) %>%
  select(Point, gridIndex:ForAR)

dat.grid <- landscape_data %>%
  select(Grid:NNdist_Opn3km, heatload:PctTrt_1kmNB) %>%
  rename(percTrt = PctTrt_1kmNB)

#___________ Plotting functions ____________#
p.stat.fn <- function(dat, ylab) {
  dt.sum <- dat %>% group_by(Trt_stat) %>%
    summarise(mn = mean(Y, na.rm = T),
              se = sd(Y, na.rm = T) / (sum(!is.na(Y)) ^ 0.5)) %>%
    mutate(lo95 = mn - se*1.96,
           hi95 = mn + se*1.96)
  p <- ggplot(dat, aes(x = Trt_stat, y = Y)) +
    geom_jitter(alpha = 0.1) +
    geom_errorbar(data = dt.sum, aes(x = Trt_stat, ymin = lo95, ymax = hi95),
                  inherit.aes = F, width = 0.3, size = 1, color = "blue") +
    geom_point(data = dt.sum, aes(x = Trt_stat, y = mn), size = 3, fill = "blue", color = "blue") +
    scale_x_continuous(breaks = c(0, 1), labels = c("Untreated", "Treated")) +
    xlab("Treatment status (Trt)") + ylab(ylab)
  return(p)
}

p.time.fn <- function(dat, ylab) {
  p <- ggplot(dat, aes(x = Trt_time, y = Y)) +
    geom_jitter(alpha = 0.1) +
    geom_smooth() +
    scale_x_continuous(breaks = 1:10) +
    xlab("Years since treatment (YST)") + ylab(ylab)
  return(p)
}

p.stitch.fn <- function(p.stat, p.time) {
  p <- ggdraw() +
    draw_plot(p.stat, x = 0, y = 0, width = .5, height = 1) +
    draw_plot(p.time, x = .5, y = 0, width = .5, height = 1)
  return(p)
}

p.ptrt.fn <- function(dat, ylab) {
  p <- ggplot(dat, aes(x = percTrt, y = Y)) +
    geom_point(alpha = 0.1) +
    geom_smooth(color = "blue") +
    xlab(NULL) + ylab(ylab)
  return(p)
}
#__________________________________#

#### Point-level relations ####
## Tabulate correlation coefficients ##
vars <- names(dat.point)[c(6:8, 10:12, 15:19, 21, 22)]
cols <- c("Trt_status", "Trt_time")
out <- matrix("", nrow = length(vars), ncol = length(cols),
              dimnames = list(vars, cols))

for(i in 1:length(vars)) {
  v1 <- dat.point %>% pull(vars[i])
  v2 <- dat.point$Trt_stat
  r <- cor(v1, v2, use = "complete") %>%
    round(digits = 3)
  n <- sum(!is.na(v1) & !is.na(v2))
  p <- cor.test(v1, v2, alternative = "two.sided", method = "pearson")$p.value
  ifelse(p < 0.05, p.sig <- "*", p.sig <- "")
  out[vars[i], "Trt_status"] <- str_c(r, " (", n, ")", p.sig)

  v1 <- dat.point %>% filter(Trt_stat == 1) %>% pull(vars[i])
  v2 <- dat.point %>% filter(Trt_stat == 1) %>% pull(Trt_time)
  r <- cor(v1, v2, use = "complete") %>%
    round(digits = 3)
  n <- sum(!is.na(v1) & !is.na(v2))
  p <- cor.test(v1, v2, alternative = "two.sided", method = "pearson")$p.value
  ifelse(p < 0.05, p.sig <- "*", p.sig <- "")
  out[vars[i], "Trt_time"] <- str_c(r, " (", n, ")", p.sig)
}

out %>% write.csv("Trt_hab_correlations_point.csv", row.names = T)

## Multiple regression models ##
vars <- names(dat.point)[c(6:8, 10:12, 15, 17:19, 21)]

mod <- glm(str_c("Trt_stat ~ ", str_c(vars, collapse = "+")), data = dat.point, family = "binomial")
summary(mod)$coefficients %>% write.csv("TrtStat_hab_relations_point.csv", row.names = T)

mod <- glm(str_c("Trt_time ~ ", str_c(vars, collapse = "+")), data = dat.point %>% filter(Trt_stat == 1), family = "poisson")
summary(mod)$coefficients %>% write.csv("TrtTime_hab_relations_point.csv", row.names = T)

## Correlations among habitat variables ##
cor(dat.point[, c(6:8, 10:12, 15:19, 21, 22)], use = "complete") %>%
  round(digits = 3) %>%
  write.csv("Correlations_all_points.csv", row.names = T)

cor((dat.point %>% filter(Trt_stat == 1))[, c(6:8, 10:12, 15:19, 21, 22)], use = "complete") %>%
  round(digits = 3) %>%
  write.csv("Correlations_treated_points.csv", row.names = T)

## Plot vegetation (mehanistic factors) VS outbreak metrics ##

# Canopy cover #
p.stat <- p.stat.fn(dat.point %>% rename(Y = CanCov), ylab = "CanCov")
p.time <- p.time.fn(dat.point %>% rename(Y = CanCov) %>% filter(Trt_stat == 1), ylab = "CanCov")
p.CanCov <- p.stitch.fn(p.stat, p.time)

# Canopy height #
p.stat <- p.stat.fn(dat.point %>% rename(Y = CanHt), ylab = "CanHt")
p.time <- p.time.fn(dat.point %>% rename(Y = CanHt) %>% filter(Trt_stat == 1), ylab = "CanHt")
p.CanHt <- p.stitch.fn(p.stat, p.time)

# Number of snags #
p.stat <- p.stat.fn(dat.point %>% rename(Y = NumSnags), ylab = "NSnag")
p.time <- p.time.fn(dat.point %>% rename(Y = NumSnags) %>% filter(Trt_stat == 1), ylab = "NSnag")
p.NumSnags <- p.stitch.fn(p.stat, p.time)

# Ponderosa pine dominance #
p.stat <- p.stat.fn(dat.point %>% rename(Y = RCOV_PP), ylab = "PIPO")
p.time <- p.time.fn(dat.point %>% rename(Y = RCOV_PP), ylab = "PIPO")
p.RCOV_PP <- p.stitch.fn(p.stat, p.time)

# Douglas fir dominance #
p.stat <- p.stat.fn(dat.point %>% rename(Y = RCOV_DF), ylab = "PSME")
p.time <- p.time.fn(dat.point %>% rename(Y = RCOV_DF), ylab = "PSME")
p.RCOV_DF <- p.stitch.fn(p.stat, p.time)

# Aspen dominance #
p.stat <- p.stat.fn(dat.point %>% rename(Y = RCOV_AS), ylab = "POTR5")
p.time <- p.time.fn(dat.point %>% rename(Y = RCOV_AS), ylab = "POTR5")
p.RCOV_AS <- p.stitch.fn(p.stat, p.time)

# Shrub volume #
p.stat <- p.stat.fn(dat.point %>% rename(Y = ShrubVol), ylab = "ShrbVol")
p.time <- p.time.fn(dat.point %>% rename(Y = ShrubVol), ylab = "ShrbVol")
p.ShrubVol <- p.stitch.fn(p.stat, p.time)

# # Shrub diversity #
# p.stat <- p.stat.fn(dat.point %>% rename(Y = ShrubDiv), ylab = "Shrub diversity")
# p.time <- p.time.fn(dat.point %>% rename(Y = ShrubDiv), ylab = "Shrub diversity")
# p.ShrubDiv <- p.stitch.fn(p.stat, p.time)

# Ladder fuel shrubs #
p.stat <- p.stat.fn(dat.point %>% rename(Y = RSCV_Ladder), ylab = "LadFuel")
p.time <- p.time.fn(dat.point %>% rename(Y = RSCV_Ladder), ylab = "LadFuel")
p.RSCV_Ladder <- p.stitch.fn(p.stat, p.time)

# # Berry shrubs #
# p.stat <- p.stat.fn(dat.point %>% rename(Y = RSCV_Ber), ylab = "Berry shrub RC")
# p.time <- p.time.fn(dat.point %>% rename(Y = RSCV_Ber), ylab = "Berry shrub RC")
# p.RSCV_Ber <- p.stitch.fn(p.stat, p.time)
# 
# # Aspen shrubs #
# p.stat <- p.stat.fn(dat.point %>% rename(Y = RSCV_AS), ylab = "Aspen shrub RC")
# p.time <- p.time.fn(dat.point %>% rename(Y = RSCV_AS), ylab = "Aspen shrub RC")
# p.RSCV_AS <- p.stitch.fn(p.stat, p.time)

# Herb-grass volume #
p.stat <- p.stat.fn(dat.point %>% rename(Y = HerbGrassVol), ylab = "Herb")
p.time <- p.time.fn(dat.point %>% rename(Y = HerbGrassVol), ylab = "Herb")
p.HerbGrassVol <- p.stitch.fn(p.stat, p.time)

# # Shrub-overstory height ratio #
# p.stat <- p.stat.fn(dat.point %>% rename(Y = SOHtRatio), ylab = "Shrub-overstory Ratio")
# p.time <- p.time.fn(dat.point %>% rename(Y = SOHtRatio), ylab = "Shrub-overstory Ratio")
# p.SOHtRatio <- p.stitch.fn(p.stat, p.time)

p.point <- ggdraw() + 
  draw_plot(p.CanCov, x = 0, y = 0.6666667, width = 0.3333333, height = 0.3333333) +
  draw_plot(p.CanHt, x = 0.3333333, y = 0.6666667, width = 0.3333333, height = 0.3333333) +
  draw_plot(p.NumSnags, x = 0.6666667, y = 0.6666667, width = 0.3333333, height = 0.3333333) +
  draw_plot(p.RCOV_PP, x = 0, y = 0.3333333, width = 0.3333333, height = 0.3333333) +
  draw_plot(p.RCOV_DF, x = 0.3333333, y = 0.3333333, width = 0.3333333, height = 0.3333333) +
  draw_plot(p.RCOV_AS, x = 0.6666667, y = 0.3333333, width = 0.3333333, height = 0.3333333) +
  draw_plot(p.ShrubVol, x = 0, y = 0, width = 0.3333333, height = 0.3333333) +
  #draw_plot(p.ShrubDiv, x = 0.5, y = 0.4285714, width = 0.5, height = 0.1428571) +
  draw_plot(p.RSCV_Ladder, x = 0.3333333, y = 0, width = 0.3333333, height = 0.3333333) +
  #draw_plot(p.RSCV_Ber, x = 0.5, y = 0.2857143, width = 0.5, height = 0.1428571) +
  #draw_plot(p.RSCV_AS, x = 0, y = 0.1428571, width = 0.5, height = 0.1428571) +
  draw_plot(p.HerbGrassVol, x = 0.6666667, y = 0, width = 0.3333333, height = 0.3333333)
  #draw_plot(p.SOHtRatio, x = 0.5, y = 0, width = 0.5, height = 0.1428571)


save_plot("figure_trt_vs_hab_point.tiff", p.point, ncol = 3, nrow = 3, dpi = 600)

#### Grid-level relations ####
## Tabulate correlation coefficients ##
vars <- names(dat.grid)[c(5, 6, 12, 14, 17:18)]

mod <- glm(str_c("percTrt ~ ", str_c(vars, collapse = "+")), data = dat.grid, family = quasibinomial())
summary(mod)$coefficients %>% write.csv("TrtPerc_hab_relations_grid.csv", row.names = T)

cols <- c("percTrt")
out <- matrix("", nrow = length(vars), ncol = length(cols),
              dimnames = list(vars, cols))

for(i in 1:length(vars)) {
  v1 <- dat.grid %>% pull(vars[i])
  v2 <- dat.grid$percTrt
  r <- cor(v1, v2, use = "complete") %>%
    round(digits = 3)
  n <- sum(!is.na(v1) & !is.na(v2))
  p <- cor.test(v1, v2, alternative = "two.sided", method = "pearson")$p.value
  ifelse(p < 0.05, p.sig <- "*", p.sig <- "")
  out[vars[i], "percTrt"] <- str_c(r, " (", n, ")", p.sig)
}

out %>% write.csv("Trt_hab_correlations_grid.csv", row.names = T)

## Plot vegetation (mehanistic factors) VS outbreak metrics ##

# Percent area #
p.PACC10 <- p.ptrt.fn(dat.grid %>% rename(Y = PACC10_3km), ylab = "PACCGap")
p.PACC40 <- p.ptrt.fn(dat.grid %>% rename(Y = PACC40_3km), ylab = "PACCOpn")
#p.mnPtchAr_Gap <- p.ptrt.fn(dat.grid %>% rename(Y = mnPtchAr_Gap3km), ylab = "Mean gap area")
#p.mnPtchAr_Opn <- p.ptrt.fn(dat.grid %>% rename(Y = mnPtchAr_Opn3km), ylab = "Mean open forest")
#p.mnPerArRatio_Gap <- p.ptrt.fn(dat.grid %>% rename(Y = mnPerArRatio_Gap3km), ylab = "Mean gap PA ratio")
p.mnPerArRatio_Opn <- p.ptrt.fn(dat.grid %>% rename(Y = mnPerArRatio_Opn3km), ylab = "PAROpn")
#p.NNdist_Gap <- p.ptrt.fn(dat.grid %>% rename(Y = NNdist_Gap3km), ylab = "Mean gap NN dist")
#p.NNdist_Opn <- p.ptrt.fn(dat.grid %>% rename(Y = NNdist_Opn3km), ylab = "Mean open NN dist")

p.grid <- ggdraw() + 
  draw_plot(p.PACC10, x = 0, y = 0.6667, width = 1, height = 0.3333) +
  draw_plot(p.PACC40, x = 0, y = 0.3333, width = 1, height = 0.3333) +
  #draw_plot(p.mnPtchAr_Gap, x = 0, y = 0.5, width = 0.5, height = 0.25) +
  #draw_plot(p.mnPtchAr_Opn, x = .5, y = 0.5, width = 0.5, height = 0.25) +
  #draw_plot(p.mnPerArRatio_Gap, x = 0, y = 0.25, width = 0.5, height = 0.25) +
  draw_plot(p.mnPerArRatio_Opn, x = 0, y = 0, width = 1, height = 0.3333)
  #draw_plot(p.NNdist_Gap, x = 0, y = 0, width = 0.5, height = 0.25) +
  #draw_plot(p.NNdist_Opn, x = 0.5, y = 0, width = 0.5, height = 0.25)
p.grid <- ggdraw() +
  draw_plot(p.grid, x = 0, y = 0.05, width = 1, height = 0.95) +
  draw_plot_label("Percent treated (percTrt)", x = 0.33, y = 0.05, hjust = 0)

save_plot("figure_trt_vs_hab_grid.tiff", p.grid, ncol = 1, nrow = 3, dpi = 600)

## Correlations among habitat variables ##
cor(dat.grid[vars], use = "complete") %>%
  round(digits = 3) %>%
  write.csv("Correlations_grid.csv", row.names = T)
