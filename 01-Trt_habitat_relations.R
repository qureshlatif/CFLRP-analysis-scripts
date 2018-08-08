require(QSLpersonal)
require(ggplot2)
require(cowplot)

setwd("C:/Users/Quresh.Latif/files/projects/FS/CFLRP")
load("Data_compiled.RData")

dat.point <- Cov %>% tbl_df %>%
  mutate(Point = row.names(Cov)) %>%
  select(Point, gridIndex:Rdens)

dat.grid <- landscape_data %>%
  left_join(
    dat.point %>% group_by(gridIndex) %>%
      summarise(percTrt = (sum(Trt_stat) / n()) * 100,
                Trt_time = mean(Trt_time, na.rm = T)),
    by = "gridIndex")

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
    xlab(NULL) + ylab(ylab)
  return(p)
}

p.time.fn <- function(dat, ylab) {
  p <- ggplot(dat, aes(x = Trt_time, y = Y)) +
    geom_jitter(alpha = 0.1) +
    geom_smooth() +
    scale_x_continuous(breaks = 1:10) +
    xlab("Years since trt") + ylab(ylab)
  return(p)
}

p.stitch.fn <- function(p.stat, p.time) {
  p <- ggdraw() +
    draw_plot(p.stat, x = 0, y = 0, width = .5, height = 1) +
    draw_plot(p.time, x = .5, y = 0, width = .5, height = 1)
  return(p)
}
#__________________________________#

#### Point-level relations ####
## Tabulate correlation coefficients ##
vars <- names(dat.point)[c(6:12, 15:19, 21, 22)]
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

  v2 <- dat.point$Trt_time
  r <- cor(v1, v2, use = "complete") %>%
    round(digits = 3)
  n <- sum(!is.na(v1) & !is.na(v2))
  p <- cor.test(v1, v2, alternative = "two.sided", method = "pearson")$p.value
  ifelse(p < 0.05, p.sig <- "*", p.sig <- "")
  out[vars[i], "Trt_time"] <- str_c(r, " (", n, ")", p.sig)
}

out %>% write.csv("Trt_hab_correlations_point.csv", row.names = T)

## Plot vegetation (mehanistic factors) VS outbreak metrics ##

# Canopy cover #
p.stat <- p.stat.fn(dat.point %>% rename(Y = CanCov), ylab = "Canopy cover")
p.time <- p.time.fn(dat.point %>% rename(Y = CanCov) %>% filter(Trt_stat == 1), ylab = "Canopy cover")
p.CanCov <- p.stitch.fn(p.stat, p.time)

# Canopy height #
p.stat <- p.stat.fn(dat.point %>% rename(Y = CanHt), ylab = "Canopy height")
p.time <- p.time.fn(dat.point %>% rename(Y = CanHt) %>% filter(Trt_stat == 1), ylab = "Canopy height")
p.CanHt <- p.stitch.fn(p.stat, p.time)

# Number of snags #
p.stat <- p.stat.fn(dat.point %>% rename(Y = NumSnags), ylab = "Number of snags")
p.time <- p.time.fn(dat.point %>% rename(Y = NumSnags) %>% filter(Trt_stat == 1), ylab = "Number of snags")
p.NumSnags <- p.stitch.fn(p.stat, p.time)

# Ponderosa pine dominance #
p.stat <- p.stat.fn(dat.point %>% rename(Y = RCOV_PP), ylab = "Ponderosa pine RC")
p.time <- p.time.fn(dat.point %>% rename(Y = RCOV_PP), ylab = "Ponderosa pine RC")
p.RCOV_PP <- p.stitch.fn(p.stat, p.time)

# Douglas fir dominance #
p.stat <- p.stat.fn(dat.point %>% rename(Y = RCOV_DF), ylab = "Douglas fir RC")
p.time <- p.time.fn(dat.point %>% rename(Y = RCOV_DF), ylab = "Douglas fir RC")
p.RCOV_DF <- p.stitch.fn(p.stat, p.time)

# Aspen dominance #
p.stat <- p.stat.fn(dat.point %>% rename(Y = RCOV_AS), ylab = "Aspen RC")
p.time <- p.time.fn(dat.point %>% rename(Y = RCOV_AS), ylab = "Aspen RC")
p.RCOV_AS <- p.stitch.fn(p.stat, p.time)

# Shrub volume #
p.stat <- p.stat.fn(dat.point %>% rename(Y = ShrubVol), ylab = "Shrub volume")
p.time <- p.time.fn(dat.point %>% rename(Y = ShrubVol), ylab = "Shrub volume")
p.ShrubVol <- p.stitch.fn(p.stat, p.time)

# Shrub diversity #
p.stat <- p.stat.fn(dat.point %>% rename(Y = ShrubDiv), ylab = "Shrub diversity")
p.time <- p.time.fn(dat.point %>% rename(Y = ShrubDiv), ylab = "Shrub diversity")
p.ShrubDiv <- p.stitch.fn(p.stat, p.time)

# Ladder fuel shrubs #
p.stat <- p.stat.fn(dat.point %>% rename(Y = RSCV_Ladder), ylab = "Ladder fuel shrub RC")
p.time <- p.time.fn(dat.point %>% rename(Y = RSCV_Ladder), ylab = "Ladder fuel shrub RC")
p.RSCV_Ladder <- p.stitch.fn(p.stat, p.time)

# Berry shrubs #
p.stat <- p.stat.fn(dat.point %>% rename(Y = RSCV_Ber), ylab = "Berry shrub RC")
p.time <- p.time.fn(dat.point %>% rename(Y = RSCV_Ber), ylab = "Berry shrub RC")
p.RSCV_Ber <- p.stitch.fn(p.stat, p.time)

# Aspen shrubs #
p.stat <- p.stat.fn(dat.point %>% rename(Y = RSCV_AS), ylab = "Aspen shrub RC")
p.time <- p.time.fn(dat.point %>% rename(Y = RSCV_AS), ylab = "Aspen shrub RC")
p.RSCV_AS <- p.stitch.fn(p.stat, p.time)

# Herb-grass volume #
p.stat <- p.stat.fn(dat.point %>% rename(Y = HerbGrassVol), ylab = "Herbaceous volume")
p.time <- p.time.fn(dat.point %>% rename(Y = HerbGrassVol), ylab = "Herbaceous volume")
p.HerbGrassVol <- p.stitch.fn(p.stat, p.time)

# Shrub-overstory height ratio #
p.stat <- p.stat.fn(dat.point %>% rename(Y = SOHtRatio), ylab = "Shrub-overstory Ratio")
p.time <- p.time.fn(dat.point %>% rename(Y = SOHtRatio), ylab = "Shrub-overstory Ratio")
p.SOHtRatio <- p.stitch.fn(p.stat, p.time)

p.point <- ggdraw() + 
  draw_plot(p.CanCov, x = 0, y = 0.8571429, width = 0.5, height = 0.1428571) +
  draw_plot(p.CanHt, x = .5, y = 0.8571429, width = 0.5, height = 0.1428571) +
  draw_plot(p.NumSnags, x = 0, y = 0.7142857, width = 0.5, height = 0.1428571) +
  draw_plot(p.RCOV_PP, x = .5, y = 0.7142857, width = 0.5, height = 0.1428571) +
  draw_plot(p.RCOV_DF, x = 0, y = 0.5714286, width = 0.5, height = 0.1428571) +
  draw_plot(p.RCOV_AS, x = 0.5, y = 0.5714286, width = 0.5, height = 0.1428571) +
  draw_plot(p.ShrubVol, x = 0, y = 0.4285714, width = 0.5, height = 0.1428571) +
  draw_plot(p.ShrubDiv, x = 0.5, y = 0.4285714, width = 0.5, height = 0.1428571) +
  draw_plot(p.RSCV_Ladder, x = 0, y = 0.2857143, width = 0.5, height = 0.1428571) +
  draw_plot(p.RSCV_Ber, x = 0.5, y = 0.2857143, width = 0.5, height = 0.1428571) +
  draw_plot(p.RSCV_AS, x = 0, y = 0.1428571, width = 0.5, height = 0.1428571) +
  draw_plot(p.HerbGrassVol, x = 0, y = 0, width = 0.5, height = 0.1428571) +
  draw_plot(p.SOHtRatio, x = 0.5, y = 0, width = 0.5, height = 0.1428571)


save_plot("figure_trt_vs_hab_point.tiff", p.point, ncol = 2, nrow = 3.5, dpi = 300)

#### Grid-level relations ####
## Tabulate correlation coefficients ##
vars <- names(dat.grid)[5:18]
cols <- c("percTrt", "Trt_time")
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
  
  v2 <- dat.grid$Trt_time
  r <- cor(v1, v2, use = "complete") %>%
    round(digits = 3)
  n <- sum(!is.na(v1) & !is.na(v2))
  p <- cor.test(v1, v2, alternative = "two.sided", method = "pearson")$p.value
  ifelse(p < 0.05, p.sig <- "*", p.sig <- "")
  out[vars[i], "Trt_time"] <- str_c(r, " (", n, ")", p.sig)
}

out %>% write.csv("Trt_hab_correlations_grid.csv", row.names = T)

## Plot vegetation (mehanistic factors) VS outbreak metrics...(not done yet) ##


save_plot("figure_trt_vs_hab_grid.tiff", p.point, ncol = 2, nrow = 3.5, dpi = 300)

