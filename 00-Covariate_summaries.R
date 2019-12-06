require(QSLpersonal)
setwd("C:/Users/Quresh.Latif/files/projects/FS/CFLRP")
load("Data_compiled.RData")

dat.grid <- landscape_data %>%
  select(-ends_with("2km"))

## Tabulate summary statistics ##
# Point level #
vars <- names(veg_data)[-which(names(veg_data) %in% c("Point_year", "Forest", "RCOV_OT", "RSCV_OT", "TWIP", "TWI", "heatload", "Rdens", "ForAR"))]
vars <- vars[c((length(vars)-1):length(vars), 1:(length(vars)-2))]
binary <- which(vars == "Trt_stat")
out <- SumStats_df(veg_data, vars, binary)

write.csv(out, "Covariate_sum_points.csv", row.names = T)

# Grid level #
vars <- names(dat.grid)[-which(names(dat.grid) %in% c("Grid", "Year", "gridIndex", "YearInd"))]
vars <- vars[c(17, 1:8, 12, 13)]
out <- SumStats_df(dat.grid, vars)

write.csv(out, "Covariate_sum_grids.csv", row.names = T)
