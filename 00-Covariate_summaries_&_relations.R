require(QSLpersonal)
setwd("C:/Users/Quresh.Latif/files/projects/FS/CFLRP")
load("Data_compiled.RData")

dat.grid <- landscape_data %>%
  full_join(veg_data %>%
              mutate(Grid = str_sub(Point_year, 1, -9),
                     Year = str_sub(Point_year, -4, -1)) %>%
              group_by(Grid, Year) %>%
              summarise(PercTrt = mean(Trt_stat) * 100,
                        TWIP = mean(TWIP),
                        Rdens = mean(Rdens),
                        Trt_time = mean(Trt_time, na.rm = T)),
            by = c("Grid", "Year")) %>%
  select(-ends_with("2km"))

## Tabulate summary statistics ##
# Point level #
vars <- names(veg_data)[-which(names(veg_data) %in% c("Point_year", "Forest", "RCOV_OT", "RSCV_OT", "TWIP", "Rdens"))]
vars <- vars[c((length(vars)-1):length(vars), 1:(length(vars)-2))]
binary <- which(vars == "Trt_stat")
out <- SumStats_df(veg_data, vars, binary)

write.csv(out, "Covariate_sum_points.csv", row.names = T)

# Grid level #
vars <- names(dat.grid)[-which(names(dat.grid) %in% c("Grid", "Year", "gridIndex", "YearInd"))]
vars <- vars[c(9, 12, 1:8, 10, 11)]
out <- SumStats_df(dat.grid, vars)

write.csv(out, "Covariate_sum_grids.csv", row.names = T)
