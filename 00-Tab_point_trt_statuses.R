library(foreign)
library(dplyr)
library(stringr)
setwd("C:/Users/Quresh.Latif/files/projects/FS/CFLRP")
load("Data_compiled.RData")

# dat <- read.dbf("Bird_survey_point_coords.dbf", as.is = T) %>%
#   mutate(Point = str_c(TransectNu, "-", Point %>% str_pad("0", width = 2, side = "left"))) %>%
#   mutate(SuperStrat = str_sub(TransectNu, 1, 8))

# Treatment status by year surveyed and strata #
sum.point <- veg_data %>%
  mutate(Year = Point_year %>% str_sub(-4, -1)) %>%
  group_by(Year) %>%
  summarize(Trt_CFLRP = sum(Trt_stat == 1 & (str_sub(Point_year, 4, 8) == "CFLRP")),
            Cnt_CFLRP = sum(Trt_stat == 0 & (str_sub(Point_year, 4, 8) == "CFLRP")),
            Trt_IMBCR = sum(Trt_stat == 1 & (str_sub(Point_year, 4, 8) == "BCR16")),
            Cnt_IMBCR = sum(Trt_stat == 0 & (str_sub(Point_year, 4, 8) == "BCR16"))) %>%
  bind_rows(data.frame(Year = "ALL", stringsAsFactors = F) %>%
              mutate(Trt_CFLRP = veg_data %>%
                       filter(str_sub(Point_year, 4, 8) == "CFLRP" &
                                Trt_stat ==1) %>%
                       pull(Point_year) %>% str_sub(1, -6) %>% unique %>% length,
                     Cnt_CFLRP = veg_data %>%
                       filter(str_sub(Point_year, 4, 8) == "CFLRP" &
                                Trt_stat == 0) %>%
                       pull(Point_year) %>% str_sub(1, -6) %>% unique %>% length,
                     Trt_IMBCR = veg_data %>%
                       filter(str_sub(Point_year, 4, 8) == "BCR16" &
                                Trt_stat == 1) %>%
                       pull(Point_year) %>% str_sub(1, -6) %>% unique %>% length,
                     Cnt_IMBCR = veg_data %>%
                       filter(str_sub(Point_year, 4, 8) == "BCR16" &
                                Trt_stat == 0) %>%
                       pull(Point_year) %>% str_sub(1, -6) %>% unique %>% length))

sum.grid <- landscape_data %>%
  group_by(Year) %>%
  summarize(Trt_CFLRP = sum(PctTrt_1kmNB > 0 & (str_sub(Grid, 4, 8) == "CFLRP")),
            Cnt_CFLRP = sum(PctTrt_1kmNB == 0 & (str_sub(Grid, 4, 8) == "CFLRP")),
            Trt_IMBCR = sum(PctTrt_1kmNB > 0 & (str_sub(Grid, 4, 8) == "BCR16")),
            Cnt_IMBCR = sum(PctTrt_1kmNB == 0 & (str_sub(Grid, 4, 8) == "BCR16"))) %>%
  bind_rows(data.frame(Year = "ALL", stringsAsFactors = F) %>%
              mutate(Trt_CFLRP = landscape_data %>%
                       filter(str_sub(Grid, 4, 8) == "CFLRP" &
                                PctTrt_1kmNB > 0) %>%
                       pull(Grid) %>% unique %>% length,
                     Cnt_CFLRP = landscape_data %>%
                       filter(str_sub(Grid, 4, 8) == "CFLRP" &
                                PctTrt_1kmNB == 0) %>%
                       pull(Grid) %>% unique %>% length,
                     Trt_IMBCR = landscape_data %>%
                       filter(str_sub(Grid, 4, 8) == "BCR16" &
                                PctTrt_1kmNB > 0) %>%
                       pull(Grid)%>% unique %>% length,
                     Cnt_IMBCR = landscape_data %>%
                       filter(str_sub(Grid, 4, 8) == "BCR16" &
                                PctTrt_1kmNB == 0) %>%
                       pull(Grid) %>% unique %>% length))

sum.table <- sum.grid %>%
  mutate(Level = "Grid") %>%
  bind_rows(sum.point %>%
              mutate(Level = "Point")) %>%
  select(Level, Year:Cnt_IMBCR)

write.csv(sum.table, "Treatment_sampling_summary.csv", row.names = F)

#### Tabulate sample sizes by treatment status and timing ####
dat <- Cov %>% tbl_df %>%
  mutate(Point = row.names(Cov)) %>%
  select(Point, gridIndex:Rdens)

rows <- c("Untreated", "Treated", str_c("Trt_yr", sort(unique(dat$Trt_time))))
cols <- c("n.point", "n.grid")
out <- matrix(NA, nrow = length(rows), ncol = length(cols),
              dimnames = list(rows, cols))

out["Untreated", ] <- c(
  sum(dat$Trt_stat == 0),
  sum(tapply(dat$Trt_stat, dat$gridIndex, function(x) any(x == 0)))
)

out["Treated", ] <- c(
  sum(dat$Trt_stat == 1),
  sum(tapply(dat$Trt_stat, dat$gridIndex, function(x) any(x == 1)))
)

for(i in sort(unique(dat$Trt_time))) 
  out[str_c("Trt_yr", i), ] <- c(
    sum(dat$Trt_time == i, na.rm = T),
    sum(tapply(dat$Trt_time, dat$gridIndex, function(x) any(x == i, na.rm = T)))
  )

out %>% tbl_df %>%
  mutate(Status = row.names(out)) %>%
  mutate(n = str_c(n.point, "(", n.grid, ")")) %>%
  select(Status, n) %>%
  write.csv("Sample_size_table.csv", row.names = F)

# Adjust data following review of this figure
dat.LP <- dat.LP %>%
  mutate(YSO = replace(YSO, which(Outbreak == 0), -1))

dat.SF <- dat.SF %>%
  mutate(YSO = replace(YSO, which(Outbreak == 0), -1))
