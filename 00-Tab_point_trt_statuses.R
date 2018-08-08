library(foreign)
library(dplyr)
library(stringr)
setwd("C:/Users/Quresh.Latif/files/projects/FS/CFLRP")

dat <- read.dbf("Bird_survey_point_coords.dbf", as.is = T) %>%
  mutate(Point = str_c(TransectNu, "-", Point %>% str_pad("0", width = 2, side = "left"))) %>%
  mutate(SuperStrat = str_sub(TransectNu, 1, 8))

# Treatment status by year surveyed and strata #
sum.point <- dat %>%
  summarize(Trt_CFLRP = sum(Trt_status <= 2014 & SuperStrat == "CO-CFLRP", na.rm = T),
            Cnt_CFLRP = sum(Trt_status > 2014 & SuperStrat == "CO-CFLRP", na.rm = T) +
              sum(SuperStrat == "CO-CFLRP" & is.na(Trt_status)),
            Trt_IMBCR = sum(Trt_status <= 2014 & SuperStrat == "CO-BCR16", na.rm = T),
            Cnt_IMBCR = sum(Trt_status > 2014 & SuperStrat == "CO-BCR16", na.rm = T) +
              sum(SuperStrat == "CO-BCR16" & is.na(Trt_status))) %>%
  mutate(Year = 2014) %>%
  bind_rows(dat %>%
              summarize(Trt_CFLRP = sum(Trt_status <= 2015 & SuperStrat == "CO-CFLRP", na.rm = T),
                        Cnt_CFLRP = sum(Trt_status > 2015 & SuperStrat == "CO-CFLRP", na.rm = T) +
                          sum(SuperStrat == "CO-CFLRP" & is.na(Trt_status)),
                        Trt_IMBCR = sum(Trt_status <= 2015 & SuperStrat == "CO-BCR16", na.rm = T),
                        Cnt_IMBCR = sum(Trt_status > 2015 & SuperStrat == "CO-BCR16", na.rm = T) +
                          sum(SuperStrat == "CO-BCR16" & is.na(Trt_status))) %>%
              mutate(Year = 2015)) %>%
  bind_rows(dat %>%
              summarize(Trt_CFLRP = sum(Trt_status <= 2016 & SuperStrat == "CO-CFLRP", na.rm = T),
                        Cnt_CFLRP = sum(Trt_status > 2016 & SuperStrat == "CO-CFLRP", na.rm = T) +
                          sum(SuperStrat == "CO-CFLRP" & is.na(Trt_status)),
                        Trt_IMBCR = sum(Trt_status <= 2016 & SuperStrat == "CO-BCR16", na.rm = T),
                        Cnt_IMBCR = sum(Trt_status > 2016 & SuperStrat == "CO-BCR16", na.rm = T) +
                          sum(SuperStrat == "CO-BCR16" & is.na(Trt_status))) %>%
              mutate(Year = 2016)) %>%
  select(Year, Trt_CFLRP:Cnt_IMBCR)

sum.grid <- dat %>%
  group_by(TransectNu) %>%
  summarize(Trt = sum(Trt_status <= 2014, na.rm = T)) %>%
  ungroup %>%
  mutate(SuperStrat = str_sub(TransectNu, 1, 8)) %>%
  summarize(Trt_CFLRP = sum(Trt > 0 & SuperStrat == "CO-CFLRP"),
            Cnt_CFLRP = sum(Trt == 0 & SuperStrat == "CO-CFLRP"),
            Trt_IMBCR = sum(Trt > 0 & SuperStrat == "CO-BCR16"),
            Cnt_IMBCR = sum(Trt == 0 & SuperStrat == "CO-BCR16")) %>%
  mutate(Year = 2014) %>%
  bind_rows(dat %>%  group_by(TransectNu) %>%
              summarize(Trt = sum(Trt_status <= 2015, na.rm = T)) %>%
              ungroup %>%
              mutate(SuperStrat = str_sub(TransectNu, 1, 8)) %>%
              summarize(Trt_CFLRP = sum(Trt > 0 & SuperStrat == "CO-CFLRP"),
                        Cnt_CFLRP = sum(Trt == 0 & SuperStrat == "CO-CFLRP"),
                        Trt_IMBCR = sum(Trt > 0 & SuperStrat == "CO-BCR16"),
                        Cnt_IMBCR = sum(Trt == 0 & SuperStrat == "CO-BCR16")) %>%
              mutate(Year = 2015)) %>%
  bind_rows(dat %>%  group_by(TransectNu) %>%
              summarize(Trt = sum(Trt_status <= 2016, na.rm = T)) %>%
              ungroup %>%
              mutate(SuperStrat = str_sub(TransectNu, 1, 8)) %>%
              summarize(Trt_CFLRP = sum(Trt > 0 & SuperStrat == "CO-CFLRP"),
                        Cnt_CFLRP = sum(Trt == 0 & SuperStrat == "CO-CFLRP"),
                        Trt_IMBCR = sum(Trt > 0 & SuperStrat == "CO-BCR16"),
                        Cnt_IMBCR = sum(Trt == 0 & SuperStrat == "CO-BCR16")) %>%
              mutate(Year = 2016)) %>%
  select(Year, Trt_CFLRP:Cnt_IMBCR)

sum.table <- sum.grid %>%
  mutate(Level = "Grid") %>%
  bind_rows(sum.point %>%
              mutate(Level = "Point")) %>%
  select(Level, Year:Cnt_IMBCR)

write.csv(sum.table, "Treatment_sampling_summary.csv", row.names = F)

#### Tabulate sample sizes by treatment status and timing ####
load("Data_compiled.RData")

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
