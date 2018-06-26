library(foreign)
library(dplyr)
library(stringr)

dat <- read.dbf("C:/Users/Quresh.Latif/files/GIS/FS/CFLRP/Bird_survey_point_coords.dbf", as.is = T) %>%
  mutate(Point = str_c(TransectNu, "-", Point %>% str_pad("0", width = 2, side = "left"))) %>%
  mutate(SuperStrat = str_sub(TransectNu, 1, 8))

# Summary tables #
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

write.csv(sum.table, "C:/Users/Quresh.Latif/files/projects/FS/CFLRP/Treatment_sampling_summary.csv", row.names = F)
