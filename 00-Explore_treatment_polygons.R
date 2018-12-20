library(foreign)
library(dplyr)

dat <- read.dbf("C:/Users/Quresh.Latif/files/GIS/FS/CFLRP/FACTS_TrtPolys_edited.dbf", as.is = T)
unique(dat$FUND_CODE)
sum(dat$Shape_Area)

dat.trt.compl <- dat %>% filter(!Year_Com_1 %in% c("Unk", "NotTreated", "Not treated", "2016", "2017", "2018"))
sum(dat.trt.compl$Shape_Area)

unique(dat.trt.compl$TREATMENT_)
sum(dat.trt.compl %>% filter(TREATMENT_ == "Thinning") %>% pull(Shape_Area))
sum(dat.trt.compl %>% filter(TREATMENT_ == "Broadcast Burn") %>% pull(Shape_Area))

sum(dat %>% filter(Year_Com_1 == "Unk") %>% pull(Shape_Area))
