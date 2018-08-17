library(stringr)
library(BCRDataAPI)
library(dplyr)

#setwd("/home/RMBO.LOCAL/quresh.latif/CFLRP")
setwd("C:/Users/Quresh.Latif/files/projects/FS/CFLRP")

################## Inputs ####################
trunc.pct <- 0.95
strata <- c("CO-CFLRP-CF", "CO-BCR16-RC", "CO-BCR16-PC")
SampDesign <- c("IMBCR", "GRTS")
nG <- 10 # number of distance categories
##############################################

# Data grab #
BCRDataAPI::set_api_server('192.168.137.180')

BCRDataAPI::reset_api()
BCRDataAPI::add_columns(c('TransectNum|str',
                          'Point|int',
                          'Stratum|str',
                          'Year|int',
                          'Date|str',
                          'PointVisitStartTime|str',
                          'easting|int',
                          'northing|int',
                          'zone|int',
                          'BirdCode|str',
                          'radialDistance|int',
                          'TimePeriod|int',
                          'CL_ID|str',
                          'CL_Count|int'))
BCRDataAPI::filter_on(c('SelectionMethod in IMBCR,GRTS',
                      str_c('Stratum in ', str_c(strata, collapse = ",")),
                      'BirdCode in RESQ',
                      #'TransectVisitExcludeAnalysis = FALSE',
                      #'TransectExcludeAnalysis = FALSE',
                      'ninetynine = 0',
                      'eightyeight = 0',
                      'How <> F',
                      'Sex <> J',
                      'Migrant = 0',
                      'radialDistance > -1',
                      'TimePeriod > -1',
                      'Year in 2014,2015,2016'))
grab <- BCRDataAPI::get_data(interpolate_effort=TRUE)

# Derive parameters for distance sampling #
cutoff <- quantile(grab$radialDistance, trunc.pct, na.rm=TRUE) # truncation distance
area.circle <- as.numeric(pi * (cutoff / 1000) ^ 2) # area of point count circle in km^2
breaks <- seq(0, cutoff, length.out = nG + 1) # breaks for distance categories
area.band <- (pi * breaks[-1]^2) - (pi * breaks[-(nG+1)]^2) # area of each distance category
area.prop <- area.band / sum(area.band)

#***Note: No clusters broken up on separate lines for RESQ, so no processing of clusters needed***
library(timeDate)
grab.proc <- grab %>%
  mutate(CL_ID = replace(CL_ID, which(radialDistance < 0 | radialDistance >= cutoff), NA)) %>%
  mutate(CL_Count = replace(CL_Count, which(radialDistance < 0 | radialDistance >= cutoff), NA)) %>%
  mutate(radialDistance = replace(radialDistance, which(radialDistance < 0 | radialDistance >= cutoff), NA)) %>%
  mutate(radialDistance = replace(radialDistance, which(radialDistance == 0), 0.01)) %>%
  mutate(dclass = ceiling(radialDistance / breaks[2])) %>%
  
  ## Trim dates, compile day of year & start time in minutes ##
  mutate(Day = Date %>% str_sub(6, 7)) %>%
  mutate(Month = Date %>% str_sub(9, 11)) %>%
  mutate(MM = "05") %>%
  mutate(MM = replace(MM, which(Month == "Jun"), "06")) %>%
  mutate(MM = replace(MM, which(Month == "Jul"), "07")) %>%
  mutate(Date = str_c(MM, Day, sep = "-")) %>%
  mutate(DOY = Year %>% str_c("-", Date) %>% timeDate() %>% dayOfYear()) %>%
  mutate(PointVisitStartTime = PointVisitStartTime %>%
           replace(which(PointVisitStartTime == "0"), NA)) %>%
  mutate(HR = PointVisitStartTime %>% str_sub(1, -3) %>% as.integer) %>%
  mutate(MIN = PointVisitStartTime %>% str_sub(-2, -1) %>% as.integer) %>%
  mutate(Time_min = HR*60 + MIN) %>%
  mutate(Point_year = str_c(TransectNum, "-", str_pad(Point, width = 2, pad = "0", side = "left"), "-", Year)) %>%
  select(Point_year, TransectNum:Date, DOY, Time_min, radialDistance, dclass, TimePeriod, PointVisitStartTime:zone, CL_Count)

# Get covariate data #
# Canopy data #
BCRDataAPI::reset_api()
BCRDataAPI::set_api_server('192.168.137.180')
BCRDataAPI::add_columns(c('TransectNum|str',
                          'Point|int',
                          'Year|int',
                          'Stratum|str',
                          'o_canopy_percent|int'
))

BCRDataAPI::filter_on(c(str_c('Stratum in ', str_c(strata, collapse = ",")),
                        str_c('SelectionMethod in ', str_c(SampDesign, collapse = ",")),
                        'Year in 2014,2015,2016',
                        'UnusableDataOverstorySpecies = FALSE'))
BCRDataAPI::group_by(c('TransectNum', 'Point', 'Year', 'Stratum', 'o_canopy_percent'))
veg_data <- BCRDataAPI::get_data() %>%
  mutate(Point_year = str_c(TransectNum, "-", str_pad(Point, width = 2, pad = "0", side = "left"), "-", Year)) %>%
  rename(CanCov = o_canopy_percent) %>%
  mutate(CanCov = replace(CanCov, which(CanCov == -1), NA)) %>%
  select(Point_year, CanCov) %>%
  filter(Point_year %in% grab.proc$Point_year)

# Shrub data #
BCRDataAPI::reset_api()
BCRDataAPI::set_api_server('192.168.137.180')
BCRDataAPI::add_columns(c('TransectNum|str',
                          'Point|int',
                          'Year|int',
                          'shrub_cover|int',
                          'shrub_mean_height|num'
))

BCRDataAPI::filter_on(c(str_c('Stratum in ', str_c(strata, collapse = ",")),
                        str_c('SelectionMethod in ', str_c(SampDesign, collapse = ",")),
                        'Year in 2014,2015,2016',
                        'ShrubSpecies = FALSE'
))
BCRDataAPI::group_by(c('TransectNum', 'Point', 'Year', 'shrub_cover', 'shrub_mean_height'))
shrub_data <- BCRDataAPI::get_data() %>%
  mutate(Point_year = str_c(TransectNum, "-", str_pad(Point, width = 2, pad = "0", side = "left"), "-", Year)) %>%
  rename(ShrubHt = shrub_mean_height) %>%
  mutate(shrub_cover = replace(shrub_cover, which(shrub_cover == -1), NA),
         ShrubHt = replace(ShrubHt, which(ShrubHt == -1), NA)) %>%
  mutate(ShrubVol = (((shrub_cover / 100) * pi * 50^2) * ShrubHt) ^ (1/3)) %>%
  select(Point_year, shrub_cover, ShrubVol)

# Impute missing shrub volume values based on shrub cover #
mod <- lm(ShrubVol ~ shrub_cover + I(shrub_cover^2), data = shrub_data)
ind.missing <- which(is.na(shrub_data$ShrubVol))
shrub_data$ShrubVol[ind.missing] <- predict(mod, shrub_data %>% filter(is.na(ShrubVol)))

shrub_data <- shrub_data %>%
  select(-shrub_cover)

veg_data <- veg_data %>% left_join(shrub_data, by = "Point_year")
rm(shrub_data, grab, ind.missing, mod)

# Spatial data #
dat.cov <- foreign::read.dbf("C:/Users/Quresh.Latif/files/GIS/FS/CFLRP/Bird_survey_point_coords.dbf", as.is = T) %>%
  tbl_df() %>%
  mutate(Trt_status = replace(Trt_status, which(Trt_status == "Not treat*"), NA) %>% as.integer())
dat.cov <- dat.cov %>%
  mutate(Year = 2014) %>%
  bind_rows(dat.cov %>%
              mutate(Year = 2015)) %>%
  bind_rows(dat.cov %>%
              mutate(Year = 2016)) %>%
  mutate(Point_year = str_c(TransectNu, "-", str_pad(Point, width = 2, pad = "0", side = "left"), "-", Year)) %>%
  mutate(Trt_status = replace(Trt_status, which(is.na(Trt_status)), 9999)) %>%
  mutate(Trt_stat = (Year > Trt_status) %>% as.integer)
dat.cov <- dat.cov %>%
  mutate(Trt_status = replace(Trt_status, which(Trt_status == 9999), NA)) %>%
  mutate(Trt_time = Year - Trt_status) %>%
  mutate(Trt_time = replace(Trt_time, which(Trt_stat == 0), NA)) %>%
  select(Point_year, Trt_stat, Trt_time, TWIP, Rdens_1km) %>%
  rename(Rdens = Rdens_1km)

dat.cov <- dat.cov %>% left_join(veg_data, by = "Point_year")
rm(veg_data)

## Compile detection data ##
cov.names <- c("gridIndex", "yearIndex", "DayOfYear", "Time", names(dat.cov)[-1])
pointXyear.list <- grab.proc$Point_year %>%
  unique %>% sort
Y.dist <- tapply(grab.proc$CL_Count, grab.proc$Point_year, sum, na.rm = T)
dclass <- grab.proc %>%
  filter(!is.na(grab.proc$CL_Count)) %>%
  mutate(PtYrInd = match(Point_year, names(Y.dist))) %>%
  select(PtYrInd, CL_Count, dclass) %>%
  as.matrix()

Cov <- matrix(NA, nrow = length(pointXyear.list), ncol = length(cov.names),
                 dimnames = list(pointXyear.list, cov.names))
Cov[, "gridIndex"] <- pointXyear.list %>% str_sub(1, -9) %>% as.factor %>% as.integer
Cov[, "yearIndex"] <- pointXyear.list %>% str_sub(-4, -1) %>% as.factor %>% as.integer
Cov[, "DayOfYear"] <- (grab.proc %>% filter(Point_year %in% pointXyear.list) %>%
                            select(Point_year, DOY) %>% unique %>% arrange(Point_year))$DOY
Cov[, "Time"] <- (grab.proc %>% filter(Point_year %in% pointXyear.list) %>%
                       select(Point_year, Time_min) %>% unique %>% arrange(Point_year))$Time_min
ind.vals <- which(pointXyear.list %in% (dat.cov$Point_year))
Cov[ind.vals, -c(1:4)] <- dat.cov %>%
  filter(Point_year %in% pointXyear.list) %>%
  arrange(Point_year) %>%
  select(-Point_year) %>%
  as.matrix

rm(ind.vals)
save.image("Data_compiled_RESQ.RData")
