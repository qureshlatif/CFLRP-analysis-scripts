library(stringr)
library(BCRDataAPI)
library(dplyr)

#setwd("/home/RMBO.LOCAL/quresh.latif/CFLRP")
setwd("C:/Users/Quresh.Latif/files/projects/FS/CFLRP")

################## Inputs ####################
# Get latest AOU checklist with tax names and order #
aou.checklist <- read.csv("C:/Users/Quresh.Latif/files/data/NACC_list_bird_species_downloaded_20180319.csv",
                          header = T, stringsAsFactors = F) %>% tbl_df %>%
  mutate(tax_ord = row_number())

spp.exclude <- c("Squirrel, Red", "Ruffed Grouse", "Turkey Vulture", "Wild Turkey",
                 "Sandhill Crane", "Bald Eagle", "American Kestrel", "Red-tailed Hawk",
                 "Great Blue Heron", "Swainson's Hawk", "Canada Goose", "Squirrel, Abert's",
                 "Northern Pygmy-Owl", "Northern Goshawk", "Sharp-shinned Hawk", "Green-winged Teal",
                 "Cooper's Hawk", "Great Horned Owl", "Pika", "Gambel's Quail", "Osprey",
                 "Common Merganser", "White-tailed Ptarmigan", "Peregrine Falcon",
                 "Boreal Owl", "Spotted Owl", "Black-crowned Night-Heron", "Ring-necked Duck",
                 "California Gull", "Northern Saw-whet Owl", "Long-eared Owl", "Flammulated Owl",
                 "Prairie Falcon", "Northern Harrier", "American White Pelican", "Western Screech-Owl",
                 "Double-crested Cormorant", "Bufflehead", "Thicket Tinamou", "Dusky Grouse", "Mallard",
                 "Golden Eagle")
strata <- c("CO-CFLRP-CF", "CO-BCR16-RC", "CO-BCR16-PC", "CO-BCR16-VO", "CO-BCR16-PO")
SampDesign <- c("IMBCR", "GRTS")
##############################################

#### Compile species list ####
BCRDataAPI::reset_api()
BCRDataAPI::set_api_server('192.168.137.180')
BCRDataAPI::add_columns(c('Year|int',
                          'primaryHabitat|str',
                          'Stratum|str',
                          'BirdCode|str',
                          'Species|str')
)
BCRDataAPI::filter_on('SelectionMethod in IMBCR,GRTS')
BCRDataAPI::filter_on('Migrant = 0')
BCRDataAPI::filter_on('BirdCode <> NOBI')
BCRDataAPI::filter_on('BCR = 16')
grab <- BCRDataAPI::get_data()

spp.list <- grab %>%
  filter(Year %in% 2008:2017) %>%
  filter(primaryHabitat %in% c("LP", "MC", "II", "PP")) %>% #II?
  select(primaryHabitat, BirdCode, Species) %>%
  unique %>%
  filter(!str_sub(BirdCode, 1, 2) == "UN") %>%
  filter(!Species %in% spp.exclude)

# Collapsing sub-species and renamed species #
ss <- BCRDataAPI::subspecies()
spp.list <- spp.list %>%
  mutate(BirdCode = ss[BirdCode] %>% as.character)

spp.out <- spp.list %>%
  select(BirdCode, Species) %>%
  unique
dup.codes <- spp.out$BirdCode[which(duplicated(spp.out$BirdCode))] %>%
  unlist %>% as.character %>% unique
spp.out <- spp.out %>%
  group_by(BirdCode) %>%
  mutate(min_length = min(nchar(Species))) %>%
  mutate(Species = str_sub(Species, 1, min_length)) %>%
  select(BirdCode, Species) %>%
  # Additional tweaks #
  mutate(Species = replace(Species, which(Species %in% c("Western Scrub-Jay", "Woodhouse's Scrub")), "Woodhouse's Scrub-Jay")) %>%
  ungroup %>%
  unique

#sum(!spp.out$Species %in% aou.checklist$common_name) # check - should be zero
spp.out <- spp.out %>%
  rename(common_name = Species) %>%
  left_join(aou.checklist, by = "common_name") %>%
  select(tax_ord, BirdCode, common_name, species) %>%
  arrange(tax_ord)

spp.excluded <- grab %>%
  filter(Year %in% 2008:2017) %>%
  filter(primaryHabitat %in% c("LP", "MC", "II", "PP")) %>% #II?
  select(primaryHabitat, BirdCode, Species) %>%
  unique %>%
  filter(!str_sub(BirdCode, 1, 2) == "UN") %>%
  filter(Species %in% spp.exclude) %>%
  select(BirdCode, Species) %>%
  unique %>%
  rename(common_name = Species) %>%
  left_join(aou.checklist, by = "common_name") %>%
  select(tax_ord, BirdCode, common_name, species) %>%
  arrange(tax_ord)

#### Detection data ####
 BCRDataAPI::reset_api()
BCRDataAPI::set_api_server('192.168.137.180')
BCRDataAPI::add_columns(c('TransectNum|str',
                          'Point|int',
                          'Year|int',
                          'Date|str',
                          'PointVisitStartTime|str',
                          'easting|int',
                          'northing|int',
                          'zone|int',
                          'Stratum|str',
                          'radialDistance|int',
                          'CL_count|int',
                          'BirdCode|str',
                          'How|str',
                          'Sex|str',
                          'TimePeriod|int'
))

BCRDataAPI::filter_on(str_c('Stratum in ', str_c(strata, collapse = ",")))
BCRDataAPI::filter_on(str_c('SelectionMethod in ', str_c(SampDesign, collapse = ",")))
BCRDataAPI::filter_on('ninetynine = 0')
BCRDataAPI::filter_on('eightyeight = 0')
BCRDataAPI::filter_on('How <> F')
BCRDataAPI::filter_on('Sex <> J')
BCRDataAPI::filter_on('Migrant = 0')
BCRDataAPI::filter_on('TimePeriod > -1')
BCRDataAPI::filter_on('radialDistance < 125')
grab <- BCRDataAPI::get_data() %>%
  mutate(BirdCode = ss[BirdCode] %>% as.character)

point.coords <- grab %>%
  select(TransectNum, Point, easting, northing, zone) %>%
  unique
point.list <- unique(str_c(point.coords$TransectNum,
                           str_pad(point.coords$Point, width = 2, side = "left", pad = "0"), sep = "-")) %>%
  sort
grid.list <- unique(point.coords$TransectNum) %>% sort

## Point X years surveyed ##
pointXyears.list <- unique(str_c(grab$TransectNum,
                                 str_pad(grab$Point, width = 2,
                                         side = "left", pad = "0"),
                                 grab$Year, sep = "-")) %>% sort

## Make sure all spp in data are in spp.list ##

## Add number of detections and count summaries to spp.out by stratum ##
smry <- grab %>% select(BirdCode, TransectNum, Point, Year) %>%
  unique %>% group_by(BirdCode) %>% count() %>%
  rename(Detections = n)
spp.out <- spp.out %>% left_join(smry, by = "BirdCode")

smry <- grab %>% select(BirdCode, CL_count) %>%
  group_by(BirdCode) %>%
  summarise(sumCount = sum(CL_count))
spp.out <- spp.out %>% left_join(smry, by = "BirdCode")

spp.out <- spp.out %>% # replace NAs with zeros
  mutate_at(vars(Detections, sumCount), (function(x) replace(x, is.na(x), 0)))

spp.out <- spp.out %>% # compile ratio of count totals to number of detections for spp with > 30 detections #
  mutate(RatioCountToDet_LP = sumCount_LP / Detections_LP) %>%
  mutate(RatioCountToDet_LP = replace(RatioCountToDet_LP, which(Detections_LP < 30), NA)) %>%
  mutate(RatioCountToDet_SF = sumCount_SF / Detections_SF) %>%
  mutate(RatioCountToDet_SF = replace(RatioCountToDet_SF, which(Detections_SF < 30), NA))

maxDetPossible <- length(pointXyears.list) # max possible by stratum
names(spp.out)[which(names(spp.out) == "Detections")] <-
  str_c("Detections (max = ", maxDetPossible, ")")

write.csv(spp.out, "Spp_list.csv", row.names = F)
rm(smry)

## Add number of detections and count summaries to excluded species by stratum ##
smry <- grab %>% filter(Stratum == "CO-CPW-LP") %>%
  select(BirdCode, TransectNum, Point) %>% unique %>%
  group_by(BirdCode) %>% count() %>%
  rename(Detections_LP = n)
spp.excluded <- spp.excluded %>% left_join(smry, by = "BirdCode")

smry <- grab %>% filter(Stratum == "CO-CPW-SF") %>%
  select(BirdCode, TransectNum, Point) %>% unique %>%
  group_by(BirdCode) %>% count() %>%
  rename(Detections_SF = n)
spp.excluded <- spp.excluded %>% left_join(smry, by = "BirdCode")

smry <- grab %>% filter(Stratum == "CO-CPW-LP") %>%
  select(BirdCode, TransectNum, Point, CL_count) %>%
  group_by(BirdCode) %>%
  summarise(sumCount_LP = sum(CL_count))
spp.excluded <- spp.excluded %>% left_join(smry, by = "BirdCode")

smry <- grab %>% filter(Stratum == "CO-CPW-SF") %>%
  select(BirdCode, TransectNum, Point, CL_count) %>%
  group_by(BirdCode) %>%
  summarise(sumCount_SF = sum(CL_count))
spp.excluded <- spp.excluded %>% left_join(smry, by = "BirdCode")

spp.excluded <- spp.excluded %>% # replace NAs with zeros
  mutate_at(vars(Detections_LP, Detections_SF, sumCount_LP, sumCount_SF),
            (function(x) replace(x, is.na(x), 0)))

write.csv(spp.excluded, "Spp_excluded.csv", row.names = F)
rm(smry)

#**** Stopped here **** Did not finish adapting the rest yet from CPW analysis - need to wait for covariates

## Trim dates, compile day of year & start time in minutes ##
library(timeDate)
grab <- grab %>%
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
  select(TransectNum:Date, DOY, Time_min, PointVisitStartTime:TimePeriod)
  
## Habitat covariate data ##
cov_tab_import <- read.csv("Covariates.csv", header = T, stringsAsFactors = F) %>%
  tbl_df() %>%
  mutate(YSI = replace(YSI, which(YSI == -1), NA)) %>%
  rename(DeadConif = DeadConif_RCov) %>%
  rename(YSO = YSI) %>%
  select(Point, TWIP, DeadConif, YSO, CanCov, RCOV_AS, RCOV_ES, RCOV_Pine, shrub_cover,
         RCShrb_UD, HerbCov)

## Compile multidimensional detection data array ##
spp.list <- spp.out$BirdCode
cov.names <- c("gridIndex", "DayOfYear", "Time", names(cov_tab_import)[-1])
grab <- grab %>%
  mutate(Point = TransectNum %>% str_sub(8, -1) %>%
           str_c("-", (Point %>% str_pad(width = 2, pad = 0, side = "left"))))

# Lodgepole pine stratum #
point.list.LP <- point.list[which(str_sub(point.list, 8, 9) == "LP")] %>%
  str_sub(8, -1) %>% sort
Y.LP <- matrix(NA, nrow = length(point.list.LP), ncol = length(spp.list),
               dimnames = list(point.list.LP, spp.list))
T.LP <- matrix(6, nrow = length(point.list.LP), ncol = length(spp.list),
               dimnames = list(point.list.LP, spp.list))
for(sp in 1:length(spp.list)) {
  obs <- grab %>% filter(BirdCode == spp.list[sp] & str_sub(Point, 1, 2) == "LP")
  if(nrow(obs) > 0) {
    Y.LP[, sp] <- (point.list.LP %in% obs$Point) %>% as.integer
    tvec <- tapply(obs$TimePeriod, obs$Point, min)
    tvec <- tvec[order(names(tvec))]
    T.LP[which(point.list.LP %in% obs$Point), sp] <- tvec
  } else {
    Y.LP[, sp] <- 0
  }
}

Cov.LP <- matrix(NA, nrow = length(point.list.LP), ncol = length(cov.names),
                 dimnames = list(point.list.LP, cov.names))
Cov.LP[, "gridIndex"] <- point.list.LP %>% str_sub(1, -4) %>% as.factor %>% as.integer
Cov.LP[, "DayOfYear"] <- (grab %>% filter(Stratum == "CO-CPW-LP") %>%
                            select(Point, DOY) %>% unique %>% arrange(Point))$DOY
Cov.LP[, "Time"] <- (grab %>% filter(Stratum == "CO-CPW-LP") %>%
                            select(Point, Time_min) %>% unique %>% arrange(Point))$Time_min
ind.vals <- which(point.list.LP %in% (cov_tab_import$Point %>% str_sub(8, -1)))
Cov.LP[ind.vals, -c(1:3)] <- cov_tab_import %>%
  mutate(Point = Point %>% str_sub(8, -1)) %>%
  filter(Point %in% point.list.LP) %>%
  arrange(Point) %>%
  select(-Point) %>%
  as.matrix

# Spruce-fir stratum #
point.list.SF <- point.list[which(str_sub(point.list, 8, 9) == "SF")] %>%
  str_sub(8, -1) %>% sort
Y.SF <- matrix(NA, nrow = length(point.list.SF), ncol = length(spp.list),
               dimnames = list(point.list.SF, spp.list))
T.SF <- matrix(6, nrow = length(point.list.SF), ncol = length(spp.list),
               dimnames = list(point.list.SF, spp.list))
for(sp in 1:length(spp.list)) {
  obs <- grab %>% filter(BirdCode == spp.list[sp] & str_sub(Point, 1, 2) == "SF")
  if(nrow(obs) > 0) {
    Y.SF[, sp] <- (point.list.SF %in% obs$Point) %>% as.integer
    tvec <- tapply(obs$TimePeriod, obs$Point, min)
    tvec <- tvec[order(names(tvec))]
    T.SF[which(point.list.SF %in% obs$Point), sp] <- tvec
  } else {
    Y.SF[, sp] <- 0
  }
}

Cov.SF <- matrix(NA, nrow = length(point.list.SF), ncol = length(cov.names),
                 dimnames = list(point.list.SF, cov.names))
Cov.SF[, "gridIndex"] <- point.list.SF %>% str_sub(1, -4) %>% as.factor %>% as.integer
Cov.SF[, "DayOfYear"] <- (grab %>% filter(Stratum == "CO-CPW-SF") %>%
                            select(Point, DOY) %>% unique %>% arrange(Point))$DOY
Cov.SF[, "Time"] <- (grab %>% filter(Stratum == "CO-CPW-SF") %>%
                            select(Point, Time_min) %>% unique %>% arrange(Point))$Time_min
ind.vals <- which(point.list.SF %in% (cov_tab_import$Point %>% str_sub(8, -1)))
Cov.SF[ind.vals, -c(1:3)] <- cov_tab_import %>%
  mutate(Point = Point %>% str_sub(8, -1)) %>%
  filter(Point %in% point.list.SF) %>%
  arrange(Point) %>%
  select(-Point) %>%
  as.matrix

rm(ind.vals, obs, maxDetPossible, sp, tvec)
save.image("Data_compiled.RData")

## Correlation matrices ##
Cov.LP[, c("TWIP", "DeadConif", "CanCov", "RCOV_AS", "RCOV_ES", "RCOV_Pine",
           "shrub_cover", "RCShrb_UD", "HerbCov")] %>% cor(use = "complete")
#plot(Cov.LP[, "RCOV_Pine"], Cov.LP[, "RCOV_AS"]) #Drop RCOV_ES from LP strata analysis.

Cov.SF[, c("TWIP", "DeadConif", "CanCov", "RCOV_AS", "RCOV_ES", "RCOV_Pine",
           "shrub_cover", "RCShrb_UD", "HerbCov")] %>% cor(use = "complete")
