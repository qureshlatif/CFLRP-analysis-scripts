library(BCRDataAPI)
library(timeDate)
library(stringr)
library(dplyr)

setwd("C:/Users/Quresh.Latif/files/projects/FS/CFLRP")

#### Variables ####
trunc.pct <- 0.95
strata <- c("CO-CFLRP-CF", "CO-BCR16-RC", "CO-BCR16-PC")
nG <- 10 # number of distance categories
Spp <- c("BTLH", "WISA", "WEWP", "WAVI", "STJA", "CLNU", "CORA", "TRES", # These are all spp with >= 60 detections and a count:detection >= 1.2
         "VGSW", "MOCH", "PYNU", "HOWR", "RCKI", "WEBL", "MOBL", "HETH",
         "AMRO", "EVGR", "CAFI", "RECR", "PISI", "GTTO", "SPTO", "CHSP",
         "LISP", "DEJU", "VIWA", "YRWA", "WETA", "RESQ")
SampDesign <- c("IMBCR", "GRTS")
dropGrids <- c("CO-BCR16-PC24", "CO-BCR16-PC27", "CO-BCR16-PC28", "CO-BCR16-PC31", "CO-BCR16-RC16")
###################

# Data grab #
BCRDataAPI::set_api_server('analysis.api.birdconservancy.org')

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
                      str_c('TransectNum not_in ', str_c(dropGrids, collapse = ",")), 
                      str_c('BirdCode in ', str_c(Spp, collapse = ",")),
                      str_c('SelectionMethod in ', str_c(SampDesign, collapse = ",")),
                      'ninetynine = 0',
                      'eightyeight = 0',
                      'How <> F',
                      'Sex <> J',
                      'Migrant = 0',
                      'radialDistance > -1',
                      'TimePeriod > -1'))
grab <- BCRDataAPI::get_data(interpolate_effort=TRUE) %>%
  mutate(Year = str_sub(Date, 13, 16) %>% as.integer) %>%
  filter(Year %in% 2014:2016)
any(!Spp %in% grab$BirdCode) # Should be false

# Derive parameters for distance sampling #
cutoff.list <- area.circle.list <- numeric(length = length(Spp))
breaks.list <- area.band.list <- area.prop.list <- rep(list(NULL), length(Spp))
for(i in 1:length(Spp)) {
  dat <- grab %>% filter(BirdCode == Spp[i])
  cutoff.list[i] <- quantile(dat$radialDistance, trunc.pct, na.rm=TRUE) # truncation distance
  area.circle.list[i] <- as.numeric(pi * (cutoff.list[i] / 1000) ^ 2) # area of point count circle in km^2
  breaks.list[[i]] <- seq(0, cutoff.list[i], length.out = nG + 1) # breaks for distance categories
  area.band.list[[i]] <- (pi * breaks.list[[i]][-1]^2) - (pi * breaks.list[[i]][-(nG+1)]^2) # area of each distance category
  area.prop.list[[i]] <- area.band.list[[i]] / sum(area.band.list[[i]])
}
rm(dat, i)

## Consolidate clusters ##
#grab %>% filter(!is.na(CL_ID)) %>% # Check within-cluster range of distances
#  dplyr::group_by(TransectNum, Point,
#                  Year, TimePeriod, CL_ID) %>%
#  summarize(DistDiff = max(radialDistance) - min(radialDistance)) %>%
#  View

grab.proc <- grab %>% filter(is.na(CL_ID)) %>%
  select(TransectNum, Point, Year, TimePeriod, BirdCode, Stratum, Date,
         PointVisitStartTime, radialDistance, CL_Count) %>%
  bind_rows(
    grab %>% filter(!is.na(CL_ID)) %>%
      dplyr::group_by(TransectNum, Point, Year,
                    TimePeriod, BirdCode, CL_ID) %>%
      summarize(Stratum = unique(Stratum),
                Date = unique(Date),
                PointVisitStartTime = unique(PointVisitStartTime),
                radialDistance = mean(radialDistance),
                CL_Count = sum(CL_Count)) %>%
      select(-CL_ID)
  )

## Additional processing ##
grab.proc <- grab.proc %>%
  mutate(radialDistance = replace(radialDistance, which(radialDistance == 0), 0.01)) %>%
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
  mutate(Point_year = str_c(TransectNum,
                            str_pad(Point, width = 2, side = "left", pad = "0"),
                            Year, sep = "-")) %>%
  select(Point_year, TransectNum:Date, DOY, Time_min, PointVisitStartTime, radialDistance, CL_Count)

for(i in 1:length(Spp)) {
  dat.spp <- grab.proc %>%
    filter(BirdCode == Spp[i] & !radialDistance >= cutoff.list[i]) %>%
    mutate(dclass = ceiling(radialDistance / breaks.list[[i]][2])) %>%
    select(Point_year:radialDistance, dclass, CL_Count)
  assign(str_c("detects.", Spp[i]), dat.spp)
}

dat.detect.all <- grab.proc
rm(grab.proc, dat.spp)

#### Covariate data ####
pointXyears.list <- dat.detect.all %>%
  pull(Point_year) %>% unique %>% sort

# Canopy data #
BCRDataAPI::reset_api()
BCRDataAPI::set_api_server('analysis.api.birdconservancy.org')
BCRDataAPI::add_columns(c('TransectNum|str',
                          'Point|int',
                          'Year|int',
                          'Stratum|str',
                          'o_canopy_percent|int',
                          'OverstorySpecies|str',
                          'OverstorySpeciesAbundance|int',
                          'OverstoryCommonName|str',
                          'OverstoryScientificName|str',
                          'NumSnags|int',
                          'primaryHabitat|str',
                          'HabitatCommonName',
                          'o_mean_height|int'
))

BCRDataAPI::filter_on(c(str_c('Stratum in ', str_c(strata, collapse = ",")),
                        str_c('TransectNum not_in ', str_c(dropGrids, collapse = ",")), 
                        str_c('SelectionMethod in ', str_c(SampDesign, collapse = ",")),
                        'Year in 2014,2015,2016',
                        'UnusableDataOverstorySpecies = FALSE'))
grab <- BCRDataAPI::get_data() %>%
  mutate(Point_year = str_c(TransectNum, "-", str_pad(Point, width = 2, pad = "0", side = "left"), "-", Year)) %>%
  rename(CanCov = o_canopy_percent, CanHt = o_mean_height,
         Species = OverstorySpecies, Abundance = OverstorySpeciesAbundance) %>%
  mutate(CanCov = replace(CanCov, which(CanCov == -1), NA),
         CanHt = replace(CanHt, which(CanHt == -1), NA),
         NumSnags = replace(NumSnags, which(NumSnags == -1), NA)) %>%
  filter(Abundance != -1) # There are two records here. 100 - sum(Abundance for other species) </= 1, so decided to delete them.

# Primary habitats table #
#grab %>% select(primaryHabitat, HabitatCommonName) %>%
#  unique %>% arrange(primaryHabitat) %>%
#  View

# Overstory species table (refer to this for developing categories). #
#grab %>% select(Species, OverstoryCommonName, OverstoryScientificName) %>%
#  rename(ComName = OverstoryCommonName, SciName = OverstoryScientificName) %>%
#  unique %>% arrange(Species) %>%
#  View

veg_data <- data.frame(Point_year = pointXyears.list, stringsAsFactors = F) %>%
  left_join(grab %>%
              select(Point_year, CanCov, CanHt, primaryHabitat, NumSnags) %>%
              unique %>%
              # Calculate forest indicator
              mutate(Forest = (primaryHabitat %in% c("AS", "BU", "II", "LO", "LP", "MC", "OA", "PP", "SF")) %>% as.integer) %>%
              select(-primaryHabitat), by = "Point_year") %>%
  # Species group relative covers #
  left_join(grab %>%
              reshape2::dcast(Point_year ~ Species, sum, value.var = "Abundance") %>%
              rename(RCOV_PP = PP, RCOV_DF = DF, RCOV_AS = AS) %>%
              mutate(RCOV_Dead = DC + DD + DA + BC + BD + DJ) %>%
              select(Point_year, RCOV_PP, RCOV_DF, RCOV_AS, RCOV_Dead), by = "Point_year") %>%
  # Everything else (Keep track of what the majority of this is.) #
  left_join(grab %>%
              filter(!Species %in% c("PP", "DF", "AS", "DC", "DD", "DA", "BC", "BD", "DJ")) %>%
              reshape2::dcast(Point_year ~ Species, sum, value.var = "Abundance") %>%
              mutate(RCOV_OT = select(., -Point_year) %>% apply(1, sum)) %>%
              select(Point_year, RCOV_OT), by = "Point_year") %>%
  mutate_at(vars(RCOV_PP:RCOV_OT), funs(replace(.,is.na(.),0))) %>%
  mutate(RCOV_TOT = RCOV_PP + RCOV_DF + RCOV_AS + RCOV_Dead + RCOV_OT) %>%
  # RCOV_TOT values != 100 (n = 83) are within 10% of 100, so rescaling them to sum to 100.
  mutate(RCOV_PP = (RCOV_PP / RCOV_TOT) *100) %>%
  mutate(RCOV_DF = (RCOV_DF / RCOV_TOT) *100) %>%
  mutate(RCOV_AS = (RCOV_AS / RCOV_TOT) *100) %>%
  mutate(RCOV_Dead = (RCOV_Dead / RCOV_TOT) *100) %>%
  mutate(RCOV_OT = (RCOV_OT / RCOV_TOT) *100) %>%
  select(-RCOV_TOT)

# Shrub data #
BCRDataAPI::reset_api()
BCRDataAPI::set_api_server('analysis.api.birdconservancy.org')
BCRDataAPI::add_columns(c('TransectNum|str',
                          'Point|int',
                          'Year|int',
                          'shrub_cover|int',
                          'shrub_mean_height|num',
                          'ShrubLayerSpecies|str',
                          'ShrubLayerCommonName|str',
                          'ShrubLayerScientificName|str',
                          'ShrubLayerSpeciesAbundance|int'
))

BCRDataAPI::filter_on(c(str_c('Stratum in ', str_c(strata, collapse = ",")),
                        str_c('TransectNum not_in ', str_c(dropGrids, collapse = ",")), 
                        str_c('SelectionMethod in ', str_c(SampDesign, collapse = ",")),
                        'Year in 2014,2015,2016',
                        'ShrubSpecies = FALSE'
))
grab <- BCRDataAPI::get_data() %>%
  mutate(Point_year = str_c(TransectNum, "-", str_pad(Point, width = 2, pad = "0", side = "left"), "-", Year)) %>%
  rename(ShrubHt = shrub_mean_height,
         Species = ShrubLayerSpecies, Abundance = ShrubLayerSpeciesAbundance) %>%
  mutate(shrub_cover = replace(shrub_cover, which(shrub_cover == -1), NA),
         ShrubHt = replace(ShrubHt, which(ShrubHt == -1), NA)) %>%
  filter(!Species %in% c("BG", "SN")) %>% # Remove stuff that isn't a shrub
  mutate(Abundance = replace(Abundance, which(Point_year == "CO-CFLRP-CF28-15-2014" & Species == "SU"), 30)) # Fill in missing value with 100 - sum(other spp)

# Major shrub spp (> 1000 records): CJ & JU (juniper), DF (Douglas fir), GB (Gooseberry / Currant / Ribes), AS (Aspen), PP (ponderosa pine)
# Minor shrub spp (300-500 records): CR (cliff rose, Purshia spp.), MM (Mountain Mahogany), LP (lodgepole pine), WX (waxflower), MA (Rocky Mountain Maple), ES (Englemann Spruce)
# Incidental spp (< 300 records): WI (Willow), WR (Wild rose), SC (Shrubby Cinquefoil), YU (Yucca), BB (blackberry), LM (Limber pine), CC (Choke cherry),
# BF (buffaloberry), SA (sagebrush), NB (ninebark), SY (snowberry), BS (blue spruce), UD (unknown deciduous), SB (service berry), SU (subalpine fir),
# MZ (manzanita), SK (skunkbrush), GO (Gambel oak), OT (other??), MS (mountain spray), HA (hawthorn), HB (huckleberry), BU (ragweed/bursage), AM (Apache plume),
# RA (rabbitbrush), DC (dead conifer), AL (alder), TW (twinberry), BR (bristlecone pine), WB (water birch), UC (unknown conifer), SW (snakeweed), WF (white fir),
# BC (burnt conifer), WN (winterfat), RD (Red-osier dogwood), DD (dead deciduous), VI (Viburnum), TA (Tamarisk), CA (Ceanothus), SL (Saltbush),
# PY (Pinyon pine), NC (narrow-leaf cottonwood), EB (elderberry), XX (not listed), DA (dead aspen), AP (American plum), PI (??), PC (Plains cottonwood),
# GW (greasewood), SI (??), SE (single-leaf ash), OB (oak bush), LU (??), LT (??), DJ (dead juniper), CE (creosote), BI (birch), AH (ash)

# Shrub categories:
# Ladder fuels - c("CJ", "JU", "DJ", "DF", "PP", "LP", "ES", "UC", "SU", "BR", "DD", "BC", "PY", "BS", "LM", "WF", "GO")
# Berry - c("GB", "SY", "BB", "CC", "SB", "TW", "HA", "EB")
# Aspen - c("AS", "DA")
# Xeric - c("CR", "MM", "WX", "YU", "BF", "SA", "MZ", "SK", "GO", "RA", "SL", "CE", "GW", "WN")
# Mesic - c("WI", "WR", "AL", "WB", "PC", "BI", "TA", "NC", "RD")
# Other (left out) - c("MA", "SC", "NB", "UD", "OT", "MS", "HA", "BU", "AM", "SW", "DD", "VI", "CA", "XX", "AP", "PI", "SI", "SE", "OB", "LU", "LT", "AH")

shrub_data <- grab %>%
  mutate(ShrubVol = (((shrub_cover / 100) * pi * 50^2) * ShrubHt) ^ (1/3)) %>%
  select(Point_year, shrub_cover, ShrubHt, ShrubVol) %>%
  unique %>%
  # Shrub diversity #
  left_join(grab %>%
              mutate(p = Abundance / 100) %>%
              dplyr::group_by(Point_year) %>%
              summarise(ShrubDiv = -1*sum(p * log(p))), by = "Point_year") %>%
  # Species group relative covers #
  left_join(grab %>%
              reshape2::dcast(Point_year ~ Species, sum, value.var = "Abundance") %>%
              rename(RSCV_AS = AS) %>%
              mutate(RSCV_Ladder = CJ + JU + DF + PP + LP + ES + UC + SU +
                       BR + DD + BC + PY + BS + LM + WF + GO + OB,
                     RSCV_Ber = GB + SY + BB + CC + SB + TW + HA + EB) %>%
              # Also considered but dropped:
              # Xeric: CR + MM + WX + YU + BF + SA + MZ + SK + GO + RA + SL + CE + GW + WN
              # Mesic: WI + WR + AL + WB + PC + BI + TA + NC + RD
              select(Point_year, RSCV_Ladder, RSCV_Ber, RSCV_AS), by = "Point_year") %>%
  # Everything else #
  left_join(grab %>%
              filter(!Species %in% c("CJ", "JU", "DJ","DF", "PP", "LP", "ES", "UC", "SU", "BR", "DD", "BC", "PY", "BS", "LM", "WF",
                                     "GB", "SY", "BB", "CC", "SB", "TW", "HB", "OG", "EB",
                                     "AS", "DA")) %>%
              reshape2::dcast(Point_year ~ Species, sum, value.var = "Abundance") %>%
              mutate(RSCV_OT = select(., -Point_year) %>% apply(1, sum)) %>%
              select(Point_year, RSCV_OT), by = "Point_year") %>%
  mutate_at(vars(RSCV_Ladder:RSCV_OT), funs(replace(.,is.na(.),0))) %>%
  mutate(RCOV_TOT = RSCV_Ladder + RSCV_Ber + RSCV_AS + RSCV_OT) %>% # RSCV_Jun + RSCV_Xer + RSCV_Mes + 
  # Rescale where TOT within 10% of 100...
  #mutate(RSCV_Jun = (RSCV_Jun / RCOV_TOT) *100) %>%
  mutate(RSCV_Ladder = (RSCV_Ladder / RCOV_TOT) *100) %>%
  mutate(RSCV_Ber = (RSCV_Ber / RCOV_TOT) *100) %>%
  mutate(RSCV_AS = (RSCV_AS / RCOV_TOT) *100) %>%
  #mutate(RSCV_Xer = (RSCV_Xer / RCOV_TOT) *100) %>%
  #mutate(RSCV_Mes = (RSCV_Mes / RCOV_TOT) *100) %>%
  mutate(RSCV_OT = (RSCV_OT / RCOV_TOT) *100) %>%
  # ...and dump the rest.
  #mutate(RSCV_Jun = replace(RSCV_Jun, which(RCOV_TOT < 90 | RCOV_TOT > 110), NA)) %>%
  mutate(RSCV_Ladder = replace(RSCV_Ladder, which(RCOV_TOT < 90 | RCOV_TOT > 110), NA)) %>%
  mutate(RSCV_Ber = replace(RSCV_Ber, which(RCOV_TOT < 90 | RCOV_TOT > 110), NA)) %>%
  mutate(RSCV_AS = replace(RSCV_AS, which(RCOV_TOT < 90 | RCOV_TOT > 110), NA)) %>%
  #mutate(RSCV_Xer = replace(RSCV_Xer, which(RCOV_TOT < 90 | RCOV_TOT > 110), NA)) %>%
  #mutate(RSCV_Mes = replace(RSCV_Mes, which(RCOV_TOT < 90 | RCOV_TOT > 110), NA)) %>%
  mutate(RSCV_OT = replace(RSCV_OT, which(RCOV_TOT < 90 | RCOV_TOT > 110), NA)) %>%
  select(-RCOV_TOT)

# Impute missing shrub volume values based on shrub cover #
mod <- lm(ShrubVol ~ shrub_cover + I(shrub_cover^2), data = shrub_data)
ind.missing <- which(is.na(shrub_data$ShrubVol))
shrub_data$ShrubVol[ind.missing] <- predict(mod, shrub_data %>% filter(is.na(ShrubVol)))

shrub_data <- shrub_data %>%
  select(-shrub_cover)

veg_data <- veg_data %>% left_join(shrub_data, by = "Point_year")
rm(shrub_data, grab, ind.missing, mod)

# Ground cover #
BCRDataAPI::reset_api()
BCRDataAPI::set_api_server('analysis.api.birdconservancy.org')
BCRDataAPI::add_columns(c('TransectNum|str',
                          'Point|int',
                          'Year|int',
                          'gc_woody|int',
                          'gc_live_grass|int',
                          'gc_grass|int',
                          'gc_herb|int',
                          'gc_bare_litter|int',
                          'gc_live_grass_height|int',
                          'gc_grass_height|int'
))
BCRDataAPI::group_by(c('TransectNum', 'Point', 'Year', 'gc_woody', 'gc_live_grass', 'gc_grass',
                       'gc_herb', 'gc_bare_litter', 'gc_live_grass_height', 'gc_grass_height'))
BCRDataAPI::filter_on(c(str_c('Stratum in ', str_c(strata, collapse = ",")),
                        str_c('TransectNum not_in ', str_c(dropGrids, collapse = ",")), 
                        str_c('SelectionMethod in ', str_c(SampDesign, collapse = ",")),
                        'Year in 2014,2015,2016',
                        'UnusableDataOverstorySpecies = FALSE'))
grab <- BCRDataAPI::get_data() %>%
  mutate(Point_year = str_c(TransectNum, "-", str_pad(Point, width = 2, pad = "0", side = "left"), "-", Year)) %>%
  rename(LiveGrass = gc_live_grass,
         DeadGrass = gc_grass,
         Herb = gc_herb,
         LGrassHt = gc_live_grass_height,
         DGrassHt = gc_grass_height) %>%
  mutate(LiveGrass = replace(LiveGrass, which(LiveGrass == -1), NA),
         DeadGrass = replace(DeadGrass, which(DeadGrass == -1), NA),
         Herb = replace(Herb, which(Herb == -1), NA)) %>%
  mutate(HerbGrass = LiveGrass + DeadGrass + Herb) %>%
  mutate(LGrassHt = replace(LGrassHt, which(LiveGrass == 0), 0),
         LGrassHt = replace(LGrassHt, which(is.na(LiveGrass)), NA),
         LGrassHt = replace(LGrassHt, which(LGrassHt == -1), NA),
         DGrassHt = replace(DGrassHt, which(DeadGrass == 0), 0),
         DGrassHt = replace(DGrassHt, which(is.na(DeadGrass)), NA),
         DGrassHt = replace(DGrassHt, which(DGrassHt == -1), NA)) %>%
  mutate(GrassHt = (LGrassHt * LiveGrass + DGrassHt * DeadGrass) / (LiveGrass + DeadGrass),
         HerbGrassVol = (((HerbGrass / 100) * (pi * 50^2)) * # coverage in m^2
                           (GrassHt / 100))^(1/3)) %>% # ht in m
  select(Point_year, HerbGrass, GrassHt, HerbGrassVol) %>%
  mutate(HerbGrassVol = replace(HerbGrassVol, which(HerbGrass == 0), 0))

# Impute missing shrub volume values based on shrub cover #
mod <- lm(HerbGrassVol ~ HerbGrass + I(HerbGrass^2), data = grab)
ind.missing <- which(is.na(grab$HerbGrassVol))
grab$HerbGrassVol[ind.missing] <- predict(mod, grab %>% filter(is.na(HerbGrassVol)))

grab <- grab %>%
  select(-one_of("HerbGrass", "GrassHt"))

veg_data <- veg_data %>% left_join(grab, by = "Point_year")
rm(grab, mod, ind.missing)

# Compile understory-overstory height ratio #
veg_data <- veg_data %>%
  mutate(SOHtRatio = ShrubHt / CanHt) %>%
  select(-ShrubHt)

# Get GIS covariates #
dat.gis <- foreign::read.dbf("C:/Users/Quresh.Latif/files/GIS/FS/CFLRP/Bird_survey_point_coords.dbf", as.is = T) %>%
  tbl_df() %>%
  mutate(Trt_status = replace(Trt_status, which(Trt_status == "Not treat*"), NA) %>% as.integer())
dat.gis <- dat.gis %>%
  mutate(Year = 2014) %>%
  bind_rows(dat.gis %>%
              mutate(Year = 2015)) %>%
  bind_rows(dat.gis %>%
              mutate(Year = 2016)) %>%
  mutate(Point_year = str_c(TransectNu, "-", str_pad(Point, width = 2, pad = "0", side = "left"), "-", Year)) %>%
  mutate(Trt_status = replace(Trt_status, which(is.na(Trt_status)), 9999)) %>%
  mutate(Trt_stat = (Year > Trt_status) %>% as.integer)
dat.gis <- dat.gis %>%
  mutate(Trt_status = replace(Trt_status, which(Trt_status == 9999), NA)) %>%
  mutate(Trt_time = Year - Trt_status) %>%
  mutate(Trt_time = replace(Trt_time, which(Trt_stat == 0), NA)) %>%
  mutate(ForAR = (Forest == "AR") %>% as.integer) %>%
  select(Point_year, Trt_stat, Trt_time, TWIP, heatload, TWI, Rdens_1km, ForAR) %>%
  rename(Rdens = Rdens_1km)

veg_data <- veg_data %>% left_join(dat.gis, by = "Point_year")

# Get grid-level landscape structure (LANDFIRE) covariates #
landscape_data <- data.frame(Grid = pointXyears.list %>% str_sub(1, -9),
                             Year = pointXyears.list %>% str_sub(-4, -1), stringsAsFactors = F) %>%
  mutate(gridIndex = Grid %>% as.factor %>% as.integer,
         YearInd = Year %>% as.factor %>% as.integer) %>%
  unique

PA_data <- read.csv("LANDFIRE_vars/CFLRP Percent Areas 2014.csv", header = T, stringsAsFactors = F) %>%
  tbl_df %>%
  mutate(Year = "2014") %>%
  bind_rows(
    read.csv("LANDFIRE_vars/CFLRP Percent Areas 2015.csv", header = T, stringsAsFactors = F) %>%
      tbl_df %>%
      mutate(Year = "2015")
  ) %>%
  bind_rows(
    read.csv("LANDFIRE_vars/CFLRP Percent Areas 2016.csv", header = T, stringsAsFactors = F) %>%
      tbl_df %>%
      mutate(Year = "2016")
  ) %>%
  rename(PACC10_3km = Gap_3km_square,
         PACC40_3km = Open_3km_square,
         PACC10_2km = Gap_2km_circle,
         PACC40_2km = Open_2km_circle)

landscape_data <- landscape_data %>%
  left_join(
    PA_data %>%
      select(TransectNum, Year, PACC10_3km, PACC40_3km),
    by = c("Grid" = "TransectNum", "Year" = "Year")
  ) %>%
  left_join(read.csv("LANDFIRE_vars/Patch_struct_&_config.csv", header = T, stringsAsFactors = F) %>%
              tbl_df %>%
              mutate(Year = as.character(Year),
                     # Rescale area and distance covariates to ha and km
                     mnPtchAr_Gap2km = mnPtchAr_Gap2km / 10000,
                     mnPtchAr_Opn2km = mnPtchAr_Opn2km / 10000,
                     mnPtchAr_Gap3km = mnPtchAr_Gap3km / 10000,
                     mnPtchAr_Opn3km = mnPtchAr_Opn3km / 10000,
                     NNdist_Gap2km = NNdist_Gap2km / 1000,
                     NNdist_Opn2km = NNdist_Opn2km / 1000,
                     NNdist_Gap3km = NNdist_Gap3km / 1000,
                     NNdist_Opn3km = NNdist_Opn3km / 1000),
            by = c("Grid" = "TransNum", "Year" = "Year"))

dat.gis <- dat.gis %>%
  mutate(Grid = str_sub(Point_year, 1, -9)) %>%
  mutate(Year = str_sub(Point_year, -4, -1)) %>%
  dplyr::group_by(Grid, Year) %>%
  summarize(PctTrt = mean(Trt_stat) * 100,
            Trt_time = mean(Trt_time, na.rm = T),
            TWIP = mean(TWIP),
            heatload = mean(heatload),
            TWI = mean(TWI),
            Rdens = mean(Rdens),
            ForAR = mean(ForAR))

landscape_data <- landscape_data %>%
  left_join(dat.gis, by = c("Grid", "Year"))

rm(PA_data)
rm(dat.gis)

# Get alternate grid-level treatment extents for 1-km radius neighborhoods #
landscape_data <- landscape_data %>% left_join(
  read.csv("Landscape_Treatment_1km-r_NB.csv", header = T, stringsAsFactors = F) %>% tbl_df %>%
    rename(Grid = TransNum, YST_1kmNB = MEAN_YST, PctTrt_1kmNB = Treatment.Area) %>%
    mutate(Year = as.character(Year)),
  by = c("Grid", "Year"))


## Compile detection data arrays ##
cov.names <- c("gridIndex", "YearInd", "DayOfYear", "Time", names(veg_data)[-1])

for(sp in 1:length(Spp)) {
  dat <- eval(as.name(str_c("detects.", Spp[sp])))
  dat <- dat %>%
    bind_rows(
      dat.detect.all %>%
        mutate(BirdCode = NA,
               TimePeriod = NA,
               radialDistance = NA,
               CL_Count = NA) %>%
        unique %>%
        filter(!Point_year %in% dat$Point_year)
    ) %>%
    filter(Point_year %in% pointXyears.list)
  Y <- tapply(dat$CL_Count, dat$Point_year, sum, na.rm = T)
  dclass <- dat %>%
    filter(!is.na(CL_Count)) %>%
    filter(Point_year %in% pointXyears.list) %>%
    select(Point_year, CL_Count, dclass) %>%
    mutate(Point_year = match(Point_year, names(Y))) %>%
    rename(Yindex = Point_year) %>%
    as.matrix()
  
  assign(str_c("Y.", Spp[sp], ".dist"), Y)
  assign(str_c("dclass.", Spp[sp]), dclass)
}

Cov <- matrix(NA, nrow = length(pointXyears.list), ncol = length(cov.names),
                 dimnames = list(pointXyears.list, cov.names))
Cov[, "gridIndex"] <- pointXyears.list %>% str_sub(1, -9) %>% as.factor %>% as.integer
Cov[, "YearInd"] <- pointXyears.list %>% str_sub(-4, -1) %>% as.factor %>% as.integer
Cov[, "DayOfYear"] <- (dat.detect.all %>% select(Point_year, DOY) %>%
                            filter(Point_year %in% pointXyears.list) %>% unique %>% arrange(Point_year))$DOY
Cov[, "Time"] <- (dat.detect.all %>% select(Point_year, Time_min) %>%
                       filter(Point_year %in% pointXyears.list) %>% unique %>% arrange(Point_year))$Time_min
Cov[, -c(1:4)] <- veg_data %>%
  filter(Point_year %in% pointXyears.list) %>%
  arrange(Point_year) %>%
  select(-Point_year) %>%
  as.matrix

rm(i, dat, dclass, Y, sp)
save.image("Data_compiled_abundance.RData")
