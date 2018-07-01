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
                 "Golden Eagle", "Gadwall", "Virginia Rail")
strata <- c("CO-CFLRP-CF", "CO-BCR16-RC", "CO-BCR16-PC")
SampDesign <- c("IMBCR", "GRTS")
##############################################

#### Compile species list ####
BCRDataAPI::reset_api()
BCRDataAPI::set_api_server('192.168.137.180')
BCRDataAPI::add_columns(c('BirdCode|str',
                          'Species|str')
)
BCRDataAPI::filter_on('SelectionMethod in IMBCR,GRTS')
BCRDataAPI::filter_on('Migrant = 0')
BCRDataAPI::filter_on('BirdCode <> NOBI')
BCRDataAPI::filter_on('BCR = 16')
BCRDataAPI::filter_on('primaryHabitat in LP,MC,II,PP')
BCRDataAPI::filter_on(str_c('Year in ', str_c(2008:2017, collapse = ",")))
BCRDataAPI::group_by(c('BirdCode', 'Species'))
grab <- BCRDataAPI::get_data()

spp.out <- grab %>%
  filter(!str_sub(BirdCode, 1, 2) == "UN") %>%
  filter(!Species %in% spp.exclude)

# Collapsing sub-species and renamed species #
ss <- BCRDataAPI::subspecies()
spp.out <- spp.out %>%
  mutate(BirdCode = ss[BirdCode] %>% as.character)%>%
  dplyr::group_by(BirdCode) %>%
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

# Remove additional implausible members of the metacommunity (based on review of BNA range maps and habitat accounts) #
spp.out <- spp.out %>%
  filter(!BirdCode %in% c("RUHU", "PAWR", "OLWA", "AMPI", "WWCR", "SAGS"))

spp.excluded <- grab %>%
  select(BirdCode, Species) %>%
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
                          'Species|str',
                          'How|str',
                          'Sex|str',
                          'TimePeriod|int'
))

BCRDataAPI::filter_on(str_c('Stratum in ', str_c(strata, collapse = ",")))
BCRDataAPI::filter_on(str_c('SelectionMethod in ', str_c(SampDesign, collapse = ",")))
BCRDataAPI::filter_on('Year in 2014,2015,2016')
BCRDataAPI::filter_on(str_c('BirdCode in ', str_c(spp.out$BirdCode, collapse = ",")))
BCRDataAPI::filter_on('ninetynine = 0')
BCRDataAPI::filter_on('eightyeight = 0')
BCRDataAPI::filter_on('How <> F')
BCRDataAPI::filter_on('Sex <> J')
BCRDataAPI::filter_on('Migrant = 0')
BCRDataAPI::filter_on('TimePeriod > -1')
BCRDataAPI::filter_on('radialDistance < 125')
grab <- BCRDataAPI::get_data(interpolate_effort = T) %>%
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

## Add number of detections and count summaries to spp.out by stratum ##
smry <- grab %>% select(BirdCode, TransectNum, Point, Year) %>%
  unique %>% dplyr::group_by(BirdCode) %>% count() %>%
  rename(Detections = n)
spp.out <- spp.out %>% left_join(smry, by = "BirdCode")

smry <- grab %>% select(BirdCode, CL_count) %>%
  dplyr::group_by(BirdCode) %>%
  summarise(sumCount = sum(CL_count))
spp.out <- spp.out %>% left_join(smry, by = "BirdCode")

spp.out <- spp.out %>% # replace NAs with zeros
  mutate_at(vars(Detections, sumCount), (function(x) replace(x, is.na(x), 0)))

spp.out <- spp.out %>% # compile ratio of count totals to number of detections for spp with > 30 detections #
  mutate(RatioCountToDet = sumCount / Detections) %>%
  mutate(RatioCountToDet = replace(RatioCountToDet, which(Detections < 60), NA))

maxDetPossible <- length(pointXyears.list) # max possible by stratum
names(spp.out)[which(names(spp.out) == "Detections")] <-
  str_c("Detections (max = ", maxDetPossible, ")")

write.csv(spp.out, "Spp_list.csv", row.names = F)
rm(smry)

## Add number of detections and count summaries to excluded species ##
smry <- grab %>% select(BirdCode, TransectNum, Point, Year) %>% unique %>%
  dplyr::group_by(BirdCode) %>% count() %>%
  rename(Detections = n)
spp.excluded <- spp.excluded %>% left_join(smry, by = "BirdCode")

smry <- grab %>% select(BirdCode, TransectNum, Point, CL_count) %>%
  dplyr::group_by(BirdCode) %>%
  summarise(sumCount = sum(CL_count))
spp.excluded <- spp.excluded %>% left_join(smry, by = "BirdCode")

spp.excluded <- spp.excluded %>% # replace NAs with zeros
  mutate_at(vars(Detections, sumCount),
            (function(x) replace(x, is.na(x), 0)))

write.csv(spp.excluded, "Spp_excluded.csv", row.names = F)
rm(smry)

bird_data <- grab %>%  # Store bird survey data for later use.
  mutate(Point_year = str_c(TransectNum, "-", str_pad(Point, width = 2, pad = "0", side = "left"), "-", Year))

#### Covariate data ####
# Canopy data #
BCRDataAPI::reset_api()
BCRDataAPI::set_api_server('192.168.137.180')
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
BCRDataAPI::set_api_server('192.168.137.180')
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
  select(Point_year, shrub_cover, ShrubHt) %>%
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
              mutate(RSCV_Ladder = CJ + JU + DJ + DF + PP + LP + ES + UC + SU +
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

veg_data <- veg_data %>% left_join(shrub_data, by = "Point_year")
rm(shrub_data, grab)

# Ground cover #
BCRDataAPI::reset_api()
BCRDataAPI::set_api_server('192.168.137.180')
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
         HerbGrassVol = ((HerbGrass / 100) * (pi * 50^2)) * # coverage in m^2
           (GrassHt / 100)) %>% # ht in m^2
  select(Point_year, HerbGrassVol)

veg_data <- veg_data %>% left_join(grab, by = "Point_year")
rm(grab)

# Compile understory-overstory height ratio #
veg_data <- veg_data %>%
  mutate(SOHtRatio = ShrubHt / CanHt)

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
  mutate(Trt_stat = (Year >= Trt_status) %>% as.integer)
dat.gis <- dat.gis %>%
  mutate(Trt_status = replace(Trt_status, which(Trt_status == 9999), NA)) %>%
  mutate(Trt_time = Year - Trt_status) %>%
  mutate(Trt_time = replace(Trt_time, which(Trt_stat == 0), NA)) %>%
  select(Point_year, Trt_stat, Trt_time, TWIP)

veg_data <- veg_data %>% left_join(dat.gis, by = "Point_year")
rm(dat.gis)

## Trim dates, compile day of year & start time in minutes ##
library(timeDate)
bird_data <- bird_data %>%
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

## Compile multidimensional detection data array ##
spp.list <- spp.out$BirdCode
cov.names <- c("gridIndex", "YearInd", "DayOfYear", "Time", names(veg_data)[-1])

bird_data <- bird_data %>%
  mutate(Point_year = str_c(TransectNum, "-", str_pad(Point, width = 2, pad = "0", side = "left"), "-", Year))
Y.mat <- matrix(NA, nrow = length(pointXyears.list), ncol = length(spp.list),
               dimnames = list(pointXyears.list, spp.list))
TR.mat <- matrix(6, nrow = length(pointXyears.list), ncol = length(spp.list),
               dimnames = list(pointXyears.list, spp.list))
for(sp in 1:length(spp.list)) {
  obs <- bird_data %>% filter(BirdCode == spp.list[sp])
  if(nrow(obs) > 0) {
    Y.mat[, sp] <- (pointXyears.list %in% obs$Point_year) %>% as.integer
    tvec <- tapply(obs$TimePeriod, obs$Point_year, min)
    tvec <- tvec[order(names(tvec))]
    TR.mat[which(pointXyears.list %in% obs$Point_year), sp] <- tvec
  } else {
    Y.mat[, sp] <- 0
  }
}

Cov <- matrix(NA, nrow = length(pointXyears.list), ncol = length(cov.names),
              dimnames = list(pointXyears.list, cov.names))
Cov[, "gridIndex"] <- pointXyears.list %>% str_sub(1, -9) %>% as.factor %>% as.integer
Cov[, "YearInd"] <- pointXyears.list %>% str_sub(-4, -1) %>% as.factor %>% as.integer
Cov[, "DayOfYear"] <- (bird_data %>%
                            select(Point_year, DOY) %>% unique %>% arrange(Point_year))$DOY
Cov[, "Time"] <- (bird_data %>%
                            select(Point_year, Time_min) %>% unique %>% arrange(Point_year))$Time_min
ind.vals <- which(pointXyears.list %in% (veg_data$Point_year))
Cov[ind.vals, -c(1:4)] <- veg_data %>%
  filter(Point_year %in% pointXyears.list) %>%
  arrange(Point_year) %>%
  select(-Point_year) %>%
  as.matrix

rm(ind.vals, obs, maxDetPossible, sp, ss, tvec)
save.image("Data_compiled.RData")

## Correlation matrix ##
Cov[, cov.names[-c(1:4, 12:13, 20, 25)]] %>% cor(use = "complete")
