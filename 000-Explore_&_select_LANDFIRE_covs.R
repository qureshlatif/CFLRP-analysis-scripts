library(stringr)
library(dplyr)

#setwd("/home/RMBO.LOCAL/quresh.latif/CFLRP")
setwd("C:/Users/Quresh.Latif/files/projects/FS/CFLRP")


Patch_structure <- read.csv("LANDFIRE_vars/Gap_2km.csv", header = T, stringsAsFactors = F) %>%
  tbl_df

# Explore correlation matrix and dump correlated variables #
#cor(Patch_structure[, c("total.edge", "mean.patch.area", "mean.perim.area.ratio", "total.core.area", "mean.patch.core.area", "patch.cohesion.index")], use = "complete")
# Dumping core area - correlated with mean patch area (r = 0.84)
# Patch cohesion index correlated with log area (r = 0.85)
#cor(Patch_structure[, c("mean.patch.area", "mean.perim.area.ratio")], use = "complete")

Patch_structure <- Patch_structure %>%
  rename(Year = Year1,
         mnPtchAr_Gap2km = mean.patch.area,
         mnPerArRatio_Gap2km = mean.perim.area.ratio) %>%
  select(TransNum, Year, mnPtchAr_Gap2km, mnPerArRatio_Gap2km) %>%
  full_join(read.csv("LANDFIRE_vars/Open_2km.csv", header = T, stringsAsFactors = F) %>%
              tbl_df %>%
              rename(Year = Year1,
                     mnPtchAr_Opn2km = mean.patch.area,
                     mnPerArRatio_Opn2km = mean.perim.area.ratio) %>%
              select(TransNum, Year, mnPtchAr_Opn2km, mnPerArRatio_Opn2km),
            by = c("TransNum", "Year")) %>%
  full_join(read.csv("LANDFIRE_vars/Gap_3km.csv", header = T, stringsAsFactors = F) %>%
              tbl_df %>%
              rename(Year = Year1,
                     mnPtchAr_Gap3km = mean.patch.area,
                     mnPerArRatio_Gap3km = mean.perim.area.ratio) %>%
              select(TransNum, Year, mnPtchAr_Gap3km, mnPerArRatio_Gap3km),
            by = c("TransNum", "Year")) %>%
  full_join(read.csv("LANDFIRE_vars/Open_3km.csv", header = T, stringsAsFactors = F) %>%
              tbl_df %>%
              rename(Year = Year1,
                     mnPtchAr_Opn3km = mean.patch.area,
                     mnPerArRatio_Opn3km = mean.perim.area.ratio) %>%
              select(TransNum, Year, mnPtchAr_Opn3km, mnPerArRatio_Opn3km),
            by = c("TransNum", "Year"))

#cor(Patch_structure[, names(Patch_structure)[-c(1:2)]], use = "complete")
#plot(Patch_structure$mnPerArRatio_Gap2km, Patch_structure$mnPerArRatio_Gap3km)
#plot(Patch_structure$mnPerArRatio_Opn2km, Patch_structure$mnPerArRatio_Opn3km)

Patch_structure <- Patch_structure %>%
#  select(-mnPtchAr_Gap3km) %>%
#  select(-mnPtchAr_Opn3km) %>%
  full_join(read.csv("LANDFIRE_vars/Distance_Gap_2km.csv", header = T, stringsAsFactors = F) %>%
              tbl_df %>%
              rename(NNdist_Gap2km = observed.mean.distance) %>%
              select(TransNum, Year, NNdist_Gap2km),
            by = c("TransNum", "Year")) %>%
  full_join(read.csv("LANDFIRE_vars/Distance_Open_2km.csv", header = T, stringsAsFactors = F) %>%
              tbl_df %>%
              rename(NNdist_Opn2km = observed.mean.distance) %>%
              select(TransNum, Year, NNdist_Opn2km),
            by = c("TransNum", "Year")) %>%
  full_join(read.csv("LANDFIRE_vars/Distance_Gap_3km.csv", header = T, stringsAsFactors = F) %>%
              tbl_df %>%
              rename(NNdist_Gap3km = observed.mean.distance) %>%
              select(TransNum, Year, NNdist_Gap3km),
            by = c("TransNum", "Year")) %>%
  full_join(read.csv("LANDFIRE_vars/Distance_Open_3km.csv", header = T, stringsAsFactors = F) %>%
              tbl_df %>%
              rename(NNdist_Opn3km = observed.mean.distance) %>%
              select(TransNum, Year, NNdist_Opn3km),
            by = c("TransNum", "Year"))

#cor(Patch_structure %>% select(starts_with("NN")), use = "complete")
#plot(Patch_structure$NNdist_Gap2km, Patch_structure$NNInd_Gap2km)
#View(Patch_structure %>% filter(NNInd_Gap2km == max(NNInd_Gap2km, na.rm = T)))

#Patch_structure <- Patch_structure %>%
#  select(-one_of(c("NNInd_Gap2km", "NNInd_Opn2km", "NNdist_Gap3km",
#                   "NNInd_Gap3km", "NNdist_Opn3km", "NNInd_Opn3km")))

#ind.2km <- which(str_detect(names(Patch_structure), "2km"))
#ind.3km <- which(str_detect(names(Patch_structure), "3km"))

write.csv(Patch_structure, "LANDFIRE_vars/Patch_struct_&_config.csv", row.names = F)

# Tabulate SDs #
rows <- names(Patch_structure)[-c(1:2)] %>% str_sub(1, -4) %>% unique
cols <- c("Nsquare_3km", "Ncircle_2km")
out <- matrix(NA, nrow = length(rows), ncol = length(cols),
              dimnames = list(rows, cols))

for(i in 1:length(rows)) {
  out[i, "Ncircle_2km"] <- sd(Patch_structure %>% pull(str_c(rows[i], "2km")), na.rm = T)
  out[i, "Nsquare_3km"] <- sd(Patch_structure %>% pull(str_c(rows[i], "3km")), na.rm = T)
}

out <- out %>%
  as.data.frame %>%
  tbl_df %>%
  mutate(Var = row.names(out)) %>%
  mutate(Ratio_2kmTo3km = Ncircle_2km / Nsquare_3km) %>%
  select(Var, Ncircle_2km, Nsquare_3km, Ratio_2kmTo3km)

write.csv(out, "LANDFIRE_vars/SDs_by_scale.csv", row.names = F)
