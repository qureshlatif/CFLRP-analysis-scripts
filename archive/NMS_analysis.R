#############################################################################################
#### Script for Multivariate ordination of habitat data
#### Author: Quresh Latif (quresh.latif@birdconservancy.org)
#### Institution: Bird Conservancy of the Rockies
#### Date Created: 11/12/2018
#### Last Modified: 11/12/2018
#############################################################################################
# The goal of this script is to create  Non-metric multi-dimensional scaling ordination for
# grid-level habitat data derived from LANDFIRE (8 variables) and point-level habitat data
# from field measurements (12 variables). This script is adapted from a script written by
# Jeff Cannon for seedling and sapling data.
#############################################################################################

#---> Load libraries and data
require(openxlsx)
require(vegan)
require(dplyr)
require(stringr)
require(scatterplot3d)
require(shape)
require(reshape)

setwd("C:/Users/Quresh.Latif/files/projects/FS/CFLRP")
load("Data_compiled.RData")

#---> Load data (cite sources)
# Lookup tables
grid_LU <- landscape_data %>%
  mutate(Grid_year = str_c(Grid, "_", Year)) %>%
  select(Grid_year, YST_1kmNB, PctTrt_1kmNB)
point_LU <- Cov %>% as.data.frame %>%
  mutate(Point_year = dimnames(Cov)[[1]]) %>%
  select(Point_year, Trt_stat, Trt_time)

# Seedling and sapling species matrices
grid_vars <- landscape_data %>%
  mutate(Grid_year = str_c(Grid, "_", Year)) %>%
  select(Grid_year, PACC10_3km, PACC40_3km, mnPtchAr_Gap3km:mnPerArRatio_Opn3km, NNdist_Gap3km, NNdist_Opn3km)
point_vars <- Cov %>% as.data.frame %>%
  mutate(Point_year = dimnames(Cov)[[1]]) %>%
  select(Point_year, CanCov:NumSnags, RCOV_PP:RCOV_AS, ShrubVol, RSCV_Ladder:RSCV_AS, HerbGrassVol)

#---> Input variables (Name all necessary input parameters to be soft-coded)
run_stress_plot = TRUE # create MDS stress plots?
run_MDS = TRUE         # run MDS 
rand_starts = 20     # Random starts for NMS (best to be 100-1000)
n_corr_vectors = 5     # number of corr. vectors for NMS plots
colors = c('#3182bd', '#9ecae1')
pchs = c(21, 25)

#---> Output variables (Name all filenames for outputs e.g., pdf or csv filenames/locations)
# Figure outputs
MDS_stress_output = 'DATA/OUTPUT/mds_stress_log.txt'
MDS_stress = 'DATA/OUTPUT/FIGS/mds_stress.pdf'
NMS_sap_fig   = 'DATA/OUTPUT/FIGS/NMS_sap_fig.pdf'
NMS_seed_fig   = 'DATA/OUTPUT/FIGS/NMS_seed_fig.pdf'

# Tabular outputs
NMS_sap_tab  = 'DATA/OUTPUT/NMS_saps.csv'
NMS_seed_tab = 'DATA/OUTPUT/NMS_seed.csv'

#---> Load custom functions with explanations (as needed)
#########################
# FUNCTION: clean_up_matrix()
# Convert species matrix from PNWR study to format for vegan::metaMDS
#########################
clean_up_matrix = function(mx){
  mx$samp = with(mx, paste(plot, dir, year, sep = '_'))
  mx$samp = with(mx, paste(plot, dir, year, sep = '_'))
  rownames(mx) = mx$samp
  mx = mx[grep(
    'plot|dir|year|samp|winch|burn|NA',
    colnames(mx),
    invert = TRUE,
    value = TRUE
  )]
  mx$dummy = min(mx[mx > 0])
  for (c in colnames(mx))
    mx[, c] = as.numeric(as.character(mx[, c]))
  mx = as.matrix(mx)
  mx = log(mx + 1)
  return(mx)
}

#########################
# FUNCTION: mds_stress_test()
# Outputs stress for a range of selected k values for vegan::metaMDS
#########################
mds_stress_test = function(k_range = 1:5, mx){
  output = data.frame(k = numeric(), stress = numeric())
  for(k in k_range){
    cat('======Running Stress: k = ', k, '======\n') 
    mds = metaMDS(mx, k = k)
    tmp_stress = data.frame(k = k, stress = mds$stress)
    output = rbind(output, tmp_stress)
  }
  cat('Stress test complete.\n')
  return(output)
}

#########################
# FUNCTION: get_mds_output()
# Runs vegan::metaMDS and cleans up data for plotting
#########################
get_mds_output = function(mx, veg, treat_LU, try = 20){
  mds = metaMDS(mx, k = 3, distance = 'bray', try = try)
  print(mds)
  stressplot(mds)
  out = veg[, c('plot', 'dir', 'year')]
  out$sub = with(out, paste(plot, dir, year, sep = '_'))
  out_mds = as.data.frame(mds$points)
  out_mds$sub = rownames(out_mds)
  out = merge(out, out_mds)
  out = merge(treat_LU, out)
  return(out)
}

#########################
# FUNCTION: get_mds_means()
# Calculates mean MDS scores from vegan::metaMDS by treatment type
#########################
get_mds_means = function(mds){
  mds = ddply(mds, .(winch, burn, year), summarize,
              MDS1 = mean(MDS1),
              MDS2 = mean(MDS2),
              MDS3 = mean(MDS3))
}
#########################
# FUNCTION: MDS_3dscatter()
# Plots 3d scatter of output from get_MDS_means() using default settings
#########################
MDS_3dscatter = function(mds, main, cols, pchs) {
  mds$col = with(mds, ifelse(burn == 'control', cols[1], cols[2]))
  mds$pch = with(mds, ifelse(winch == 'control', pchs[1], pchs[2]))
  par(fig = c(0, 1, 0, 1))
  s = with(mds,
           scatterplot3d(
             x = MDS1, y = MDS2, z = MDS3,
             xlab = 'NMS1', ylab = 'NMS2', zlab = 'NMS3',
             main = main,
             mar = c(3,3,1,2),
             bg = mds$col, pch = pch, type = 'h', cex.symbols = 2))
  return(s)
}

#########################
# FUNCTION: MDS_add_arrows()
# Adds arrows to time series data from get_MDS_means() to 3dscatter
#########################
MDS_add_arrows = function(mds, s, treats = 4, years = 4, arr.length = 0.15){
  for(i in 0:(treats-1)) {
    for(j in 1:(years-1)){
      row = ((years) * i) + j
      p1 = with(mds, s$xyz.convert(MDS1[row], MDS2[row], MDS3[row]))
      p2 = with(mds, s$xyz.convert(MDS1[row + 1], MDS2[row + 1], MDS3[row +1]))
      Arrows(p1$x,p1$y,p2$x,p2$y,arr.length = arr.length,arr.adj = 1)
    }
  }
}
#########################
# FUNCTION: add_mds_legend()
# Default legend settings
#########################
add_mds_legend = function(location = 'topleft', cols, pchs) {
  legend(location,c('Intact unburned','Intact burned','Winched unburned','Winched burned'),
         bg = 'white', pch = c(pchs[1], pchs[1], pchs[2], pchs[2]),
         pt.bg =  c(cols[1], cols[2]),
         bty = 'o', box.col = NA,
         pt.cex = 2)
}
#########################
# FUNCTION: get_spp_vect()
# Calculates correlation of species from metaMDS() output
#########################
get_spp_vect = function(mx, mds, head = 5) {
  vegDat = as.data.frame(mx)
  vegDat$dummy = NULL
  vegDat$sub = rownames(vegDat)
  vegDat = merge(mds, vegDat)
  species = grep(
    'sub|plot|dir|year|winch|burn|MDS',
    colnames(vegDat),
    invert = TRUE,
    value = TRUE
  )
  corr_mat = data.frame(
    species = character(0),
    MDS1 = numeric(),
    MDS2 = numeric(),
    MDS3 = numeric()
  )
  for (spp in species) {
    corr_tmp = vegDat[, c(grep('MDS', colnames(vegDat), value = TRUE), spp)]
    corr_tmp = cor(corr_tmp)[spp, grep('MDS', colnames(corr_tmp), value = TRUE)]
    corr_tmp = as.data.frame(t(corr_tmp))
    corr_tmp = cbind(data.frame(species = spp), corr_tmp)
    corr_mat = rbind(corr_mat, corr_tmp)
  }
  corr_mat$SS = with(corr_mat, sqrt(MDS1 ^ 2 + MDS2 ^ 2 + MDS3 ^ 2))
  corr_mat = corr_mat[order(-corr_mat$SS), ]
  corr_mat = corr_mat[1:head, ]
  corr_mat$SS <- NULL
  return(corr_mat)
}

#########################
# FUNCTION: plot_corr_vects()
# Creates a 3d inset of correlation vectors from get_spp_vects()
#########################
plot_corr_vects = function(corr_vect, quad, length.arr = 0.15, adj = 1.25, cex = 0.8) {
  par(fig = quad, new = TRUE)
  q = with(corr_vect, 
           scatterplot3d(
             x = c(0, MDS1),
             y = c(0, MDS2),
             z = c(0, MDS3),
             box = FALSE, grid = FALSE,
             axis = FALSE, pch = NA,
             color = 'grey', type = 'h'))
  for (i in 1:nrow(corr_vect)) {
    p1 = q$xyz.convert(0, 0, 0)
    p2 = with(corr_vect, q$xyz.convert(MDS1[i], MDS2[i], MDS3[i]))
    Arrows(p1$x,p1$y,p2$x,p2$y,arr.length = length.arr,arr.adj = 1)
    text(p2$x * adj, p2$y * adj, label = corr_vect[i, 'species'], cex = cex)
  }
}

########################################END SET UP###########################################

######################################START ANALYSIS#########################################
#---> Clean up sapling and seedling data
saps_mx = clean_up_matrix(saps_veg)
seed_mx = clean_up_matrix(seeds_veg)

#---> Run stress plot to determine optimal k value for MDS
set.seed(6000) # set randomly generated seed for repeateable results
if(run_stress_plot == TRUE) {
  mds_stress_saps = mds_stress_test(k_range = 1:5, mx = saps_mx)
  mds_stress_seed = mds_stress_test(k_range = 1:5, mx = seed_mx)
}

#---> Run MDS test with 10,000 starts, k = 3, bray-curtis distances
if(run_MDS == TRUE){
  sink(MDS_stress_output)
  cat('************ SAPLING NMS ************\n\n')
  mds_sap = get_mds_output(saps_mx, saps_veg, treat_LU, try = rand_starts)
  cat('\n\n************ SEEDLING NMS ************\n\n')
  mds_seed = get_mds_output(seed_mx, seeds_veg, treat_LU, try = rand_starts)
  sink()
}
save.image('DATA/OUTPUT/MDS_progress.RData')

#---> Calculate mean scores for each treatment level x year
mds_sap_means = get_mds_means(mds_sap)
mds_seed_means = get_mds_means(mds_seed)

goodness(mds_sap)

#######################################END ANALYSIS##########################################

######################################START OUTPUTS##########################################
#---> Output seedling and sapling NMS scores
write.csv(mds_sap, file = NMS_sap_tab, row.names = FALSE)
write.csv(mds_seed, file = NMS_seed_tab, row.names = FALSE)

#######################################END OUTPUTS###########################################

#####################################START GRAPHICS##########################################
#---> Sapling NMS figure
for(i in 1:0){
  if(i) {pdf(NMS_sap_fig, width = 8, height = 6)} else {
    png(gsub('\\.pdf','.\\png',NMS_sap_fig,ignore.case=TRUE), width =8, height=6, units='in', res=300)}
  par(mfrow = 1:2)
  s = MDS_3dscatter(mds_sap_means, main = 'Sapling compositional change (NMS)', cols = colors, pchs = pchs)
  MDS_add_arrows(mds_sap_means, s)
  add_mds_legend('topleft', colors, pchs)
  corr_vect = get_spp_vect(saps_mx, mds_sap, head = n_corr_vectors)
  quad = c(0.1, 0.5, 0, 0.5)
  plot_corr_vects(corr_vect, quad)
  dev.off()
}

for(i in 1:0){
  if(i) {pdf(NMS_seed_fig, width = 8, height = 6)} else {
    png(gsub('\\.pdf','.\\png',NMS_seed_fig,ignore.case=TRUE), width =8, height=6, units='in', res=300)}
  s = MDS_3dscatter(mds_seed_means, main = 'Seedling composition change (NMS)', cols = colors, pchs = pchs)
  MDS_add_arrows(mds_seed_means, s, treats = 4, years = 3)
  add_mds_legend('topleft', colors, pchs)
  corr_vect = get_spp_vect(seed_mx, mds_seed, head = n_corr_vectors)
  quad = c(0.4, 0.8, 0, 0.5)
  plot_corr_vects(corr_vect, quad)
dev.off()
}

#---> Output stress plots
if(run_stress_plot == TRUE) {
  pdf('DATA/OUTPUT/FIGS/mds_stress.pdf', width = 4, height = 4)
  plot(stress ~ k, data = mds_stress_saps, type = 'b', pch = 20, main = 'Sapling MDS stress')
  plot(stress ~ k, data = mds_stress_seed, type = 'b', pch = 20, main = 'Seedling MDS stress')
  dev.off()
}
#######################################END GRAPHICS##########################################
