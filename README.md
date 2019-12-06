# CFLRP-analysis-scripts

### Basic function of each model script in this repository ###
# Note: This readme does not include descriptions of script files in the "archive" folder.

##________ Scripts for peer-reviewed effectiveness monitoring manuscript________##
# Data and tabulation compilation # 
# 000-Explore_&_select_LANDFIRE_covs.R - Processes and calculates correlations among several potential landscape metrics based on LANDFIRE. After reviewing metrics with this script, a subset was selected for use in analysis.
# 00-CFLRP_point_coords.R - Queries and compiles spatial coordinates for survey points for mapping.
# 00-Covariate_summaries.R - Tabulates summary (descriptive) statistics for covariates.
# 00-Data_wrangle.R - Queries raw data from Bird Conservancy database using API. Compiled data are saved in a 'Data_compiled.RDATA' workspace.
# 00-Explore_treatment_polygons.R - Queries hectares of CFLRP units treated with only mechanical thinning versus also treated with prescribed fire.
# 00-Tab_point_trt_statuses.R - Tabulates sample sizes (number of points and grids) by treatment status (Table 1 of manuscript).

# Analysis implementation and compilation of results #
# 01-Analyze_habitat_relations.R - Implements community occupancy model with habitat covariates ("habitat model" in manuscript).
# 01-Analyze_treatment_relations.R - Implements community occupancy model with treatment covariates ("treatment model" in manuscript).
# 01-Trt_habitat_relations.R - Calculates correlations between treatment and habitat covariates (Appendix S2).
# 02-Analyze_spec_richness_habitat.R - Implements and plots post hoc derivation and analysis of species richness habitat relationships (Figure 5, Table 5).
# 02-Analyze_spec_richness_treatment.R - Implements and plots post hoc derivation and analysis of species richness treatment relationships (Figure 4, Table 5).
# 02-Plot_habitat_effects.R - Plots parameter estimates for species habitat relationships (Appendix S9).
# 02-Plot_spp_trt_relations.R - Plots species occupancy against treatment conditions where statistically supported (Appendix S8).
# 02-Plot_treatment_effects.R - Plots parameter estimates for species treatment relationships (Figure 3).
# 02-Tabulate_detection_habitat_model.R - Tabulates detection parameter estimates from habitat model (Appendix S7: Table S2).
# 02-Tabulate_detection_treatment_model.R - Tabulates detection parameter estimates from treatment model (Appendix S7: Table S1).
# 02-Tabulate_hyper-parameters.R - Tabulates community hyper-parameter estimates (Table 5, Appendix S3).
# 02-Tabulate_parameters.R - Tabulates species parameters from community occupancy models for review.

# JAGS model code #
# model_habitat.jags - Habitat model
# model_treatment.jags - Treatment model

##________ Additional scripts presentations (not relevant to peer-reviewed manuscript)________##
# 02-Plot_spec_richness_trt_presentations.R
# 02-Plots_lightening_talk.R
# 02-Plots_webinar.R
