
## -----------------------------------------------------------------------------------------------------------------
# assess_cm_burnin
# 
# Amelia Bertozzi-Villa, Institute for Disease Modeling
# March 2022 
# 
# analyze the impact of running CM-only sims from burnin, followed by added ITNs
## -----------------------------------------------------------------------------------------------------------------------


library(ggplot2)
library(data.table)

analysis_subdir <- "20220310_cm_burnin_test"
analysis_metric <- "inc"

main_dir <- file.path(Sys.getenv("HOME"), 
                      "Dropbox (IDM)/Malaria Team Folder/projects/map_intervention_impact/intervention_impact",
                      analysis_subdir, "results/raw", analysis_metric)


burnin <- fread(file.path(main_dir, "MAP_20210511_itn_burnin_Burnin_incidence.csv"))
int <- fread(file.path(main_dir, "20220310_cm_burnin_test_Int_TEST_incidence.csv"))
int[, final_inc_mean:= mean(final_inc), by=.(Site_Name, day, CM_Coverage, ITN_Start, ITN_Coverage, ITN_Retention_Halflife, ITN_Initial_Block, ITN_Initial_Kill, x_Temporary_Larval_Habitat)]



ggplot(int[ITN_Start==3650], aes(x=day, y=final_inc_mean, color=factor(CM_Coverage))) +
  geom_line(aes(linetype=factor(ITN_Start))) +
  geom_line(aes(y=final_inc, group=interaction(Run_Number, CM_Coverage)), alpha=0.5) +
  # scale_linetype_manual(values= c("dashed", "solid")) +
  geom_vline(aes(xintercept=ITN_Start, linetype=factor(ITN_Start))) +
  facet_wrap(. ~ Site_Name)
