###############################################################################################################
## alt_01_calculate_intervention_impact.r
## Amelia Bertozzi-Villa
## December 2019
## 
## 
##############################################################################################################


library(data.table)
library(ggplot2)

rm(list=ls())

analysis_subdir <- "20210325_itn_suite"
analysis_metric <- "inc"
final_day <- 1095
suffix <- ""

main_dir <- file.path(Sys.getenv("HOME"), 
                      "Dropbox (IDM)/Malaria Team Folder/projects/map_intervention_impact/intervention_impact",
                      analysis_subdir)
out_dir <- file.path(main_dir,"results", "clean")
dir.create(out_dir, recursive = T, showWarnings = F)

# read in data
in_dir <- file.path(main_dir, "results", "raw", analysis_metric)
fnames <- list.files(in_dir)
initial <- fread(file.path(in_dir, fnames[fnames %like% "Burnin"]))
final <- fread(file.path(in_dir, fnames[fnames %like% "Int"]))

# in some cases, the final year will be duplicated, with "0" in the prevalence column. Flag and remove these.
# (should only happen for the older prevalence datasets)
remove_duplicates <- function(this_dt){
  dt_names <- names(this_dt)[!names(this_dt) %like% "prev"]
  this_dt[, count:=seq_len(.N), by=dt_names]
  this_dt <- this_dt[count==1]
  this_dt[, count:=NULL]
  return(this_dt)
}

if (analysis_metric=="prev"){
  initial <- remove_duplicates(initial)
  final <- remove_duplicates(final)
}


# for now, ignore initial, compare results to no-ITN scenario
control <- final[ITN_Coverage==0, list(Site_Name, Run_Number, x_Temporary_Larval_Habitat, day, control_inc=final_inc, control_severe_inc=final_severe_inc)]
int <- final[ITN_Coverage>0]
compare <- merge(control, int)

compare[, x_Temporary_Larval_Habitat:= round(x_Temporary_Larval_Habitat, 4)]

compare_summary <- compare[, lapply(.SD, mean), by = .(Site_Name, day, ITN_Coverage, ITN_Retention_Halflife, ITN_Blocking_Halflife, ITN_Initial_Kill, x_Temporary_Larval_Habitat), 
                           .SDcols = c("control_inc", "control_severe_inc", "final_inc", "final_severe_inc")]

day_to_plot <- 365
sites_to_plot <-  c(3, 4, 10)
kill_to_plot <- c(0.2)
ret_to_plot <- unique(compare$ITN_Retention_Halflife) # c(730)
block_to_plot <- c(730)

subset <- compare[day==day_to_plot & Site_Name %in% sites_to_plot & ITN_Initial_Kill %in% kill_to_plot & ITN_Retention_Halflife %in% ret_to_plot & ITN_Blocking_Halflife %in% block_to_plot]
subset_summary <- compare_summary[day==day_to_plot & Site_Name %in% sites_to_plot & ITN_Initial_Kill %in% kill_to_plot & ITN_Retention_Halflife %in% ret_to_plot & ITN_Blocking_Halflife %in% block_to_plot] 


ggplot(subset, aes(x=control_inc, y=final_inc, color=factor(ITN_Retention_Halflife))) + 
  geom_abline() + 
  geom_point(alpha=0.25) +
  geom_smooth(se=F, span=0.3) + 
  facet_grid(ITN_Coverage ~ Site_Name)


ggplot(control[x_Temporary_Larval_Habitat<10], aes(x=x_Temporary_Larval_Habitat, y=control_inc)) + 
  geom_point() +
  facet_wrap(~Site_Name)

x_temps_to_plot <- c(0.1995, 0.5012, 5.0119)
subset_summary[, pct_change:= 100*(final_inc-control_inc)/control_inc]
subset_summary <- subset_summary[x_Temporary_Larval_Habitat %in% x_temps_to_plot]
subset_summary[, transmission:= factor(x_Temporary_Larval_Habitat, levels=x_temps_to_plot, labels=c("Low", "Medium", "High"))]


ggplot(subset_summary, aes(x=transmission, y=pct_change, fill=factor(ITN_Retention_Halflife))) +
  geom_bar(stat="identity", position = "dodge", alpha=0.75) +
  facet_grid(ITN_Coverage ~ Site_Name)








