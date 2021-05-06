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

analysis_subdir <- "20210331_itn_suite2"
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

control[, x_Temporary_Larval_Habitat:= round(x_Temporary_Larval_Habitat, 4)]
ggplot(control[Site_Name<11 & Site_Name>1], aes(x=day, y=control_inc, color=factor(x_Temporary_Larval_Habitat))) +
      geom_point() + 
      geom_smooth() + 
      facet_grid(Site_Name~ x_Temporary_Larval_Habitat)+
      theme(legend.position = "none")

compare[, x_Temporary_Larval_Habitat:= round(x_Temporary_Larval_Habitat, 4)]

# assign intervention ids 
ints <- unique(compare[, list(ITN_Coverage, ITN_Retention_Halflife, ITN_Blocking_Halflife, ITN_Initial_Kill, ITN_Initial_Block, ITN_Start)])[order(ITN_Coverage, ITN_Retention_Halflife, ITN_Blocking_Halflife, ITN_Initial_Kill, ITN_Initial_Block, ITN_Start)]
ints[, int_id:= as.integer(rownames(ints))]
compare <- merge(compare, ints, by=c("ITN_Coverage", "ITN_Retention_Halflife", "ITN_Blocking_Halflife", "ITN_Initial_Kill", "ITN_Initial_Block", "ITN_Start"))

compare_summary <- compare[, lapply(.SD, mean), by = .(Site_Name, day, int_id, ITN_Coverage, ITN_Retention_Halflife, ITN_Initial_Block, ITN_Initial_Kill, x_Temporary_Larval_Habitat), 
                           .SDcols = c("control_inc", "control_severe_inc", "final_inc", "final_severe_inc")]


get_smooth <- function(x, y){
  if (max(y)<0.05){
    return(y)
  }else{
    lo <- loess(y[y>0]~x[y>0])
    predictions <- c(y[y==0], predict(lo))
    return(pmax(predictions, rep(0, length(predictions))))
  }
}

smooth_compare <- lapply(unique(compare_summary$Site_Name), function(site_name){
  sub_list <- lapply(unique(compare_summary$int_id), function(int_name){
    sub_sub_list <- lapply(unique(compare_summary$day), function(this_day){
      subset <- compare_summary[Site_Name==site_name & int_id==int_name & day==this_day]
      # subset[, smooth_min:= get_smooth(mean_initial, min_final)]
      # subset[, smooth_max:= get_smooth(mean_initial, max_final)]
      subset[, smooth_mean:= get_smooth(control_inc, final_inc)]
      return(subset)
    })
    sub_sub_list <- rbindlist(sub_sub_list)
  })
  sub_list <- rbindlist(sub_list)
})
smooth_compare <- rbindlist(smooth_compare)

day_to_plot <- unique(compare$day)  # 365
sites_to_plot <-  c(4)
kill_to_plot <- c(0.4)
ret_to_plot <- unique(compare$ITN_Retention_Halflife) # c(730)
block_to_plot <- c(0.3)
cov_to_plot <- unique(compare$ITN_Coverage)
x_temps_to_plot <- c(0.1259, 0.5012, 5.0119)

subset <- compare[day %in% day_to_plot & Site_Name %in% sites_to_plot & ITN_Initial_Kill %in% kill_to_plot & ITN_Retention_Halflife %in% ret_to_plot & ITN_Initial_Block %in% block_to_plot & ITN_Coverage %in% cov_to_plot]
subset_summary <- smooth_compare[day %in% day_to_plot & Site_Name %in% sites_to_plot & ITN_Initial_Kill %in% kill_to_plot & ITN_Retention_Halflife %in% ret_to_plot & ITN_Initial_Block %in% block_to_plot & ITN_Coverage %in% cov_to_plot]

subset_summary[, pct_change:= 100*(smooth_mean-control_inc)/control_inc]
subset_summary[, pct_reduction:= -pct_change]

subset_for_bars <- subset_summary[x_Temporary_Larval_Habitat %in% x_temps_to_plot]
subset_for_bars[, transmission:= factor(x_Temporary_Larval_Habitat, levels=x_temps_to_plot, labels=c("Low", "Medium", "High"))]


ggplot(subset_for_bars, aes(x=transmission, y=pct_reduction, fill=factor(ITN_Retention_Halflife))) +
  geom_bar(stat="identity", position = "dodge", color="black") +
  scale_fill_brewer(type="seq", palette = "YlGnBu", name="Retention\nHalf-Life\n(Days)") + 
  theme_minimal() + 
  facet_grid(ITN_Coverage ~ day) +
  labs(x="Transmission Level",
       y="% Reduction in Incidence",
       title="Reduction in Incidence by Transmission Intensity, ITN Retention Half-Life, Time in Days (columns), and ITN Coverage (rows)")



line_plots <- ggplot(subset_summary, aes(x=control_inc, color=factor(ITN_Retention_Halflife))) + 
                  geom_vline(data=unique(subset_summary[x_Temporary_Larval_Habitat %in% x_temps_to_plot & Site_Name %in% sites_to_plot,
                                                        list(Site_Name, day, x_Temporary_Larval_Habitat, control_inc)]),
                             aes(xintercept=control_inc)) +
                  geom_abline() + 
                  geom_point(aes(y=final_inc), alpha=0.5) +
                  scale_color_brewer(type="seq", palette = "YlGnBu", name="Retention\nHalf-Life") + 
                  # geom_smooth(se=F, span=0.3) + 
                  geom_line(aes(y=smooth_mean))+ 
                  facet_grid(ITN_Coverage ~ day)


# ggplot(control[x_Temporary_Larval_Habitat<10], aes(x=x_Temporary_Larval_Habitat, y=control_inc)) + 
#   geom_point() +
#   facet_wrap(~Site_Name)









