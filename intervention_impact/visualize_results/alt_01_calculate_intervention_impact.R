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

analysis_subdir <- "20210511_itn_burnin"
analysis_metric <- "inc"
final_day <- 1095
suffix <- ""

main_dir <- file.path(Sys.getenv("HOME"), 
                      "Dropbox (IDM)/Malaria Team Folder/projects/map_intervention_impact/intervention_impact",
                      analysis_subdir)
# out_dir <- file.path(main_dir,"results", "figs")
# dir.create(out_dir, recursive = T, showWarnings = F)

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

final[, x_Temporary_Larval_Habitat:= round(x_Temporary_Larval_Habitat, 4)]
initial[, x_Temporary_Larval_Habitat:= round(x_Temporary_Larval_Habitat, 4)]

setnames(final, paste0("final_", analysis_metric), "final_val")
setnames(initial, paste0("initial_", analysis_metric), "initial_val")

if (analysis_metric=="inc"){
  final[, final_severe_inc:=NULL]
  initial[, initial_severe_inc:= NULL]
}

# for now, ignore initial, compare results to no-ITN scenario
control <- final[ITN_Coverage==0, list(Site_Name, Run_Number, x_Temporary_Larval_Habitat, day, control_val=final_val)]
if (nrow(control)==0){
  backup_control_dir <- file.path(main_dir, "../20210511_itn_burnin/results/raw/", analysis_metric)
  fnames <- list.files(backup_control_dir)
  control <- fread(file.path(backup_control_dir, fnames[fnames %like% "Int"]))
  control <- control[ITN_Coverage==0, list(Site_Name, Run_Number, x_Temporary_Larval_Habitat=round(x_Temporary_Larval_Habitat, 4), day, control_val=final_val)]
}

int <- final[ITN_Coverage>0]



# compare to burnin, not no-intervention

initial_lastyear <- initial[day==max(day), list(Site_Name, Run_Number, x_Temporary_Larval_Habitat, control_val=initial_val)]
compare <- merge(control, int)



# ggplot(initial[Site_Name %in% c(3,4) & day>16790], aes(x=day, y=initial_val, color=factor(x_Temporary_Larval_Habitat))) +
#   geom_point() +
#   geom_smooth(color="black") +
#   facet_grid(Site_Name ~ x_Temporary_Larval_Habitat)+
#   theme(legend.position = "none")
# 
# 
# ggplot(initial[Site_Name %in% c(3,4) & day>10000 ], aes(x=day, y=initial_val, color=factor(x_Temporary_Larval_Habitat))) +
#   geom_point() +
#   geom_smooth(color="black") +
#   facet_grid(Site_Name ~ x_Temporary_Larval_Habitat)+
#   theme(legend.position = "none")

# ggplot(initial[Site_Name<11 & Site_Name>1 & day>12409], aes(x=day, y=initial_val, color=factor(x_Temporary_Larval_Habitat))) +
#   geom_point() +
#   geom_smooth() +
#   facet_grid(Site_Name~ x_Temporary_Larval_Habitat)+
#   theme(legend.position = "none")

# ggplot(control[Site_Name<5 & Site_Name>2], aes(x=day, y=control_val, color=factor(x_Temporary_Larval_Habitat))) +
#   geom_point() +
#   geom_smooth(color="black") +
#   facet_grid(Site_Name~ x_Temporary_Larval_Habitat)+
#   theme(legend.position = "none")

# assign intervention ids 
ints <- unique(compare[, list(ITN_Coverage, ITN_Retention_Halflife, ITN_Blocking_Halflife, ITN_Initial_Kill, ITN_Initial_Block, ITN_Start)])[order(ITN_Coverage, ITN_Retention_Halflife, ITN_Blocking_Halflife, ITN_Initial_Kill, ITN_Initial_Block, ITN_Start)]
ints[, int_id:= as.integer(rownames(ints))]
compare <- merge(compare, ints, by=c("ITN_Coverage", "ITN_Retention_Halflife", "ITN_Blocking_Halflife", "ITN_Initial_Kill", "ITN_Initial_Block", "ITN_Start"))

compare_summary <- compare[, lapply(.SD, mean), by = .(Site_Name, day, int_id, ITN_Coverage, ITN_Retention_Halflife, ITN_Initial_Block, ITN_Initial_Kill, x_Temporary_Larval_Habitat), 
                           .SDcols = c("control_val", "final_val")]


compare[, final_val_mean:= mean(final_val), by=.(Site_Name, day, int_id, ITN_Coverage, ITN_Retention_Halflife, ITN_Initial_Block, ITN_Initial_Kill, x_Temporary_Larval_Habitat)]

day_to_plot <- unique(compare$day)  # 365
kill_to_plot <- c(0.6)
ret_to_plot <- unique(compare$ITN_Retention_Halflife) # c(730)
block_to_plot <- c(0.6)
cov_to_plot <- unique(compare$ITN_Coverage)
# sites_to_plot <-  c(4)
# x_temps_to_plot <- c(0.1, 0.5012, 5.0119)
sites_to_plot <-  c(4)
x_temps_to_plot <- c(0.1, 0.5012, 5.0119)

subset <- compare[day %in% day_to_plot & Site_Name %in% sites_to_plot & ITN_Initial_Kill %in% kill_to_plot & ITN_Retention_Halflife %in% ret_to_plot & ITN_Initial_Block %in% block_to_plot & ITN_Coverage %in% cov_to_plot & x_Temporary_Larval_Habitat %in% x_temps_to_plot]
subset_summary <- compare_summary[day %in% day_to_plot & Site_Name %in% sites_to_plot & ITN_Initial_Kill %in% kill_to_plot & ITN_Retention_Halflife %in% ret_to_plot & ITN_Initial_Block %in% block_to_plot & ITN_Coverage %in% cov_to_plot & x_Temporary_Larval_Habitat %in% x_temps_to_plot]

subset[, transmission:= factor(x_Temporary_Larval_Habitat, levels=x_temps_to_plot, labels=c("Low", "Medium", "High"))]
subset_summary[, transmission:= factor(x_Temporary_Larval_Habitat, levels=x_temps_to_plot, labels=c("Low", "Medium", "High"))]

subset_summary[, change:= final_val-control_val]
subset_summary[, reduction:= -change]
subset_summary[, pct_change:= 100*(final_val-control_val)/control_val]
subset_summary[, pct_reduction:= -pct_change]

metric_label <- ifelse(analysis_metric=="inc", "Incidence", "Prevalence")

subset[, ITN_Retention_Halflife:=factor(ITN_Retention_Halflife, labels=c("6 Months", "1 Year", "2 Years", "3 Years"))]
subset_summary[, ITN_Retention_Halflife:=factor(ITN_Retention_Halflife, labels=c("6 Months", "1 Year", "2 Years", "3 Years"))]

subset[, year:=day/365]
subset_summary[, year:=day/365]

# mean metric values in control sims
subset[, mean(control_val), by="transmission"]

pdf(file.path(out_dir, "time_lines.pdf"), height=6, width=8)

ggplot(subset_summary, aes(x=year, color=ITN_Retention_Halflife)) +
  geom_line(aes(y=final_val)) +
  geom_line(aes(y=control_val), color="black") +
  theme_minimal() + 
  geom_point(data=subset, aes(y=final_val), size=0.5) + 
  geom_point(data=subset, aes(y=control_val), color="black", size=0.5) +
  facet_grid(ITN_Coverage~transmission) +
  theme(text=element_text(size=6)) +
  labs(x="Years since distribution",
       y=metric_label,
       color="Mean Retention\nTime",
       title=paste(metric_label, "over Time for No-Intervention Scenario (Black) vs Intervention Scenarios"))

graphics.off()

subset_summary[, year_label:=paste(year, "Years Since Dist.")]
subset_summary[, year_label_short:=paste0("Yr", year)]

pdf(file.path(out_dir, "bars_bytime.pdf"), height=6, width=8)
    ggplot(subset_summary, aes(x=ITN_Retention_Halflife, y=pct_reduction, group=year_label_short, fill=ITN_Retention_Halflife)) +
      geom_bar(stat="identity", position = "dodge", color="black") +
      geom_text(aes(label=year_label_short, y=3), position = position_dodge(width=1), size=1.5) +
      scale_fill_brewer(type="seq", palette = "YlGnBu", name="Mean Retention\nTime") + 
      theme_minimal() + 
      theme(text=element_text(size=6)) +
      facet_grid(ITN_Coverage ~ transmission) +
      labs(x="Mean Retention Time",
           y=paste("% Reduction in", metric_label),
           title=paste("Reduction in", metric_label, "by Transmission Intensity (Columns), Mean Retention Time,\nYears Since Distribution (Grouped Bars), and ITN Coverage (Rows)"))
graphics.off()

pdf(file.path(out_dir, "bars_byret.pdf"), height=6, width=8)
ggplot(subset_summary, aes(x=transmission, y=pct_reduction, fill=ITN_Retention_Halflife)) +
  geom_bar(stat="identity", position = "dodge", color="black") +
  scale_fill_brewer(type="seq", palette = "YlGnBu", name="Mean Retention\nTime") + 
  theme_minimal() + 
  theme(text=element_text(size=6)) +
  facet_grid(ITN_Coverage ~ year_label) +
  labs(x="Transmission Level",
       y=paste("% Reduction in", metric_label),
       title=paste("Reduction in", metric_label, "by Transmission Intensity, Mean Retention Time,\nYears Since Distribution (Columns), and ITN Coverage (Rows)"))
graphics.off()

write.csv(control, file=file.path(out_dir, "../clean/control_results.csv"), row.names = F)

# Exponential decay curves
find_decay <- function(init, mean_decay, time){
  return(init*exp(-(1/mean_decay)*time))
}


kill_mean_dur <- 1460
block_mean_dur <- 730
days <- 0:1095
decay_dt <- data.table(day=days,
                       kill=find_decay(kill_to_plot, kill_mean_dur, days),
                       block=find_decay(block_to_plot, block_mean_dur, days))

decay_dt <- melt(decay_dt, id.vars = "day")

pdf(file.path(out_dir, "killblock_decay.pdf"), height=3, width=5)
ggplot(decay_dt, aes(x=day/365, linetype=variable, y=value)) +
  ylim(0,1) +
  geom_line() +
  theme_classic() +
  theme(legend.position = "none",
        text=element_text(size=6)) +
  labs(x="Year",
       y="Effectiveness")
graphics.off()

retention_decay <- data.table(expand.grid(ret_to_plot, days))
names(retention_decay) <- c("mean_retention", "day")
retention_decay[, decay:= find_decay(1, mean_retention, day)]

pdf(file.path(out_dir, "retention_decay.pdf"), height=3, width=5)
ggplot(retention_decay, aes(x=day/365, color=factor(mean_retention), y=decay)) +
  ylim(0,1) +
  geom_line() +
  theme_classic() +
  theme(legend.position = "none",
        text=element_text(size=6)) +
  labs(x="Year",
       y="Prop. of Nets Retained")
graphics.off()

# smoothing for line plots
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
      subset[, smooth_mean:= get_smooth(control_val, final_val)]
      return(subset)
    })
    sub_sub_list <- rbindlist(sub_sub_list)
  })
  sub_list <- rbindlist(sub_list)
})
smooth_compare <- rbindlist(smooth_compare)




# line_plots <- ggplot(subset_summary, aes(x=control_val, color=factor(ITN_Retention_Halflife))) + 
#                   geom_vline(data=unique(subset_summary[x_Temporary_Larval_Habitat %in% x_temps_to_plot & Site_Name %in% sites_to_plot,
#                                                         list(Site_Name, day, x_Temporary_Larval_Habitat, control_val)]),
#                              aes(xintercept=control_val)) +
#                   geom_abline() + 
#                   geom_point(aes(y=final_val), alpha=0.5) +
#                   scale_color_brewer(type="seq", palette = "YlGnBu", name="Retention\nHalf-Life") + 
#                   # geom_smooth(se=F, span=0.3) + 
#                   geom_line(aes(y=smooth_mean))+ 
#                   facet_grid(ITN_Coverage ~ day)
# 

# ggplot(control[x_Temporary_Larval_Habitat<10], aes(x=x_Temporary_Larval_Habitat, y=control_val)) + 
#   geom_point() +
#   facet_wrap(~Site_Name)









