###############################################################################################################
## alt_01_calculate_intervention_impact.r
## Amelia Bertozzi-Villa
## December 2019
## 
## 
##############################################################################################################

library(grid)
library(rgeos)
library(rgdal)
library(maptools)
library(stringr)
library(raster)
library(rasterVis)
library(data.table)
library(ggplot2)
library(gridExtra)
library(latticeExtra)

rm(list=ls())

analysis_subdir <- "20210520_rerun_blocktime"
analysis_metric <- "inc"
final_day <- 1095
suffix <- ""

main_dir <- file.path(Sys.getenv("HOME"), 
                      "Dropbox (IDM)/Malaria Team Folder/projects/map_intervention_impact/intervention_impact",
                      analysis_subdir)
out_dir <- file.path(main_dir,"results", "figs")

# for mapping
arch_dir <- file.path(main_dir, "../../archetypes/results/v4_era5_bounded_transmission/africa")
africa_shp_dir <- "/Volumes/GoogleDrive/My Drive/itn_cube/input_data/general/shapefiles/Africa.shp"
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

int <- final[ITN_Coverage>0]
compare <- merge(control, int)

# ggplot(control[Site_Name<11 & Site_Name>1], aes(x=day, y=control_val, color=factor(x_Temporary_Larval_Habitat))) +
#   geom_point() +
#   geom_smooth(color="black") +
#   facet_grid(Site_Name~ x_Temporary_Larval_Habitat)+
#   theme(legend.position = "none")

# assign intervention ids 
ints <- unique(compare[, list(ITN_Coverage, ITN_Retention_Halflife, ITN_Blocking_Halflife, ITN_Initial_Kill, ITN_Initial_Block, ITN_Start)])[order(ITN_Coverage, ITN_Retention_Halflife, ITN_Blocking_Halflife, ITN_Initial_Kill, ITN_Initial_Block, ITN_Start)]
ints[, int_id:= as.integer(rownames(ints))]
compare <- merge(compare, ints, by=c("ITN_Coverage", "ITN_Retention_Halflife", "ITN_Blocking_Halflife", "ITN_Initial_Kill", "ITN_Initial_Block", "ITN_Start"))

compare[, change:= final_val-control_val]
compare[, reduction:= -change]
compare[, pct_change:= 100*(final_val-control_val)/control_val]
compare[, pct_reduction:= -pct_change]

compare_summary <- compare[, unlist(recursive=FALSE, lapply(
                                                            .(mean = mean, mid=median, iqr = IQR),
                                                            function(f) lapply(.SD, f)
                                                          )),
                           by = .(Site_Name, day, int_id, ITN_Coverage, ITN_Retention_Halflife, ITN_Initial_Block, ITN_Blocking_Halflife, ITN_Initial_Kill, x_Temporary_Larval_Habitat), 
                           .SDcols = c("pct_reduction")]

compare[, final_val_mean:= mean(final_val), by=.(Site_Name, day, int_id, ITN_Coverage, ITN_Retention_Halflife, ITN_Initial_Block, ITN_Initial_Kill, x_Temporary_Larval_Habitat)]
compare[, control_val_mean:= mean(control_val), by=.(Site_Name, day, int_id, ITN_Coverage, ITN_Retention_Halflife, ITN_Initial_Block, ITN_Initial_Kill, x_Temporary_Larval_Habitat)]

# Find which x_temps correspond to the desired EIRs
eir_dir <- file.path(in_dir, "../eir")
eir_fnames <- list.files(eir_dir)
eirs <- fread(file.path(eir_dir, eir_fnames[eir_fnames %like% "Burnin"]))

setnames(eirs, "initial_eir", "eir")
eirs[, x_Temporary_Larval_Habitat:= round(x_Temporary_Larval_Habitat, 4)]
eirs <- eirs[, lapply(.SD, mean), by = .(Site_Name, day, x_Temporary_Larval_Habitat), .SDcols = c("eir")]
eirs <- eirs[, lapply(.SD, mean), by = .(Site_Name,  x_Temporary_Larval_Habitat), .SDcols = c("eir")]
eirs <- eirs[order(Site_Name, x_Temporary_Larval_Habitat)]
eirs[, eir:= round(eir, 2)]

target_eirs <- c(0.5, 5, 50)

closest_eirs <- data.table(expand.grid(Site_Name=1:12, match_eir=target_eirs))
eirs[, match_eir:=eir]

setkeyv(eirs, c("Site_Name", 'match_eir'))
setkeyv(closest_eirs, c("Site_Name", 'match_eir'))

closest_eirs <- eirs[closest_eirs, roll='nearest']
closest_eirs[, transmission:= factor(match_eir, levels=target_eirs, labels=c("Low", "Medium", "High"))]

for_output_eirs <- dcast.data.table(closest_eirs[transmission!="Low"], Site_Name  ~ transmission, value.var = "eir")
write.csv(for_output_eirs, file.path(out_dir, "output_eirs.csv"), row.names = F)
write.csv(closest_eirs, file.path(out_dir, "output_eirs_and_larval_habs.csv"), row.names = F)

day_to_plot <- unique(compare$day)  # 365
kill_to_plot <- c(0.4)
ret_to_plot <- unique(compare$ITN_Retention_Halflife) # c(730)
block_to_plot <- c(0.3)
cov_to_plot <- c(0.4, 0.6, 0.8) # unique(compare$ITN_Coverage)
sites_to_plot <-  2:10
eirs_for_subsetting <- closest_eirs[Site_Name %in% sites_to_plot, list(Site_Name, x_Temporary_Larval_Habitat, eir, transmission)]

# x_temps_to_plot <- c(0.1, 0.5012, 5.0119)
# subset <- compare[day %in% day_to_plot & Site_Name %in% sites_to_plot & ITN_Initial_Kill %in% kill_to_plot & ITN_Retention_Halflife %in% ret_to_plot & ITN_Initial_Block %in% block_to_plot & ITN_Coverage %in% cov_to_plot & x_Temporary_Larval_Habitat %in% x_temps_to_plot]
# subset_summary <- compare_summary[day %in% day_to_plot & Site_Name %in% sites_to_plot & ITN_Initial_Kill %in% kill_to_plot & ITN_Retention_Halflife %in% ret_to_plot & ITN_Initial_Block %in% block_to_plot & ITN_Coverage %in% cov_to_plot & x_Temporary_Larval_Habitat %in% x_temps_to_plot]
# subset[, transmission:= factor(x_Temporary_Larval_Habitat, levels=x_temps_to_plot, labels=c("Low", "Medium", "High"))]
# subset_summary[, transmission:= factor(x_Temporary_Larval_Habitat, levels=x_temps_to_plot, labels=c("Low", "Medium", "High"))]

subset <- compare[day %in% day_to_plot & Site_Name %in% sites_to_plot & ITN_Initial_Kill %in% kill_to_plot & ITN_Retention_Halflife %in% ret_to_plot & ITN_Initial_Block %in% block_to_plot & ITN_Coverage %in% cov_to_plot]
subset_summary <- compare_summary[day %in% day_to_plot & Site_Name %in% sites_to_plot & ITN_Initial_Kill %in% kill_to_plot & ITN_Retention_Halflife %in% ret_to_plot & ITN_Initial_Block %in% block_to_plot & ITN_Coverage %in% cov_to_plot]

subset <- merge(eirs_for_subsetting, subset, by=c("Site_Name", "x_Temporary_Larval_Habitat"), all.x=T)
subset_summary <- merge(eirs_for_subsetting, subset_summary, by=c("Site_Name", "x_Temporary_Larval_Habitat"), all.x=T)

metric_label <- ifelse(analysis_metric=="inc", "Incidence", "Prevalence")

subset[, ITN_Retention_Halflife:=factor(ITN_Retention_Halflife, labels=c("6 Months", "1 Year", "2 Years", "3 Years"))]
subset_summary[, ITN_Retention_Halflife:=factor(ITN_Retention_Halflife, labels=c("6 Months", "1 Year", "2 Years", "3 Years"))]

subset[, year:=day/365]
subset_summary[, year:=day/365]

# correlate retention and blocking:
subset_summary <- subset_summary[(ITN_Retention_Halflife %in% c("6 Months") & ITN_Blocking_Halflife==180) |
                                 (ITN_Retention_Halflife %in% c("1 Year") & ITN_Blocking_Halflife==365) |
                                 (ITN_Retention_Halflife %in% c("2 Years") & ITN_Blocking_Halflife==730) |
                                 (ITN_Retention_Halflife %in% c("3 Years") & ITN_Blocking_Halflife==1095)
                                   ]
subset <- subset[(ITN_Retention_Halflife %in% c("6 Months") & ITN_Blocking_Halflife==180) |
                   (ITN_Retention_Halflife %in% c("1 Year") & ITN_Blocking_Halflife==365) |
                   (ITN_Retention_Halflife %in% c("2 Years") & ITN_Blocking_Halflife==730) |
                   (ITN_Retention_Halflife %in% c("3 Years") & ITN_Blocking_Halflife==1095)
]


subset_summary[, year_label:=paste(year, "Years Since Dist.")]
subset_summary[, year_label_short:=paste0("Yr", year)]
subset[, year_label:=paste(year, "Years Since Dist.")]
subset[, year_label_short:=paste0("Yr", year)]



subset[, mean(control_val), by="transmission"]


##### cluster maps and time series
nclust <- 10
palette <- c("#98B548", "#00A08A", "#8971B3", "#F2AD00", "#5392C2", "#D71B5A", "#902E57", "#F98400", "#B33539", "#367A40")

# shape
africa_shp <- readOGR(africa_shp_dir)
africa_dt <- data.table(fortify(africa_shp, region = "COUNTRY_ID"))
africa_shp <- gSimplify(africa_shp, tol=0.1, topologyPreserve=TRUE)

# raster
cluster_in_dir <- file.path(arch_dir, "02_kmeans")
cluster_raster <- raster(file.path(cluster_in_dir, paste0("map_", nclust, "_cluster", ".tif")))
cluster_dt <- data.table(rasterToPoints(cluster_raster))
names(cluster_dt) <- c("long", "lat", "value")
cluster_dt[, value:= as.factor(value)]

# time
summary_vals <- fread(file.path(cluster_in_dir,  paste0("summary_", nclust, "_cluster", ".csv")))
time_series <- summary_vals[variable_name=="month"]
time_series[, cluster:=as.factor(cluster)]

# quick bit of work for kedar & aysu coverage question-- keep 1yr retetnion time only, plot impact over year by coverage level

this_site <- 4
# main data subsetting
this_subset <- subset[Site_Name==this_site & ITN_Retention_Halflife=="1 Year"]
this_subset_summary <- subset_summary[Site_Name==this_site & ITN_Retention_Halflife=="1 Year"]

this_subset_summary[, cov_label:= factor(ITN_Coverage, labels=c("40%", "60%", "80%"))]
this_subset_summary[, transmission:=factor(transmission, labels=c("EIR ~0.5", "EIR ~5", "EIR ~50"))]

barplot_time <- ggplot(this_subset_summary[transmission!="EIR ~0.5"], aes(x=cov_label, y=mid.pct_reduction, group=year_label_short, fill=cov_label, ymin=mid.pct_reduction-iqr.pct_reduction, ymax=mid.pct_reduction+iqr.pct_reduction)) +
  geom_bar(stat="identity", position = position_dodge(width=0.9), color="black") +
  geom_errorbar (position=position_dodge(width=0.9), colour="black", size=0.25) +
  geom_text(aes(label=year_label_short, y=2), position = position_dodge(width=1), size=1.5) +
  scale_fill_brewer(type="seq", palette = "YlGnBu", name="ITN Use") + 
  theme_minimal() + 
  theme(text=element_text(size=6),
        legend.position = "bottom") +
  facet_grid(. ~ transmission) +
  labs(x="ITN Use",
       y=paste("% Reduction in", metric_label),
       title=paste("Site", this_site, ": Reduction in", metric_label, "by Transmission Intensity,\nYears Since Distribution (Grouped Bars), and ITN Use.\n Nets have a median retention time of one year."))

# mapping plots 
this_cluster_dt <- copy(cluster_dt)
this_cluster_dt[, value:=factor(value==this_site)]
this_color <-  palette[this_site]

cluster_plot <- ggplot() +
  geom_raster(data = this_cluster_dt, aes(fill = value, y = lat, x = long)) +
  geom_path(data = africa_dt, aes(x = long, y = lat, group = group), color = "black", size = 0.3) + 
  scale_fill_manual(values= c("lightgrey", this_color)) +
  coord_equal(xlim = c(-18, 52), ylim = c(-35, 38)) +
  labs(x = NULL, y = NULL, title = NULL) +
  theme_classic(base_size = 12) +
  theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "in"), legend.position = "none")

timeseries_plot <- ggplot(time_series[cluster==this_site & cov=="precip_era5"], aes(x=as.integer(variable_val), y=median, color=cluster, fill=cluster)) +
  geom_ribbon(aes(ymin=perc_25, ymax=perc_75), alpha=0.5, color=NA) +
  geom_line(size=1) +
  geom_line(aes(y=perc_05), size=0.75, linetype=2) +
  geom_line(aes(y=perc_95), size=0.75, linetype=2) +
  scale_color_manual(values = this_color) +
  scale_fill_manual(values = this_color) + 
  scale_x_continuous(breaks=seq(2,12,2), labels=c("F","A","J","A","O","D"), minor_breaks=seq(1,12,2)) +
  theme_minimal() +
  theme(legend.position="none",
        plot.title = element_text(size=8),
        strip.background = element_blank(),
        strip.text.y = element_blank()) +
  labs(title="Rainfall",
       x="",
       y="")

# save all
main_vp <- viewport(width = 0.7, height = 0.75, x = 0.35, y = 0.5)
map_vp <- viewport(width = 0.25, height = 0.5, x = 0.85, y = 0.3)
timeseries_vp <- viewport(width = 0.3, height = 0.4, x = 0.85, y = 0.75)

pdf("~/Desktop/example_impact.pdf", height=4, width=6)
print(barplot_time, vp=main_vp)
print(cluster_plot, vp=map_vp)
print(timeseries_plot, vp=timeseries_vp)
graphics.off()


## end aysu/kedar subsetting


for (this_site in unique(subset$Site_Name)){
  print(this_site)
  
  # main data subsetting
  this_subset <- subset[Site_Name==this_site]
  this_subset_summary <- subset_summary[Site_Name==this_site]
  
  # main data plots 
  lineplot <- ggplot(this_subset[Site_Name==this_site], aes(x=year, color=ITN_Retention_Halflife)) +
    geom_line(aes(y=final_val_mean)) +
    geom_line(aes(y=control_val_mean), color="black") +
    theme_minimal() + 
    geom_point( aes(y=final_val), size=0.5) + 
    geom_point(aes(y=control_val), color="black", size=0.5) +
    facet_grid(ITN_Coverage~transmission) +
    theme(text=element_text(size=6)) +
    labs(x="Years since distribution",
         y=metric_label,
         color="Median Retention\nTime",
         title=paste(metric_label, "over Time for No-Intervention Scenario (Black) vs Intervention Scenarios"))
  
  barplot_time <- ggplot(this_subset_summary[transmission!="Low"], aes(x=ITN_Retention_Halflife, y=mid.pct_reduction, group=year_label_short, fill=ITN_Retention_Halflife, ymin=mid.pct_reduction-iqr.pct_reduction, ymax=mid.pct_reduction+iqr.pct_reduction)) +
    geom_bar(stat="identity", position = position_dodge(width=0.9), color="black") +
    geom_errorbar (position=position_dodge(width=0.9), colour="black", size=0.25) +
    geom_text(aes(label=year_label_short, y=3), position = position_dodge(width=1), size=1.5) +
    scale_fill_brewer(type="seq", palette = "YlGnBu", name="Median Retention\nTime") + 
    theme_minimal() + 
    theme(text=element_text(size=6)) +
    facet_grid(ITN_Coverage ~ transmission) +
    labs(x="Median Retention Time",
         y=paste("% Reduction in", metric_label),
         title=paste("Site", this_site, ": Reduction in", metric_label, "by Transmission Intensity (Columns), Median Retention Time,\nYears Since Distribution (Grouped Bars), and ITN Coverage (Rows)"))
  
  barplot_ret <- ggplot(this_subset_summary[transmission!="Low"], aes(x=transmission, y=mid.pct_reduction, fill=ITN_Retention_Halflife, ymin=mid.pct_reduction-iqr.pct_reduction, ymax=mid.pct_reduction+iqr.pct_reduction)) +
    geom_bar(stat="identity", position = "dodge", color="black") +
    geom_errorbar (position=position_dodge(width=0.9), colour="black", size=0.25) +
    scale_fill_brewer(type="seq", palette = "YlGnBu", name="Median Retention\nTime") + 
    theme_minimal() + 
    theme(text=element_text(size=6)) +
    facet_grid(ITN_Coverage ~ year_label) +
    labs(x="Transmission Level",
         y=paste("% Reduction in", metric_label),
         title=paste("Site", this_site, ": Reduction in", metric_label, "by Transmission Intensity, Median Retention Time,\nYears Since Distribution (Columns), and ITN Coverage (Rows)"))
  
  # mapping plots 
  this_cluster_dt <- copy(cluster_dt)
  this_cluster_dt[, value:=factor(value==this_site)]
  this_color <-  palette[this_site]
  
  cluster_plot <- ggplot() +
    geom_raster(data = this_cluster_dt, aes(fill = value, y = lat, x = long)) +
    geom_path(data = africa_dt, aes(x = long, y = lat, group = group), color = "black", size = 0.3) + 
    scale_fill_manual(values= c("lightgrey", this_color)) +
    coord_equal(xlim = c(-18, 52), ylim = c(-35, 38)) +
    labs(x = NULL, y = NULL, title = NULL) +
    theme_classic(base_size = 12) +
    theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
          plot.margin = unit(c(0, 0, 0, 0), "in"), legend.position = "none")
  
  timeseries_plot <- ggplot(time_series[cluster==this_site & cov=="precip_era5"], aes(x=as.integer(variable_val), y=median, color=cluster, fill=cluster)) +
                geom_ribbon(aes(ymin=perc_25, ymax=perc_75), alpha=0.5, color=NA) +
                geom_line(size=1) +
                geom_line(aes(y=perc_05), size=0.75, linetype=2) +
                geom_line(aes(y=perc_95), size=0.75, linetype=2) +
                scale_color_manual(values = this_color) +
                scale_fill_manual(values = this_color) + 
                scale_x_continuous(breaks=seq(2,12,2), labels=c("F","A","J","A","O","D"), minor_breaks=seq(1,12,2)) +
                theme_minimal() +
                theme(legend.position="none",
                      plot.title = element_text(size=8),
                      strip.background = element_blank(),
                      strip.text.y = element_blank()) +
                labs(title="Rainfall",
                     x="",
                     y="")
  
  # save all
  main_vp <- viewport(width = 0.8, height = 1, x = 0.4, y = 0.5)
  map_vp <- viewport(width = 0.3, height = 0.3, x = 0.8, y = 0.2)
  timeseries_vp <- viewport(width = 0.3, height = 0.3, x = 0.8, y = 0.75)
  
  
  pdf(file.path(out_dir, paste0("time_lines_site_", this_site, ".pdf")), height=6, width=8)
    print(lineplot, vp=main_vp)
    print(cluster_plot, vp=map_vp)
    print(timeseries_plot, vp=timeseries_vp)
  graphics.off()
  
  pdf(file.path(out_dir, paste0("bars_bytime_site_", this_site, ".pdf")), height=5, width=6)
    print(barplot_time, vp=main_vp)
    print(cluster_plot, vp=map_vp)
    print(timeseries_plot, vp=timeseries_vp)
  graphics.off()
  
  pdf(file.path(out_dir, paste0("bars_byret_site_", this_site, ".pdf")), height=5, width=6)
    print(barplot_ret, vp=main_vp)
    print(cluster_plot, vp=map_vp)
    print(timeseries_plot, vp=timeseries_vp)
  graphics.off()
  
  
}




# Exponential decay curves
find_decay <- function(init, mean_decay, time){
  return(init*exp(-(1/mean_decay)*time))
}


kill_median_dur <- 1460
block_median_dur <- 730
days <- 0:1095
decay_dt <- data.table(day=days,
                       kill=find_decay(kill_to_plot, kill_median_dur, days),
                       block=find_decay(block_to_plot, block_median_dur, days))

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
names(retention_decay) <- c("median_retention", "day")
retention_decay[, decay:= find_decay(1, median_retention, day)]

pdf(file.path(out_dir, "retention_decay.pdf"), height=3, width=5)
ggplot(retention_decay, aes(x=day/365, color=factor(median_retention), y=decay)) +
  ylim(0,1) +
  geom_line() +
  theme_classic() +
  theme(legend.position = "none",
        text=element_text(size=6)) +
  labs(x="Year",
       y="Prop. of Nets Retained")
graphics.off()

# # smoothing for line plots
# get_smooth <- function(x, y){
#   if (max(y)<0.05){
#     return(y)
#   }else{
#     lo <- loess(y[y>0]~x[y>0])
#     predictions <- c(y[y==0], predict(lo))
#     return(pmax(predictions, rep(0, length(predictions))))
#   }
# }
# 
# smooth_compare <- lapply(unique(compare_summary$Site_Name), function(site_name){
#   sub_list <- lapply(unique(compare_summary$int_id), function(int_name){
#     sub_sub_list <- lapply(unique(compare_summary$day), function(this_day){
#       subset <- compare_summary[Site_Name==site_name & int_id==int_name & day==this_day]
#       # subset[, smooth_min:= get_smooth(mean_initial, min_final)]
#       # subset[, smooth_max:= get_smooth(mean_initial, max_final)]
#       subset[, smooth_mean:= get_smooth(control_val, final_val)]
#       return(subset)
#     })
#     sub_sub_list <- rbindlist(sub_sub_list)
#   })
#   sub_list <- rbindlist(sub_list)
# })
# smooth_compare <- rbindlist(smooth_compare)
# 



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









