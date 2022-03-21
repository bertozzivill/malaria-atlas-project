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
library(RColorBrewer)

rm(list=ls())

analysis_subdir <- "20220218_higher_cm"
analysis_metric <- "inc"
final_day <- 1095
suffix <- ""

main_dir <- file.path(Sys.getenv("HOME"), 
                      "Dropbox (IDM)/Malaria Team Folder/projects/map_intervention_impact/intervention_impact",
                      analysis_subdir)
out_dir <- file.path(main_dir,"results", "figs")
dir.create(out_dir, recursive = T, showWarnings = F)

# for mapping
arch_dir <- file.path(main_dir, "../../archetypes/results/v4_era5_bounded_transmission/africa")
africa_shp_dir <- "/Volumes/GoogleDrive/My Drive/itn_cube/input_data/general/shapefiles/Africa.shp"


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

# compare results to no-ITN scenario with different AM levels 
no_int <- final[ITN_Coverage==0 & CM_Coverage==0, list(Site_Name, Run_Number, x_Temporary_Larval_Habitat, day, no_int_val=final_val)]
control <- final[ITN_Coverage==0 & CM_Coverage>0, list(Site_Name, Run_Number, x_Temporary_Larval_Habitat, CM_Coverage, day, control_val=final_val)]
control <- merge(no_int, control)
keycols <- c(key(control), "CM_Coverage")
setkeyv(control, keycols)

int <- final[ITN_Coverage>0]
compare <- merge(control, int)

# plot control against end of initial time series to see differences

control_test <- data.table::copy(control)
control_test[, day:=day+ max(initial$day)]
last_initial <- control_test[day==max(day), list(Site_Name, Run_Number, x_Temporary_Larval_Habitat, CM_Coverage, day=max(initial$day))]
setkey(last_initial, NULL)
last_initial <- merge(last_initial, initial[day==max(day)])
setnames(last_initial, "initial_val", "control_val")
control_test <- rbind(control_test, last_initial, fill=T)
control_test[, year:=day/365]
control_test[, mean_control:=mean(control_val), by=list(Site_Name, x_Temporary_Larval_Habitat, CM_Coverage, day, year)]
control_test[, mean_noint:=mean(no_int_val), by=list(Site_Name, x_Temporary_Larval_Habitat, CM_Coverage, day, year)]

initial[, year:=day/365]
initial[, mean_initial:=mean(initial_val), by=list(Site_Name, x_Temporary_Larval_Habitat, day, year)]

initial_subset <- initial[day>16000  & x_Temporary_Larval_Habitat>0.03 & x_Temporary_Larval_Habitat<50]
control_subset <- control_test[x_Temporary_Larval_Habitat>0.03]

ggplot() +
  geom_point(data=control_subset, aes(x=year, y=no_int_val), alpha=0.5) +
  geom_point(data=control_subset, aes(x=year, y=control_val, color=factor(CM_Coverage)), alpha=0.5) +
  geom_point(data=initial_subset, aes(x=year, y=initial_val), alpha=0.5) + 
  geom_line(data=control_subset, aes(x=year, y=mean_noint)) +
  geom_line(data=control_subset, aes(x=year, y=mean_control, color=factor(CM_Coverage))) +
  geom_line(data=initial_subset, aes(x=year, y=mean_initial)) + 
  theme_minimal() + 
  facet_grid(Site_Name ~ x_Temporary_Larval_Habitat)+
  labs(x="Year", y="Incidence")

# compare different blocking scenarios
different_blocking <- T
if (different_blocking){
  #compare[ITN_Blocking_Halflife==36500, blocking_label:= "Constant"]
  compare[ITN_Blocking_Halflife==ITN_Retention_Halflife, blocking_label:= "Matched Decay"]
  compare[ITN_Blocking_Halflife==730, blocking_label := "2yr Decay"]
  compare[duplicated(compare), blocking_label:="Matched Decay"]
  compare[, blocking_label:= factor(blocking_label, levels=c("2yr Decay", "Matched Decay" #,
                                                             #"Constant"
                                                             ))]
}

# assign intervention ids 
ints <- unique(compare[, list(blocking_label, 
                              CM_Coverage,
                              ITN_Coverage, 
                              ITN_Retention_Halflife, 
                              ITN_Blocking_Halflife, 
                              ITN_Initial_Kill, 
                              ITN_Initial_Block, 
                              ITN_Start)])[order(blocking_label, 
                                                 CM_Coverage,
                                                 ITN_Coverage, 
                                                 ITN_Retention_Halflife, 
                                                 ITN_Blocking_Halflife, 
                                                 ITN_Initial_Kill, 
                                                 ITN_Initial_Block, 
                                                 ITN_Start)]
ints[, int_id:= as.integer(rownames(ints))]
compare <- merge(compare, ints, by=c("blocking_label", "CM_Coverage", "ITN_Coverage", "ITN_Retention_Halflife", "ITN_Blocking_Halflife", "ITN_Initial_Kill", "ITN_Initial_Block", "ITN_Start"))


# for incidence, change units from cases/person/yr to cases/1000 people/yr
if (analysis_metric=="inc"){
  compare[, final_val:=final_val*1000]
  compare[, control_val:=control_val*1000]
  compare[, no_int_val:= no_int_val*1000]
}

compare[, change:= final_val-control_val]
compare[, reduction:= -change]
compare[, pct_change:= 100*(final_val-control_val)/control_val]
compare[, pct_reduction:= -pct_change]

# also calculate reduction compared to the current "standard"
current_standard <- compare[blocking_label== "2yr Decay", list(Site_Name, x_Temporary_Larval_Habitat, day, CM_Coverage, ITN_Coverage, ITN_Initial_Kill, ITN_Initial_Block, Run_Number, standard_pct_reduction=pct_reduction)]

# compare <- merge(compare, current_standard, by=intersect(colnames(compare), colnames(current_standard)), allow.cartesian=T)
# compare[, diff_pct_reduction:= pct_reduction-standard_pct_reduction]

compare_summary <- compare[, unlist(recursive=FALSE, lapply(
  .(mean = mean, mid=median, iqr = IQR),
  function(f) lapply(.SD, f)
)),
by = .(Site_Name, day, int_id, CM_Coverage, ITN_Coverage, ITN_Retention_Halflife, ITN_Initial_Block, ITN_Blocking_Halflife, ITN_Initial_Kill, x_Temporary_Larval_Habitat, blocking_label), 
.SDcols = c("reduction", "pct_reduction")]

compare[, final_val_mean:= mean(final_val), by=.(Site_Name, day, int_id, ITN_Coverage, ITN_Retention_Halflife, ITN_Initial_Block, ITN_Initial_Kill, x_Temporary_Larval_Habitat)]
compare[, control_val_mean:= mean(control_val), by=.(Site_Name, day, int_id, ITN_Coverage, ITN_Retention_Halflife, ITN_Initial_Block, ITN_Initial_Kill, x_Temporary_Larval_Habitat)]
compare[, noint_val_mean:= mean(no_int_val), by=.(Site_Name, day, int_id, ITN_Coverage, ITN_Retention_Halflife, ITN_Initial_Block, ITN_Initial_Kill, x_Temporary_Larval_Habitat)]



# Find which x_temps correspond to the desired EIRs
eir_dir <- file.path(in_dir, "../eir")
eir_fnames <- list.files(eir_dir)
eirs <- fread(file.path(eir_dir, eir_fnames[eir_fnames %like% "Burnin"]))

setnames(eirs, "initial_eir", "eir")
eirs[, x_Temporary_Larval_Habitat:= round(x_Temporary_Larval_Habitat, 4)]
# eirs <- eirs[, lapply(.SD, mean), by = .(Site_Name, day, x_Temporary_Larval_Habitat), .SDcols = c("eir")]
eirs <- eirs[day==max(day)]
eirs <- eirs[, lapply(.SD, mean), by = .(Site_Name,  x_Temporary_Larval_Habitat), .SDcols = c("eir")]
eirs <- eirs[order(Site_Name, x_Temporary_Larval_Habitat)]
eirs[, eir:= round(eir, 2)]

target_eirs <- c(0.5, 5, 50)

closest_eirs <- data.table(expand.grid(Site_Name=1:10, match_eir=target_eirs))
eirs[, match_eir:=eir]

setkeyv(eirs, c("Site_Name", 'match_eir'))

closest_eirs <- eirs[closest_eirs, roll='nearest']
closest_eirs[, transmission:= factor(match_eir, levels=target_eirs, labels=c("Low Transmission", "Medium Transmission", "High Transmission"))]

for_output_eirs <- dcast.data.table(closest_eirs[transmission!="Low"], Site_Name  ~ transmission, value.var = "eir")

write.csv(for_output_eirs, file.path(out_dir, "output_eirs.csv"), row.names = F)
write.csv(closest_eirs, file.path(out_dir, "output_eirs_and_larval_habs.csv"), row.names = F)

day_to_plot <- unique(compare$day)  # 365
kill_to_plot <- c(0.4)
ret_to_plot <- unique(compare$ITN_Retention_Halflife) # c(730)
block_to_plot <- c(0.9)
cov_to_plot <- c(0.4) # unique(compare$ITN_Coverage)
sites_to_plot <-  2:10
eirs_for_subsetting <- closest_eirs[Site_Name %in% sites_to_plot, list(Site_Name, x_Temporary_Larval_Habitat, eir, transmission)]

subset <- compare[day %in% day_to_plot & Site_Name %in% sites_to_plot & ITN_Initial_Kill %in% kill_to_plot & ITN_Retention_Halflife %in% ret_to_plot & ITN_Initial_Block %in% block_to_plot & ITN_Coverage %in% cov_to_plot]
subset_summary <- compare_summary[day %in% day_to_plot & Site_Name %in% sites_to_plot & ITN_Initial_Kill %in% kill_to_plot & ITN_Retention_Halflife %in% ret_to_plot & ITN_Initial_Block %in% block_to_plot & ITN_Coverage %in% cov_to_plot]

# temp to diagnose absolute reduction differences 

test_subset <- merge(subset, eirs, by=c("Site_Name", "x_Temporary_Larval_Habitat"))

ggplot(test_subset[CM_Coverage==0.4 & Site_Name==4], aes(x=day, color=factor(ITN_Retention_Halflife), linetype=blocking_label)) +
  geom_point( aes(y=final_val), size=0.5) + 
  geom_point(aes(y=control_val), color="black", size=0.5) +
  geom_line(aes(y=final_val_mean)) +
  geom_line(aes(y=control_val_mean), color="black") +
  theme_minimal() + 
  facet_wrap(~eir) +
  # theme(text=element_text(size=6)) +
  labs(x="Years since distribution",
       y="Incidence",
       color="Median Retention\nTime",
       title=paste("Incidence over Time for No-Intervention Scenario (Black) vs Intervention Scenarios"))


# --end temp


subset <- merge(eirs_for_subsetting, subset, by=c("Site_Name", "x_Temporary_Larval_Habitat"), all.x=T)
subset_summary <- merge(eirs_for_subsetting, subset_summary, by=c("Site_Name", "x_Temporary_Larval_Habitat"), all.x=T)

subset <- subset[transmission!="Low Transmission"]
subset_summary <- subset_summary[transmission!="Low Transmission"]

metric_label <- ifelse(analysis_metric=="inc", "Incidence", "Prevalence")

subset[, ITN_Initial_Kill:= factor(ITN_Initial_Kill)]
subset_summary[, ITN_Initial_Kill:= factor(ITN_Initial_Kill)]

subset[, ITN_Retention_Halflife:=factor(ITN_Retention_Halflife, labels=c("1 Year", "2 Years", "3 Years", "4 Years"))]
subset_summary[, ITN_Retention_Halflife:=factor(ITN_Retention_Halflife, labels=c("1 Year", "2 Years", "3 Years", "4 Years"))]

subset[, year:=day/365]
subset_summary[, year:=day/365]

subset_summary[, year_label:=paste(year, "Yrs Since Dist.")]
subset_summary[, year_label_short:=paste0("Yr", year)]
subset[, year_label:=paste(year, "Yrs Since Dist.")]
subset[, year_label_short:=paste0("Yr", year)]

# for color plotting
subset_summary[, for_colors:= factor(paste(blocking_label, ITN_Retention_Halflife))]
fill_palette <- c(brewer.pal(5, "YlOrBr")[2:5],  
                  #brewer.pal(5, "YlGn")[2:5], 
                  brewer.pal(4, "YlGnBu"))

subset[, mean(control_val), by="transmission"]

# Save dataset used for plots
write.csv(subset, file.path(out_dir, "subset_data_for_plots.csv"), row.names = F)
write.csv(subset_summary, file.path(out_dir, "subset_summary_data_for_plots.csv"), row.names = F)


compare_datasets <- F

if (compare_datasets){
  # load some prior results for comparison plots
  additional_subset_summary <- fread(file.path(gsub(analysis_subdir, "20210520_rerun_blocktime", out_dir), "subset_summary_data_for_plots.csv"), stringsAsFactors = T)
  subset_summary[, c("year_label", "year_label_short")] <- lapply(subset_summary[, c("year_label", "year_label_short")], factor)
  subset_summary <- rbind(subset_summary, additional_subset_summary)
  subset_summary[, ITN_Retention_Halflife := factor(ITN_Retention_Halflife, levels=c("6 Months", "1 Year", "2 Years", "3 Years", "4 Years"))]
}


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


ggplot(subset, aes(x=year, y=control_val, color=factor(CM_Coverage))) +
  geom_point() +
  geom_line(aes(y=control_val_mean)) +
  facet_grid( transmission ~ Site_Name, scales="free")

metric_type <- "pct_reduction"; this_site <- 4

for (metric_type in c("pct_reduction")){
  
  this_out_dir <- file.path(out_dir, metric_type)
  dir.create(this_out_dir, recursive = T, showWarnings = F)
  
  for (subdir in c("lineplots", "byret", "byblock")){
    dir.create( file.path(this_out_dir, subdir), recursive = T, showWarnings = F)
  }
  
  print(metric_type)
  
  for (this_site in unique(subset$Site_Name)){
    print(this_site)
    
    # main data subsetting
    this_subset <- subset[Site_Name==this_site]
    this_subset[, cm_cov_label:= factor(CM_Coverage, labels=paste0(unique(this_subset$CM_Coverage)*100, "% CM"))  ]
    this_subset_summary <- subset_summary[Site_Name==this_site & transmission=="Medium Transmission"] # TEMP: only look at medium transmission
    this_subset_summary[, cm_cov_label:= factor(CM_Coverage, labels=paste0(unique(this_subset$CM_Coverage)*100, "% CM"))  ]
    
    # rename according to outcome metric
    if (metric_type=="reduction"){
      this_subset_summary[, plotting_median:=mid.reduction]
      this_subset_summary[, plotting_iqr:=iqr.reduction]
      metric_type_label <- "Reduction"
      axis_suffix <- "(Cases/1000/Yr)"
    }else if (metric_type=="pct_reduction"){
      this_subset_summary[, plotting_median:=mid.pct_reduction]
      this_subset_summary[, plotting_iqr:=iqr.pct_reduction]
      metric_type_label <- "% Reduction"
      axis_suffix <- ""
    }
    
    # main data plots 
    lineplot <- ggplot(this_subset, aes(x=year, color=ITN_Retention_Halflife, linetype=blocking_label)) +
      geom_line(aes(y=noint_val_mean), color="black", size=1, alpha=0.5, linetype="solid") +
      geom_line(aes(y=final_val_mean)) +
      geom_line(aes(y=control_val_mean), color="black") +
      theme_minimal() + 
      geom_point( aes(y=final_val), size=0.5) + 
      geom_point(aes(y=control_val), color="black", size=0.5) +
      guides(linetype=guide_legend(title="Blocking Type")) +
      facet_grid(transmission~cm_cov_label) +
      theme(text=element_text(size=11)) +
      labs(x="Years since distribution",
           y=metric_label,
           color="Median Retention\nTime",
           title=paste(metric_label, "over Time for No-Intervention Scenario (Gray),\nNo-Net Scenario (Black), and Net Scenarios (Color)"))
    
    barplot_ret <- ggplot(this_subset_summary, aes(x=year_label, y=plotting_median, group=ITN_Retention_Halflife, fill=for_colors, ymin=plotting_median-plotting_iqr, ymax=plotting_median+plotting_iqr)) +
                    geom_bar(stat="identity", position = "dodge", color="black") +
                    geom_errorbar (position=position_dodge(width=0.9), colour="black", size=0.25) +
                    geom_hline(yintercept=100, size=0.5) + 
                    scale_fill_manual(values=fill_palette, name="Median Retention\nTime") +
                    theme_minimal() + 
                    ylim(-25, 125) +
                    theme(text=element_text(size=12),
                          axis.text.x = element_text(angle = 15),
                          legend.position = "none") +
                    facet_grid(cm_cov_label ~ blocking_label) +
                    labs(x="Decay Type",
                         y=paste(metric_type_label, "in", metric_label, axis_suffix),
                         #title=paste("Site", this_site, ":", metric_type_label, "in", metric_label, "by Transmission Intensity, Median Retention Time,\nYears Since Distribution, and ITN Blocking")
                         title=NULL
                         )
    
    barplot_block <- ggplot(this_subset_summary, aes(x=blocking_label, y=plotting_median, group=ITN_Retention_Halflife, fill=for_colors, ymin=plotting_median-plotting_iqr, ymax=plotting_median+plotting_iqr)) +
                      geom_bar(stat="identity", position = "dodge", color="black") +
                      geom_errorbar (position=position_dodge(width=0.9), colour="black", size=0.25) +
                      geom_hline(yintercept=100, size=0.5) + 
                      scale_fill_manual(values=fill_palette, name="Median Retention\nTime") +
                      theme_minimal() + 
                      ylim(-25, 125) +
                      theme(text=element_text(size=12),
                            legend.position = "none") +
                      facet_grid(cm_cov_label ~ year_label) +
                      labs(x="Transmission Level",
                           y=paste(metric_type_label, "in", metric_label, axis_suffix),
                           # title=paste("Site", this_site, ":", metric_type_label, "in", metric_label, "by Transmission Intensity, Median Retention Time,\nYears Since Distribution, and ITN Blocking")
                           title=NULL
                           )
                    
    compare_to_std <- F
    if (metric_type=="pct_reduction" & compare_to_std){
      # plot as a relative change from "current" nets
      
      compare_to_std_byret <- ggplot(this_subset_summary, aes(x=year_label, y=mid.diff_pct_reduction, group=ITN_Retention_Halflife, fill=for_colors,  ymin=mid.diff_pct_reduction-iqr.diff_pct_reduction, ymax=mid.diff_pct_reduction+iqr.diff_pct_reduction)) +
                                geom_bar(stat="identity", position = "dodge", color="black") +
                                geom_errorbar (position=position_dodge(width=0.9), colour="black", size=0.25) +
                                scale_fill_manual(values=fill_palette, name="Median Retention\nTime") +
                                theme_minimal() + 
                                theme(text=element_text(size=12),
                                      axis.text.x = element_text(angle = 15),
                                      legend.position = "none") +
                                facet_grid(blocking_label ~ transmission) +
                                labs(x="Years Since Distribution",
                                     y="Percentage point reduction (relative to 'Standard' net)",
                                     # title=paste("Site", this_site, ":", metric_type_label, "in", metric_label, "by Transmission Intensity, Median Retention Time,\nYears Since Distribution (Columns), and ITN Coverage (Rows)")
                                     )
      
      compare_to_std_byblock <- ggplot(this_subset_summary, aes(x=blocking_label, y=mid.diff_pct_reduction, group=ITN_Retention_Halflife, fill=for_colors,  ymin=mid.diff_pct_reduction-iqr.diff_pct_reduction, ymax=mid.diff_pct_reduction+iqr.diff_pct_reduction)) +
                                  geom_bar(stat="identity", position = "dodge", color="black") +
                                  geom_errorbar (position=position_dodge(width=0.9), colour="black", size=0.25) +
                                  scale_fill_manual(values=fill_palette, name="Median Retention\nTime") +
                                  theme_minimal() + 
                                  theme(text=element_text(size=12),
                                        legend.position = "none") +
                                  facet_grid(year_label ~ transmission) +
                                  labs(x="Years Since Distribution",
                                       y="Percentage point reduction (relative to 'Standard' net)",
                                       # title=paste("Site", this_site, ":", metric_type_label, "in", metric_label, "by Transmission Intensity, Median Retention Time,\nYears Since Distribution (Columns), and ITN Coverage (Rows)")
                                  )
                                
      
      
    }
    
    
  
    # mapping plots 
    this_cluster_dt <- data.table::copy(cluster_dt)
    this_cluster_dt[, value:=factor(value==this_site)]
    this_color <-  palette[this_site]
    
    cluster_plot <- ggplot() +
      geom_raster(data = this_cluster_dt, aes(fill = value, y = lat , x = long)) +
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
            plot.title = element_text(size=12, hjust=0.5),
            strip.background = element_blank(),
            strip.text.y = element_blank(),
            axis.text.y = element_blank(),
            axis.title.x = element_text(size=10)) +
      labs(title="Rainfall",
           x="",
           y="")
    
    for_legend <- unique(this_subset_summary[, list(ITN_Retention_Halflife, blocking_label, for_colors)])
    
    legend <- ggplot(for_legend) +
      geom_raster(aes(x = ITN_Retention_Halflife, y = blocking_label, fill=for_colors), show.legend = F) +
      scale_fill_manual(values = fill_palette) +
      scale_y_discrete(limits = rev(levels(for_legend$blocking_label))) + 
      coord_equal() +
      labs(x = "ITN Retention Halflife", y = NULL, title = NULL) +
      theme_minimal() +
      theme(text=element_text(size=12),
            axis.line = element_blank(), axis.ticks = element_blank(),
            panel.grid = element_blank(), axis.text.y = element_text(angle = 45, hjust = 0.5),
            axis.text.x = element_text(angle = 45, hjust = 1),
            plot.margin = unit(c(0, 0, 0, 0), "in"), panel.border = element_rect(fill = NA, color = "black"))
    
    
    # save all 
    main_vp <- viewport(width = 0.7, height = 1, x = 0.35, y = 0.5)
    legend_vp <- viewport(width = 0.3, height = 0.4, x = 0.85, y = 0.2)
    timeseries_vp <- viewport(width = 0.3, height = 0.3, x = 0.85, y = 0.75)
    map_vp <- viewport(width = 0.3, height = 0.2, x=0.87, y=0.5)
    
    # different viewports for line plots
    pdf(file.path(this_out_dir, "lineplots", paste0("time_lines_site_", this_site, ".pdf")), height=6, width=8)
    print(lineplot, vp=viewport(width = 1, height = 1, x = 0.5, y = 0.5))
    print(cluster_plot, vp=viewport(width = 0.125, height = 0.125, x=0.875, y=0.9))
    print(timeseries_plot, vp=viewport(width = 0.175, height = 0.175, x = 0.875, y = 0.75)) 
    graphics.off()
    
    pdf(file.path(this_out_dir, "byblock", paste0("bars_byblock_", metric_type, "_site_", this_site, ".pdf")),  height=6.25, width=10)
    print(barplot_block, vp=main_vp)
    print(cluster_plot, vp=map_vp)
    print(timeseries_plot, vp=timeseries_vp)
    print(legend, vp=legend_vp)
    graphics.off()
    
    pdf(file.path(this_out_dir, "byret", paste0("bars_byret_", metric_type, "_site_", this_site, ".pdf")),  height=6.25, width=10)
    print(barplot_ret, vp=main_vp)
    print(cluster_plot, vp=map_vp)
    print(timeseries_plot, vp=timeseries_vp)
    print(legend, vp=legend_vp)
    graphics.off()
    
    if (metric_type=="pct_reduction" & compare_to_std){
      
      for (subdir in c("compare_byret", "compare_byblock")){
        dir.create( file.path(this_out_dir, subdir), recursive = T, showWarnings = F)
      }
      
      
      pdf(file.path(this_out_dir, "compare_byret", paste0("compare_to_std_byret_", metric_type, "_site_", this_site, ".pdf")), height=6.25, width=10)
      print(compare_to_std_byret, vp=main_vp)
      print(cluster_plot, vp=map_vp)
      print(timeseries_plot, vp=timeseries_vp)
      print(legend, vp=legend_vp)
      graphics.off()
      
      pdf(file.path(this_out_dir, "compare_byblock", paste0("compare_to_std_bytime_", metric_type, "_site_", this_site, ".pdf")),  height=6.25, width=10)
      print(compare_to_std_byblock, vp=main_vp)
      print(cluster_plot, vp=map_vp)
      print(timeseries_plot, vp=timeseries_vp)
      print(legend, vp=legend_vp)
      graphics.off()
      

      
    }
    
    
  }
  
  
  
  
}






# Exponential decay curves
find_decay <- function(init, half_life, time){
  return(init*2^(-(1/half_life)*time))
}


kill_median_dur <- 1460
ret_median_durs <- c(365, 730, 1095, 1460)

days <- 0:1095

this_out_dir <- file.path(out_dir, "decay_plots")
dir.create(this_out_dir, recursive = T, showWarnings = F)

palette <- brewer.pal(4, "YlGnBu")
palette[1] <- "#f8e963"

retention_decay <- data.table(expand.grid(ret_median_durs, days))
names(retention_decay) <- c("median_retention", "day")
retention_decay[, decay:= find_decay(1, median_retention, day)]
retention_decay[, decay_block:=find_decay(0.9, median_retention, day)]

pdf(file.path(this_out_dir, "retention_decay.pdf"), height=3, width=5)
ggplot(retention_decay, aes(x=day/365, color=factor(median_retention), y=decay)) +
  ylim(0,1) +
  geom_line() +
  theme_classic() +
  scale_color_manual(values = palette, name="Median Retention\nTime") +
  theme(legend.position = "none",
        text=element_text(size=6)) +
  labs(x="Year",
       y="Prop. of Nets Retained")
graphics.off()


# add different blocking decays
block_decay <- data.table(expand.grid(c(730, 36500), days))
names(block_decay) <- c("median_blocking", "day")
block_decay[, decay:= find_decay(0.9, median_blocking, day)]

pdf(file.path(this_out_dir, "block_decay_current.pdf"), height=3, width=5)
ggplot(retention_decay, aes(x=day/365, color=factor(median_retention), y=decay)) +
  ylim(0,1) +
  geom_line() +
  geom_line(data=block_decay[median_blocking==730], color="black", linetype="dashed") +
  theme_classic() +
  scale_color_manual(values = palette, name="Median Retention\nTime") +
  theme(legend.position = "none",
        text=element_text(size=6)) +
  labs(x="Year",
       y="Prop. of Nets Retained")
graphics.off()

pdf(file.path(this_out_dir, "block_decay_constant.pdf"), height=3, width=5)
ggplot(retention_decay, aes(x=day/365, color=factor(median_retention), y=decay)) +
  ylim(0,1) +
  geom_line() +
  geom_line(data=block_decay[median_blocking==36500], color="black", linetype="dashed") +
  theme_classic() +
  scale_color_manual(values = palette, name="Median Retention\nTime") +
  theme(legend.position = "none",
        text=element_text(size=6)) +
  labs(x="Year",
       y="Prop. of Nets Retained")
graphics.off()

pdf(file.path(this_out_dir, "block_decay_variable.pdf"), height=3, width=5)
ggplot(retention_decay, aes(x=day/365, color=factor(median_retention), y=decay)) +
  ylim(0,1) +
  geom_line() +
  geom_line(aes(group=factor(median_retention), y=decay_block), color="black", linetype="dashed")+
  theme_classic() +
  scale_color_manual(values = palette, name="Median Retention\nTime") +
  theme(legend.position = "none",
        text=element_text(size=6)) +
  labs(x="Year",
       y="Prop. of Nets Retained")
graphics.off()

kill_decay <- data.table(day=days,
                       kill=find_decay(0.4, kill_median_dur, days)
                       )

pdf(file.path(this_out_dir, "kill_decay.pdf"), height=3, width=5)
ggplot(kill_decay, aes(x=day/365,  y=kill)) +
  ylim(0,1) +
  geom_line() +
  theme_classic() +
  theme(legend.position = "none",
        text=element_text(size=6)) +
  labs(x="Year",
       y="Effectiveness")
graphics.off()











