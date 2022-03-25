
## -----------------------------------------------------------------------------------------------------------------
# itn_rain_seasonality_comparison.r
# 
# Amelia Bertozzi-Villa, Institute for Disease Modeling
# March 2022 
# 
# Compare seasonality of net access and rainfall by country
## -----------------------------------------------------------------------------------------------------------------------

rm(list=ls())

library(rgdal)
library(rgeos)
library(raster)
library(rasterVis)
library(data.table)
library(ggplot2)

rainfall_dir <- file.path("~/Dropbox (IDM)/Malaria Team Folder/projects/map_intervention_impact/archetypes/covariates/no_transmission_limits/africa/precip_era5")
access_dir <- "~/repos/map-itn-cube/paper_figures/figure_data/fig_2_access_use_timeseries.csv"
africa_shp_dir <- "/Volumes/GoogleDrive/My Drive/itn_cube/input_data/general/shapefiles/Africa.shp"

africa_shp <- readOGR(africa_shp_dir)
africa_dt <- data.table(fortify(africa_shp, region = "COUNTRY_ID"))
africa_shp_test <- gSimplify(africa_shp, tol=0.1, topologyPreserve=TRUE)

# load and aggregate rainfall rasters
rastlist <- list.files(path = rainfall_dir, pattern='.tif$', all.files=TRUE, full.names=T)
allrasters <- stack(rastlist)
aggregated <- extract(allrasters, africa_shp, fun=mean, na.rm=T, sp=T)
save(aggregated, file=file.path(rainfall_dir, "means_by_country.RData"))

rainfall <- data.table(aggregated@data)
rainfall <- melt.data.table(rainfall, id.vars = c("COUNTRY_ID", "name"), measure.vars = paste0("precip_era5_month_", str_pad(1:12, 2, "left", "0")), value.name = "rainfall")
setnames(rainfall, c("COUNTRY_ID", "variable"), c("iso3", "month"))
rainfall[, month:=as.integer(gsub(".*_([0-9]{2})", "\\1", month))]

#load access values
access <- fread(access_dir)
access <- access[variable=="access"]

#merge and plot
access <- merge(access, rainfall, allow.cartesian=T)
access[, mean_among_atrisk:=mean_among_atrisk*100]
access[, lower_among_atrisk:=lower_among_atrisk*100]
access[, upper_among_atrisk:=upper_among_atrisk*100]
access[, rainfall := rainfall/5]
access <- access[order(iso3, time)]
access[, increasing:= (mean_among_atrisk-data.table::shift(mean_among_atrisk, 1))>0, by=c("iso3", "variable")]

out_dir <- "~/Desktop"
pdf(file.path(out_dir, "rainfall_access_correlation.pdf"), width=5, height=3)

for (this_country in unique(access$iso3)){
  print(this_country)
  this_plot <- ggplot(access[iso3==this_country], aes(x=time)) +
    geom_vline(data=access[iso3==this_country & increasing==T], aes(xintercept=time), color="chartreuse4", alpha=0.5) +
    geom_line(aes(y=rainfall), color="deepskyblue4") +
    #geom_ribbon(aes(ymin=lower_among_atrisk, ymax=upper_among_atrisk), alpha=0.35) +
    geom_line(aes(y=mean_among_atrisk)) +
    facet_wrap(~country_name) +
    labs(x="Year",
         y="Access (%)")
  print(this_plot)
  
}
graphics.off()

ggplot(access, aes(x=time)) +
  geom_vline(data=access[increasing==T], aes(xintercept=time), color="chartreuse4", alpha=0.35) +
  geom_line(aes(y=rainfall), color="deepskyblue4") +
  geom_ribbon(aes(ymin=lower_among_atrisk, ymax=upper_among_atrisk), alpha=0.35) +
  geom_line(aes(y=mean_among_atrisk)) +
  facet_wrap(~country_name) +
  labs(x="Year",
       y="Access (%)")



