
## -----------------------------------------------------------------------------------------------------------------
# compare_eirs.r
# 
# Amelia Bertozzi-Villa, Institute for Disease Modeling
# April 2022 
# 
## -----------------------------------------------------------------------------------------------------------------------

rm(list = ls())

library(rjson)
library(data.table)
library(ggplot2)

fname <- "20220413_emodpy_replicate_itns"

in_dir <- file.path("~/Dropbox (IDM)/Malaria Team Folder/projects/map_intervention_impact/intervention_impact/", fname, "results/raw")

all_data <- lapply(c("eir", "inc", "prev"), function(varname){
  this_data <- fread(file.path(in_dir, varname, paste0(fname, "_Burnin_", varname, ".csv")))
  this_data[, metric:=varname]
  setnames(this_data, paste0("initial_", varname), "val")
  
  if (varname=="inc"){
    this_data[, initial_severe_inc:=NULL]
  }
  
  return(this_data)
})


all_data <- rbindlist(all_data)
all_data_wide <- dcast.data.table(all_data, day + Site_Name + Run_Number + x_Temporary_Larval_Habitat ~ metric, value.var ="val")

ggplot(all_data[day==max(day)], aes(x=x_Temporary_Larval_Habitat, y=val)) +
  geom_point() +
  facet_grid(metric~Site_Name, scales="free") 
  #geom_hline(yintercept=500)
