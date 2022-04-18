
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

varname <- "eir"

in_dir <- "~/Dropbox (IDM)/Malaria Team Folder/projects/map_intervention_impact/intervention_impact/20220412_burnin_bugfix_test_v2/results/raw"

# vectors
dtk_vecs <- fread(file.path(in_dir, "ReportVectorStats_dtktools.csv"))
idm_vecs <- fread(file.path(in_dir, "ReportVectorStats_idmtools.csv"))
dtk_vecs[, type:="dtktools"]
idm_vecs[, type:="idmtools"]

all_vecs <- rbind(dtk_vecs, idm_vecs, fill=T)

# ggplot(all_vecs, aes(x=Time, y=VectorPopulation, color=type, linetype=Species)) +
#   geom_line() +
#   facet_wrap(~NodeID)

all_data <- NULL

idx <- 1

for (varname in c("eir", "prev")){
  
  buggy_data <- fread(file=file.path(in_dir, varname, "MAP_20210511_itn_burnin_Burnin.csv"))
  bugfix_v1 <- fread(file=file.path(in_dir, varname, "20220406_burnin_bugfix_test_Burnin.csv"))
  bugfix_v2 <- fread(file=file.path(in_dir, varname, "20220412_burnin_bugfix_test_v2_Burnin.csv"))
  bugfix_idmtools <- fread(file=file.path(in_dir, varname, paste0("20220412_burnin_bugfix_test_v2_Burnin_",varname, ".csv")))
  bugfix_idmtools_loweracquire <- fread(file=file.path(in_dir, varname, paste0("20220412_burnin_bugfix_test_v2_Burnin_",varname, "_fix1.csv")))
  bugfix_idmtools_v3 <- fread(file=file.path(in_dir, varname, paste0("20220412_burnin_bugfix_test_v2_Burnin_",varname, "_fix2.csv")))
  
  buggy_data[, type:="Old_Bugs"]
  bugfix_v1[, type:="Bugfix_1"]
  bugfix_v2[, type:="Bugfix_2"]
  bugfix_idmtools[, type:="Bugfix_2_idmtools"]
  bugfix_idmtools_loweracquire[, type:="Bugfix_2_idmtools_lower_acquire"]
  bugfix_idmtools_v3[, type:="Bugfix_3_idmtools"]
  
  all_var_data <- rbind(buggy_data[day==1825 & Run_Number==0 & Site_Name<11],
                    bugfix_v1[day==1825],
                    bugfix_v2[day==1825],
                    bugfix_idmtools[day==1825],
                    bugfix_idmtools_loweracquire[day==1825],
                    bugfix_idmtools_v3[day==1825]
  )
  
  all_var_data[, Variable:=varname]
  setnames(all_var_data, paste0("initial_", varname), "val")
  
  all_data[[idx]] <- all_var_data
  idx <- idx + 1
  
}

all_data <- rbindlist(all_data)
all_data_wide <- dcast.data.table(all_data, day + Site_Name + x_Temporary_Larval_Habitat + type ~ Variable, value.var ="val")


ggplot(all_data_wide[x_Temporary_Larval_Habitat<=30 & type %like% "idmtools"], aes(x=x_Temporary_Larval_Habitat, y=eir)) +
  geom_point(aes(color=type)) +
  facet_wrap(~Site_Name, scales="free") +
  #geom_hline(yintercept=500) +
  labs(x="X_temp",
       y=varname)
