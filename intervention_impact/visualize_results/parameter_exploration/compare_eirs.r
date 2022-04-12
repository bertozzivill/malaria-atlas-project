
## -----------------------------------------------------------------------------------------------------------------
# compare_eirs.r
# 
# Amelia Bertozzi-Villa, Institute for Disease Modeling
# April 2022 
# 
## -----------------------------------------------------------------------------------------------------------------------

library(rjson)
library(data.table)
library(ggplot2)

varname <- "prev"

in_dir <- "~/Dropbox (IDM)/Malaria Team Folder/projects/map_intervention_impact/intervention_impact/20220406_burnin_bugfix_test/results/raw"
buggy_data <- fread(file=file.path(in_dir, varname, "MAP_20210511_itn_burnin_Burnin.csv"))
new_data <- fread(file=file.path(in_dir, varname, "20220406_burnin_bugfix_test_Burnin.csv"))

setnames(buggy_data, paste0("initial_", varname), paste0("old_", varname))
setnames(new_data, paste0("initial_", varname), paste0("new_", varname))

all_data <- merge(new_data, buggy_data)

end_year <- all_data[day==1825]
end_year_long <- melt(end_year, id.vars = c("day", "Site_Name", "Run_Number", "x_Temporary_Larval_Habitat"),
                      value.name="Param",
                      variable.name = "Model")


ggplot(end_year_long, aes(x=x_Temporary_Larval_Habitat, y=Param)) +
  geom_point(aes(color=Model)) +
  facet_wrap(~Site_Name) +
  #geom_hline(yintercept=500) +
  labs(x="X_temp",
       y=varname)
