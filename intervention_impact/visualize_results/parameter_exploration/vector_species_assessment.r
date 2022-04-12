
## -----------------------------------------------------------------------------------------------------------------
# vector_species_assessment.r
# 
# Amelia Bertozzi-Villa, Institute for Disease Modeling
# March 2022 
# 
## -----------------------------------------------------------------------------------------------------------------------

library(rjson)
library(data.table)
library(ggplot2)

in_dir <- "~/Dropbox (IDM)/Malaria Team Folder/projects/map_intervention_impact/intervention_impact/20220310_cm_burnin_test/results/raw"
data <- fread(file=file.path(in_dir, "ReportVectorStats.csv"))
desired_props <- fread("~/Dropbox (IDM)/Malaria Team Folder/projects/map_intervention_impact/intervention_impact/20220310_cm_burnin_test/input/vector/vector_proportions.csv")

setnames(data, "NodeID", "id")
data[, TotalVectorPopulation:=sum(VectorPopulation), by=.(Time, id)]
data[, VectorProportion:=VectorPopulation/TotalVectorPopulation]

desired_props <- melt(desired_props, id.vars = c("id", "continent", "lat", "lon"), variable.name = "Species", value.name = "TargetVectorProportion")
data <- merge(data, desired_props, by=c("id", "Species"))

data[, summed:=sum(AdultCount), by=.(Time)]
ggplot(data[id==4 & Species=="gambiae"], aes(x=Time, y=summed)) +
  geom_line()


subset <- data[continent=="Africa" & id>1 & Time<2500 & Species %in% c("arabiensis", "funestus", "gambiae")]

ggplot(subset, aes(x=Time, y=VectorPopulation)) + 
  geom_line(aes(color=Species)) +
  facet_wrap(~id)               

ggplot(subset, aes(x=Time, color=Species)) + 
  geom_line(aes(y=TargetVectorProportion), linetype="dashed") +
  geom_line(aes(y=VectorProportion)) +
  facet_wrap(~id) 


#try to unpack what the water veg and temp rainfall habs are doing

# water vegetation ("semi"):
update_water_veg <- function(water_veg_current, # current larval capacity
                             rain_meters, 
                             K_semi=4e8, # max larval capacity, set by user
                             decay_rate_semi=0.01, # param set by user, Semipermanent_Habitat_Decay_Rate in config
                             D_cell=1, # cell diameter
                             dt=1 # time step: one day
                             ){
  water_veg_new <- rain_meters * K_semi * D_cell^2 - water_veg_current * dt * decay_rate_semi
  return(water_veg_new)
}

# temporary rainfall ("temp")

update_temp_rainfall <- function(temp_rainfall_current, # current larval capacity
                                 rain_meters,
                                 temperature_kelvin, 
                                 rel_humidity, # 0-1
                                 K_temp=4e8, # max larval capacity, set by user
                                 decay_scalar_temp=0.05, # param set by user, Temporary_Habitat_Decay_Factor in config
                                 C_1=-5628.1, # multiplicative factor divided by temperature in exponential in C-C calculation of saturated vapor pressure
                                 C_2=5.1e11, # integration constant in C-C calculation of saturated vapor pressure
                                 C_3=0.00034457, # (0.018 / (2 * pi * R)) = 0.00034457 where R is the ideal gas constant. terms based on molecular weight of water to get evaporation rate
                                 D_cell=1, # cell diameter
                                 dt=1 # time step: one day
                                 ){
  decay_rate_temp <- C_2*exp(C_1/temperature_kelvin) * decay_scalar_temp * sqrt(C_3/temperature_kelvin) * (1-rel_humidity)
  temp_rainfall_new <- rain_meters * K_temp * D_cell^2 - temp_rainfall_current * dt * decay_rate_temp
  return(temp_rainfall_new)
}


rainfall_vals <- c(0.05, 0.1, 0.15, 0.2, 0.25, 0.3)
temperature_vals <- c(20,30,40) + 273.15
rel_hum_vals <- c(0.25, 0.5, 0.75, 1)
max_day <- 30
init_capacity <- 4e8
idx <- 1
dt_list <- NULL

for (rainfall in rainfall_vals){
  for (temp_kelvin in temperature_vals){
    for (rel_hum in rel_hum_vals){
    
      this_dt <- data.table(day=1:max_day, 
                            rainfall=rainfall,
                            temperature=temp_kelvin-273.15,
                            rel_humidity=rel_hum)
        
      veg_list <- c(init_capacity)
      temp_list <- c(init_capacity)
      day <- 2
      while (day<=max_day){
        veg_list[[day]] <- update_water_veg(veg_list[[day-1]], rainfall)
        temp_list[[day]] <- update_temp_rainfall(temp_list[[day-1]], rainfall, temp_kelvin, rel_hum)
        day <- day+1
      }
      
      this_dt[, water_veg:= veg_list]
      this_dt[, temp_rainfall:= temp_list]
      
      dt_list[[idx]] <- this_dt
      idx <- idx +1
    }
  }
}

full_dt <- rbindlist(dt_list)
# full_dt <- melt(full_dt, id.vars = c("day", "rainfall", "temperature", "rel_humidity"),
#                 value.name = "larval_capacity",
#                 variable.name = "habitat")

full_dt[, ratio := water_veg/temp_rainfall]

ggplot(full_dt[day==30], aes(x=rel_humidity, y=temp_rainfall, color=factor(temperature), shape=factor(temperature))) +
  geom_hline(aes(yintercept=water_veg)) +
  geom_point() +
  facet_grid(.~rainfall)
  
ggplot(full_dt[day==30], aes(x=rel_humidity, y=ratio, color=factor(temperature), shape=factor(rainfall))) +
 #geom_hline(aes(yintercept=water_veg)) +
  geom_point() 
