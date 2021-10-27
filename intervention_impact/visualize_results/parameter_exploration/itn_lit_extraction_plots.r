###############################################################################################################
## itn_lit_extraction_plots.r
## Amelia Bertozzi-Villa
## October 2021
## 
##############################################################################################################


library(data.table)
library(ggplot2)

rm(list=ls())

main_dir <- file.path(Sys.getenv("HOME"), 
                      "Dropbox (IDM)/Malaria Team Folder/projects/map_intervention_impact/writing_and_presentations/tropmed_2021")
out_dir <- file.path(main_dir, "lit_review_figs")
dir.create(out_dir, recursive = T, showWarnings = F)

main_lit <- fread(file.path(main_dir, "itn_impact_lit.csv"))
meta_lit <- fread(file.path(main_dir, "itn_impact_meta.csv"))

main_lit[, metric:=ifelse(outcome %like% "Parasitemia", "Prevalence", "Incidence") ]
main_lit <- main_lit[order(yrs_elapsed)]
main_lit[, row_id:=as.integer(rownames(main_lit))]
meta_lit[, metric:= capitalize(metric)]


pdf(file.path(out_dir, "int_impact.pdf"), width=6, height=4)
ggplot(main_lit) +
  geom_bar(aes(x=1, fill=yrs_elapsed, y=pct_red_morb, group=row_id), stat = "identity", position = position_dodge2(width = 0.3, preserve = "single")) +
  geom_hline(data=meta_lit[metric!="Mortality"], aes(yintercept=protective_efficacy_pct), color="red") +
  scale_fill_viridis_c("Years\nElapsed", direction=-1) +
  facet_grid(metric ~ .) +
  labs(x="Study",
       y="% Reduction") +
  theme(axis.text.x = element_blank())
graphics.off()
  
