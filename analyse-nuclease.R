
# This script is for analysing nuclease assay data
# Sequence the nuclease assay products
# Export the sizing data from genemapper
# this script will calculate percent cut
# then export the data to graph in graphpad (my preference)


#install packages and load libraries

install.packages("tidyverse")
library(dplyr)
library(readr)

#set working directory
#add your working directory here, use " " and use / not \
# eg setwd("C:/Users/documents/mydata")
setwd("")


#open file, add your own file name
peaks_uncut <- read_tsv("uncut.txt")
peaks_cut <- read_tsv("cut.txt")


#filter uncut peak
uncut <- peaks_uncut %>%
  filter(Size >= 45 & Size <= 55)

#filter cut peak
cut_endo <- peaks_cut %>%
  filter(Size >=20 & Size <= 37)

#filter cut peak
cut_exo <- peaks_cut %>%
  filter(Size <=13)

#sum area of uncut peak for each trace
uncut_areasum <- uncut %>% 
  group_by(`Sample File Name`) %>% 
  summarise(peak_area_uncut = sum(Area))


#sum area of cut peak for each trace
cutendo_areasum <- cut_endo %>% 
  group_by(`Sample File Name`) %>% 
  summarise(peak_area_cut_endo = sum(Area))

#sum area of cut peak for each trace
cutexo_areasum <- cut_exo %>% 
  group_by(`Sample File Name`) %>% 
  summarise(peak_area_cut_exo = sum(Area))

#join cut and uncut tables
percent_cut <- full_join(uncut_areasum, cutendo_areasum, by = "Sample File Name")
percent_cut <- full_join(percent_cut, cutexo_areasum, by = "Sample File Name")

# change na to 0
percent_cut[is.na(percent_cut)] <- 0

#add percent cut column FOR ENDONUCLEASE
percent_cut$peak_endo_percent_cut <- (percent_cut$peak_area_cut_endo / (percent_cut$peak_area_cut_endo+percent_cut$peak_area_uncut))*100

#add percent cut column FOR EXONUCLEASE
percent_cut$peak_exo_percent_cut <- (percent_cut$peak_area_cut_exo / (percent_cut$peak_area_cut_exo+percent_cut$peak_area_uncut))*100

#export as csv
write.csv(percent_cut, "percent_cut.csv", row.names = FALSE)




