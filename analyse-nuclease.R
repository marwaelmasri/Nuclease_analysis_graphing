
# This script is for analysing nuclease assay data
# Sequence the nuclease assay products
# Export the sizing data from genemapper
# this script will calculate percent cut
# then export the data to graph in graphpad (my preference)


#install packages and load libraries
install.packages("readxl")
install.packages("tidyverse")
library(readxl)
library(dplyr)


#set working directory
#add your working directory here, use " " and use / not \
# eg setwd("C:/Users/documents/mydata")
setwd("")


#open file, add your own file name
peaks <- read_excel("peaks.xlsx")

#filter uncut peak
uncut <- peaks %>%
  filter(Size >= 45 & Size <= 55)

#filter cut peak
cut_endo <- peaks %>%
  filter(Size >=20 & Size <= 37)

#filter cut peak
cut_exo <- peaks %>%
  filter(Size <=13)

#sum area of uncut peak for each trace
uncut_areasum <- uncut %>% 
  group_by(`Sample File Name`) %>% 
  summarise(peak_area_uncut = sum(Area))

#if you have to delete any rows
#uncut_areasum <- uncut_areasum[-60,]

#sum area of cut peak for each trace
cutendo_areasum <- cut_endo %>% 
  group_by(`Sample File Name`) %>% 
  summarise(peak_area_cut_endo = sum(Area))

#sum area of cut peak for each trace
cutexo_areasum <- cut_exo %>% 
  group_by(`Sample File Name`) %>% 
  summarise(peak_area_cut_exo = sum(Area))

#join cut and uncut tables
percent_cut <- cbind(cutendo_areasum, cutexo_areasum, uncut_areasum)

#add percent cut column FOR ENDONUCLEASE
percent_cut$peak_endo_percent_cut <- (percent_cut$peak_area_cut_endo / (percent_cut$peak_area_cut_endo/percent_cut$peak_area_uncut))*100

#add percent cut column FOR EXONUCLEASE
percent_cut$peak_exo_percent_cut <- (percent_cut$peak_area_cut_exo / (percent_cut$peak_area_cut_exo/percent_cut$peak_area_uncut))*100

#export as csv
write.csv(percent_cut, "percent_cut.csv", row.names = FALSE)



