
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
setwd( #add your working directory here, use " " and use / not \
  # eg setwd("C:/Users/documents/mydata")
  )


#open file, add your own file name
peaks <- read_excel("peaks.xlsx")

#filter uncut peak
uncut <- peaks %>%
  filter(Size >= 45 & Size <= 55)

#filter cut peak
cut <- peaks %>%
  filter(Size >=27 & Size <= 37)

#sum area of uncut peak for each trace
uncut_areasum <- uncut %>% 
  group_by(`Sample File Name`) %>% 
  summarise(peak_area_uncut = sum(Area))

#if you have to delete any rows
#uncut_areasum <- uncut_areasum[-60,]

#sum area of cut peak for each trace
cut_areasum <- cut %>% 
  group_by(`Sample File Name`) %>% 
  summarise(peak_area_cut = sum(Area))

#join cut and uncut tables
percent_cut <- cbind(cut_areasum, uncut_areasum)

#add percent cut column
percent_cut$peak_percent_cut <- (percent_cut$peak_area_cut / percent_cut$peak_area_uncut)*100

#export as csv
write.csv(percent_cut, "percent_cut.csv", row.names = FALSE)

