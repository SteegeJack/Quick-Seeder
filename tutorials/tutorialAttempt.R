library(tidyverse)
library(QuICSeedR)
library(knitr)
library(readxl)

print("We have read the packages!")

#Load data from sample 96-well plate
plate_96 = read_xlsx('./tutorials/data/20240716_s1_plate.xlsx')
raw_96 = read_xlsx('./tutorials/data/20240716_s1_raw.xlsx')

#Load data from sample 384-well plate
plate_384 = read_xlsx('./tutorials/data/20230808_M12_plate.xlsx')
raw_384 = read_xlsx('./tutorials/data/20230808_M12_raw.xlsx')

plate_time_96 = ConvertTime(raw_96)
plate_time_384 = ConvertTime(raw_384)

print(PlotPlate(raw_96, plate_time_96))

replicate_96 = GetReplicate(plate_96)
replicate_384 = GetReplicate(plate_384)

meta_384 = CleanMeta(raw_384, plate_384, replicate_384)
