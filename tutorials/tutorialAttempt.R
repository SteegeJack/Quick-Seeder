library(tidyverse)
library(QuICSeedR)
library(knitr)
library(readxl)

#Load data from sample 96-well plate
plate_96 = read_xlsx('./tutorials/data/20240716_s1_plate.xlsx')
raw_96 = read_xlsx('./tutorials/data/20240716_s1_raw.xlsx')

#Load data from sample 384-well plate
plate_384 = read_xlsx('./tutorials/data/20230808_M12_plate.xlsx')
raw_384 = read_xlsx('./tutorials/data/20230808_M12_raw.xlsx')

plate_time_96 = ConvertTime(raw_96)
plate_time_384 = ConvertTime(raw_384)

PlotPlate(raw_96, plate_time_96)

replicate_96 = GetReplicate(plate_96)
replicate_384 = GetReplicate(plate_384)

meta_96 = CleanMeta(raw_96, plate_96, replicate_96)
meta_384 = CleanMeta(raw_384, plate_384, replicate_384)

cleanraw_96 = CleanRaw(meta_96, raw_96, plate_time_96)
cleanraw_384 = CleanRaw(meta_384, raw_384, plate_time_384)

print(class(cleanraw_384))

PlotRawSingle(cleanraw_96, "pos") 