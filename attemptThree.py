import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects import r, pandas2ri, StrVector

# Activate automatic conversion between R and pandas DataFrames
pandas2ri.activate()

tidyverse = importr('tidyverse')
quicseedr = importr('QuICSeedR')
knitr = importr('knitr')
readxl = importr('readxl')
plotly = importr('plotly')

plate_96 = readxl.read_xlsx('./tutorials/data/20240716_s1_plate.xlsx')
raw_96 = readxl.read_xlsx('./tutorials/data/20240716_s1_raw.xlsx')

plate_384 = readxl.read_xlsx('./tutorials/data/20230808_M12_plate.xlsx')
raw_384 = readxl.read_xlsx('./tutorials/data/20230808_M12_raw.xlsx')

plate_time_96 = quicseedr.ConvertTime(raw_96)

plate_time_384 = quicseedr.ConvertTime(raw_384)

print(quicseedr.PlotPlate(raw_96, plate_time_96))

replicate_96 = quicseedr.GetReplicate(plate_96)
replicate_384 = quicseedr.GetReplicate(plate_384)

meta_96 = quicseedr.CleanMeta(raw_96, plate_96, replicate_96)
meta_384 = quicseedr.CleanMeta(raw_384, plate_384, replicate_384)

cleanraw_96 = quicseedr.CleanRaw(meta_96, raw_96, plate_time_96)
cleanraw_384 = quicseedr.CleanRaw(meta_384, raw_384, plate_time_384)

print(quicseedr.PlotRawMulti(cleanraw_96, StrVector(["Neg", "Pos"])))
# print(quicseedr.PlotRawSingle(cleanraw_384, "pos"))

# calculation_384 = quicseedr.GetCalculation(cleanraw_384, meta_384, norm = True, norm_ct = 'pos')

# print(quicseedr.PlotMetric(calculation_384, y = "MS"))





# plate = readxl.read_xlsx('./tutorials/data/20230808_M12_plate.xlsx')
# plate_df = pandas2ri.rpy2py(plate)
# print(plate_df.head())

# # Load your R script
# ro.r['source']('./R/CleanMeta.R')

# # Access an R function
# CleanMeta = ro.globalenv['CleanMeta']

