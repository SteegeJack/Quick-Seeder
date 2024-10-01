import rpy2.robjects as ro
import pandas as pd
from rpy2.robjects.packages import importr
from rpy2.robjects import r, pandas2ri, StrVector, default_converter, conversion

# Activate automatic conversion between R and pandas DataFrames
pandas2ri.activate()

tidyverse = importr('tidyverse')
quicseedr = importr('QuICSeedR')
knitr = importr('knitr')
readxl = importr('readxl')
plotly = importr('plotly')

#Load data from sample 96-well plate
plate_96 = readxl.read_xlsx('./tutorials/data/20240716_s1_plate.xlsx')
raw_96 = readxl.read_xlsx('./tutorials/data/20240716_s1_raw.xlsx')

#Load data from sample 384-well plate
plate_384 = readxl.read_xlsx('./tutorials/data/20230808_M12_plate.xlsx')
raw_384 = readxl.read_xlsx('./tutorials/data/20230808_M12_raw.xlsx')

plate_time_96 = quicseedr.ConvertTime(raw_96)
plate_time_384 = quicseedr.ConvertTime(raw_384)

print(type(raw_96))
print(type(plate_time_96))

print(quicseedr.PlotPlate(raw_96, plate_time_96))


ro.r('''
CleanRaw = function (meta, raw, plate_time, cycle_total) 
{
    if (missing(cycle_total) || is.null(cycle_total) || length(cycle_total) == 
        0) {
        cycle_total = nrow(raw) - 1
    }
    else {
        cycle_total <- cycle_total
    }
    raw = raw[-1, -c(1:2)]
    raw = raw[1:cycle_total, meta$well]
    raw = as.numeric(unlist(raw))
    raw = matrix(raw, nrow = cycle_total)
    rownames(raw) = unlist(plate_time[1:cycle_total, ])
    colnames(raw) = meta$content_replicate
    return(raw)
}
''')

CleanRaw = ro.r['CleanRaw']


replicate_96 = quicseedr.GetReplicate(plate_96)
replicate_384 = quicseedr.GetReplicate(plate_384)

meta_96 = quicseedr.CleanMeta(raw_96, plate_96, replicate_96)
meta_384 = quicseedr.CleanMeta(raw_384, plate_384, replicate_384)

print(type(meta_384))
print(type(replicate_96))


with conversion.localconverter(default_converter):
    result = quicseedr.CleanRaw(meta_384, raw_384, plate_time_384)
    print(quicseedr.PlotRawSingle(result, "pos"))

    exitInput = ""

    while(exitInput != "Stop"):
        exitInput = input("Type 'Stop' to end the program: ")


    


cleanraw_96 = quicseedr.CleanRaw(meta_96, raw_96, plate_time_96)
cleanraw_384 = quicseedr.CleanRaw(meta_384, raw_384, plate_time_384)
