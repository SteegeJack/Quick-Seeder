import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects import r, pandas2ri

# Activate automatic conversion between R and pandas DataFrames
pandas2ri.activate()

# Load your R script
ro.r['source']('./R/TestFunction.R')

# Access an R function
firstFunction = ro.globalenv['TestFunction']

# Load your R script
ro.r['source']('./R/TestFunction2.R')

# Access an R function
secondFunction = ro.globalenv['TestFunction']

firstFunction()
secondFunction()