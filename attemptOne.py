import rpy2.robjects as ro

# Source the R script
ro.r['source']('./tutorials/tutorialAttempt.R')

exitInput = ""

while(exitInput != "Stop"):
    exitInput = input("Type 'Stop' to end the program: ")


