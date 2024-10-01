import rpy2.robjects as ro
import pandas as pd
from rpy2.robjects.packages import importr
from rpy2.robjects import r, pandas2ri, StrVector, default_converter, conversion

# Activate automatic conversion between R and pandas DataFrames
pandas2ri.activate()

def CleanRaw(meta, raw, plate_time):
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

    with conversion.localconverter(default_converter):
        return ro.r['CleanRaw'](meta, raw, plate_time)

# Still doesn't have the split_into parameter
def CleanMeta(raw, plate, replicate, split_content = False, split_by = "_", del_na = True):
    ro.r('''
    CleanMeta <- function (raw, plate, replicate, split_content = FALSE, split_by = "_", split_into = c(),
                       del_na = TRUE) {
  
        if (split_content && length(split_into) == 0) {
            stop("If split_content is TRUE, split_into must be provided and cannot be empty.")
        }
        
        n_platecol <- ifelse(ncol(plate) == 13, 13, 25)
        plate_format <- ifelse(n_platecol == 13, 96, 384)
        
        replicate <- replicate[, -1]
        replicate <- c(t(replicate))
        
        generate_wells <- function(rows, cols) {
            wells <- character(length(rows) * length(cols))
            index <- 1
            for (r in rows) {
            for (c in cols) {
                wells[index] <- paste0(r, c)
                index <- index + 1
            }
            }
            return(wells)
        }
        
        if (plate_format == 96) {
            rows <- LETTERS[1:8]
            cols <- sprintf("%02d", 1:12)
            n_row <- 8
            n_col <- 12
        } else if (plate_format == 384) {
            rows <- LETTERS[1:16]
            cols <- sprintf("%02d", 1:24)
            n_row <- 16
            n_col <- 24
        } else {
            stop("Invalid format. Must be either 96 or 384.")
        }
        
        well <- generate_wells(rows, cols)  
        
        content <- plate[, 2:n_platecol]
        content <- c(t(content))
        
        if (del_na) {
            valid_well <- which(!is.na(replicate))
            well <- well[valid_well]
            content <- content[valid_well]
            replicate <- replicate[valid_well]
        }
        
        content_replicate <- paste(content, replicate, sep = "_")
        content_replicate = gsub ("NA_NA", NA, content_replicate)
        
        meta <- data.frame(
            well = well,
            content = content,
            replicate = replicate,
            content_replicate = content_replicate,
            format = plate_format,
            stringsAsFactors = FALSE
        )
        
        if (split_content) {
            split_df <- do.call(rbind, strsplit(as.character(meta$content), split_by))
            
            if (ncol(split_df) != length(split_into)) {
            stop(paste("Number of split columns (", ncol(split_df), 
                        ") does not match the length of 'split_into' (", 
                        length(split_into), ").", sep=""))
            colnames(split_df) <- paste0("split_", seq_len(ncol(split_df)))
            }
            
            colnames(split_df) <- split_into
            meta <- cbind(meta, split_df)
        }
        
        return(meta)
    }
    ''')

    return ro.r['CleanMeta'](raw, plate, replicate, split_content, split_by, del_na)

def ConvertTime(raw):
    ro.r('''
        ConvertTime = function (raw) {
            time = c(raw[-1, 2])
            time = unlist(time) %>% data.frame(stringsAsFactors = FALSE)
            if (grepl("min", time) == TRUE) {
                hours <- as.numeric(gsub(" h.*$", "", time$.))
                suppressWarnings({
                minutes <- ifelse(grepl("min", time$.), as.numeric(gsub("^.*?([0-9]+) min$", 
                                                                        "\\1", time$.)), 0)
                })
                decimal_hours <- hours + (minutes/60)
                time = (data.frame(. = decimal_hours))
            }
            else {
                time$. = as.numeric(as.character(time$.))
            }
            return(time)
        }
    ''')

    return ro.r['ConvertTime'](raw)

def BulkProcessing(data):
    ro.r('''
        BulkProcessing = function(data, do_analysis = TRUE, params = list(), verbose = FALSE) {
            subcalculation <- list()
            subcleanraw <- list()
            subresult <- list()
            
            log <- function(...) {
                if (verbose) cat(...)
            }
            
            for (j in 1:length(data)) {
                plate <- data[[j]]$plate
                raw <- data[[j]]$raw
                replicate <- data[[j]]$replicate
                
                log("Processing plate", j, "\n")
                log("Dimensions of raw:", dim(raw), "\n")
                
                plate_time <- do.call(ConvertTime, c(list(raw), params$ConvertTime %||% list()))
                meta <- do.call(CleanMeta, c(list(raw = raw, plate = plate, 
                                                replicate = replicate), params$CleanMeta %||% list()))
                
                log("Dimensions of meta:", dim(meta), "\n")
                log("Dimensions of plate_time:", dim(plate_time), "\n")
                
                clean_raw_params <- params$CleanRaw %||% list()
                
                raw <- tryCatch({
                do.call(CleanRaw, c(list(meta = meta, raw = raw, 
                                        plate_time = plate_time), clean_raw_params))
                }, error = function(e) {
                log("Error in CleanRaw for plate", j, ":", conditionMessage(e), "\n")
                return(NULL)
                })
                
                if (is.null(raw)) {
                log("Skipping further processing for plate", j, "\n")
                next
                }
                
                log("Dimensions of cleaned raw:", dim(raw), "\n")
                
                calculation <- tryCatch({
                do.call(GetCalculation, c(list(raw = raw, 
                                                meta = meta), params$GetCalculation %||% list()))
                }, error = function(e) {
                log("Error in GetCalculation for plate", j, ":", conditionMessage(e), "\n")
                return(NULL)
                })
                
                if (is.null(calculation)) {
                log("Skipping further processing for plate", j, "\n")
                next
                }
                
                if (do_analysis) {
                calculation_spread <- do.call(SpreadCalculation, 
                                                c(list(calculation), params$SpreadCalculation %||% list()))
                analysis <- do.call(GetAnalysis, c(list(calculation_spread), 
                                                    params$GetAnalysis %||% list()))
                }
                
                result <- tryCatch({
                if (do_analysis) {
                    do.call(SummarizeResult, c(list(analysis = analysis, 
                                                    calculation = calculation), params$SummarizeResult %||% list()))
                } else {
                    do.call(SummarizeResult, c(list(calculation = calculation), 
                                            params$SummarizeResult %||% list()))
                }
                }, error = function(e) {
                log("Error in SummarizeResult for plate", j, ":", conditionMessage(e), "\n")
                return(NULL)
                })
                
                if (!is.null(result)) {
                subcalculation[[j]] <- calculation
                subcleanraw[[j]] <- raw
                subresult[[j]] <- result
                }
            }
            
            subcalculation <- Filter(Negate(is.null), subcalculation)
            subcleanraw <- Filter(Negate(is.null), subcleanraw)
            subresult <- Filter(Negate(is.null), subresult)
            
            if (length(subcalculation) == 0) {
                warning("No plates were successfully processed.")
                return(NULL)
            }
            
            names(subcalculation) = names(data)[1:length(subcalculation)]
            names(subcleanraw) = names(data)[1:length(subcleanraw)]
            names(subresult) = names(data)[1:length(subresult)]
            
            subresult <- lapply(names(subresult), function(name) {
                data <- subresult[[name]]
                data$plate_name <- name
                data
            })
            fullresult = do.call(rbind, subresult)
            
            subcalculation <- lapply(names(subcalculation), function(name) {
                data <- subcalculation[[name]]
                data$plate_name <- name
                data
            })
            fullcalc = do.call(rbind, subcalculation)
            
            return(list(combined_calculation = fullcalc, combined_cleanraw = subcleanraw, 
                        combined_result = fullresult))
        }
    ''')

    return ro.r['BulkProcessing'](data)

def BulkReadMARS(path, plate_subfix, raw_subfix, helper_func = None):
    ro.r('''
        BulkReadMARS <- function(path, plate_subfix, raw_subfix, helper_func = NULL) {
  
            folders <- list.dirs(path = path, recursive = FALSE)
            
            mylist <- vector(mode = 'list', length = length(folders))
            
            listnames <- basename(folders)  
            names(mylist) <- listnames
            
            for (i in seq_along(folders)) { 
                folder <- folders[i]
                
                files <- list.files(path = folder, pattern = "\\.xlsx$", full.names = TRUE)
                
                plate_path <- files[grepl(plate_subfix, files, fixed = TRUE)]
                raw_path <- files[grepl(raw_subfix, files, fixed = TRUE)]
                
                if (length(plate_path) == 0 || length(raw_path) == 0) {
                warning(paste("Skipping folder", folder, "due to missing files."))
                next
                }
                
                plate_data <- read_xlsx(plate_path)
                raw_data <- read_xlsx(raw_path)
                replicate_data <- GetReplicate(plate_data)
                
                mylist[[i]] <- list(
                plate = if (is.null(helper_func)) plate_data else data.frame(lapply(plate_data, helper_func)),
                raw = raw_data,
                replicate = replicate_data
                )
            }
            mylist <- Filter(Negate(is.null), mylist)
            
            return(mylist)
        }
    ''')

    return ro.r['BulkReadMARS'](path, plate_subfix, raw_subfix, helper_func)

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
print(plate_time_96)

print(quicseedr.PlotPlate(raw_96, plate_time_96))

replicate_96 = quicseedr.GetReplicate(plate_96)
replicate_384 = quicseedr.GetReplicate(plate_384)

meta_96 = CleanMeta(raw_96, plate_96, replicate_96)
meta_384 = CleanMeta(raw_384, plate_384, replicate_384)

print(type(meta_384))
print(type(replicate_96))

result = CleanRaw(meta_384, raw_384, plate_time_384)
print(quicseedr.PlotRawSingle(result, "pos"))
exitInput = ""

while(exitInput != "Stop"):
    exitInput = input("Type 'Stop' to end the program: ")



grinder_data = quicseedr.BulkReadMARS(path = './data/grinder', plate_subfix = 'plate', raw_subfix = 'raw')

print(type(grinder_data))