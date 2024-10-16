import numpy as np
import rpy2.robjects as ro
import pandas as pd
from rpy2.robjects.packages import importr
from rpy2.robjects import r, pandas2ri, StrVector, default_converter, conversion

def CleanMeta(raw, plate, replicate, split_content=False, split_by="_", split_into=None, del_na=True):
    # Error checking
    if split_content and (split_into is None or len(split_into) == 0):
        raise ValueError("If split_content is True, split_into must be provided and cannot be empty.")
    
    # Figure out plate format
    n_platecol = 13 if plate.shape[1] == 13 else 25
    plate_format = 96 if n_platecol == 13 else 384

    # Flatten replicate data, ignoring the first column
    replicate = replicate.iloc[:, 1:].values.flatten()
    
    # Function to generate wells
    def generate_wells(rows, cols):
        return [f"{r}{c}" for r in rows for c in cols]
    
    # Set rows and columns based on plate format
    if plate_format == 96:
        rows = [chr(i) for i in range(65, 73)]  # A-H
        cols = [f"{i:02d}" for i in range(1, 13)]
        n_row, n_col = 8, 12
    elif plate_format == 384:
        rows = [chr(i) for i in range(65, 81)]  # A-P
        cols = [f"{i:02d}" for i in range(1, 25)]
        n_row, n_col = 16, 24
    else:
        raise ValueError("Invalid format. Must be either 96 or 384.")
    
    # Generate well IDs
    well = generate_wells(rows, cols)
    
    # Flatten content from plate, ignoring the first column
    content = plate.iloc[:, 1:n_platecol].values.flatten()

    # Remove NA values if del_na is True
    if del_na:
        valid_indices = ~np.isnan(replicate)
        well = np.array(well)[valid_indices]
        content = content[valid_indices]
        replicate = replicate[valid_indices]
    
    # Combine content and replicate into a single string
    content_replicate = [f"{c}_{r}" if not pd.isna(c) and not pd.isna(r) else np.nan 
                         for c, r in zip(content, replicate)]
    
    # Create a DataFrame for meta data
    meta = pd.DataFrame({
        'well': well,
        'content': content,
        'replicate': replicate,
        'content_replicate': content_replicate,
        'format': plate_format
    })
    
    # Split the content into additional columns if requested
    if split_content:
        split_df = meta['content'].str.split(split_by, expand=True)
        
        if split_df.shape[1] != len(split_into):
            raise ValueError(f"Number of split columns ({split_df.shape[1]}) does not match the length of 'split_into' ({len(split_into)}).")
        
        split_df.columns = split_into
        meta = pd.concat([meta, split_df], axis=1)
    
    return meta


def CleanRaw(meta, raw, plate_time, cycle_total=None):
    # Determine cycle_total if not provided
    if cycle_total is None or len(cycle_total) == 0:
        cycle_total = raw.shape[0] - 1
    
    # Remove the first row and first two columns from raw
    raw = raw.iloc[1:, 2:]
    
    # Select only columns corresponding to wells in meta
    raw = raw.loc[:cycle_total - 1, meta['well']]
    
    # Convert the raw data to a numeric numpy array
    raw = pd.to_numeric(raw.values.flatten(), errors='coerce')
    
    # Reshape the array into a matrix with the specified number of cycles (rows)
    raw = raw.reshape(cycle_total, -1)
    
    # Set row names using the first cycle_total rows of plate_time
    raw_df = pd.DataFrame(raw, index=plate_time.iloc[:cycle_total, 0])
    raw_df.columns = meta['content_replicate']
    
    return raw_df


def ConvertTime(raw):
    # Extract the time column, skipping the first row
    time = raw.iloc[1:, 1].copy()
    
    # Convert the time to a DataFrame for consistency
    time = pd.DataFrame(time, columns=['.'])

    # Check if the time format contains "min"
    if time['.'].str.contains('min').any():
        # Extract hours and minutes using regex
        hours = time['.'].str.extract(r'(\d+)\s*h')[0].astype(float)
        minutes = time['.'].str.extract(r'(\d+)\s*min')[0].fillna(0).astype(float)
        
        # Calculate decimal hours
        decimal_hours = hours + (minutes / 60)
        time['.'] = decimal_hours
    else:
        # Convert directly to numeric if in decimal format
        time['.'] = pd.to_numeric(time['.'], errors='coerce')
    
    return time


def BulkProcessing(data, do_analysis=True, params=None, verbose=False):
    if params is None:
        params = {}
    
    subcalculation = {}
    subcleanraw = {}
    subresult = {}
    
    def log(*args):
        if verbose:
            print(*args)
    
    # Loop through each dataset in data
    for j, experiment in enumerate(data):
        plate = experiment['plate']
        raw = experiment['raw']
        replicate = experiment['replicate']
        
        log(f"Processing plate {j + 1}")
        log(f"Dimensions of raw: {raw.shape}")
        
        # Time conversion
        plate_time = convert_time(raw, **(params.get('ConvertTime', {})))
        
        # Metadata cleaning
        meta = clean_meta(raw=raw, plate=plate, replicate=replicate, **(params.get('CleanMeta', {})))
        
        log(f"Dimensions of meta: {meta.shape}")
        log(f"Dimensions of plate_time: {plate_time.shape}")
        
        # Raw data cleaning
        clean_raw_params = params.get('CleanRaw', {})
        try:
            raw = clean_raw(meta=meta, raw=raw, plate_time=plate_time, **clean_raw_params)
        except Exception as e:
            log(f"Error in CleanRaw for plate {j + 1}: {str(e)}")
            continue
        
        if raw is None:
            log(f"Skipping further processing for plate {j + 1}")
            continue
        
        log(f"Dimensions of cleaned raw: {raw.shape}")
        
        # Calculation
        try:
            calculation = get_calculation(raw=raw, meta=meta, **(params.get('GetCalculation', {})))
        except Exception as e:
            log(f"Error in GetCalculation for plate {j + 1}: {str(e)}")
            continue
        
        if calculation is None:
            log(f"Skipping further processing for plate {j + 1}")
            continue
        
        # Analysis (if enabled)
        if do_analysis:
            calculation_spread = spread_calculation(calculation, **(params.get('SpreadCalculation', {})))
            analysis = get_analysis(calculation_spread, **(params.get('GetAnalysis', {})))
        
        # Summarize results
        try:
            result = summarize_result(
                analysis=analysis if do_analysis else None, 
                calculation=calculation, 
                **(params.get('SummarizeResult', {}))
            )
        except Exception as e:
            log(f"Error in SummarizeResult for plate {j + 1}: {str(e)}")
            continue
        
        # Store the results if no errors occurred
        if result is not None:
            subcalculation[j] = calculation
            subcleanraw[j] = raw
            subresult[j] = result
    
    # Filter out any None entries from the results
    subcalculation = {k: v for k, v in subcalculation.items() if v is not None}
    subcleanraw = {k: v for k, v in subcleanraw.items() if v is not None}
    subresult = {k: v for k, v in subresult.items() if v is not None}
    
    # Check if any plates were processed successfully
    if len(subcalculation) == 0:
        print("Warning: No plates were successfully processed.")
        return None
    
    # Combine results
    combined_calculation = pd.concat([v.assign(plate_name=k) for k, v in subcalculation.items()], ignore_index=True)
    combined_result = pd.concat([v.assign(plate_name=k) for k, v in subresult.items()], ignore_index=True)
    
    return {
        'combined_calculation': combined_calculation,
        'combined_cleanraw': subcleanraw,
        'combined_result': combined_result
    }

import os

def BulkReadMARS(path, plate_subfix, raw_subfix, helper_func=None):
    """
    Bulk Read MARS data.
    
    Args:
        path (str): Path to the parent directory containing the data folders.
        plate_subfix (str): Suffix used to identify plate data files.
        raw_subfix (str): Suffix used to identify raw data files.
        helper_func (callable, optional): Function to apply to each column of the plate data (default: None).
        
    Returns:
        list: A list containing data from each folder, including plate, raw, and replicate data.
    """
    # Get all folders in the specified path
    folders = [f for f in os.listdir(path) if os.path.isdir(os.path.join(path, f))]
    
    mylist = {}
    
    for folder in folders:
        folder_path = os.path.join(path, folder)
        
        # Get all Excel files in the folder
        files = [f for f in os.listdir(folder_path) if f.endswith('.xlsx')]
        
        # Identify plate and raw data files
        plate_files = [f for f in files if plate_subfix in f]
        raw_files = [f for f in files if raw_subfix in f]
        
        if not plate_files or not raw_files:
            print(f"Warning: Skipping folder {folder_path} due to missing files.")
            continue
        
        plate_path = os.path.join(folder_path, plate_files[0])
        raw_path = os.path.join(folder_path, raw_files[0])
        
        # Read the Excel files
        plate_data = pd.read_excel(plate_path)
        raw_data = pd.read_excel(raw_path)
        
        # Assume GetReplicate is defined elsewhere; you need to implement it
        replicate_data = get_replicate(plate_data)  # Replace with actual implementation
        
        # Apply helper function if provided
        if helper_func:
            plate_data = plate_data.apply(helper_func)
        
        mylist[folder] = {
            'plate': plate_data,
            'raw': raw_data,
            'replicate': replicate_data
        }
    
    return mylist


# Activate automatic conversion between R and pandas DataFrames
pandas2ri.activate()

tidyverse = importr('tidyverse')
quicseedr = importr('QuICSeedR')
knitr = importr('knitr')
readxl = importr('readxl')
plotly = importr('plotly')

#Load data from sample 96-well plate
raw_96_r = readxl.read_xlsx('./tutorials/data/20240716_s1_raw.xlsx')
plate_96_r = readxl.read_xlsx('./tutorials/data/20240716_s1_plate.xlsx')

#Load data from sample 384-well plate
plate_384_r = readxl.read_xlsx('./tutorials/data/20230808_M12_plate.xlsx')
raw_384_r = readxl.read_xlsx('./tutorials/data/20230808_M12_raw.xlsx')

raw_96 = pandas2ri.rpy2py(raw_96_r)
plate_96 = pandas2ri.rpy2py(plate_96_r)
# raw_384 = pandas2ri.rpy2py(raw_384)
# plate_384 = pandas2ri.rpy2py(plate_384)

print(type(raw_96))

plate_time_96 = ConvertTime(raw_96)
# plate_time_384 = ConvertTime(raw_384)

plate_time_96_alt = quicseedr.ConvertTime(raw_96_r)
# plate_time_384 = quicseedr.ConvertTime(raw_384)

print("hello")

plate_time_96_r = pandas2ri.py2rpy(plate_time_96)

print(type(raw_96_r))
print(type(plate_time_96_r))


print(quicseedr.PlotPlate(raw_96_r, plate_time_96_r))

