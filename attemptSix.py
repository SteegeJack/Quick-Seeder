import numpy as np
import rpy2.robjects as ro
import pandas as pd
import os
from rpy2.robjects.packages import importr
from rpy2.robjects import r, pandas2ri, StrVector, default_converter, conversion

def CleanMeta(rawData, plateData, replicateData, split_content=False, split_by="_", split_into=None, del_na=True):
    raw = pandas2ri.rpy2py(rawData)
    plate = pandas2ri.rpy2py(plateData)
    replicate = pandas2ri.rpy2py(replicateData)
    
    # Error checking
    if split_content and (split_into is None or len(split_into) == 0):
        raise ValueError("If split_content is True, split_into must be provided and cannot be empty.")
    
    n_platecol = 13 if plate.shape[1] == 13 else 25
    plate_format = 96 if n_platecol == 13 else 384

    replicate = replicate.iloc[:, 1:].values.flatten()
    
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
    
    well = generate_wells(rows, cols)
    
    content = plate.iloc[:, 1:n_platecol].values.flatten()

    if del_na:
        valid_indices = ~np.isnan(replicate)
        well = np.array(well)[valid_indices]
        content = content[valid_indices]
        replicate = replicate[valid_indices]

    # Combine content and replicate into a single string, making values treated as integers
    content_replicate = [
    f"{int(c) if isinstance(c, float) and c.is_integer() else c}_{int(r) if isinstance(r, float) and r.is_integer() else r}" if not pd.isna(c) and not pd.isna(r) else np.nan for c, r in zip(content, replicate)]

    
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
    
    return pandas2ri.py2rpy(meta)


def CleanRaw(metaData, rawData, plate_timeData, cycle_total=None):
    meta = pandas2ri.rpy2py(metaData)
    raw = pandas2ri.rpy2py(rawData)
    plate_time = pandas2ri.rpy2py(plate_timeData)

    if cycle_total is None or cycle_total == 0:
        cycle_total = raw.shape[0] - 1

    # Remove the first row and the first two columns
    raw = raw.iloc[1:, 2:]

    # Select only the cycles specified by cycle_total and the wells from the metadata
    raw = raw.iloc[:cycle_total].loc[:, meta['well']]

    # Flatten the raw data into a 1D array and convert it back to a matrix with the specified number of rows
    raw = np.array(raw).flatten().astype(float)
    raw = raw.reshape(cycle_total, -1)

    row_names = plate_time.iloc[:cycle_total, 0].values
    cleaned_raw = pd.DataFrame(raw, index=row_names, columns=meta['content_replicate'])

    return cleaned_raw


def ConvertTime(rawData):
    raw = pandas2ri.rpy2py(rawData)

    # Extract the time column, skipping the first row
    time = raw.iloc[1:, 1].copy()

    if time.empty:
        raise ValueError("Extracted time data is empty. Check the input DataFrame.")

    time_df = pd.DataFrame({'time': time})

    # Convert the column to strings
    time_df['time'] = time_df['time'].astype(str)
    
    # Handle any NaN values in str.contains
    if time_df['time'].str.contains('min', na=False).any():

        # Extract hours and minutes using regex
        hours = time_df['time'].str.extract(r'(\d+)\s*h')[0].astype(float).fillna(0)
        minutes = time_df['time'].str.extract(r'(\d+)\s*min')[0].fillna(0).astype(float)

        # Calculate decimal hours
        decimal_hours = hours + (minutes / 60)
        time_df['time'] = decimal_hours
    else:
        time_df['time'] = pd.to_numeric(time_df['time'], errors='coerce')

    return pandas2ri.py2rpy(time_df)


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

        # Maybe need to convert raw, plate, and replicate to pandas,
        # but I don't think so
        
        log(f"Processing plate {j + 1}")
        
        plate_time = ConvertTime(raw, **(params.get('ConvertTime', {})))
        
        meta = CleanMeta(raw=raw, plate=plate, replicate=replicate, **(params.get('CleanMeta', {})))
        
        clean_raw_params = params.get('CleanRaw', {})

        # CleanRaw
        try:
            raw = CleanRaw(meta=meta, raw=raw, plate_time=plate_time, **clean_raw_params)
        except Exception as e:
            log(f"Error in CleanRaw for plate {j + 1}: {str(e)}")
            continue
        
        if raw is None:
            log(f"Skipping further processing for plate {j + 1}")
            continue
        
        log(f"Dimensions of cleaned raw: {raw.shape}")
        
        # GetCalculation
        try:
            calculation = GetCalculation(raw=raw, meta=meta, **(params.get('GetCalculation', {})))
        except Exception as e:
            log(f"Error in GetCalculation for plate {j + 1}: {str(e)}")
            continue
        
        if calculation is None:
            log(f"Skipping further processing for plate {j + 1}")
            continue
        
        # SpreadCalculation and GetAnalysis
        if do_analysis:
            calculation_spread = SpreadCalculation(calculation, **(params.get('SpreadCalculation', {})))
            analysis = GetAnalysis(calculation_spread, **(params.get('GetAnalysis', {})))
        
        # SummarizeResult
        try:
            result = SummarizeResult(
                analysis=analysis if do_analysis else None, 
                calculation=calculation, 
                **(params.get('SummarizeResult', {}))
            )
        except Exception as e:
            log(f"Error in SummarizeResult for plate {j + 1}: {str(e)}")
            continue
        
        if result is not None:
            subcalculation[j] = calculation
            subcleanraw[j] = raw
            subresult[j] = result
    
    # Filter out any None entries from the results
    subcalculation = {k: v for k, v in subcalculation.items() if v is not None}
    subcleanraw = {k: v for k, v in subcleanraw.items() if v is not None}
    subresult = {k: v for k, v in subresult.items() if v is not None}
    
    # Check if no plates were processed successfully
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

def BulkReadMARS(path, plate_subfix, raw_subfix, helper_func=None):

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
            # Get all folders in the specified path

        
        plate_path = os.path.join(folder_path, plate_files[0])
        raw_path = os.path.join(folder_path, raw_files[0])
        
        plate_data = pd.read_excel(plate_path)
        raw_data = pd.read_excel(raw_path)
        replicate_data = quicseedr.GetReplicate(plate_data) 
        
        # Apply helper function if provided
        if helper_func:
            plate_data = plate_data.apply(helper_func)
        
        mylist[folder] = {
            'plate': plate_data,
            'raw': raw_data,
            'replicate': replicate_data
        }
    
    return mylist


