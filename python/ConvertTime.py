import pandas as pd
import re

def ConvertTime(raw):
    """
    Extract and Convert Time Data to Decimal Hours.
    
    This function extracts and converts run time information from MARS output.
    
    :param raw: A DataFrame containing the MARS output.
    :return: A DataFrame containing the time information in decimal hours.
    """
    # Extract time column, skipping the header row
    time = raw.iloc[1:, 1].values.flatten()
    time = pd.DataFrame(time, columns=['.'])
    
    # Check if any values contain "min" (indicating "hours and minutes" format)
    if time['.'].str.contains("min", regex=False).any():
        # Extract hours and convert to numeric
        time['hours'] = time['.'].apply(lambda x: float(re.search(r"(\d+)\s*h", x).group(1)) if 'h' in x else 0)
        
        # Extract minutes and convert to numeric
        time['minutes'] = time['.'].apply(lambda x: float(re.search(r"(\d+)\s*min", x).group(1)) if 'min' in x else 0)
        
        # Calculate decimal hours
        time['.'] = time['hours'] + (time['minutes'] / 60)
        
        # Drop auxiliary columns used for calculation
        time = time[['.']]
    else:
        # Convert to numeric if it's already in decimal hours format
        time['.'] = pd.to_numeric(time['.'], errors='coerce')
    
    return time


df = pd.read_excel('./tutorials/data/20240716_s1_raw.xlsx')

print(df.head())