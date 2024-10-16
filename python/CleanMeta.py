import pandas as pd
import numpy as np

def CleanMeta(raw, plate, replicate, split_content=False, split_by="_", split_into=None, del_na=True):
    """
    Get Clean Metadata.
    
    This function processes raw data, plate layout, and replicate information 
    to create a clean metadata dataframe. It can optionally split the content column into additional columns.
    
    :param raw: A DataFrame containing the raw data.
    :param plate: A DataFrame containing the plate layout information.
    :param replicate: A DataFrame containing the replicate information.
    :param split_content: Boolean, whether to split the content. Default is False.
    :param split_by: A string to split the content by. Default is "_".
    :param split_into: A list specifying names for the split columns. Required if split_content is True.
    :param del_na: Boolean, whether to drop rows containing NA. Default is True.
    :return: A DataFrame containing the cleaned metadata.
    """
    if split_content and (split_into is None or len(split_into) == 0):
        raise ValueError("If split_content is True, split_into must be provided and cannot be empty.")
    
    n_platecol = 13 if plate.shape[1] == 13 else 25
    plate_format = 96 if n_platecol == 13 else 384
    
    replicate = replicate.iloc[:, 1:].values.flatten()
    
    def generate_wells(rows, cols):
        wells = [f"{r}{c}" for r in rows for c in cols]
        return wells
    
    if plate_format == 96:
        rows = list("ABCDEFGH")
        cols = [f"{i:02}" for i in range(1, 13)]
        n_row, n_col = 8, 12
    elif plate_format == 384:
        rows = list("ABCDEFGHIJKLMNOP")
        cols = [f"{i:02}" for i in range(1, 25)]
        n_row, n_col = 16, 24
    else:
        raise ValueError("Invalid format. Must be either 96 or 384.")
    
    well = generate_wells(rows, cols)
    
    content = plate.iloc[:, 1:n_platecol].values.flatten()
    
    if del_na:
        valid_well = ~pd.isna(replicate)
        well = np.array(well)[valid_well]
        content = content[valid_well]
        replicate = replicate[valid_well]
    
    content_replicate = np.where(
        (pd.isna(content)) | (pd.isna(replicate)),
        np.nan,
        content + "_" + replicate
    )
    
    meta = pd.DataFrame({
        'well': well,
        'content': content,
        'replicate': replicate,
        'content_replicate': content_replicate,
        'format': plate_format
    })
    
    if split_content:
        split_df = meta['content'].str.split(split_by, expand=True)
        
        if split_df.shape[1] != len(split_into):
            raise ValueError(f"Number of split columns ({split_df.shape[1]}) does not match the length of 'split_into' ({len(split_into)}).")
        
        split_df.columns = split_into
        meta = pd.concat([meta, split_df], axis=1)
    
    return meta


