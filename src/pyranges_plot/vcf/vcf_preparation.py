import pyranges as pr
import pandas as pd

def split_fields(data, target_cols: str | list, field_sep: str, col_name_sep: str | None=None, col_names: list[str] | None=None, col_types: list[str] | None=None, keep_col: bool=False):
    """
    Splits a column or columns into multiple columns based on specified separators.

    Parameters
    ----------
    data (pd.DataFrame): 
            The input DataFrame containing the columns to be split.
    target_cols (str | list): 
            Column name(s) in the DataFrame to be split. Can be a single column (str) 
            or a list of column names.
    field_sep (str): 
            Separator used to split the fields in the target column(s).
    col_name_sep (str | None, optional): 
            If provided, this separator is used to split each field into a column name 
            and value. For example, `"key=value"` will generate a column named `key` 
            with the corresponding `value`. Defaults to None.
    col_names (list[str] | None, optional): 
            A list of names for the new columns. If not provided, column names are 
            generated automatically based on the target column name and field index. 
            If `col_name_sep` is specified, the column names can be inferred from the field keys. 
            Defaults to None.
    col_types (list[str] | None, optional): 
            A list of data types for the new columns. If not provided, columns will 
            retain their default inferred types. Defaults to None.
    keep_col (bool, optional): 
            Whether to retain the original target column(s) in the output DataFrame. 
            Defaults to False (the original column(s) will be removed).
    
    Returns
    -------
        pd.DataFrame: 
            A DataFrame with the new columns added (and the target columns removed if `keep_col` is False).

    Raises
    ------
        ValueError: 
            If any specified `target_cols` are not present in the DataFrame.
        ValueError: 
            If the number of provided `col_names` does not match the number of new columns generated.
        ValueError: 
            If the number of provided `col_types` does not match the number of new columns generated.

    
    """
    result_data = data.copy()
    if isinstance(target_cols, str):
        target_cols = [target_cols]

    for target_col in target_cols:
        if target_col not in data.columns:
            raise ValueError(f"Target column '{target_col}' does not exist in the DataFrame.")
        target_series = result_data[target_col]
    
        # Initialize an empty dictionary to hold new columns
        new_columns = {}
    
        # Iterate over each row in the target column
        for idx, row in target_series.items():
            if pd.isna(row):
                continue  # Skip NaN values
    
            # Split the row by the field separator
            fields = row.split(field_sep)
    
            for field_idx, field in enumerate(fields):
                if col_name_sep and col_name_sep in field:
                    # Split by col_name_sep if present
                    column_name, value = field.split(col_name_sep, 1)
                else:
                    # Generate column name dynamically if no separator
                    column_name, value = f"{target_col}_{field_idx}", field
    
                # Ensure the column exists in the dictionary
                if column_name not in new_columns:
                    new_columns[column_name] = [None] * len(result_data)
    
                # Assign the value to the correct row
                new_columns[column_name][idx] = value

        # If col_names is provided, validate and apply it
        if col_names:
            if len(col_names) != len(new_columns):
                raise ValueError(
                    "The number of provided column names does not match the number of resulting columns."
                )
            new_columns = dict(zip(col_names, new_columns.values()))
    
        # Add the new columns to the DataFrame
        for col_name, values in new_columns.items():
            result_data[col_name] = values

        # Remove the original target column if keep_col is False
        if not keep_col:
            result_data.drop(columns=target_col, inplace=True)

    return result_data