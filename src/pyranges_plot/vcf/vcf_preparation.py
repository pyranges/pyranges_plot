import pandas as pd


def split_fields(
    data,
    target_cols: str | list,
    field_sep: str,
    col_name_sep: str | None = None,
    col_names: list[str] | None = None,
    col_types: list[str] | None = None,
    keep_col: bool = False,
):
    """
    Splits a column or columns into multiple columns based on specified separators.

    Parameters
    ----------
    data: pd.DataFrame
            The input DataFrame containing the columns to be split.
    target_cols: {str or list of strings}
            Column name(s) in the DataFrame to be split. Can be a single column (str)
            or a list of column names.
    field_sep: str
            Separator used to split the fields in the target column(s).
    col_name_sep: str, default None
            If provided, this separator is used to split each field into a column name
            and value. For example, `"key=value"` will generate a column named `key`
            with the corresponding `value`. Defaults to None.
    col_names: list[str], default None
            A list of names for the new columns. If not provided, column names are
            generated automatically based on the target column name and field index.
            If `col_name_sep` is specified, the column names can be inferred from the field keys.
            Defaults to None.
    col_types: list[str], default None
            A list of data types for the new columns. If not provided, columns will
            retain their default inferred types. Defaults to None.
    keep_col: bool , default False
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

    Example
    -------
    >>> vcf = prp.example_data.ncbi_vcf()
    >>> vcf
    index    |    Chromosome    Start     ID            REF       ALT       QUAL      FILTER      ...
    int64    |    object        int32     object        object    object    object    category    ...
    -------  ---  ------------  --------  ------------  --------  --------  --------  ----------  -----
    0        |    1             943995    rs761448939   C         G,T       nan       .           ...
    1        |    1             964512    rs756054473   C         A,T       nan       .           ...
    2        |    1             976215    rs7417106     A         C,G,T     nan       .           ...
    3        |    1             1013983   rs1644247121  G         A         nan       .           ...
    ...      |    ...           ...       ...           ...       ...       ...       ...         ...
    242182   |    Y             2787592   rs104894975   A         T         nan       .           ...
    242183   |    Y             2787600   rs104894977   G         A         nan       .           ...
    242184   |    Y             7063898   rs199659121   A         T         nan       .           ...
    242185   |    Y             12735725  rs778145751   TAAGT     T         nan       .           ...
    PyRanges with 242186 rows, 9 columns, and 1 index columns. (2 columns not shown: "INFO", "End").
    Contains 25 chromosomes.
    >>> prp.vcf.split_fields(vcf,target_cols="INFO",field_sep=";",col_name_sep="=")
    index    |    Chromosome    Start     ID            REF       ALT       QUAL      FILTER      End       INFO_0     TSA       INFO_2                  INFO_3                  ...
    int64    |    object        int32     object        object    object    object    category    int32     object     object    object                  object                  ...
    -------  ---  ------------  --------  ------------  --------  --------  --------  ----------  --------  ---------  --------  ----------------------  ----------------------  -----
    0        |    1             943995    rs761448939   C         G,T       nan       .           943996    dbSNP_156  SNV       E_Freq                  E_Cited                 ...
    1        |    1             964512    rs756054473   C         A,T       nan       .           964513    dbSNP_156  SNV       E_Freq                  E_Cited                 ...
    2        |    1             976215    rs7417106     A         C,G,T     nan       .           976216    dbSNP_156  SNV       E_Freq                  E_1000G                 ...
    3        |    1             1013983   rs1644247121  G         A         nan       .           1013984   dbSNP_156  SNV       E_Phenotype_or_Disease  CLIN_pathogenic         ...
    ...      |    ...           ...       ...           ...       ...       ...       ...         ...       ...        ...       ...                     ...                     ...
    242182   |    Y             2787592   rs104894975   A         T         nan       .           2787593   dbSNP_156  SNV       E_Cited                 E_Phenotype_or_Disease  ...
    242183   |    Y             2787600   rs104894977   G         A         nan       .           2787601   dbSNP_156  SNV       E_Cited                 E_Phenotype_or_Disease  ...
    242184   |    Y             7063898   rs199659121   A         T         nan       .           7063899   dbSNP_156  SNV       E_Freq                  E_Cited                 ...
    242185   |    Y             12735725  rs778145751   TAAGT     T         nan       .           12735726  dbSNP_156  indel     E_Freq                  E_Cited                 ...
    PyRanges with 242186 rows, 31 columns, and 1 index columns. (19 columns not shown: "INFO_4", "INFO_5", "INFO_6", ...).
    Contains 25 chromosomes.
    """
    result_data = data.copy()
    if isinstance(target_cols, str):
        target_cols = [target_cols]

    for target_col in target_cols:
        if target_col not in data.columns:
            raise ValueError(
                f"Target column '{target_col}' does not exist in the DataFrame."
            )
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
