import pandas as pd
import pyranges as pr
from pathlib import Path
from io import StringIO


def read_vcf(f: str | Path, nrows: bool | None = None):
    """
    Read a VCF (Variant Call Format) file and convert it into a PyRanges object.

    This function processes a VCF file by reading the data, extracting the header and
    data lines, and creating a PyRanges object for genomic analysis. The metadata
    lines (lines starting with '##') are ignored, and the column names are extracted
    from the header line (starting with '#CHROM').

    Parameters
    ----------
    f : str | Path
        The file path to the VCF file to be read.
    nrows : bool | None, optional
        The number of rows to read from the file. If None, reads the entire file.

    Returns
    -------
    pr.PyRanges
        A PyRanges object containing the VCF data, adding the following columns:
        - Chromosome: Chromosome names (from 'CHROM' in the VCF).
        - Start: Start positions of variants (from 'POS' in the VCF).
        - End: End positions of variants (calculated as Start + 1).

    Raises
    ------
    FileNotFoundError
        If the provided file path does not exist.
    ValueError
        If the VCF file is malformed or missing essential fields.

    Notes
    -----
    - Missing quality scores ('.') are replaced with pandas.NA.
    - The function reads the file in chunks for large VCF files to handle memory usage.
    - Columns 'CHROM' and 'POS' are renamed to 'Chromosome' and 'Start' respectively,
      to align with PyRanges conventions.

    Examples
    --------
    >>> vcf_pyranges = prp.vcf.read_vcf("example.vcf")
    >>> vcf_ranges
    index    |    Chromosome    Start    ID          REF       ALT       QUAL      FILTER      INFO                       End
    int64    |    category      int32    category    object    object    object    category    object                     int32
    -------  ---  ------------  -------  ----------  --------  --------  --------  ----------  -------------------------  -------
    0        |    1             500      .           A         T         <NA>      PASS        TRANSCRIPT=t1;SECOND_ID=a  501
    1        |    1             3500     .           A         T         <NA>      PASS        TRANSCRIPT=t1;SECOND_ID=a  3501
    2        |    1             300      .           A         T         <NA>      PASS        TRANSCRIPT=t2;SECOND_ID=a  301
    3        |    1             1300     .           A         T         <NA>      PASS        TRANSCRIPT=t2;SECOND_ID=a  1301
    ...      |    ...           ...      ...         ...       ...       ...       ...         ...                        ...
    5        |    1             4500     .           A         T         <NA>      PASS        TRANSCRIPT=t3;SECOND_ID=b  4501
    6        |    1             4900     .           A         T         <NA>      PASS        TRANSCRIPT=t3;SECOND_ID=b  4901
    7        |    1             5600     .           A         T         <NA>      PASS        TRANSCRIPT=t3;SECOND_ID=b  5601
    8        |    1             6000     .           A         T         <NA>      PASS        TRANSCRIPT=t4;SECOND_ID=b  6001
    PyRanges with 9 rows, 9 columns, and 1 index columns.
    Contains 1 chromosomes.
    """

    path = Path(f)
    dtypes = {
        "CHROM": "category",
        "POS": "int32",
        "ID": "category",
        "REF": "str",
        "ALT": "str",
        "QUAL": "category",
        "FILTER": "category",
    }

    # Read the file to find the header line and ignore metadata
    with open(path, "r") as file:
        lines = []
        for line in file:
            if line.startswith("##"):
                continue  # Skip other metadata lines
            elif line.startswith("#CHROM"):
                # This is the header line; extract the column names
                column_names = line.strip("#").strip().split("\t")
            else:
                # Append data lines
                lines.append(line.strip())

    # Convert the list of lines into a DataFrame
    data = "\n".join(lines)

    df_iter = pd.read_csv(
        StringIO(data),
        sep="\t",
        comment="#",  # Ignore lines that start with '#'
        usecols=list(range(len(column_names))),  # Only read essential columns
        header=None,
        names=column_names,
        dtype=dtypes,
        chunksize=int(1e5),  # Use chunks for large files
        nrows=nrows,
    )

    # Process the chunks
    dfs = []
    for chunk in df_iter:
        # Replace missing quality scores represented by "."
        chunk["QUAL"] = chunk["QUAL"].astype(object).replace(".", pd.NA)
        dfs.append(chunk)

    # Concatenate all chunks into a single DataFrame
    df = pd.concat(dfs, sort=False)

    # Renaming the CHROM and POS columns
    vcf_df = df.rename(columns={"CHROM": "Chromosome", "POS": "Start"})
    vcf_df["End"] = vcf_df["Start"] + 1

    # Converting it to Pyranges object
    pyranges_vcf = pr.PyRanges(vcf_df)

    return pyranges_vcf
