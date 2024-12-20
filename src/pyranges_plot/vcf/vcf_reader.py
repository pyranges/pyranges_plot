import pandas as pd
import pyranges as pr
from pathlib import Path
from io import StringIO

def read_vcf(f: str | Path, nrows: bool | None = None):
    path = Path(f)
    dtypes = {
        "CHROM": "category",
        "POS": "int32",
        "ID": "category",
        "REF": "str",
        "ALT": "str",
        "QUAL": "category",
        "FILTER": "category"
    }

    # Read the file to find the header line and ignore metadata
    with open(path, 'r') as file:
        lines = []
        for line in file:
            if line.startswith("##"):
                continue  # Skip other metadata lines
            elif line.startswith("#CHROM"):
                # This is the header line; extract the column names
                column_names = line.strip('#').strip().split('\t')
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
        chunk['QUAL'] = chunk['QUAL'].astype(object).replace('.', pd.NA)
        dfs.append(chunk)
    
    # Concatenate all chunks into a single DataFrame
    df = pd.concat(dfs, sort=False)

    # Renaming the CHROM and POS columns
    vcf_df = df.rename(columns={"CHROM": "Chromosome", "POS": "Start"})
    vcf_df["End"] = vcf_df["Start"] + 1

    # Converting it to Pyranges object
    pyranges_vcf = pr.PyRanges(vcf_df)

    return pyranges_vcf