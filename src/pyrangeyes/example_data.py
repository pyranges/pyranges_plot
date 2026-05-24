import pyranges1 as pr
import os
from .vcf.vcf_reader import read_vcf


p1 = pr.PyRanges(
    {
        "Chromosome": [1, 1, 2, 2, 2, 2, 2, 3],
        "Strand": ["+", "+", "-", "-", "-", "-", "-", "+"],
        "Start": [1, 40, 10, 70, 85, 110, 150, 140],
        "End": [11, 60, 25, 80, 100, 115, 180, 152],
        "transcript_id": ["t1", "t1", "t2", "t2", "t3", "t3", "t3", "t4"],
        "feature1": ["a", "a", "b", "b", "c", "c", "c", "d"],
        "feature2": ["A", "A", "B", "B", "C", "C", "C", "D"],
    }
)

p2 = pr.PyRanges(
    {
        "Chromosome": [1, 1, 2, 2, 2, 2, 2, 3, 4, 4, 4, 4, 4, 4],
        "Strand": [
            "+",
            "+",
            "-",
            "-",
            "+",
            "+",
            "+",
            "+",
            "-",
            "-",
            "-",
            "-",
            "+",
            "+",
        ],
        "Start": [
            1,
            40,
            10,
            70,
            85,
            110,
            150,
            140,
            30100,
            30150,
            30500,
            30647,
            29850,
            29970,
        ],
        "End": [
            11,
            60,
            25,
            80,
            100,
            115,
            180,
            152,
            30300,
            30300,
            30700,
            30700,
            29900,
            30000,
        ],
        "transcript_id": [
            "t1",
            "t1",
            "t2",
            "t2",
            "t3",
            "t3",
            "t3",
            "t4",
            "t5",
            "t5",
            "t5",
            "t5",
            "t6",
            "t6",
        ],
        "feature1": [
            "1",
            "1",
            "1",
            "1",
            "1",
            "2",
            "2",
            "2",
            "2",
            "2",
            "2",
            "2",
            "2",
            "2",
        ],
        "feature2": [
            "A",
            "A",
            "B",
            "B",
            "C",
            "C",
            "C",
            "D",
            "E",
            "E",
            "E",
            "E",
            "F",
            "F",
        ],
        "Feature": [
            "exon",
            "exon",
            "CDS",
            "CDS",
            "CDS",
            "CDS",
            "CDS",
            "exon",
            "exon",
            "CDS",
            "CDS",
            "exon",
            "CDS",
            "CDS",
        ],
    }
)

p3 = pr.PyRanges(
    {
        "Chromosome": ["1"] * 10 + ["2"] * 10,
        "Strand": ["+", "+", "+", "+", "-", "-", "-", "-", "+", "+"]
        + ["+", "+", "+", "+", "-", "-", "-", "-", "+", "+"],
        "Start": [90, 61, 104, 228, 9, 142, 52, 149, 218, 151]
        + [6, 27, 37, 47, 1, 7, 42, 37, 60, 80],
        "End": [92, 64, 113, 229, 12, 147, 57, 155, 224, 153]
        + [8, 32, 40, 50, 5, 10, 46, 40, 70, 90],
        "transcript_id": ["t1", "t1", "t1", "t1", "t2", "t2", "t2", "t2", "t3", "t3"]
        + ["t4", "t4", "t4", "t4", "t5", "t5", "t5", "t5", "t6", "t6"],
    }
)


mrna1 = pr.PyRanges(
    {
        "Chromosome": ["chr1"] * 14,
        "Strand": ["+"] * 14,
        "Feature": [
            "exon",
            "CDS",
            "exon",
            "CDS",
            "exon",
            "CDS",
            "exon",
            "CDS",
            "exon",
            "CDS",
            "exon",
            "CDS",
            "exon",
            "exon",
        ],
        "Start": [10, 25, 80, 80, 140, 140, 210, 225, 275, 275, 340, 340, 430, 470],
        "End": [55, 55, 120, 120, 190, 175, 255, 255, 315, 315, 390, 365, 455, 500],
        "transcript_id": ["tx1"] * 6 + ["tx2"] * 6 + ["lnc1"] * 2,
    }
)

mrna2 = pr.PyRanges(
    {
        "Chromosome": ["chr1"] * 18,
        "Strand": ["+"] * 18,
        "Feature": [
            "exon",
            "CDS",
            "exon",
            "CDS",
            "exon",
            "CDS",
            "exon",
            "exon",
            "CDS",
            "exon",
            "CDS",
            "exon",
            "CDS",
            "exon",
            "CDS",
            "exon",
            "exon",
            "exon",
        ],
        "Start": [
            35,
            50,
            115,
            115,
            170,
            170,
            245,
            45,
            62,
            130,
            130,
            205,
            205,
            285,
            285,
            365,
            420,
            465,
        ],
        "End": [
            90,
            90,
            150,
            150,
            225,
            210,
            295,
            100,
            100,
            175,
            175,
            245,
            245,
            330,
            315,
            405,
            445,
            505,
        ],
        "transcript_id": ["alt_tx1"] * 7 + ["alt_tx2"] * 9 + ["alt_lnc"] * 2,
    }
)


p_ala = pr.PyRanges(
    {
        "Start": [10, 50, 90] + [13, 60, 72, 120],
        "End": [20, 75, 130] + [16, 63, 75, 123],
        "Chromosome": [1] * 7,
        "id": ["gene1"] * 7,
        "trait1": ["exon"] * 3 + ["aa"] * 4,
        "trait2": ["gene_1"] * 3 + ["Ala"] * 4,
        "depth": [0] * 3 + [1] * 4,
        "thick": [0.3] * 3 + [0.6] * 4,
    }
)

p_cys = pr.PyRanges(
    {
        "Start": [10, 50, 90] + [15, 55, 62, 100, 110],
        "End": [20, 75, 130] + [18, 58, 65, 103, 113],
        "Chromosome": [1] * 8,
        "id": ["gene1"] * 8,
        "trait1": ["exon"] * 3 + ["aa"] * 5,
        "trait2": ["gene_1"] * 3 + ["Cys"] * 5,
        "depth": [0] * 3 + [1] * 5,
        "thick": [0.3] * 3 + [0.6] * 5,
    }
)

# Define the path to the data folder
DATA_DIR = os.path.join(os.path.dirname(__file__), "data")


def ncbi_gff():
    """
    Load the example NCBI GFF3 file as a PyRanges object.

    Returns:
        PyRanges: A PyRanges object containing the GFF3 data.
    """
    file_path = os.path.join(DATA_DIR, "ncbi.gff3")
    if not os.path.exists(file_path):
        raise FileNotFoundError(
            "The file 'ncbi.gff3' was not found in the data folder."
        )
    return pr.read_gff3(file_path)


def ncbi_vcf():
    """
    Load the example VCF file as a PyRanges object.

    Returns:
        PyRanges: A PyRanges object containing the VCF data.
    """
    file_path = os.path.join(DATA_DIR, "homo_sapiens_clinically_associated.vcf")
    if not os.path.exists(file_path):
        raise FileNotFoundError(
            "The file 'homo_sapiens_clinically_associated.vcf' was not found in the data folder."
        )
    return read_vcf(file_path)
