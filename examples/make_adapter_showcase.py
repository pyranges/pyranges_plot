"""Generate adapter tutorial figures and a PR showcase PDF."""

from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import pyranges1 as pr

import pyrangeyes as pre

matplotlib.use("Agg")

ROOT = Path(__file__).resolve().parents[1]
DOCS_IMAGES = ROOT / "docs" / "images"
ARTIFACTS = ROOT / "artifacts"
DOCS_IMAGES.mkdir(parents=True, exist_ok=True)
ARTIFACTS.mkdir(parents=True, exist_ok=True)


def tutorial_mrna():
    return pre.example_data.p2


def tutorial_snps():
    return pr.PyRanges(
        {
            "Chromosome": [1, 1, 2, 2, 4],
            "Start": [7, 45, 13, 73, 29870],
            "End": [8, 46, 14, 74, 29871],
            "ID": ["rs1", "rs2", "rs3", "rs4", "rs5"],
            "REF": ["A", "G", "C", "T", "A"],
            "ALT": ["T", "A", "G", "C", "G"],
        }
    )


def tutorial_unpacked():
    return pr.PyRanges(
        {
            "Chromosome": ["1", "1", "1", "1"],
            "Start": [5, 13, 18, 4],
            "End": [10, 15, 22, 7],
            "transcript_id": ["first", "second", "second", "third"],
        }
    )


def save(fig, path):
    fig.savefig(path, dpi=120, bbox_inches="tight")
    plt.close(fig)


def make_figures():
    pre.set_engine("matplotlib")

    pp = tutorial_mrna()
    snps = tutorial_snps()

    fig = pre.plot(pp, "mRNA", return_plot="fig", warnings=False)
    save(fig, DOCS_IMAGES / "prp_rtd_35.png")

    fig = pre.plot(
        snps,
        "SNP",
        color_col="ALT",
        width=40,
        return_plot="fig",
        warnings=False,
    )
    save(fig, DOCS_IMAGES / "prp_rtd_36.png")

    fig = pre.plot(
        [pp, snps],
        adapter=["mRNA", "SNP"],
        width=40,
        return_plot="fig",
        warnings=False,
    )
    save(fig, DOCS_IMAGES / "prp_rtd_37.png")

    fig = pre.plot(
        tutorial_unpacked(),
        id_col="transcript_id",
        packed=False,
        return_plot="fig",
        warnings=False,
        title_chr=" ",
    )
    save(fig, DOCS_IMAGES / "prp_rtd_38.png")


def make_pdf():
    from matplotlib.backends.backend_pdf import PdfPages

    make_figures()
    pdf_path = ARTIFACTS / "adapter_showcase.pdf"
    with PdfPages(pdf_path) as pdf:
        for image_name, title in [
            ("prp_rtd_35.png", "mRNA adapter: thin UTR/exon, thick CDS"),
            ("prp_rtd_36.png", "SNP adapter: diamond variant markers"),
            ("prp_rtd_37.png", "Per-track adapters: mRNA + SNP"),
            ("prp_rtd_38.png", "text=None: labels hidden for unpacked plots"),
        ]:
            img = plt.imread(DOCS_IMAGES / image_name)
            fig, ax = plt.subplots(figsize=(11, 6))
            ax.imshow(img)
            ax.set_title(title)
            ax.axis("off")
            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)
    return pdf_path


if __name__ == "__main__":
    print(make_pdf())
