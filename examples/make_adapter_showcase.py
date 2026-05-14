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


def save(fig, path):
    fig.savefig(path, dpi=120, bbox_inches="tight")
    plt.close(fig)


def make_figures():
    """Generate the exact figures referenced by docs/tutorial.rst."""

    pre.set_engine("matplotlib")

    pp = tutorial_mrna()
    snps = tutorial_snps()

    fig = pre.plot(pp, "mRNA", return_plot="fig", warnings=False)
    save(fig, DOCS_IMAGES / "prp_rtd_35.png")

    fig = pre.plot(
        snps,
        "SNP",
        color_col="ALT",
        shape="diamond",
        return_plot="fig",
        warnings=False,
    )
    save(fig, DOCS_IMAGES / "prp_rtd_36.png")

    fig = pre.plot(
        [pp, snps],
        adapter=["mRNA", "SNP"],
        shape="triangle-up",
        return_plot="fig",
        warnings=False,
    )
    save(fig, DOCS_IMAGES / "prp_rtd_37.png")


def _showcase_page(pdf, fig, title):
    fig.suptitle(title)
    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)


def make_pdf():
    from matplotlib.backends.backend_pdf import PdfPages

    make_figures()
    pdf_path = ARTIFACTS / "adapter_showcase.pdf"
    pp = tutorial_mrna()
    snps = tutorial_snps()

    with PdfPages(pdf_path) as pdf:
        pre.set_engine("matplotlib")
        _showcase_page(
            pdf,
            pre.plot(pp, "mRNA", return_plot="fig", warnings=False),
            "mRNA adapter: thin UTR/exon, thick CDS",
        )
        for shape in ["diamond", "triangle-up", "triangle-down", "circle"]:
            _showcase_page(
                pdf,
                pre.plot(
                    snps,
                    "SNP",
                    shape=shape,
                    color_col="ALT",
                    return_plot="fig",
                    warnings=False,
                ),
                f"SNP adapter: fixed-size {shape} markers",
            )
        _showcase_page(
            pdf,
            pre.plot(
                [pp, snps],
                adapter=["mRNA", "SNP"],
                shape="triangle-up",
                return_plot="fig",
                warnings=False,
            ),
            "Per-track adapters: mRNA + SNP triangle-up markers",
        )
    return pdf_path


if __name__ == "__main__":
    print(make_pdf())
