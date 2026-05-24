import pyranges1 as pr
import pyrangeyes as pre

pre.set_engine("plotly")

## Sars cov 2
# Download data: http://ftp.ensemblgenomes.org/pub/viruses/gtf/sars_cov_2/Sars_cov_2.ASM985889v3.101.gtf.gz
p = pr.read_gtf(
    "../../performance/sars_cov_2/Sars_cov_2.ASM985889v3.101.gtf"
)  ## modify if needed

p_cds = p[p["Feature"] == "CDS"]

pre.plot(
    p_cds,
    id_col="gene_name",
    label=False,
    legend=True,
    to_file=("fig_4.png", (600, 500)),
    arrow_size=0.05,
)
