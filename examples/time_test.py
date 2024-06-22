import pyranges as pr
import pyranges_plot as prp
import time
from memory_profiler import memory_usage


prp.set_engine("plt")
prp.set_warnings(False)

# Download data: https://ftp.ensembl.org/pub/release-112/gtf/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.46.112.gtf.gz

# Define chromosome and feature to test
chrom = "X"
feat = "CDS"

start_loading_t = time.time()
start_loading_m = memory_usage()[0]

path = "drosophila/Drosophila_melanogaster.BDGP6.46.112.gtf"  ## modify if needed
d = pr.read_gtf(path)

end_loading_t = time.time()
end_loading_m = memory_usage()[0]

# Get CDS
d_cds = d[d["Feature"] == feat]

# Check number of genes in each chrom
# d_cds.groupby("Chromosome").apply(
#     lambda x: print(
#         "CHROM:"
#         + str(x["Chromosome"].iloc[0])
#         + ":   "
#         + str(len(x.groupby("gene_id")))
#     )
# )

# Chromosome subset
d_cds_chrom = d_cds[d_cds["Chromosome"] == chrom]


end_subset_t = time.time()
end_subset_m = memory_usage()[0]

# Perform plot
prp.plot(d_cds_chrom, id_col="gene_id", to_file="_tmp.png", max_shown=3000)

end_plot_t = time.time()
end_plot_m = memory_usage()[0]

# Print results
print("==============================================================")
print(f"Drosophila melanogaster | Chromosome {chrom} | Feature {feat}\n")

print("Loading time: \t" + str(end_loading_t - start_loading_t) + " (s)")
print("Loading memory:\t" + str(end_loading_m - start_loading_m) + " (MiB)")
print()
print("Subsetting time: \t" + str(end_subset_t - end_loading_t) + " (s)")
print("Subsetting memory: \t" + str(end_subset_m - end_loading_m) + " (MiB)")
print()
print("Plotting and exporting time: \t" + str(end_plot_t - end_subset_t) + " (s)")
print("Plotting and exporting memory: \t" + str(end_plot_m - end_subset_m) + " (MiB)")
print("\n\n\n")
