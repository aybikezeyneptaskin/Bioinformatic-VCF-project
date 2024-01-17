import numpy as np
import pysam
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster.hierarchy import linkage
import matplotlib.colors as mcolors

# Define the function to read variants from a VCF file
def read_variants(vcf_file):
    with pysam.VariantFile(vcf_file) as vcf:
        variants = {(rec.chrom, rec.pos, rec.ref, tuple(rec.alts)) for rec in vcf}
    return variants

def read_and_filter_variants(vcf_file):
    with pysam.VariantFile(vcf_file) as vcf:
        # Filtering out records where FILTER is not 'PASS' or '.'
        variants = {
            (rec.chrom, rec.pos, rec.ref, tuple(rec.alts))
            for rec in vcf
            if 'PASS' in rec.filter.keys() or '.' in rec.filter.keys()
        }
    return variants

# List of VCF files
vcf_files = [
    "mutect-bwa-recal/all_normal_final.recode.vcf",  
    "mutect-bwa-norecal/all_normal_final.recode.vcf",
    "mutect-bowtie-recal/all_normal_final.recode.vcf",
    "mutect-bowtie-norecal/all_normal_final.recode.vcf",
    "strelka-bwa-recal/all_normal_final.recode.vcf",
    "strelka-bwa-norecal/all_normal_final.recode.vcf",
    "strelka-bowtie-recal/all_normal_final.recode.vcf",
    "strelka-bowtie-norecal/all_normal_final.recode.vcf"

]

somatic_vcf_files = [
    "somaticsniper-bwa-recal/output-somaticsniper_finall.vcf",
    "somaticsniper-bwa-norecal/output-somaticsniper_finall.vcf",
    "somaticsniper-bowtie-recal/output-somaticsniper_finall.vcf",
    "somaticsniper-bowtie-norecal/output-somaticsniper_finall.vcf"
]

names = []
for file in vcf_files:
    names.append(file.split("/")[0])

for file in somatic_vcf_files:
    names.append(file.split("/")[0])

# Read variants from each VCF file
variant_sets = [read_and_filter_variants(vcf_file) for vcf_file in vcf_files]

# Read variants from each VCF file
somatic_variant_sets = [read_variants(vcf_file) for vcf_file in somatic_vcf_files]

variant_sets = variant_sets + somatic_variant_sets

# Create a comprehensive index of all variants
variant_to_index = {}
for variant_set in variant_sets:
    for variant in variant_set:
        if variant not in variant_to_index:
            variant_to_index[variant] = len(variant_to_index)

# Create a binary matrix
binary_matrix = []
for variant_set in variant_sets:
    binary_vector = [1 if variant in variant_set else 0 for variant in variant_to_index]
    binary_matrix.append(binary_vector)

binary_matrix = np.array(binary_matrix)

# Perform hierarchical clustering
row_linkage = linkage(binary_matrix, method='average')
col_linkage = linkage(binary_matrix.T, method='average')

# Create a clustermap with dendrograms
colors = ["pink", "crimson"]
cmap = mcolors.LinearSegmentedColormap.from_list("", colors)

g = sns.clustermap(binary_matrix, row_linkage=row_linkage, col_linkage=col_linkage,
                   cmap=cmap, figsize=(32, 16), yticklabels=names, xticklabels=False)

g.fig.suptitle("Heatmap with Dendrograms for Pipelines and Variants")
plt.savefig("heatmap_variants.png")
plt.show()