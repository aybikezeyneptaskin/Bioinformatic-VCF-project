import os
import numpy as np
import pysam
import matplotlib.pyplot as plt
import seaborn as sns

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

# Define the function to calculate Jaccard distance
def jaccard_distance(set1, set2):
    intersection = len(set1 & set2)
    union = len(set1 | set2)
    return 1 - intersection / union

# List of VCF files
vcf_files = [
    "/mnt/c/Users/aybik/VCF/mutect-bwa-recal/all_normal_final.recode.vcf",  
    "/mnt/c/Users/aybik/VCF/mutect-bwa-norecal/all_normal_final.recode.vcf",
    "/mnt/c/Users/aybik/VCF/mutect-bowtie-recal/all_normal_final.recode.vcf",
    "/mnt/c/Users/aybik/VCF/mutect-bowtie-norecal/all_normal_final.recode.vcf",
    "/mnt/c/Users/aybik/VCF/strelka-bwa-recal/all_normal_final.recode.vcf",
    "/mnt/c/Users/aybik/VCF/strelka-bwa-norecal/all_normal_final.recode.vcf",
    "/mnt/c/Users/aybik/VCF/strelka-bowtie-recal/all_normal_final.recode.vcf",
    "/mnt/c/Users/aybik/VCF/strelka-bowtie-norecal/all_normal_final.recode.vcf"

]

somatic_vcf_files = [
    "/mnt/c/Users/aybik/VCF/somaticsniper-bwa-recal/output-somaticsniper_finall.vcf",
    "/mnt/c/Users/aybik/VCF/somaticsniper-bwa-norecal/output-somaticsniper_finall.vcf",
    "/mnt/c/Users/aybik/VCF/somaticsniper-bowtie-recal/output-somaticsniper_finall.vcf",
    "/mnt/c/Users/aybik/VCF/somaticsniper-bowtie-norecal/output-somaticsniper_finall.vcf"
]

names = []
for file in vcf_files:
    names.append(file.split("/")[6])

for file in somatic_vcf_files:
    names.append(file.split("/")[6])

# Read variants from each VCF file
variant_sets = [read_and_filter_variants(vcf_file) for vcf_file in vcf_files]

# Read variants from each VCF file
somatic_variant_sets = [read_variants(vcf_file) for vcf_file in somatic_vcf_files]

variant_sets = variant_sets + somatic_variant_sets

# Calculate Jaccard distances
jaccard_matrix = np.zeros((len(variant_sets), len(variant_sets)))

reversed_set = variant_sets.copy()
reversed_set.reverse()

for i, variants_i in enumerate(variant_sets):
    print(f"number of variants for {names[i]} : {len(variant_sets[i])}")
    for j, variants_j in enumerate(reversed_set):
        #print(type(variants_i),type(variants_j))
        jaccard_matrix[i, j] = jaccard_distance(variants_i, variants_j)

reverse_names = names.copy()
reverse_names.reverse()
# Generate heatmap
plt.figure(figsize=(16, 16))
sns.heatmap(jaccard_matrix, annot=True, fmt=".2f", cmap='coolwarm', xticklabels=names, yticklabels=reverse_names)
plt.title("Jaccard Distance Heatmap")
plt.show()
plt.savefig('heatmap.png')