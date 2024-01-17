import numpy as np
import pysam
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

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

# Function to read and filter variants from a VCF file
def read_and_filter_variants(vcf_file):
    with pysam.VariantFile(vcf_file) as vcf:
        # Only include variants that have 'PASS' or '.' in the filter field
        return {
            (record.chrom, record.pos, record.ref, str(record.alts)): record
            for record in vcf
            if 'PASS' in record.filter.keys() or '.' in record.filter.keys()
        }

# Function to read variants without filtering for somatic_vcf_files
def read_variants(somatic_vcf_file):
    with pysam.VariantFile(somatic_vcf_file) as vcf:
        return {(rec.chrom, rec.pos, rec.ref, tuple(rec.alts)): rec for rec in vcf}

# Create the variant index
variant_to_index = {}

# Create a comprehensive index of all variants
for file in vcf_files + somatic_vcf_files:
    if 'somaticsniper' in file.lower():
        file_variants = read_variants(file)
    else:
        file_variants = read_and_filter_variants(file)
    for variant in file_variants:
        if variant not in variant_to_index:
            variant_to_index[variant] = len(variant_to_index)

# Create the binary vectors for each file
variant_matrix = []
for file in vcf_files + somatic_vcf_files:
    if 'somaticsniper' in file.lower():
        file_variants = read_variants(file)
    else:
        file_variants = read_and_filter_variants(file)
    file_vector = [1 if variant in file_variants else 0 for variant in variant_to_index]
    variant_matrix.append(file_vector)

variant_matrix = np.array(variant_matrix)

# Perform PCA
pca = PCA(n_components=3)  # reduce to three components for 3D visualization
reduced_data = pca.fit_transform(variant_matrix)

# Explained variance ratios
explained_variances = pca.explained_variance_ratio_ * 100  # Convert to percentages

# Generate names for the legend
names = []
for file in vcf_files + somatic_vcf_files:
    names.append(file.split("/")[6])

# Define colors for each variant caller type
colors = {
    'mutect': 'red',
    'somaticsniper': 'blue',
    'strelka': 'green'
}

# Visualization
fig = plt.figure(figsize=(12, 10))
ax = fig.add_subplot(111, projection='3d')

# Scatter plot for the first three principal components with color coding
for i, file in enumerate(names):
    variant_caller = file.split('-')[0].lower()
    color = colors.get(variant_caller, 'black')  # Default to black if caller not found
    ax.text(reduced_data[i, 0], reduced_data[i, 1], reduced_data[i, 2], file, color="black")
    ax.scatter(reduced_data[i, 0], reduced_data[i, 1], reduced_data[i, 2], color=color, label=variant_caller)

handles, labels = plt.gca().get_legend_handles_labels()
by_label = dict(zip(labels, handles))
ax.legend(by_label.values(), by_label.keys())

# Update axis labels to include explained variance
ax.set_xlabel(f'Principal Component 1 ({explained_variances[0]:.2f}%)')
ax.set_ylabel(f'Principal Component 2 ({explained_variances[1]:.2f}%)')
ax.set_zlabel(f'Principal Component 3 ({explained_variances[2]:.2f}%)')
ax.set_title('3D PCA of Filtered VCF Files')

plt.show()
plt.savefig("PCA3_v2t.png")