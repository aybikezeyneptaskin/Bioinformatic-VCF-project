import pysam

def read_vcf(file):
    with pysam.VariantFile(file) as vcf:
        # Store variants by chromosome and position
        return {(record.chrom, record.pos): record for record in vcf}

def read_and_filter_vcf(file):
    with pysam.VariantFile(file) as vcf:
        filtered_records = {
            (record.chrom, record.pos): record
            for record in vcf
            if 'PASS' in record.filter.keys() or '.' in record.filter.keys() # Filter out records where FILTER is not 'PASS' or '.'
        }
    return filtered_records

def compare_variants(truth, predictions):
    truth_set = set(truth.keys())
    predictions_set = set(predictions.keys())

    tp = len(truth_set & predictions_set)
    fp = len(predictions_set - truth_set)
    fn = len(truth_set - predictions_set)

    return tp, fp, fn

def calculate_metrics(tp, fp, fn):
    precision = tp / (tp + fp) if tp + fp > 0 else 0
    recall = tp / (tp + fn) if tp + fn > 0 else 0
    if precision + recall > 0:
        f1 = 2 * (precision * recall) / (precision + recall)
    else:
        f1 = 0
    accuracy = tp / (tp + fp + fn) if tp + fp + fn > 0 else 0

    return round(precision,2), round(recall,2), round(f1,2), round(accuracy,2)

# File paths for the VCF files
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
# Read ground truth VCF
ground_truth = read_vcf("/mnt/c/Users/aybik/VCF/hc_bed_filtered.recode.vcf")

# Process each VCF file and calculate metrics
for file in vcf_files:
    dir_name = file.split("/")[6]
    vcf_data = read_and_filter_vcf(file)
    tp, fp, fn = compare_variants(ground_truth, vcf_data)
    precision, recall, f1, accuracy = calculate_metrics(tp, fp, fn)
    print(f"Metrics for {dir_name}: Precision={precision}, Recall={recall}, F1-Score={f1}, Accuracy={accuracy}")

for file in somatic_vcf_files:
    dir_name = file.split("/")[6]
    vcf_data = read_vcf(file)
    tp, fp, fn = compare_variants(ground_truth, vcf_data)
    precision, recall, f1, accuracy = calculate_metrics(tp, fp, fn)
    print(f"Metrics for {dir_name}: Precision={precision}, Recall={recall}, F1-Score={f1}, Accuracy={accuracy}")

