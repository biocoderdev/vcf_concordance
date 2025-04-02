import os
import sys
import matplotlib
matplotlib.use('Agg')  # Set non-interactive backend
from matplotlib import pyplot as plt
from matplotlib_venn import venn2, venn3
from upsetplot import from_contents, plot
import pandas as pd
from cyvcf2 import VCF

def read_vcf_variants(vcf_file):
    """Read variants from a VCF file using cyvcf2"""
    variants = set()
    try:
        reader = VCF(vcf_file)
        for variant in reader:
            # Use CHROM, POS, REF, ALT as unique identifier
            alt = variant.ALT[0] if variant.ALT else ""
            variant_id = (variant.CHROM, variant.POS, variant.REF, alt)
            variants.add(variant_id)
    except Exception as e:
        print(f"Error reading {vcf_file}: {str(e)}", file=sys.stderr)
    return variants

def create_venn_diagram(vcf_files, output_file=None):
    """Create a Venn diagram for 2-3 VCF files"""
    variant_sets = {}
    labels = [os.path.basename(f).replace('.vcf', '') for f in vcf_files]
    
    for vcf_file, label in zip(vcf_files, labels):
        variant_sets[label] = read_vcf_variants(vcf_file)
    
    plt.figure(figsize=(10, 8))
    
    if len(vcf_files) == 2:
        venn2([variant_sets[labels[0]], variant_sets[labels[1]]], set_labels=labels)
    elif len(vcf_files) == 3:
        venn3([variant_sets[labels[0]], variant_sets[labels[1]], variant_sets[labels[2]]], 
              set_labels=labels)
    else:
        print("Venn diagrams are not recommended for more than 3 sets. Using UpSet plot instead.")
        create_upset_plot(vcf_files, output_file)
        return
    
    plt.title("Variant Concordance Between VCF Files")
    
    if output_file:
        plt.savefig(output_file, bbox_inches='tight', dpi=300)
        print(f"Venn diagram saved to {output_file}")
    else:
        plt.show()

def create_upset_plot(vcf_files, output_file=None):
    """Create an UpSet plot for multiple VCF files"""
    variant_sets = {}
    labels = [os.path.basename(f).replace('.vcf', '') for f in vcf_files]
    
    for vcf_file, label in zip(vcf_files, labels):
        variant_sets[label] = read_vcf_variants(vcf_file)
    
    upset_data = from_contents(variant_sets)
    plt.figure(figsize=(12, 8))
    plot(upset_data, show_counts=True, sort_by='cardinality')
    plt.suptitle("Variant Concordance Between VCF Files")
    
    if output_file:
        plt.savefig(output_file, bbox_inches='tight', dpi=300)
        print(f"UpSet plot saved to {output_file}")
    else:
        plt.show()

def analyze_vcf_concordance(vcf_files, output_file=None):
    """Analyze concordance between VCF files"""
    if len(vcf_files) < 2:
        print("Please provide at least 2 VCF files for comparison.")
        return
    
    if len(vcf_files) <= 3:
        create_venn_diagram(vcf_files, output_file)
    else:
        create_upset_plot(vcf_files, output_file)

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python vcf_concordance.py file1.vcf file2.vcf [file3.vcf ...] [--output output.png]")
        sys.exit(1)
    
    vcf_files = []
    output_file = None
    i = 1
    while i < len(sys.argv):
        arg = sys.argv[i]
        if arg.endswith('.vcf') or arg.endswith('.vcf.gz'):
            vcf_files.append(arg)
            i += 1
        elif arg == '--output':
            output_file = sys.argv[i+1]
            i += 2
        else:
            i += 1
    
    if len(vcf_files) < 2:
        print("Please provide at least 2 VCF files for comparison.")
        sys.exit(1)
    
    # If output is directory, create default filename
    if output_file and os.path.isdir(output_file):
        output_file = os.path.join(output_file, "vcf_concordance.png")
    
    analyze_vcf_concordance(vcf_files, output_file)
