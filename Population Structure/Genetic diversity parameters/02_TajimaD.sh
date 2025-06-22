#!/bin/bash
# Script for Tajima's D analysis
# This script subdivides a VCF by population, filters variants,
# collapses chromosomes into a single coordinate system,
# and calculates Tajima's D across the entire genome.

# Input VCF file
VCF_IN="hqvariants_yayaout_snps.recode.vcf"

# Step 1: Subdivide VCF by population
# For each sample list, create a filtered VCF containing only those samples
for LIST in Eur1.txt Eur2.txt SoHC.txt IndNA.txt Mosaic.txt W29.txt
do
    BASENAME=$(basename "$LIST" .txt)
    vcftools --vcf "$VCF_IN" \
             --keep "$LIST" \
             --recode --recode-INFO-all \
             --out "${BASENAME}_filtered" \
             --max-missing 1 \
             --non-ref-ac-any 1
done

# Step 2: Filter for only biallelic SNPs
for VCF in *_filtered.recode.vcf
do
    BASENAME=$(basename "$VCF" .recode.vcf)
    vcftools --vcf "$VCF" \
             --max-missing 1 \
             --max-alleles 2 \
             --remove-indels \
             --recode --recode-INFO-all \
             --out "${BASENAME}_biallelic" \
             --non-ref-ac-any 1
done

# Step 3: Collapse chromosomes into a single coordinate system (YALI0A)
# This helps in treating the genome as one continuous sequence for sliding window analyses
for VCF in *_biallelic.recode.vcf
do
    BASENAME=$(basename "$VCF" .recode.vcf)
    OUT="${BASENAME}_recoded_cumulative.vcf"

    awk -v OFS='\t' '
    BEGIN {
        # Define cumulative offsets for each chromosome based on their lengths
        offset["YALI0A"] = 0
        offset["YALI0B"] = 2303260
        offset["YALI0C"] = 2303260 + 3066374
        offset["YALI0D"] = 2303260 + 3066374 + 3272609
        offset["YALI0E"] = 2303260 + 3066374 + 3272609 + 3661369
        offset["YALI0F"] = 2303260 + 3066374 + 3272609 + 3661369 + 4224102
    }
    /^#/ { print; next }  # Print header lines unchanged
    {
        chr = $1
        pos = $2
        if (chr in offset) {
            $1 = "YALI0A"           # Assign all variants to chromosome YALI0A
            $2 = pos + offset[chr]  # Add the cumulative offset to position
        }
        print
    }
    ' "$VCF" > "$OUT"
done

# Step 4: Repeat chromosome collapsing for the global VCF (biallelicvariants.recode.vcf)
awk '
BEGIN {
    offset["YALI0A"] = 0
    offset["YALI0B"] = 2303260
    offset["YALI0C"] = 2303260 + 3066374
    offset["YALI0D"] = 2303260 + 3066374 + 3272609
    offset["YALI0E"] = 2303260 + 3066374 + 3272609 + 3661369
    offset["YALI0F"] = 2303260 + 3066374 + 3272609 + 3661369 + 4224102
}
/^##/ { print; next }        # Print meta-information lines unchanged
/^#CHROM/ { print; next }    # Print header line unchanged
{
    chr = $1
    pos = $2
    if (chr == "YALI0B" || chr == "YALI0C" || chr == "YALI0D" || chr == "YALI0E" || chr == "YALI0F") {
        $1 = "YALI0A"
        $2 = pos + offset[chr]
    }
    print
}
' biallelicvariants.recode.vcf > biallelicvariants_recoded_cumulative.vcf

# Step 5: Calculate Tajima's D in sliding windows across the genome
WINDOW=20531075

for VCF in *_recoded_cumulative.vcf
do
    BASENAME=$(basename "$VCF" .vcf)
    vcftools --vcf "$VCF" \
             --TajimaD "$WINDOW" \
             --out "tajima_${BASENAME}"
done