#!/bin/bash
# Annotation of variants with SnpEff
# ----------------------------------
# This pipeline prepares a SnpEff database for Yarrowia lipolytica using Genolevure annotations
# and performs annotation of SNPs identified by Fast-LMM.

# -------------------------------------------------------
# STEP 1: Install and prepare AGAT for EMBL → GFF3 conversion
# -------------------------------------------------------

# Activate the Fast-LMM environment
conda activate fastlmm

# Attempt to install AGAT with conda (note: it might fail)
conda install -c bioconda agat  # This does not work reliably

# Clone AGAT from GitHub and install dependencies manually
git clone https://github.com/NBISweden/AGAT.git
cd AGAT

# Install required Perl modules
cpanm File::Basename Bio::SeqIO Getopt::Long Pod::Usage
sudo cpan Bio::SeqIO

# Set AGAT paths in the Perl environment
export PERL5LIB=$PERL5LIB:/mnt/c/Users/USUARIO/Desktop/Acetate/Fast-LMM/snpEff/AGAT/lib
export PERL5LIB=/mnt/c/Users/USUARIO/Desktop/Acetate/Fast-LMM/snpEff/AGAT/share:$PERL5LIB
export PERL5LIB=/mnt/c/Users/USUARIO/Desktop/Acetate/Fast-LMM/snpEff/AGAT:$PERL5LIB
export AGAT_CONFIG_DIR=/mnt/c/Users/USUARIO/Desktop/Acetate/Fast-LMM/snpEff/AGAT/share/

# Check if PERL5LIB was set correctly
perl -V

# Install additional Perl modules required by AGAT
cpan File::ShareDir
cpan Sort::Naturally

# -------------------------------------------------------
# STEP 2: Convert EMBL to GFF3 (do not concatenate beforehand)
# -------------------------------------------------------

# Convert the reference chromosome (example with YALI0A)
#perl bin/agat_convert_embl2gff.pl \
#  --embl /mnt/c/Users/USUARIO/Desktop/Acetate/Fast-LMM/snpEff/YALI0A.embl \
#  -o YALI0A.gff3 \
#  -config /mnt/c/Users/USUARIO/Desktop/Acetate/Fast-LMM/snpEff/AGAT/share/agat_config.yaml

# Convert all chromosomes A-F
for letter in A B C D E F; do
  perl bin/agat_convert_embl2gff.pl \
    --embl /mnt/c/Users/USUARIO/Desktop/Acetate/Fast-LMM/snpEff/YALI0${letter}.embl \
    -o /mnt/c/Users/USUARIO/Desktop/Acetate/Fast-LMM/snpEff/YALI0${letter}.gff3 \
    -config /mnt/c/Users/USUARIO/Desktop/Acetate/Fast-LMM/snpEff/AGAT/share/agat_config.yaml
done

# -------------------------------------------------------
# STEP 3: Clean unwanted FASTA sections from GFF3 files
# -------------------------------------------------------

# Remove ##FASTA and following lines from all GFF3 files
for f in *.gff3; do
  awk '/##FASTA/ { exit } { print }' "$f" > temp && mv temp "$f"
done

# -------------------------------------------------------
# STEP 4: Concatenate chromosome GFF3 files
# -------------------------------------------------------

# Start the combined GFF3 file
echo "##gff-version 3" > YALI0.gff3

# Concatenate in chromosome order (A to F)
for chr in YALI0A YALI0B YALI0C YALI0D YALI0E YALI0F; do
  grep -v '^##' "${chr}.gff3" >> YALI0.gff3
done

# -------------------------------------------------------
# STEP 5: Prepare SnpEff database directory
# -------------------------------------------------------

# Create SnpEff data directory for Yarrowia genome
mkdir -p ./YALI0

# Move your reference genome FASTA and annotations to this folder
# You must ensure that chromosome names match exactly between the FASTA and the GFF3

# Add the following line to your snpEff.config file:
# YALI0.genome : Yarrowia_lipolytica
# (YALI0 must match the folder name and annotation label)

# -------------------------------------------------------
# STEP 6: Build the SnpEff database
# -------------------------------------------------------

# Run database build using the GFF3 annotation
java -jar /mnt/c/Users/USUARIO/Desktop/Acetate/Fast-LMM/snpEff/snpEff/snpEff.jar build -gff3 -v YALI0

# -------------------------------------------------------
# STEP 7: If Genolevure lacks gene/exon features
# -------------------------------------------------------

# Genolevure annotations may lack complete gene/exon definitions.
# Reference: Neuvéglise et al. 2024 https://microbialcellfactories.biomedcentral.com/articles/10.1186/s12934-024-02354-9

# Use genome assembly ASM252v1 (fna and gtt files).
# Ensure to change chromosome identifiers as needed.

# -------------------------------------------------------
# STEP 8: Additional input formatting for database
# -------------------------------------------------------

# Required files and locations:
# protein.fa  → .data/YALI0/
# cds.fa      → .data/YALI0/ 
#   - Requires formatting: change identifiers to `TRANSCRIPT_gene-[locus_tag]`
#   - Use the custom script `cdsformatting.py` to do this
# sequence.fa or YALI0.fa → .data/genomes/
# genes.gtf   → .data/YALI0/

# Rebuild the database after adjusting files
java -jar /mnt/c/Users/USUARIO/Desktop/Acetate/Fast-LMM/snpEff/snpEff/snpEff.jar build -gff3 -v YALI0

# -------------------------------------------------------
# STEP 9: Annotate variants using SnpEff
# -------------------------------------------------------

cd /mnt/c/Users/USUARIO/Desktop/Acetate

# Annotate the VCF file (3000 bp upstream/downstream included)
java -Xmx8g -jar /mnt/c/Users/USUARIO/Desktop/Acetate/Fast-LMM/snpEff/snpEff/snpEff.jar -ud 3000 -v YALI0 phenosamples.recode.vcf > phenosamples1.ann.vcf

# -------------------------------------------------------
# OUTPUTS:
# -------------------------------------------------------
# phenosamples1.ann.vcf → Annotated VCF with SnpEff results
# snpEff_summary.html   → Summary report of the annotation
# snpEff/               → Gene-level effect annotations

