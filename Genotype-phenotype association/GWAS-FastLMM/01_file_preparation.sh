#!/bin/bash
#Pipeline to carry out genome-wide association study with Fast-LMM.

##Genrate SNPs vcfs adapted with samples for the analysis
vcftools --vcf all_raw_variants_yaya_outgroup.vcf --remove remove.txt --minQ 30 --min-meanDP 20 --non-ref-ac-any 1 --recode --recode-INFO-all --out hqvariantswithyayanomutants.vcf


#First vcftools is used to filter SNPs, individuals, etc. For the first analyisis we will use SNPs that were called in all individuals, and we only retain biallelic sites
vcftools --max-missing 1 --max-alleles 2 --vcf hqvariantswithyayanomutants.vcf.recode.vcf --recode --recode-INFO-all --out filtered_snps_max

#Create a text file with the name of the individual to remove (remove)
#Important to filter according to minimum allele frequency
vcftools --keep Acetatesamplescont.txt --vcf filtered_snps_max.recode.vcf --recode --recode-INFO-all --remove-indels --non-ref-ac-any 1 --maf 0.05 --out ./Fast-LMM/phenosamples



#Create an absence/presence matrix of the SNPs
#If you use plink:
#Remove underscores from sample names
sed -E 's/(_)(\t|$)/\2/g' phenosamples.recode.vcf > phenosamples.nounderscore.vcf
#Change chromosomes names
sed -E 's/^YALI0A/1/;
         s/^YALI0B/2/;
         s/^YALI0C/3/;
         s/^YALI0D/4/;
         s/^YALI0E/5/;
         s/^YALI0F/6/' phenosamples.nounderscore.vcf > phenosamples.renamed.vcf





##Fast-LMM
#This tool is used in "Extensive simulations assess theperformance of genome-wide associationmapping in various Saccharomycescerevisiae subpopulations", Schacherer
#It is a statistical method based on linear mixed model, so that it can correct for population structure effect.
#Input: SNPs to be tested, SNP matrix to correct for pop structure (can be the same) and phenotype data.

#SNP data must be in plink format. Since we can input continous phenotype variables, we will use the larger vcf,

#Transforming SNP data into plink FORMAT
#Annotate (id field is empty
#Annotate ID column
bcftools annotate --set-id '%CHROM\_%POS' \
  -o phenosamples.renamedids.vcf \
  -O v \
  phenosamples.recode.vcf
 

#Change chromosomes to numbers
awk 'BEGIN{FS=OFS="\t"}
     /^#/ {print; next}
     $1 ~ /^YALI0A/ { $1="1" }
     $1 ~ /^YALI0B/ { $1="2" }
     $1 ~ /^YALI0C/ { $1="3" }
     $1 ~ /^YALI0D/ { $1="4" }
     $1 ~ /^YALI0E/ { $1="5" }
     $1 ~ /^YALI0F/ { $1="6" }

     $3 ~ /^YALI0A/ { $3="1_"$2 }
     $3 ~ /^YALI0B/ { $3="2_"$2 }
     $3 ~ /^YALI0C/ { $3="3_"$2 }
     $3 ~ /^YALI0D/ { $3="4_"$2 }
     $3 ~ /^YALI0E/ { $3="5_"$2 }
     $3 ~ /^YALI0F/ { $3="6_"$2 }

     {print}' phenosamples.renamedids.vcf > phenosamples.renamed.vcf
		 
#Change
plink --vcf ./phenosamples.renamed.vcf \
      --make-bed \
      --allow-extra-chr --double-id \
      --out fastlmm_input

##Create input files for the analysis
	  
#Create the correct input phenotype file: two columns of sample names and one of phenotypes
#Ordering samples on the vcf file order, which is now stored on the .fam file.
#Save phenotype xlsx as a txt delimited by tabs. We need to retrieve phenotype data from this file and order it accordingly
#two columns of ids
#cut -d' ' -f1,2 fastlmm_input.fam > ids.txt
awk '{print $1"\t"$2}' fastlmm_input.fam | sort > ids_sorted.tsv
#Extract phenotype data and order
awk -F'\t' 'NR>1 {print $9}' ./Acetate2.txt > col1.txt
awk -F'\t' 'NR>1 {print $2}' ./Acetate2.txt > col2.txt
dos2unix col1.txt
dos2unix col2.txt

paste col1.txt col2.txt | sort > growth_data.tsv
cut -f2 growth_data.tsv > col2_growth_data.txt


#joining:
paste ids_sorted.tsv col2_growth_data.txt > pheno.txt
#Numbers must be coded with ".", not ","
sed 's/,/./g' pheno.txt > pheno1.txt

#We are going to prepare covar file, which should have same structure as pheno1.txt, with principal component from pop structure analysis
#We take the information from file yali.pca.eve (pca coordinates for PC1 and PC2).
awk '!/^#/ { print $1, $1, $2, $3 }' yali.pca.evec > covariables.txt
sort covariables.txt > covariables_sorted.txt
#Eliminate samples not to be used (those not phenotyped).
awk 'FNR==NR {ids[$1]; next} $1 in ids' col1.txt covariables_sorted.txt > covariables_final.txt





#Install Fast_LMM in a new environment
conda activate fastlmm
pip install fastlmm

#Run fastlmm1.py
PYTHONPATH=/home/andres/.local/lib/python3.10/site-packages python3.10 fastlmm1.py
