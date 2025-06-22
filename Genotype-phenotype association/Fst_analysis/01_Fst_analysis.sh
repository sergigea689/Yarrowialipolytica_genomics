#Fst analysis of best and worst growers in acetate
#goal is to compare top 25 and bottom 25 growers in acetate with the Fst index, which is related to allele frequencies, to detect regions of divergence between this groups.
vcftools --vcf ../hqvariantswithyayanomutants.vcf.recode.vcf --keep acetatestrains.txt --minQ 30 --min-meanDP 20 --non-ref-ac-any 1 --recode --recode-INFO-all --out strainsacetate


#First vcftools is used to filter SNPs, individuals, etc. For the first analyisis we will use SNPs that were called in all individuals, and we only retain biallelic sites
vcftools --max-missing 1 --max-alleles 2 --vcf strainsacetate.recode.vcf --remove-indels --recode --recode-INFO-all --out filtered_snps_max 

#Fst analysis
vcftools --vcf biallelicvariants.recode.vcf   --weir-fst-pop SoHC.txt   --weir-fst-pop  popnoshc.txt   --out fst_results