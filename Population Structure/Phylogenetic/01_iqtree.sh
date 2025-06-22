#Tree making


conda install -c bioconda iqtree
git clone https://github.com/edgardomortiz/vcf2phylip.git

#Filtering the vcf to take out the unwanted strains
vcftools --vcf hqvariants_yayaout.recode.vcf --remove remove.txt --recode --recode-INFO-all --out small-remove.vcf --max-missing 1 --non-ref-ac-any 1

vcftools --max-missing 1 --max-alleles 2 --vcf small-remove.vcf.recode.vcf --remove-indels --recode --recode-INFO-all --out biallelicvariants --non-ref-ac-any 1

#First we transform the VCF file to phylip format using a custom python script developed by edgardomortiz, and then we use it for IQTREE
#These are steps to follow, not a script
python3 vcf2phylip.py -i biallelicvariants.recode.vcf
iqtree -s filtered_outgroup_map50_snp.recode.min4.phy -st DNA -o Chis -m TEST -nt 8 #outgroup is CBS707 as it is a different species
iqtree -s biallelicvariants.recode.min4.phy -st DNA -o YL_124_concatenado -m GTR+F+G4 -nt 8 -bb 1000 -redo
iqtree -s filtered_outgroup_map50_snp.recode.min4.phy -st DNA -o Chis -m TEST -nt 8 -bb 1000
