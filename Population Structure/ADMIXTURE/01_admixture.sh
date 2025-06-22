#Admixture
#Requires preparation of the vcf file to take out unwanted samples, filter by quality, remove indels and sites with less than 2 alleles and keep only biallellic variants.
vcftools --vcf hqvariants_yayaout.recode.vcf --remove remove.txt --recode --recode-INFO-all --out small-remove.vcf --max-missing 1 --non-ref-ac-any 1

vcftools --max-missing 1 --max-alleles 2 --vcf small-remove.vcf.recode.vcf --remove-indels --recode --recode-INFO-all --out biallelicvariants --non-ref-ac-any 1
#LD pruning
 ./plink_pruning_prep.sh ../biallelicvariants.recode.vcf
 plink --vcf biallelicvariants.recode_annot.vcf --double-id --allow-extra-chr --indep-pairwise 50 10 0.2 --maf 0.05 --out ld_pruned --make-bed --threads 2
 
 awk '{$1=0; print $0}' ld_pruned.bim > ld_pruned.bim.tmp
 mv ld_pruned.bim.tmp ld_pruned.bim
 
 #Admixture

 for k in {1..6}
	do
		./admixture32 --cv ld_pruned.bed ${k} > admixture.${k}.log
	done	
	
#In a line:
for k in {1..6}; do ./admixture32 --cv ld_pruned.bed ${k} > admixture.${k}.log; done

#We can view the results with pong
pip3 install pong
conda deactivate
#Run pong
pong -m pong_filemap -i ind2pop.txt -n listlabelstreeorderR.txt -l colours.txt