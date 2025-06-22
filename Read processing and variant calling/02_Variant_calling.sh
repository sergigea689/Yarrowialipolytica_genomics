#Variant calling pipeline based on GATK

################################
#1. Build dictionary for reference fasta
################################

gatk CreateSequenceDictionary -R YALI0-1.fasta

################################
#2. Variant callign with HaplotypeCaller
################################

#ChrA
cd /media/usach/matadatitos1/Sergio/BWA/MD
for SAMPLE in $(< archivos_bam.txt)
do
/home/usach/miniconda3/pkgs/gatk4-4.3.0.0-py36hdfd78af_0/bin/gatk  HaplotypeCaller --intervals YALI0A --emit-ref-confidence GVCF -R refgenome/YALI0-1.fasta -I ${SAMPLE}.RG.MD.bam -O YALI0A_chr/${SAMPLE}_YALI0A.vcf
done 
sleep 30
#ChrB
cd /media/usach/matadatitos1/Sergio/BWA/MD
for SAMPLE in $(< archivos_bam.txt)
do
/home/usach/miniconda3/pkgs/gatk4-4.3.0.0-py36hdfd78af_0/bin/gatk  HaplotypeCaller --intervals YALI0B --emit-ref-confidence GVCF -R refgenome/YALI0-1.fasta -I ${SAMPLE}.RG.MD.bam -O YALI0B_chr/${SAMPLE}_YALI0B.vcf
done 
sleep 30
#ChrC
cd /media/usach/matadatitos1/Sergio/BWA/MD
for SAMPLE in $(< archivos_bam.txt)
do
/home/usach/miniconda3/pkgs/gatk4-4.3.0.0-py36hdfd78af_0/bin/gatk  HaplotypeCaller --intervals YALI0C --emit-ref-confidence GVCF -R refgenome/YALI0-1.fasta -I ${SAMPLE}.RG.MD.bam -O YALI0C_chr/${SAMPLE}_YALI0C.vcf
done 
sleep 30
#ChrD
cd /media/usach/matadatitos1/Sergio/BWA/MD
for SAMPLE in $(< archivos_bam.txt)
do
/home/usach/miniconda3/pkgs/gatk4-4.3.0.0-py36hdfd78af_0/bin/gatk  HaplotypeCaller --intervals YALI0D --emit-ref-confidence GVCF -R refgenome/YALI0-1.fasta -I ${SAMPLE}.RG.MD.bam -O YALI0D_chr/${SAMPLE}_YALI0D.vcf
done
sleep 30
#CHrE
cd /media/usach/matadatitos1/Sergio/BWA/MD
for SAMPLE in $(< archivos_bam.txt)
do
/home/usach/miniconda3/pkgs/gatk4-4.3.0.0-py36hdfd78af_0/bin/gatk  HaplotypeCaller --intervals YALI0E --emit-ref-confidence GVCF -R refgenome/YALI0-1.fasta -I ${SAMPLE}.RG.MD.bam -O YALI0E_chr/${SAMPLE}_YALI0E.vcf
done
sleep 30
#ChrF
cd /media/usach/matadatitos1/Sergio/BWA/MD
for SAMPLE in $(< archivos_bam.txt)
do
/home/usach/miniconda3/pkgs/gatk4-4.3.0.0-py36hdfd78af_0/bin/gatk  HaplotypeCaller --intervals YALI0F --emit-ref-confidence GVCF -R refgenome/YALI0-1.fasta -I ${SAMPLE}.RG.MD.bam -O YALI0F_chr/${SAMPLE}_YALI0F.vcf
done

################################
#3. Building database with GenomicsDBImport
################################

cd /media/usach/matadatitos1/Sergio/BWA/MD
/home/usach/miniconda3/pkgs/gatk4-4.3.0.0-py36hdfd78af_0/bin/gatk GenomicsDBImport --genomicsdb-workspace-path /media/usach/matadatitos1/Sergio/BWA/MD/secall_the_variants_A/ --intervals YALI0A   --sample-name-map YALI0A_chr/listYALI0A.txt &
/home/usach/miniconda3/pkgs/gatk4-4.3.0.0-py36hdfd78af_0/bin/gatk  GenomicsDBImport --genomicsdb-workspace-path /media/usach/matadatitos1/Sergio/BWA/MD/secall_the_variants_B/ --intervals YALI0B   --sample-name-map  YALI0B_chr/listYALI0B.txt &
/home/usach/miniconda3/pkgs/gatk4-4.3.0.0-py36hdfd78af_0/bin/gatk GenomicsDBImport --genomicsdb-workspace-path /media/usach/matadatitos1/Sergio/BWA/MD/secall_the_variants_C/ --intervals YALI0C   --sample-name-map  YALI0C_chr/listYALI0C.txt &
/home/usach/miniconda3/pkgs/gatk4-4.3.0.0-py36hdfd78af_0/bin/gatk GenomicsDBImport --genomicsdb-workspace-path /media/usach/matadatitos1/Sergio/BWA/MD/secall_the_variants_D/ --intervals YALI0D   --sample-name-map  YALI0D_chr/listYALI0D.txt &
/home/usach/miniconda3/pkgs/gatk4-4.3.0.0-py36hdfd78af_0/bin/gatk  GenomicsDBImport --genomicsdb-workspace-path /media/usach/matadatitos1/Sergio/BWA/MD/secall_the_variants_E/ --intervals YALI0E   --sample-name-map YALI0E_chr/listYALI0E.txt &
/home/usach/miniconda3/pkgs/gatk4-4.3.0.0-py36hdfd78af_0/bin/gatk GenomicsDBImport --genomicsdb-workspace-path /media/usach/matadatitos1/Sergio/BWA/MD/secall_the_variants_F/ --intervals YALI0F   --sample-name-map  YALI0F_chr/listYALI0F.txt &

################################
#4. Calling genotypes with GenotypeGVCFs
################################

cd /media/usach/matadatitos1/Sergio/BWA/MD
/home/usach/miniconda3/pkgs/gatk4-4.3.0.0-py36hdfd78af_0/bin/gatk GenotypeGVCFs -R ./refgenome/YALI0-1.fasta -V gendb://secall_the_variants_A/  -G StandardAnnotation -O /media/usach/matadatitos1/Sergio/BWA/MD/GVCF/secvariants_chrA.vcf 
/home/usach/miniconda3/pkgs/gatk4-4.3.0.0-py36hdfd78af_0/bin/gatk GenotypeGVCFs -R ./refgenome/YALI0-1.fasta  -V gendb://secall_the_variants_B/  -G StandardAnnotation -O /media/usach/matadatitos1/Sergio/BWA/MD/GVCF/secvariants_chrB.vcf 
/home/usach/miniconda3/pkgs/gatk4-4.3.0.0-py36hdfd78af_0/bin/gatk  GenotypeGVCFs -R ./refgenome/YALI0-1.fasta  -V gendb://secall_the_variants_C/  -G StandardAnnotation -O /media/usach/matadatitos1/Sergio/BWA/MD/GVCF/secvariants_chrC.vcf 
/home/usach/miniconda3/pkgs/gatk4-4.3.0.0-py36hdfd78af_0/bin/gatk  GenotypeGVCFs -R ./refgenome/YALI0-1.fasta -V gendb://secall_the_variants_D/  -G StandardAnnotation -O /media/usach/matadatitos1/Sergio/BWA/MD/GVCF/secvariants_chrD.vcf 
/home/usach/miniconda3/pkgs/gatk4-4.3.0.0-py36hdfd78af_0/bin/gatk  GenotypeGVCFs -R ./refgenome/YALI0-1.fasta  -V gendb://secall_the_variants_E/  -G StandardAnnotation -O /media/usach/matadatitos1/Sergio/BWA/MD/GVCF/secvariants_chrE.vcf 
/home/usach/miniconda3/pkgs/gatk4-4.3.0.0-py36hdfd78af_0/bin/gatk  GenotypeGVCFs -R ./refgenome/YALI0-1.fasta  -V gendb://secall_the_variants_F/  -G StandardAnnotation -O /media/usach/matadatitos1/Sergio/BWA/MD/GVCF/secvariants_chrF.vcf 

################################
#5. Merging variants with MergeVcfs
################################
cd /media/usach/matadatitos1/Sergio/BWA/MD/GVCF

/home/usach/miniconda3/pkgs/gatk4-4.3.0.0-py36hdfd78af_0/bin/gatk MergeVcfs -I /media/usach/matadatitos1/Sergio/BWA/MD/GVCF/secvariants_chrA.vcf -I /media/usach/matadatitos1/Sergio/BWA/MD/GVCF/secvariants_chrB.vcf -I /media/usach/matadatitos1/Sergio/BWA/MD/GVCF/secvariants_chrC.vcf -I /media/usach/matadatitos1/Sergio/BWA/MD/GVCF/secvariants_chrD.vcf -I /media/usach/matadatitos1/Sergio/BWA/MD/GVCF/secvariants_chrE.vcf -I /media/usach/matadatitos1/Sergio/BWA/MD/GVCF/secvariants_chrF.vcf -O all_raw_secvariants_yayanochis.vcf

################################
#5. Quality filtering with vcftools
################################

vcftools --vcf all_raw_variants_outgroup.vcf --recode --recode-INFO-all --minQ 30 --min-meanDP 20 --out FILE_outgroup
vcftools --max-missing 1 --max-alleles 2 --vcf remove_noyaya_noChis_map50.recode.vcf --recode --recode-INFO-all --out biallelicvariants_yayaout