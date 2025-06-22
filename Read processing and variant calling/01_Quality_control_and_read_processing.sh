#Quality control and read processing of short-read Illumina paired-end reads before variant-calling

################################
#1. Quality control of raw reads with fastqc
################################

#Directory
directorio_entrada="/mnt/lustre/home/sergiz01/genomes/Ylipgenomes/fastq"
#Iterate over samples
# Output directroy
directorio_salida="/mnt/lustre/home/sergiz01/genomes/Ylipgenomes/fastq/fastqc"
# Iterate over .fq files
for archivo in "$directorio_entrada"/*.fq.gz; do

    # Fastqc
    fastqc "$archivo" -o "$directorio_salida"
done

#Inspect the output manually or use Multiqc to generate a summary report
conda install -c bioconda multiqc
multiqc fastqc/

################################
#2. Trimming with fastp
################################

cd /mnt/lustre/home/sergiz01/genomes/Ylipgenomes/fastq/rawreads
#Iterate over samples
for Variable in $(< listagenomes2.txt) #list with the name of the files without read number.fq.gz
do
fastp -i ./${Variable}1.fq.gz -o ./trimmed_reads/${Variable}1.fq.gz -I ./${Variable}2.fq.gz -O ./trimmed_reads/${Variable}2.fq.gz -q 20 --trim_poly_g --adapter_sequence=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT --adapter_sequence_r2=GATCGGAAGAGCACACGTCTGAACTCCAGTCACGGATGACTATCTCGTATGCCGTCTTCTGCTTG 
done

################################
#3. Alignement to referece genome with BWA-mem
################################

#We use the reference genome of E150 strain, the revised version found at https://doi.org/10.57745/UZFAEN.
cd /mnt/lustre/home/sergiz01/genomes/Ylipgenomes/fastq/rawreads/trimmed_reads
mkdir bwa
#Iterate over samples
for Variable in $(< ../listagenomes2.txt)
do
bwa mem ./ref/YALI0-1.fasta ./${Variable}1.fastq.gz ./${Variable}2.fastq.gz | samtools sort -o ./bwa/"$Variable".bam
done

################################
#4. Quality control of alignment with bamqc
################################

#multibammultibamqcinputT.txt is a file with file names on the first columns and path to .bam or .bamqc files
qualimap multi-bamqc -d multibamqcinputT.txt -outdir /media/usach/matadatitos1/Sergio/BWA/multibamqc

#Inspect results and act accordingly.

################################
#5. Indexing with samtools
################################

cd /media/usach/matadatitos1/Sergio/BWA/
mkdir MD #in preparation for next steps
cd /media/usach/matadatitos1/Sergio/BWA/MD

#In this case, the txt file must contain the .bam extension
#Iterate over samples
for Variable in $(cat tarchivos_bam.txt)
do
	samtools index ./${Variable}
done

################################
#6. Processing for variant calling with Picard
################################

#Add or Replace Read Groups
cd /media/usach/matadatitos1/Sergio/BWA #your directory with alignment files
for Variable in $(cat archivos_bam.txt) #List of file names without .bam extension
do
java -jar /home/usach/miniconda3/pkgs/picard-2.27.4-hdfd78af_0/share/picard-2.27.4-0/picard.jar AddOrReplaceReadGroups \
I= ./${Variable}.bam \
O= ./RG/${Variable}.RG.bam \
SORT_ORDER= coordinate \
RGID= DontKnow \
RGLB= 01 \
RGPL= illumina \
RGSM= ${Variable} \
RGPU= 001 \
CREATE_INDEX= True 
done

#Mark Duplicates
for Variable in $(cat archivos_bam.txt) #List of file names without .bam extension
do
java -jar /home/usach/miniconda3/pkgs/picard-2.27.4-hdfd78af_0/share/picard-2.27.4-0/picard.jar MarkDuplicates \
I= ./RG/${Variable}.RG.bam \
O=  ./MD/${Variable}.RG.MD.bam \
M= ./MD/${Variable}.RG.MD.rmdup.txt
done
