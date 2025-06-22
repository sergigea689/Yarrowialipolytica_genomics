## FINESTRUCTURE ANALYSIS PIPELINE

# STEP 1: Filter VCF to retain only strains of interest (including outgroup)
vcftools --vcf hqvariants_yayaout.recode.vcf --keep poppurek6.txt --recode --recode-INFO-all --non-ref-ac-any 1 --out purek6.vcf

# STEP 2: Filter SNPs - retain only biallelic sites called in all individuals
vcftools --max-missing 1 --max-alleles 2 --vcf purek6.vcf.recode.vcf --recode --recode-INFO-all --out biallelicvariants

# STEP 3: Optional - rename chromosome names (if needed)
sed -i 's/YALI0A/1/g' biallelicvariants.recode.vcf

# STEP 4: Filter for Linkage Disequilibrium (LD) using custom script
wget https://github.com/joanam/scripts/raw/master/ldPruning.sh
chmod +x ldPruning.sh
./ldPruning.sh biallelicvariants.recode.vcf
gzip -d biallelicvariants.recode.LDpruned.vcf.gz

# STEP 5: Install required tools
conda install snpsift plink
conda install -c compbiocore perl-switch

# STEP 6: Download and prepare FineSTRUCTURE and phasing tools
wget https://people.maths.bris.ac.uk/~madjl/finestructure/fs_4.1.1.zip
wget https://github.com/gusevlab/germline/raw/master/phasing_pipeline.tar.gz
tar -xf phasing_pipeline.tar.gz
wget https://faculty.washington.edu/browning/beagle/recent.versions/beagle_3.0.4_05May09.zip
# → Copy beagle.jar into the phasing pipeline folder #Not needed if working with haploids

# STEP 7: Split VCF file by chromosome using SnpSift
SnpSift split biallelicvariants.recode.LDpruned.vcf

# STEP 8: Convert per-chromosome VCFs to PLINK format
for i in {1..6}; do
    chr=$(printf "%X" $i)
    plink --vcf biallelicvariants.recode.LDpruned.$i.vcf \
          --recode12 --allow-extra-chr --double-id --geno 1 --out YALI0${chr}
done

# STEP 9: Convert PLINK files to ChromoPainter format (outside conda env)
# Ensure Switch.pm is in the working directory
for i in A B C D E F; do
    perl ./fs_4.1.1/plink2chromopainter.pl -p=YALI0${i}.ped -m=YALI0${i}.map -o=YALI0${i}.chromopainter -f
done

# STEP 10: Convert ChromoPainter input to haploid format (ONLY for haploid organisms)
for i in A B C D E F; do
    (awk 'NR == 1 { print $1 / 2 }' YALI0${i}.chromopainter;
     sed '2,3!d' YALI0${i}.chromopainter;
     sed '1,3d' YALI0${i}.chromopainter | sed '0~2d') > YALI0${i}.chromopainter.haploid
done

# STEP 11: Generate recombination rate files (use constant recombination rate of 0.000004)
# Edit line 47 of makeuniformrecfile.pl to set this constant
for i in A B C D E F; do
    perl ./fs_4.1.1/makeuniformrecfile.pl YALI0${i}.chromopainter.haploid YALI0${i}_rec.chromopainter.haploid
done

# STEP 12: Run ChromoPainter (v2) in all-vs-all mode (per chromosome)
sudo apt install libgomp1
for i in A B C D E F; do
    ./fs_4.1.1/fs_linux_glibc2.3 chromopainter \
        -g YALI0${i}.chromopainter.haploid \
        -r YALI0${i}_rec.chromopainter.haploid \
        -t cut.txt \
        -o cp_YALI0${i} -a 0 0 -j 1
done

# STEP 13: Combine chromosome-specific ChromoPainter outputs
mkdir datos
mv cp_* datos
./fs_4.1.1/fs_linux_glibc2.3 chromocombine -d datos

# STEP 14: Run FineSTRUCTURE MCMC analysis on the combined output
./fs_4.1.1/fs_linux_glibc2.3 finestructure -x 100000 -y 100000 -z 1000 output.chunkcounts.out out.YALI.mcmc.xml

# STEP 15: Run FineSTRUCTURE tree building (if needed)
./fs_4.1.1/fs_linux_glibc2.3 finestructure \
    -m T -k 2 -T 1 -e popidfile \
    -x 100000 -y 100000 -z 1000 \
    output.chunkcounts.out \
    out.YALI.mcmc.xml \
    out.YALI.tree.xml






 
