import sys
import os
import glob
import argparse  # Nueva importación para gestionar argumentos

import pandas as pd
pd.set_option('display.width', 1000)

sys.path.insert(0, "/mnt/c/Users/USUARIO/Desktop/genomes/MappingTools-master/mappingTools")
#sys.path.insert(0, "/mnt/lustre/home/javievic/batch_scripts/variantcall/cnv/")

#import files
from wrappers.freec import launch_bam
from wrappers.processing import multiprocess_this, subprocess_this
from freec_analysis import process_fname

from Bio import SeqIO
from collections import Counter

dname = os.path.dirname
bname = os.path.basename

freec_p = "/mnt/c/Users/USUARIO/Desktop/genomes/FREEC-11.6/src/freec"
samtools_p = "samtools"
sample_funname = lambda fname : bname(fname).split(".")[0]

def get_file(name) :
	fname = os.path.realpath(__file__)
	return os.path.join(dname(fname), name)

def launch_freec_args(fname, freec_p, samtools_p, reference, outdir, ignore_exist, config, annos, defaultcn) :
	print(f"Parámetros recibidos: {fname}, {freec_p}, {samtools_p}, {reference}, {ignore_exist}, {config}, {outdir}")
	# launch freec
	launch_bam(fname, freec_p, samtools_p, reference, ignore_exist=ignore_exist, config=config, outdir= outdir)

	# analyse result
	outfile = os.path.join(outdir, bname(fname) + "_CNVs")
	if not os.path.isfile(outfile) :
		print ("Fail for : %s" %(fname))
		return None

	return process_fname(outfile, annos, defaultcn, ignore_exist=False)

def get_min_max_gc(reference, window=1000) :
	fdata = SeqIO.parse(reference, "fasta")
	fdata = [str(sequence.seq) for sequence in fdata]

	values = []
	for sequence in fdata :
		for i in range(0, len(sequence), window) :
			subseq = sequence[i:i+window].lower()
			if len(subseq) < window : continue
			count = Counter(subseq)
			values.append((count.get("c", 0) + count.get("g", 0)) / len(subseq))

	return min(values), max(values)

def launch(bamfiles, ncore=20, ignore_exist=True):
	reference = get_file("/mnt/c/Users/USUARIO/Desktop/genomes/ref/YALI0-1.fasta")
	mingc, maxgc = get_min_max_gc(reference)

	ploidy = get_file("yaliploidy2.csv")
	ploidy = pd.read_csv(ploidy, sep=";", index_col=0)
	ploidy = ploidy.set_index("name")["Ploidy"].to_dict()
	ploidy = {sample: int(value) for sample, value in ploidy.items()}

	bamfiles = glob.glob(bamfiles) if isinstance(bamfiles, str) else bamfiles
	args = []

	annos = get_file("./ref/anno_yali0.tsv")
	annos = pd.read_csv(annos, sep="\t", index_col=0)

	# we remove the gff contig feature and remove the chr part due to Freec BUG
	# WARNING: IN THIS SCRIPT WE ONLY LOOK AT CDS. CHANGE IF NEED

	annos = annos[annos["Kind"] == "CDS"]
	annos["End"] = pd.to_numeric(annos["End"], errors='coerce')
	annos["Start"] = pd.to_numeric(annos["Start"], errors='coerce')
	annos = annos.dropna(subset=["End", "Start"])
	annos["Contig"] = annos["Contig"].str.slice(start=0)
	annos = annos[annos["End"] - annos["Start"] > 0]


	for fname in bamfiles :

		sample = sample_funname(fname)
		sample_ploidy = ploidy.get(sample)
		if not sample_ploidy : print ("Ploidy not found for : %s (%s)" %(sample, fname)); continue

		#config = {"sample" : {"mateOrientation" : "FR"}, "general" : {"ploidy" : sample_ploidy,
		#"minExpectedGC" : "%.1f" %(mingc), "maxExpectedGC" : "%.1f" %(maxgc)}}
		config = {"sample" : {"mateOrientation" : "0"}, "general" : {"ploidy" : sample_ploidy}}

		soutdir = fname + "6.freec"
		args.append((fname, freec_p, samtools_p, reference, soutdir, ignore_exist, config, annos, sample_ploidy))

	args = sorted(args, key = lambda x : x[0])
	res = multiprocess_this(launch_freec_args, args, ncore=ncore)

if __name__ == "__main__":
    # Configuración de argparse para manejar argumentos de línea de comandos
    parser = argparse.ArgumentParser(description="Run Control-FREEC for CNV analysis on BAM files.")
    parser.add_argument("-b", "--bam", required=True, help="Pattern or list of BAM files to process")
    parser.add_argument("-n", "--ncore", type=int, default=1, help="Number of cores to use for multiprocessing")
    
    # Parseo de argumentos
    args = parser.parse_args()
    
    # Obtener valores de los argumentos
    bam = args.bam
    ncore = args.ncore

    # Lanzar el análisis
    resfiles = launch(bam, ncore=ncore)