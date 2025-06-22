import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import os
import math

# Configura rutas y archivos
ORTHOLOGS_FILE = "orthologs.tsv"  # Archivo con ortólogos
FASTA_DIR = "/mnt/lustre/home/sergiz01/genomes/Ylipgenomes/fastq/rawreads/trimmed_reads/assemblies/annotation/protfasta/weird/format"
OUTPUT_PREFIX = "part_"  # Prefijo para los archivos de salida

# Leer el archivo orthologs.tsv
df = pd.read_csv(ORTHOLOGS_FILE, sep="\t")

# Obtener las cepas desde la cuarta columna
cepas = df.columns[3:]

# Dividir en bloques de 100 filas
num_blocks = math.ceil(len(df) / 100)  # Número total de archivos
print(f"Dividiendo en {num_blocks} partes...")

for block_idx in range(num_blocks):
    start_row = block_idx * 100
    end_row = (block_idx + 1) * 100
    df_block = df.iloc[start_row:end_row]  # Sub-dataframe de 100 filas

    output_file = f"{OUTPUT_PREFIX}{block_idx + 1}.fasta"  # Nombre del archivo
    if os.path.exists(output_file):
        os.remove(output_file)  # Eliminar si existe

    # Procesar cada cepa dentro de este bloque
    for cepa in cepas:
        concatenated_sequence = ""
        print(f"Procesando cepa: {cepa} en {output_file}")

        seq_ids = df_block[cepa].dropna().tolist()  # IDs de secuencias en este bloque

        for seq_id in seq_ids:
            fasta_file = os.path.join(FASTA_DIR, f"{cepa}")

            if not os.path.exists(fasta_file):
                print(f"¡Advertencia! El archivo {fasta_file} no existe.")
                continue

            found_sequence = False
            for record in SeqIO.parse(fasta_file, "fasta"):
                if record.id == seq_id:
                    concatenated_sequence += str(record.seq)
                    found_sequence = True
                    break

            if not found_sequence:
                print(f"¡Advertencia! No se encontró {seq_id} en {fasta_file}. Rellenando con gaps.")
                avg_length = sum(len(rec.seq) for rec in SeqIO.parse(fasta_file, "fasta")) // max(1, len(list(SeqIO.parse(fasta_file, "fasta"))))
                concatenated_sequence += "-" * avg_length

        if concatenated_sequence:
            with open(output_file, "a") as output_handle:
                record = SeqRecord(concatenated_sequence, id=cepa, description="")
                SeqIO.write(record, output_handle, "fasta")

print("Proceso completado. Archivos generados:")
print("\n".join([f"{OUTPUT_PREFIX}{i+1}.fasta" for i in range(num_blocks)]))

