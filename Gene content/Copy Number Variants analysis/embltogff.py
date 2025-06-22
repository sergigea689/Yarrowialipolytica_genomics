from BCBio import GFF
from Bio import SeqIO

input_file = "YALI0F.embl"
output_file = "YALI0F.gff"

with open(input_file) as embl_handle:
    records = SeqIO.parse(embl_handle, "embl")
    with open(output_file, "w") as gff_handle:
        # Escribir sin valores de /translation
        GFF.write(records, gff_handle)