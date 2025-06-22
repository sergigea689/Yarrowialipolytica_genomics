import re

# Input and output files
fasta_file = "cds.fa"              # Original FASTA file
new_fasta_file = "cds_renamed.fa"  # New FASTA file with renamed headers

# Read the FASTA file and rename headers
with open(fasta_file, "r") as fasta, open(new_fasta_file, "w") as new_fasta:
    for line in fasta:
        if line.startswith(">"):  # If it's a header line
            header = line.strip()
            # Extract the locus_tag from the header
            match = re.search(r'locus_tag=([^\s\]]+)', header)
            if match:
                locus_tag = match.group(1)
                # Create the new header with the desired format
                new_header = f">TRANSCRIPT_gene-{locus_tag} " + header[1:]  # Remove '>' and add the new format
                new_fasta.write(new_header + "\n")
            else:
                new_fasta.write(header + "\n")  # If locus_tag is not found, keep the original header
        else:
            new_fasta.write(line)  # Write sequences without changes