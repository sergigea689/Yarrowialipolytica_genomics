#!/bin/bash
#SBATCH -p short
#SBATCH -n 32
#SBATCH -J mafft
#SBATCH --mem=32072
#SBATCH -t 0-23:30:15   
#SBATCH -o slurm/braker.%j.out
#SBATCH -e slurm/braker.%j.err
#Loop to iterate over concatenate files to align with maaft
#!/bin/bash

# Directorio de salida para los alineamientos
OUTPUT_DIR="aligned_parts"
mkdir -p "$OUTPUT_DIR"

# Rango de archivos a procesar (modificar según el trabajo)
START=6   # Cambia para cada trabajo (Ej: 6, 11)
END=10     # Cambia para cada trabajo (Ej: 10, 16)

# Ejecutar MAFFT en cada archivo
for i in $(seq $START $END); do
    INPUT_FILE="part_${i}.fasta"
    OUTPUT_FILE="${OUTPUT_DIR}/aligned_${i}.fasta"

    if [[ -f "$INPUT_FILE" ]]; then
        echo "Alineando $INPUT_FILE..."
        mafft --thread 20 --auto "$INPUT_FILE" > "$OUTPUT_FILE"
    else
        echo "¡Advertencia! No se encontró $INPUT_FILE"
    fi
done

echo "Alineamientos completados para archivos part_${START}.fasta a part_${END}.fasta"
