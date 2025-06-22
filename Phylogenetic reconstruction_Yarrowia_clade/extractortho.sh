#!/bin/bash

# Configura rutas y archivos
ORTHOLOGS_FILE="orthologs.tsv"   # Archivo con la tabla de ortólogos
FASTA_DIR="/mnt/lustre/home/sergiz01/genomes/Ylipgenomes/fastq/rawreads/trimmed_reads/assemblies/annotation/protfasta/weird/format"          # Directorio con las secuencias individuales
OUTPUT_FILE="concatenated_sequences.fasta"  # Archivo de salida concatenado
: > "$OUTPUT_FILE"  # Limpiar el archivo de salida si ya existe

# Leer encabezado del archivo orthologs.tsv
header=$(head -n 1 "$ORTHOLOGS_FILE")
column_names=($(echo "$header" | cut -f4-)) # Cepas desde la cuarta columna

# Iterar sobre las filas del archivo orthologs.tsv (excluyendo el encabezado)
tail -n +2 "$ORTHOLOGS_FILE" | while read -r line; do
    # Extraer la fila de genes (desde la 4ª columna)
    genes=$(echo "$line" | cut -f4-)
    
    # Convertir los genes de la fila en un array
    read -ra gene_array <<< "$genes"
    
    # Crear un identificador para esta entrada (genérico, puede ajustarse)
    entry_id="entry_$(echo "$line" | cut -f1)"  # Usar el nombre de la primera columna para identificar la entrada
    
    # Imprimir el encabezado en formato FASTA
    echo ">$entry_id" >> "$OUTPUT_FILE"
    
    # Concatenar las secuencias de cada ortólogo para esta entrada
    concatenated_sequence=""
    for i in "${!column_names[@]}"; do
        aa_file="${FASTA_DIR}/${column_names[$i]}"  # Archivo FASTA para la cepa
        gene_id="${gene_array[$i]}"                 # ID del gen correspondiente
        
        # Verificar si el archivo FASTA existe
        if [[ -f "$aa_file" ]]; then
            # Buscar la secuencia correspondiente al ID del gen
            sequence=$(grep -A 1 -w ">${gene_id}" "$aa_file" | tail -n 1)
            
            # Si no se encuentra la secuencia, añadir gaps del tamaño promedio
            if [[ -z "$sequence" ]]; then
                echo "¡Advertencia! No se encontró $gene_id en $aa_file. Rellenando con gaps."
                avg_length=$(grep -v "^>" "$aa_file" | awk '{total+=length($0); count++} END {if (count > 0) print int(total/count); else print 0}')
                sequence=$(printf '%*s' "$avg_length" | tr ' ' '-')
            fi
            
            # Concatenar la secuencia
            concatenated_sequence+="$sequence"
        else
            echo "¡Error! No se encontró el archivo $aa_file."
        fi
    done

    # Escribir la secuencia concatenada para esta entrada en el archivo de salida
    echo "$concatenated_sequence" >> "$OUTPUT_FILE"
done

echo "Proceso completado. Archivo concatenado guardado en $OUTPUT_FILE."
