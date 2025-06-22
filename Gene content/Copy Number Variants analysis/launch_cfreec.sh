#!/bin/bash

# Lee cada línea del archivo nombres_carpetas.txt
while IFS= read -r carpeta; do
    echo "Procesando carpeta: $carpeta"
    
    # Busca un archivo .bam en la carpeta
    bam_file=$(find "$carpeta" -maxdepth 1 -name "*.bam" | head -n 1)
    
    # Verifica si se encontró un archivo .bam
    if [[ -n "$bam_file" ]]; then
        echo "Archivo .bam encontrado: $bam_file"
        
        # Ejecuta el script de Python con el archivo .bam en segundo plano
        python3 cnfreec_yali.py -b "$bam_file" 
    else
        echo "No se encontró ningún archivo .bam en la carpeta: $carpeta"
    fi
    
done < nombres_carpetas.txt

# Espera a que todos los procesos en segundo plano terminen
wait