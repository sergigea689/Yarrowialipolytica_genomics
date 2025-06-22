#!/bin/bash

# Nombre del archivo de anotación
input_file="ann_yali0D.57745%2FXQ6GVR"
output_file="yali0D.tsv"

# Inicializamos variables
chr="YALI0A"  # Nombre del cromosoma, en este caso es constante
count=1  # Contador para Nº gen/chr



in_note=false  # Variable para saber si estamos dentro de una nota multilinea
note="None"  # Inicializar note como None

# Creamos el archivo de salida con encabezado
echo -e "Nº gen/chr\tchr\tStart\tEnd\tKind\tName\tNote" > "$output_file"

# Leemos el archivo línea por línea
while IFS= read -r line; do
    # Verificamos si es un tipo de característica (source, misc_RNA, CDS, mRNA, rRNA, etc.)
    if [[ $line =~ ^FT\ +(source|CDS|mRNA|misc_RNA|rRNA) ]]; then
        # Guardamos cualquier anotación anterior antes de comenzar una nueva
        if [[ $kind != "None" ]]; then
            echo -e "$count\t$chr\t${start:-None}\t${end:-None}\t${kind:-None}\t${name:-None}\t${note:-None}" >> "$output_file"
            count=$((count + 1))
        fi
        
        # Capturamos el tipo de elemento
        kind=$(echo "$line" | awk '{print $2}')
        
        # Capturamos el Start y End, manejando tanto formatos normales como complementarios
        location=$(echo "$line" | awk '{print $3}')
        
        # Verificamos si es un formato de 'join(...)' para manejar múltiple posiciones
        if [[ $location =~ join ]]; then
            start=$(echo "$location" | grep -oP '\d+' | head -n 1 || echo "None")
            end=$(echo "$location" | grep -oP '\d+' | tail -n 1 || echo "None")
        else
            start=$(echo "$location" | grep -oP '(\d+)(?=\.\.)' || echo "None")
            end=$(echo "$location" | grep -oP '(?<=\.\.)(\d+)' || echo "None")
        fi

        # Reseteamos campos
        name="None"
        note="None"
        in_note=false

    # Capturamos el Name (locus_tag)
    elif [[ $line =~ /locus_tag= ]]; then
        name=$(echo "$line" | sed 's/.*locus_tag="//;s/".*//')

      # Capturamos el Note, considerando el caso de nota multilinea
    elif [[ $line =~ /note=\" ]]; then
        in_note=true
        note=$(echo "$line" | sed 's/.*note="//;s/".*//')

        # Verificamos si el note termina en esta línea
        if [[ ! $line =~ \"[^\"]*$ ]]; then
            note=${note%\"*}  # Removemos cualquier cierre de comillas final
            in_note=false
        fi

    # Continuamos acumulando el Note hasta encontrar el cierre de comillas
    elif $in_note; then
        # Agregamos la línea actual al contenido de la nota
        note="$note $(echo "$line" | sed 's/".*//')"
        
        # Verificamos si el Note termina en esta línea
        if [[ $line =~ \"[^\"]*$ ]]; then
            in_note=false
        fi
    # Cuando encontramos la línea de traducción o final del bloque de una anotación
    elif [[ $line =~ ^FT\ +/translation || $line =~ ^XX ]]; then
        # Guardamos la información en el archivo de salida, incluso si hay campos faltantes
        echo -e "$count\t$chr\t${start:-None}\t${end:-None}\t${kind:-None}\t${name:-None}\t${note:-None}" >> "$output_file"
        
        # Incrementamos el contador y reseteamos las variables específicas de la anotación
        count=$((count + 1))
        kind="None"
        start="None"
        end="None"
        name="None"
        note="None"
        in_note=false
    fi
done < "$input_file"

# Guardamos la última anotación si existe
if [[ $kind != "None" ]]; then
    echo -e "$count\t$chr\t${start:-None}\t${end:-None}\t${kind:-None}\t${name:-None}\t${note:-None}" >> "$output_file"
fi