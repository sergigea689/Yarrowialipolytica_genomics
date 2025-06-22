#This script aims at analazying copy number variants (CNV) of BAM files as an input

#Se ha de tener el archivo de referencia, en este caso YALI0-1.fasta, por cromosomas
samtools faidx YALI0-1.fasta YALI0A > YALI0A.fasta
#Se repite por cada cromosoma.
#Se necesita el archivo de anotación de la referencia en formato TSV
#Hay que adaptarlo ya que los archivos de anotación que tenemos están en formato embl. https://entrepot.recherche.data.gouv.fr/dataset.xhtml?persistentId=doi:10.57745/UZFAEN
#Lo siguiente es un script para scrapear las partes de interés de la anotación y hacer un tsv

#!/bin/bash

# Nombre del archivo de anotación
input_file=".57745%2F5X3X9R"
output_file="yali0A.tsv"

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

#Por si se olvida cambiar el cromosoma
sed 's/YALI0A/YALI0F/g' yali0F_g.tsv > yali0F_2.tsv

#Después de este script, es necesario realizar una transformación de las anotaciones para eliminar información sobrante
awk -F'\t' '{gsub(/ FT[^\t ]*/, "", $7); gsub(/[ \t]+/, " ", $7); gsub(/^ +| +$/, "", $7); print}' yali0E.tsv > yali0E_g.tsv

#Generar el archivo de anotaciones completo
cat yali0A_g.tsv yali0B_2.tsv yali0C_2.tsv yali0D_2.tsv yali0E_2.tsv yali0F_2.tsv > anno_yali0.tsv

#Se necesita un archivo de plodías con una columna con los nombres de las cepas y otra columna con las plodías. El separador debe ser ;. importante que tenga las columnas que tiene (4) aunque solo se vaya a usar dos.
#Parece que poner Standardized name da problema, así que lo cambié por name y añadí las columnas restantes de forma manual.
sed 's/,/;/g' yaliploidy.csv > archivo_convertido.csv

#Detener todos los procesos en segundo plano relacionados con tu script:
ps aux | grep '[c]nfreec_yali.py' | awk '{print $2}' | xargs kill

#Se tiene que modificar la función split_fasta, que forma parte de launch desde el archivo freec.py en el wrapper de mappingtools
#con el objetivo de nombrar los archivos de los cromosomas de forma correcta
def split_fasta(fasta_file, output_dir):
    """
    Divide un archivo FASTA en varios archivos, uno por cada cromosoma,
    y guarda cada archivo en el directorio de salida especificado.
    El nombre del archivo se basa en la primera palabra del encabezado de cada cromosoma.
    """
    fdata = read_fasta(fasta_file)
    for contig, sequence in fdata.items():
        # Extrae la primera palabra del encabezado
        chrom_name = contig.split()[0]  # Usa la primera palabra del encabezado como nombre del archivo
        fasta_file_name = f"{chrom_name}.fa"
        fasta_file_path = os.path.join(output_dir, fasta_file_name)
        
        # Escribe la secuencia del cromosoma en un archivo
        sfdata = {chrom_name: sequence}
        write_fasta(sfdata, fasta_file_path)

    os.remove(fasta_file)
#La original
def split_fasta(fname) :
    fdata = read_fasta(fname)
    for contig, sequence in fdata.items() :
        sfdata = {contig : sequence}
        outfile = os.path.join(dname(fname), contig + ".fa")
        write_fasta(sfdata, outfile)

    os.remove(fname)
#implica modificar la función launch en la llamada
split_fasta(fasta_file) pasa a ser split_fasta(fasta_file, fasta_dir)