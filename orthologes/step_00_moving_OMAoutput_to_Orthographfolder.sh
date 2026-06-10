#!/bin/bash

# Directorios de entrada y salida
OG_DIR="OrthologousGroupsFasta"
OUT_TABLE="ogs_table.tsv"
OUT_FASTA="split_by_species_fna"

# Crear directorio de salida y vaciar tabla si existe
mkdir -p "$OUT_FASTA"
> "$OUT_TABLE"

# Procesar cada OG
for ogfile in "$OG_DIR"/OG*.fa; do
    OG=$(basename "$ogfile" .fa)

    awk -v og="$OG" -v out_table="$OUT_TABLE" -v out_fa="$OUT_FASTA" '
    BEGIN { 
        # Array para guardar genes ya escritos por especie
        PROCINFO["sorted_in"] = "@ind_str_asc"
    }
    /^>/ {
        header=$0
        sub(/^>/,"",header)
        split(header,a," ")
        gene=a[1]
        match(header,/\[(.*)\]/,m)
        taxon=m[1]

        # Tomar el código corto: penúltimo segmento después de "_"
        split(taxon,b,"_")
        if(length(b) > 1){
            short_tax=b[length(b)-1]
        } else {
            short_tax=b[length(b)]
        }

        # Saltar si este gen ya fue agregado para esta especie
        key=short_tax":"gene
        if(seen[key]++) next

        # Escribir fila en tabla
        print og, gene, short_tax >> out_table

        # Definir archivo FASTA de salida (.fna para CDS)
        outfile=out_fa "/" short_tax ".fna"
        print ">"gene >> outfile
        next
    }
    {
        print >> outfile
    }' "$ogfile"
done

echo "Procesamiento completo. Archivos de especies en '$OUT_FASTA' y tabla en '$OUT_TABLE'."
