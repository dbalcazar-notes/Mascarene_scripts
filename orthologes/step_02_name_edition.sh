#!/bin/bash
# Directorio donde están los resultados de TransDecoder
TRANSDECODER_DIR="."  # apunta al directorio actual
OUT_DIR="split_by_species_fna/fixed_for_orthograph"
mkdir -p "$OUT_DIR"

# Recorre todos los subdirectorios que terminan en .transdecoder_dir
for dir in "$TRANSDECODER_DIR"/*.transdecoder_dir; do
    # Buscar el CDS y el PEP dentro del subdirectorio
    cds_file=$(find "$dir" -maxdepth 1 -name "*.cds" | head -n1)
    pep_file=$(find "$dir" -maxdepth 1 -name "*.pep" | head -n1)

    # Revisar que existan
    if [[ -z "$cds_file" ]]; then
        echo "No se encontró CDS en $dir, salteando..."
        continue
    fi
    if [[ -z "$pep_file" ]]; then
        echo "No se encontró PEP en $dir, salteando..."
        continue
    fi

    # Tomar nombre base para salida (FORM, MASC, joint, etc.)
    base=$(basename "$dir" .fna.transdecoder_dir)

    # Corregir headers y generar .fna y .pep
    sed -E 's/^>(\S+)\.p1.*/>\1/' "$cds_file" > "${OUT_DIR}/${base}.fna"
    sed -E 's/^>(\S+)\.p1.*/>\1/' "$pep_file" > "${OUT_DIR}/${base}.pep"

    echo "Procesado $base"
done

echo "Todos los archivos listos en $OUT_DIR"
