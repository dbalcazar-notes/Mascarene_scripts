#!/bin/bash
#SBATCH --job-name=TransDecoder
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=10G
#SBATCH --partition=day
#SBATCH --time=10:00:00


# Directorio con archivos .fna por especie
IN_DIR="split_by_species_fna"

# Ejecutar TransDecoder para cada archivo .fna
for fna in "$IN_DIR"/*.fna; do
    base=$(basename "$fna" .fna)

    echo "Procesando $base"

    # Identificar ORFs largos
    TransDecoder.LongOrfs -t "$fna"

    # Predecir CDS/proteínas
    TransDecoder.Predict -t "$fna"

    echo "Terminado $base"
done

echo "Todos los archivos fueron procesados con TransDecoder."
