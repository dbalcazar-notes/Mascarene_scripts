import csv
import gzip
import os
import argparse

def procesar_vcf(input_path, output_path):
    with gzip.open(input_path, 'rt') as vcf_file, open(output_path, 'w', newline='') as csv_file:
        vcf_reader = csv.reader(vcf_file, delimiter='\t')

        # Leer el encabezado y encontrar los nombres de los individuos
        header = next(vcf_reader)
        while header[0].startswith('##'):
            header = next(vcf_reader)

        individuals = header[9:]

        # Escribir encabezados en el archivo CSV
        csv_writer = csv.writer(csv_file, delimiter='\t')
        csv_writer.writerow(['Chromosome', 'Position', 'Reference Allele'] + individuals)

        # Iterar sobre las líneas del archivo VCF
        for line in vcf_reader:
            chromosome = line[0]
            position = line[1]
            ref_allele = line[3]

            # Obtener información de cada individuo en las columnas restantes
            individuals_data = [1 if './.' not in gt else 0 for gt in line[9:]]

            # Escribir la fila en el archivo CSV
            csv_writer.writerow([chromosome, position, ref_allele] + individuals_data)

def main():
    parser = argparse.ArgumentParser(description='Convertir archivos VCF a CSV.')
    parser.add_argument('input_file', help='Ruta del archivo VCF de entrada')
    parser.add_argument('output_file', help='Ruta del archivo CSV de salida')
    args = parser.parse_args()

    input_path = args.input_file
    output_path = args.output_file

    procesar_vcf(input_path, output_path)

if __name__ == "__main__":
    main()
