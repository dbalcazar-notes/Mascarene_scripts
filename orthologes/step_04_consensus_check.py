import os
import argparse
from Bio import SeqIO
import shutil
import csv

def check_fasta(fasta_file, max_n_ratio=0.2):
    """Revisa si el archivo FASTA tiene problemas: vacío, solo encabezado, o alta proporción de Ns."""
    
    if not os.path.exists(fasta_file):
        return False, "Archivo no encontrado", None, None
    
    with open(fasta_file, "r") as file:
        records = list(SeqIO.parse(file, "fasta"))
        
        if len(records) == 0:
            return False, "Archivo vacío o solo encabezado", None, None
        
        for record in records:
            sequence = str(record.seq)
            seq_length = len(sequence)
            n_count = sequence.count('N')
            n_ratio = n_count / seq_length if seq_length > 0 else 1
            
            if n_ratio > max_n_ratio:
                return False, f"Alta proporción de N ({n_ratio:.2%})", seq_length, n_ratio
    
    return True, "OK", seq_length, n_ratio

def main(fasta_dir, max_n_ratio):
    fail_dir = os.path.join(fasta_dir, "fail")
    os.makedirs(fail_dir, exist_ok=True)

    report_file = os.path.join(fasta_dir, "fasta_check_report.csv")

    with open(report_file, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Archivo", "Estado", "Detalle", "Longitud", "Proporción_N"])

        for fasta_file in os.listdir(fasta_dir):
            if fasta_file.endswith(".fasta"):
                fasta_path = os.path.join(fasta_dir, fasta_file)
                status, detail, length, n_ratio = check_fasta(fasta_path, max_n_ratio=max_n_ratio)

                if not status:
                    shutil.move(fasta_path, os.path.join(fail_dir, fasta_file))
                    writer.writerow([fasta_file, "FAIL", detail, length, f"{n_ratio:.4f}" if n_ratio is not None else "NA"])
                else:
                    writer.writerow([fasta_file, "PASS", detail, length, f"{n_ratio:.4f}"])

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Verifica archivos FASTA en un directorio.")
    parser.add_argument("--input_dir", required=True, help="Ruta al directorio que contiene los archivos FASTA")
    parser.add_argument("--max_n_ratio", type=float, default=0.2,
                        help="Proporción máxima permitida de Ns (ejemplo: 0.2 = 20%)")

    args = parser.parse_args()
    main(args.input_dir, args.max_n_ratio)
