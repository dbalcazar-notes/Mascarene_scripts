#python step_03_sbatch_samples_bam.py --bed_dir genes_to_keep_bed/fixed_bed --bam sample.bam --out_dir genes_to_keep
import subprocess
import argparse
import os
import glob

def process_bed(bed_file, bam_file, output_dir):
    bam_base_name = os.path.basename(bam_file)
    temp_bam_file = f"{bam_base_name}.temp.bam"
    
    # Extraer regiones del BAM
    subprocess.run(["samtools", "view", "-b", "-M", "-L", bed_file, bam_file, "-o", temp_bam_file])
    subprocess.run(["samtools", "index", temp_bam_file])
    
    # Leer BED
    with open(bed_file, 'r') as bed:
        bed_lines = bed.readlines()
    
    consensus_sequences = []
    
    for line in bed_lines:
        if not line.strip():
            continue
        fields = line.strip().split()
        if len(fields) < 3:
            print(f"Skipping malformed line: {line}")
            continue
        chrom = fields[0]
        start = str(int(fields[1]) + 1)
        end = fields[2]
        region = f"{chrom}:{start}-{end}"
        
        # Generar consenso
        result = subprocess.run(
            ["samtools", "consensus", "-d", "5", "--mode", "simple", "-A", "-q",
             "--show-ins", "no", "--show-del", "no", "-c", "0.75", "-H", "0.5",
             temp_bam_file, "-r", region],
            capture_output=True, text=True
        )
        
        consensus_output = result.stdout.strip().split('\n')
        if len(consensus_output) > 1:
            sequence = ''.join(consensus_output[1:])
            consensus_sequences.append(sequence)
    
    combined_sequence = ''.join(consensus_sequences)
    header = os.path.splitext(os.path.basename(bed_file))[0]
    
    output_file = os.path.join(output_dir, f"{header}_{bam_base_name}.fasta")
    with open(output_file, 'w') as f:
        f.write(f">{header}\n")
        for i in range(0, len(combined_sequence), 60):
            f.write(combined_sequence[i:i+60] + '\n')
    
    subprocess.run(["rm", temp_bam_file, f"{temp_bam_file}.bai"])
    print(f"Processed {bed_file} → {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generar consensos para múltiples BEDs")
    parser.add_argument("--bed_dir", required=True, help="Directorio con archivos BED")
    parser.add_argument("--bam", required=True, help="Archivo BAM")
    parser.add_argument("--out_dir", required=True, help="Directorio de salida FASTA")
    args = parser.parse_args()
    
    os.makedirs(args.out_dir, exist_ok=True)
    bed_files = glob.glob(os.path.join(args.bed_dir, "*.bed"))
    
    for bed_file in bed_files:
        process_bed(bed_file, args.bam, args.out_dir)
