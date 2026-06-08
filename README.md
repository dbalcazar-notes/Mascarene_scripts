## Step 1: VCF missing genotype recoding

This script converts a gzipped VCF file into a tab-delimited table summarizing missing genotype information for each sample.

For each variant site, genotype calls are recoded as:

* `1`: genotype is present
* `0`: genotype is missing, identified as `./.`

The output table includes the chromosome, position, reference allele, and one column per individual.

### Usage

```bash
python step_01_missing_data_convert.py input.vcf.gz output_missing.tsv
```

### Arguments

| Argument      | Description                                                                |
| ------------- | -------------------------------------------------------------------------- |
| `input_file`  | Input gzipped VCF file (`.vcf.gz`).                                        |
| `output_file` | Output tab-delimited file containing recoded missing genotype information. |

### Output columns

| Column             | Description                                                                                  |
| ------------------ | -------------------------------------------------------------------------------------------- |
| `Chromosome`       | Chromosome or scaffold identifier.                                                           |
| `Position`         | Variant position.                                                                            |
| `Reference Allele` | Reference allele at the variant site.                                                        |
| Sample columns     | Recoded genotype presence/absence values for each individual (`1` = present, `0` = missing). |



## Window-based Missing Data Summary

This script summarizes missing genotype calls across genomic windows. It reads a tab-delimited file containing chromosome, position, and sample-level genotype information, then groups variants into fixed-size genomic bins. For each window, the script counts the number of sites coded as `0` for each individual and reports the total number of sites contained in the bin.

Empty windows (i.e., windows without any variants) are omitted from the output.

### Usage

```bash
python step_02_missing_data_windows.py input_file.tsv output_file.csv --window_size 100000 --step_size 50000
```

### Arguments

| Argument        | Description                                                                         |
| --------------- | ----------------------------------------------------------------------------------- |
| `input_file`    | Input tab-delimited file containing chromosome, position, and sample genotype data. |
| `output_file`   | Output CSV file containing window summaries.                                        |
| `--window_size` | Size of the genomic window in base pairs (default: `100`).                          |
| `--step_size`   | Step size between consecutive windows in base pairs (default: `100`).               |

### Input format

The input file must contain at least the following columns:

```text
Chromosome    Position    ...
```

Columns from the fourth column onward are assumed to contain sample-specific genotype information.

### Output

The output CSV contains:

* `Chromosome`: chromosome/scaffold identifier.
* `Start`: window start position.
* `End`: window end position.
* `Sites_in_bin`: number of variant sites within the window.
* One column per sample containing the count of sites coded as `0` within the window.

### Notes

* Windows are generated starting at position 1 for each chromosome.
* Variants are assigned to windows according to their genomic position.
* Empty windows are skipped and reported to the console.
