# Assignment: RNA-seq Differential Expression Analysis

## Description

You are a bioinformatician tasked with analyzing RNA-seq data from a bacterial gene expression experiment. You will work with alignment files (BAM format) to extract gene expression counts and perform differential expression analysis.

This experiment examines gene expression changes in *E. coli* under two conditions:

- **Control**: Multiple replicate samples under standard growth conditions
- **Treatment**: Multiple replicate samples under stress conditions

**Note**: All required files for this assignment are provided in the attached zip folder here: 

## Data Files Provided

You will be given the following files in (found in the `data` directory):

1. **BAM files**: Alignment files named `control_XX.sorted.bam` and `treatment_XX.sorted.bam`
   - Control samples: `control_01.sorted.bam`, `control_02.sorted.bam`, etc.
   - Treatment samples: `treatment_01.sorted.bam`, `treatment_02.sorted.bam`, etc.
   - Each BAM file contains paired-end reads aligned to the *E. coli* genome

2. **Gene annotations**: `ecoli_genes.gff` - Gene feature annotations in GFF3 format
3. **Sample info file**: `samples.tsv` - A TSV file that maps bam files to their sample type.


## Assignment Tasks

### Part 1: Read Counting (50 points)

Extract gene expression counts from the BAM files:

1. **Read BAM files** using `pysam` library
2. **Count reads per gene** that overlap gene coordinates from the GFF file
3. **Generate count matrix** with genes as rows and samples as columns
4. **Output**: `gene_counts.tsv` - tab-delimited file with raw read counts

**Required format for `gene_counts.tsv`:**
```
gene_id	control_01	control_02	treatment_01	treatment_02	...
gene-b0001	245	198	267	189	...
gene-b0002	1089	1156	1034	1201	...
gene-b0003	0	0	0	0	...
...
```
- First column: Gene IDs (from GFF file)
- Remaining columns: Raw read counts for each sample
- Include all genes from GFF file (even those with zero counts)

**Requirements:**
- Use only reads with mapping quality ≥ 10
- Count reads that overlap gene features by ≥ 50% of read length
- Handle paired-end reads correctly (count each pair once)

### Part 2: Normalization and Statistics (50 points)

Perform differential expression analysis:

1. **Filter low-expression genes**: Remove genes with fewer than 10 total reads across all samples
2. **Normalize counts** using TPM (Transcripts Per Million) method
3. **Calculate statistics** for each gene:
   - Mean and median normalized expression (control vs treatment)
   - Log₂ fold change: log₂(treatment/control)
   - Statistical significance using Mann-Whitney U test

4. **Output**: `differential_expression.tsv` with columns:
   ```
   gene_id	mean_control	median_control	mean_treatment	median_treatment	log2_fold_change	p_value
   ```

**Note**: The simulated dataset contains reads from only a subset of E. coli genes, so most genes will have zero or very low read counts.

## Technical Requirements

### Command Line Interface

Your script must accept command line arguments:

```bash
python firstname_lastname_rnaseq.py <sample_info_tsv> <annotation_gff> <output_directory>
```

The script will be run from the assignment directory that contains the data files. A sample info TSV file will be provided.

Example:
```bash
python firstname_lastname_rnaseq.py samples.tsv data/ecoli_genes.gff results/
```

### Sample Information File

A TSV file will be provided that specifies the BAM files and their sample types with the format:

```
bam_path	sample_type
data/control_01.sorted.bam	control
data/control_02.sorted.bam	control
data/treatment_01.sorted.bam	treatment
data/treatment_02.sorted.bam	treatment
...
```

- **bam_path**: Relative path to each BAM file
- **sample_type**: Either "control" or "treatment"

### Required Environment

Your script must run in the provided conda environment: `rna_seq.yaml`

You may use the following Python libraries (all included in the environment):

- **pysam** - for BAM file reading
- **pandas** - for data manipulation  
- **numpy** - for numerical operations
- **scipy.stats** - for statistical testing
- Standard library modules (os, sys, argparse, etc.)

### Output Files

Your script must generate:
1. `gene_counts.tsv` - Raw read counts matrix
2. `differential_expression.tsv` - DE analysis results  


## Submission Requirements

Submit only your Python script named: `firstname_lastname_rnaseq.py`

Your script must:
- Run in the provided `rna_seq.yaml` conda environment
- Accept the required command line arguments
- Generate the specified output files


## Technical Notes

### BAM File Processing
- BAM files contain paired-end reads aligned to *E. coli* genome
- Use `pysam.AlignmentFile()` to read BAM files
- Filter reads by mapping quality and proper pairing flags

### Gene Counting Strategy
- Count reads overlapping gene coordinates (start to end)
- Handle overlapping genes by proportional assignment

### Normalization - TPM (Transcripts Per Million)

TPM normalizes for both gene length and sequencing depth. **Calculate TPM separately for each sample:**

1. **For each gene in the sample**: `RPK = raw_count / gene_length_kb`
2. **Sum all RPK values in that sample**: `total_RPK = sum(all RPK values for sample)`  
3. **Calculate per-million scaling factor**: `scaling_factor = total_RPK / 1,000,000`
4. **For each gene**: `TPM = RPK / scaling_factor`

**Key points:**
- TPM is calculated per sample (each sample has its own set of TPM values)
- Within each sample, all TPM values sum to 1,000,000
- TPM allows comparison across samples and genes
