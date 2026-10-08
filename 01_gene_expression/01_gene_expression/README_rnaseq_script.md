# evan_fitzhugh_rnaseq.py

This document describes the script `evan_fitzhugh_rnaseq.py`.
The text follows ASD-STE100 Simplified Technical English rules.
The file `README.md` is the assignment specification. This file does not replace it.

## Contents

1. Purpose
2. Requirements and environment
3. Command line
4. Input files
5. Output files
6. How the script works (Part 1: read counting)
7. How the script works (Part 2: statistics)
8. Table of constants
9. Design decisions and assumptions
10. Known limitations
11. Worked example
12. Validation
13. Troubleshooting
14. Function reference

## 1. Purpose

The script counts RNA-seq read pairs for each E. coli gene in a set of BAM files.
Then it compares gene expression between control samples and treatment samples.

The script does two tasks:

- Part 1 makes a count matrix. Each row is a gene. Each column is a sample.
- Part 2 removes genes with low counts. It converts the counts to TPM. It then calculates statistics for each gene.

The script writes two tab-separated files. The names are `gene_counts.tsv` and `differential_expression.tsv`.

## 2. Requirements and environment

The script runs in the conda environment that the file `rna_seq.yaml` defines.

Create and activate the environment:

```bash
conda env create -f rna_seq.yaml
conda activate rna_seq
```

If the environment already exists, use only the `conda activate` command.
Use the environment name that is written in the first lines of `rna_seq.yaml`. The name `rna_seq` is the expected name.

The script imports these libraries:

| Library | Use in the script |
|---|---|
| `argparse`, `os`, `sys` (standard library) | Command line, file paths, messages on standard error |
| `numpy` | Arrays, sorted search, running maximum, median, log2 |
| `pandas` | Tables for genes, counts, TPM, and results |
| `pysam` | Read the BAM files |
| `scipy.stats` | Mann-Whitney U test |

The BAM files must be sorted by coordinate. An index file (`.bai`) is not necessary.
The script reads each BAM file from start to end (`until_eof=True`).

## 3. Command line

Syntax:

```bash
python evan_fitzhugh_rnaseq.py <sample_info_tsv> <annotation_gff> <output_directory>
```

| Argument | Meaning |
|---|---|
| `sample_info_tsv` | Path to the sample table (`samples.tsv`) |
| `annotation_gff` | Path to the GFF3 gene annotation file |
| `output_directory` | Directory for the output files. The script makes it if it does not exist |

Example. Run the command from the assignment directory that contains `data/`:

```bash
python evan_fitzhugh_rnaseq.py data/samples.tsv data/ecoli_genes.gff results/
```

The assignment text shows `samples.tsv` in the current directory. The BAM paths in the file are relative to the directory where you run the script.
So you must run the script from the directory that makes those paths correct.

Messages on standard error:

- `counting <sample name> ...` appears before each BAM file is read.
- `<N> / <M> genes passed the count filter` appears at the end. N is the number of genes kept. M is the number of genes in the GFF file.

The script prints nothing on standard output.

## 4. Input files

### 4.1 samples.tsv

The file is tab-separated. It has a header line.
The script uses two columns by name:

| Column | Meaning |
|---|---|
| `bam_path` | Path to one BAM file, relative to the current directory |
| `sample_type` | `control` or `treatment` |

Example:

```
bam_path	sample_type
data/control_01.sorted.bam	control
data/control_02.sorted.bam	control
data/treatment_01.sorted.bam	treatment
data/treatment_02.sorted.bam	treatment
```

Rules that the script applies:

- The sample name is the file name of the BAM path. The script removes the folder part first. Then it removes `.sorted.bam`. If that text is not there, it removes `.bam`. The path `data/control_01.sorted.bam` gives the sample name `control_01`.
- The script removes spaces around `sample_type` and changes it to lowercase. Therefore `Control` and ` control ` are valid.
- The order of the rows sets the order of the columns in the output.
- Each sample name must be unique. Two samples with the same name overwrite each other in the count table.
- Other column names are ignored. Other values in `sample_type` are counted in Part 1. Part 2 ignores them.
- Part 2 needs at least one `control` sample and at least one `treatment` sample.

### 4.2 GFF3 gene annotation

The file is a tab-separated GFF3 file with nine columns on each feature line.

Rules that the script applies:

- A line that starts with `#` is a comment. The script skips it.
- A line with fewer than nine tab-separated fields is skipped.
- Only lines with `gene` in column 3 (type) are used. All other types (CDS, exon, and so on) are skipped.
- Column 9 (attributes) is split at `;`. Each part is split at the first `=`. The script reads the `ID` attribute as the gene ID. A gene line with no `ID` attribute causes a `KeyError` and the script stops.
- Column 1 (sequence name) is read and stored. The script does not use it later (see section 9).
- Column 4 (start) and column 5 (end) are 1-based and inclusive.

Example gene line:

```
NC_000913.3	RefSeq	gene	190	255	.	+	.	ID=gene-b0001;Name=thrL
```

The script converts this line to the 0-based, half-open interval `[189, 255)`. The length is `255 - 189 = 66` bp.

Gene IDs should be unique. The script keeps the order of the genes in the file.

### 4.3 BAM files

The script expects paired-end reads. Each BAM file must have these properties:

- The file is sorted by coordinate. (Name order is not necessary.)
- Each pair has two records with the same query name.
- The aligner set the flags. The script uses flags 0x2 (proper pair), 0x4 (unmapped), 0x100 (secondary), 0x200 (QC fail), 0x400 (duplicate), and 0x800 (supplementary).
- The MAPQ field is set.

The BAM index is not necessary. Single-end data does not work. The script finds no pairs in single-end data and returns zero counts.

## 5. Output files

The script writes two files in the output directory. An existing file with the same name is overwritten.

### 5.1 gene_counts.tsv

Tab-separated. Has a header line. Has one row for each gene in the GFF file, in GFF order. Genes with zero counts are included.

| Column | Definition |
|---|---|
| `gene_id` | The `ID` attribute from the GFF file |
| one column for each sample | Integer number of read pairs (fragments) that the script assigned to the gene in that sample. The column name is the sample name |

The counts are integers. The unit is read pairs, not single reads.

Example:

```
gene_id	control_01	control_02	treatment_01	treatment_02
gene-b0001	245	198	267	189
gene-b0002	1089	1156	1034	1201
gene-b0003	0	0	0	0
```

### 5.2 differential_expression.tsv

Tab-separated. Has a header line. Has one row for each gene that passes the count filter (section 7.1). The rows keep the GFF order. The file has no index column.
The script writes the numbers with the format `%.10g` (up to 10 significant digits).

| Column | Definition |
|---|---|
| `gene_id` | The gene ID |
| `mean_control` | Mean TPM of the control samples |
| `median_control` | Median TPM of the control samples |
| `mean_treatment` | Mean TPM of the treatment samples |
| `median_treatment` | Median TPM of the treatment samples |
| `log2_fold_change` | `log2((mean_treatment + 0.01) / (mean_control + 0.01))` |
| `p_value` | Two-sided Mann-Whitney U p-value, treatment against control. The value is 1.0 if the test cannot run |

Example (the numbers are for illustration only):

```
gene_id	mean_control	median_control	mean_treatment	median_treatment	log2_fold_change	p_value
gene-b0001	120.5	118.2	310.7	305.1	1.36	0.0285714286
gene-b0002	540.1	538	260.3	259.8	-1.05	0.0285714286
```

A positive `log2_fold_change` means higher expression in the treatment samples.

## 6. How the script works (Part 1: read counting)

The `main` function calls the steps below in this order.

### Step 1. Parse the GFF file (`parse_gff`)

1. Open the file and read it line by line.
2. Skip comment lines and lines with fewer than nine fields.
3. Skip lines where column 3 is not `gene`.
4. Read the attributes into a dictionary. Get the `ID`.
5. Convert the coordinates: `start = GFF_start - 1` and `end = GFF_end`. The result is a 0-based, half-open interval. This is the same coordinate system as pysam.
6. Calculate `length = end - start`. This is the gene length in bp.
7. Return a table with the columns `gene_id`, `chrom`, `start`, `end`, and `length`.

### Step 2. Read the sample table (`read_samples`)

The function returns three lists in file order: sample names, BAM paths, and sample types. Section 4.1 gives the rules.

### Step 3. Build the gene index (`GeneIndex`)

The index finds the genes that overlap an interval without a scan of all genes.

1. Sort the genes by start position. The sort is stable, so genes with the same start keep the GFF order.
2. Store the sorted starts (`starts`) and the sorted ends (`ends`).
3. Calculate `max_end`. It is the running maximum of `ends`. It does not decrease.
4. Keep `order`. It maps a sorted position back to the original gene row.

The method `overlaps(start, end)` works as follows:

1. Find `hi`. It is the first sorted position whose gene starts at or after `end`. Use a binary search on `starts`. Genes at positions before `hi` start before the interval ends.
2. Find `lo`. It is the first sorted position where `max_end` is greater than `start`. Use a binary search on `max_end` with `side="right"`. No gene before `lo` can reach the interval.
3. For each position from `lo` to `hi - 1`, calculate `ov = min(end, gene_end) - max(start, gene_start)`.
4. If `ov > 0`, return the pair (original gene index, `ov`).

The overlap is the number of bases that the interval and the gene share. An interval that only touches a gene edge has an overlap of zero and does not count.

### Step 4. Count one BAM file (`count_bam`)

The function makes an integer array `counts` with one value for each gene. The values start at zero. It also makes an empty dictionary `pending`.

For each record in the BAM file, the function does the steps below.

**4a. Read filters.** The function skips the record (it uses `continue`) if any of these conditions is true:

| Filter | Test | Constant or flag |
|---|---|---|
| Unmapped | `read.is_unmapped` | flag 0x4 |
| Secondary alignment | `read.is_secondary` | flag 0x100 |
| Supplementary alignment | `read.is_supplementary` | flag 0x800 |
| QC fail | `read.is_qcfail` | flag 0x200 |
| Duplicate | `read.is_duplicate` | flag 0x400 |
| Not a proper pair | `not read.is_proper_pair` | flag 0x2 is not set |
| Low mapping quality | `read.mapping_quality < MIN_MAPQ` | `MIN_MAPQ = 10` |

A record with MAPQ of exactly 10 passes. The filter uses "less than 10".

**4b. Mate pairing.** The function uses the query name to find the mate.

1. Remove the query name from `pending`. If `pending` has a record with that name, this record is the second mate.
2. If `pending` has no such record, store this record in `pending` and go to the next record. This record is the first mate.
3. If a mate was found, the pair is complete. Continue to step 4c.

The mate is removed from `pending` when the pair is complete. A name that appears again later starts a new pair.
A record whose mate was filtered out stays in `pending` until the end of the file. The function never counts it.

**4c. Template span.** The template is the region from the leftmost aligned base to the rightmost aligned base of both mates.

- `t_start = min(mate.reference_start, read.reference_start)`
- `t_end = max(mate.reference_end, read.reference_end)`
- `template length = t_end - t_start`

`reference_start` is 0-based. `reference_end` is the position after the last aligned base. The template includes the gap between the two mates (the unsequenced part of the fragment).
Soft-clipped bases are not part of the aligned span.

**4d. Overlap and assignment.**

1. Call `index.overlaps(t_start, t_end)`. Sort the results by overlap in descending order. The result is the list `hits`.
2. If `hits` is empty, no gene overlaps the template. Skip the pair.
3. If `hits` has two or more genes and the two largest overlaps are equal, the pair is ambiguous. Skip the pair. The pair is not counted for any gene.
4. Apply the 50% gate. If `hits[0]` overlap is less than `MIN_OVERLAP_FRAC * template length`, skip the pair. `MIN_OVERLAP_FRAC = 0.5`. An overlap of exactly 50% passes.
5. Otherwise, add 1 to the count of the best gene. The pair goes wholly to this gene. The script does not split a pair between genes.

The tie check (step 3) runs before the 50% gate (step 4). The script compares only the two largest overlaps. This is enough to find a tie at the top.

### Step 5. Build the count matrix

For each sample, `main` calls `count_bam` and stores the result as a column. After all samples, it writes `gene_counts.tsv`. The table keeps all genes, including genes with zero counts.

## 7. How the script works (Part 2: statistics)

### 7.1 Count filter

`keep = counts.sum(axis=1) >= MIN_TOTAL_COUNT`. `MIN_TOTAL_COUNT = 10`.

The script adds the counts of one gene across all samples. It keeps the gene if the total is 10 or more. The script removes genes with a total of 9 or less.
The filter uses all samples, both control and treatment.

### 7.2 TPM (`compute_tpm`)

The script calculates TPM separately for each sample (each column), only for the kept genes:

1. `RPK = count / (gene length in bp / 1000)`
2. `total_RPK = sum of RPK over all kept genes in the sample`
3. `scaling_factor = total_RPK / 1,000,000`
4. `TPM = RPK / scaling_factor`

If `total_RPK` is zero for a sample, the scaling factor is replaced with NaN. The function then changes the NaN results to 0.0. Such a sample has a TPM of 0 for all genes.
Because the script filters before it calculates TPM, the TPM values of each sample sum to 1,000,000 over the kept genes only.

### 7.3 Statistics (`differential_expression`)

1. Select the TPM columns where the sample type is `control`. Select the TPM columns where the type is `treatment`.
2. For each gene, calculate the mean and the median of each group (`numpy`).
3. Calculate `log2_fold_change = log2((mean_treatment + PSEUDOCOUNT) / (mean_control + PSEUDOCOUNT))`. `PSEUDOCOUNT = 0.01`. The pseudocount is added to the group means, not to the counts. It prevents division by zero and `log2(0)`.
4. Calculate `p_value` with `scipy.stats.mannwhitneyu(treatment_values, control_values, alternative="two-sided")`. The function uses its default method. SciPy chooses the exact or the asymptotic method.
5. If SciPy raises a `ValueError` (for example, all values are identical), set the p-value to 1.0. If the p-value is NaN, set it to 1.0.
6. Write the table to `differential_expression.tsv`.

The script does not correct the p-values for multiple tests. The output has no adjusted p-value column.

## 8. Table of constants

All constants are at the top of the script.

| Constant | Value | Use |
|---|---|---|
| `MIN_MAPQ` | `10` | Minimum mapping quality. Records with a lower MAPQ are skipped |
| `MIN_TOTAL_COUNT` | `10` | Minimum total count over all samples for a gene to enter Part 2 |
| `MIN_OVERLAP_FRAC` | `0.5` | Minimum fraction of the template length that must lie in the assigned gene |
| `PSEUDOCOUNT` | `0.01` | Added to each group mean TPM before the log2 ratio |

Other fixed values in the code: `1000.0` (bp to kb), `1e6` (per million), and `float_format="%.10g"` (output precision).

## 9. Design decisions and assumptions

- **Proper-pair filter is kept.** The assignment says to filter reads by proper pairing flags. The script obeys this. Pairs without the 0x2 flag are dropped. About 5% of the reads have no 0x2 flag, so the script does not count them.
- **Duplicates, QC-fail, secondary, and supplementary records are skipped.** The script does not count these records, to prevent double counting of one fragment.
- **A pair is one fragment.** The script counts each pair once. It does not count each mate on its own.
- **Overlap is measured on the template.** This follows the TA hints. The assignment text says "50% of read length". The script uses the template length of the pair. The template includes the gap between the mates.
- **Winner takes the whole pair.** The pair goes wholly to the gene with the greater overlap. The assignment text mentions "proportional assignment" for overlapping genes. The script does not split counts. This keeps the counts as integers.
- **Equal overlap means ambiguous.** The pair is not counted.
- **Counts are integers.** The array type is `int64`.
- **The chromosome name is not checked.** The script compares only positions. It does not compare the sequence name in the GFF with the reference name in the BAM file. The script assumes that the genome has one sequence and that the GFF and BAM use the same coordinates. This is true for E. coli K-12 in this data set.
- **Reads are not fetched by region.** The script reads each file in full (`until_eof=True`). This is simple and needs no index.
- **Genes only.** Only `gene` features are used. Gene length is the full gene interval, not the sum of exons.
- **TPM after filtering.** The script keeps genes with a total of 10 or more counts, then calculates TPM on those genes.
- **Mann-Whitney direction.** The test is two-sided, so the order of the groups does not change the p-value.

## 10. Known limitations

- The script does not check that the sample types include both groups. A missing group causes an error or empty statistics.
- The script keeps unpaired records in memory (`pending`) until the end of the file. A very large BAM file with many unpaired records uses more memory.
- The pairing uses only the query name. It does not check that the two records are read 1 and read 2.
- Pairs in which one mate is filtered out (for example, a low MAPQ) are not counted.
- The 50% gate uses the template length. A pair with a long gap between the mates can fail the gate even when both mates lie in the gene.
- The script does not correct for multiple tests. With few samples per group, the smallest possible Mann-Whitney p-value is large. For example, 4 against 4 samples gives a smallest two-sided p-value of about 0.029.
- The pseudocount of 0.01 is small. A gene with a mean of zero in one group can have a very large fold change.
- Duplicate gene IDs in the GFF file can cause errors in Part 2.
- The script does not handle strandedness. Counts do not depend on the strand.

## 11. Worked example

Use these genes. The table shows the GFF coordinates and the converted coordinates.

| Gene | GFF start-end (1-based) | Script interval (0-based, half-open) |
|---|---|---|
| A | 1001-1500 | [1000, 1500) |
| B | 1401-1900 | [1400, 1900) |
| C | 2001-2100 | [2000, 2100) |
| D | 2101-2200 | [2100, 2200) |

### Case 1. Clear winner

Mate 1 aligns at `reference_start = 1100` and `reference_end = 1200`.
Mate 2 aligns at `reference_start = 1350` and `reference_end = 1450`.
Both records pass all filters. They have the same query name.

1. Mate 1 enters `pending`. Mate 2 finds it, so the pair is complete.
2. `t_start = min(1100, 1350) = 1100`. `t_end = max(1200, 1450) = 1450`. Template length = 350.
3. Overlap with A = `min(1450, 1500) - max(1100, 1000) = 350`.
4. Overlap with B = `min(1450, 1900) - max(1100, 1400) = 50`.
5. `hits` sorted: A (350), B (50). The two largest values differ, so there is no tie.
6. Gate: `350 >= 0.5 * 350 = 175`. The pair passes.
7. The script adds 1 to the count of gene A. Gene B gets nothing.

### Case 2. Tie

Mate 1 aligns at 2050-2080. Mate 2 aligns at 2120-2150.

1. Template: `t_start = 2050`, `t_end = 2150`. Length = 100.
2. Overlap with C = `min(2150, 2100) - max(2050, 2000) = 50`.
3. Overlap with D = `min(2150, 2200) - max(2050, 2100) = 50`.
4. `hits` has two genes with overlap 50. The two largest overlaps are equal.
5. The pair is ambiguous. The script does not count it for C or D.

### Case 3. Gate failure

Mate 1 aligns at 1850-1900. Mate 2 aligns at 1950-1990.

1. Template: 1850 to 1990. Length = 140.
2. Only gene B overlaps the template. Overlap with B = `min(1990, 1900) - max(1850, 1400) = 50`.
3. There is one hit, so there is no tie.
4. Gate: `50 < 0.5 * 140 = 70`. The pair fails the gate.
5. The script does not count the pair.

## 12. Validation

A comparison with featureCounts used the same BAM files and the same GFF file.

- With the proper-pair filter on, 99.48% of the cells in the count matrix matched featureCounts.
- With the proper-pair filter relaxed, 99.975% of the cells matched.

The difference comes from pairs without the 0x2 flag. These pairs are about 5% of the reads. The script drops them, as the assignment specification says.
The submitted script keeps the proper-pair filter. The relaxed run was a test only.

## 13. Troubleshooting

| Problem | Cause | Action |
|---|---|---|
| `ModuleNotFoundError` for `pysam`, `pandas`, `numpy`, or `scipy` | The conda environment is not active | Run `conda activate rna_seq`. Run the script again |
| `FileNotFoundError` for a BAM file | The BAM path in `samples.tsv` is relative to a different directory | Run the script from the assignment directory. Make sure the paths are correct |
| `KeyError: 'bam_path'` or `KeyError: 'sample_type'` | The header of `samples.tsv` is wrong, or the file uses spaces instead of tabs | Use the tab character. Use the exact header `bam_path` and `sample_type` |
| `KeyError: 'ID'` | A `gene` line in the GFF file has no `ID` attribute | Add the `ID` attribute to the line, or remove the line |
| All counts are zero | The BAM file has no proper pairs, or the GFF and BAM coordinates do not match | Check the flags with `samtools flagstat`. Make sure the GFF file is for the same genome |
| Counts are lower than expected | The proper-pair filter, the MAPQ filter, or the 50% gate removed pairs | Read section 6. Count how many pairs each filter removes |
| Zero rows in `differential_expression.tsv` | No gene has a total count of 10 or more | Check `gene_counts.tsv` |
| Error in Part 2 about empty arrays | The sample table has no `control` or no `treatment` rows | Correct the `sample_type` values |
| `p_value` is 1.0 for many genes | The values are identical, or the sample size is small | This is expected behavior (section 7.3) |
| `sample_type` values are not used | The values have a different spelling | Use `control` and `treatment` (case does not matter) |
| Warning or error about a missing BAM index | The code reads with `until_eof=True`, so the index is not necessary | Look for a different cause, such as a damaged file. Run `samtools quickcheck` |

## 14. Function reference

| Function or class | Input | Output |
|---|---|---|
| `parse_gff(path)` | GFF file path | Table of genes: `gene_id`, `chrom`, `start`, `end`, `length` |
| `read_samples(path)` | Sample table path | Sample names, BAM paths, lowercase sample types |
| `GeneIndex(genes)` | Gene table | Interval lookup object |
| `GeneIndex.overlaps(start, end)` | Interval (0-based, half-open) | Generator of (original gene index, overlap in bp) |
| `count_bam(bam_path, index, n_genes)` | BAM path, index, number of genes | Integer array of pair counts for each gene |
| `compute_tpm(counts, lengths)` | Count table (genes by samples), gene lengths in bp | TPM table of the same shape |
| `differential_expression(tpm, types)` | TPM table, sample types | Result table with the seven output columns |
| `main()` | Command line arguments | Writes the two output files |
