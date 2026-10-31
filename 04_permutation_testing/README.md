# Assignment: Genomic Region Overlap Analysis and Permutation Testing

## Description

You are a computational biologist investigating the relationship between transcription factor binding sites and regulatory chromatin regions. You will implement a permutation test to determine whether the overlap between two sets of genomic regions is statistically significant, a fundamental technique in genomics and regulatory biology.

This analysis examines potential co-localization between:
- **Set A**: Transcription factor binding sites identified by ChIP-seq
- **Set B**: Active chromatin regions identified by histone modification mapping

The question you are answering is: **do Set A and Set B overlap more than we would expect if Set B's regions were scattered at random along their chromosomes?**

## Getting Started

Your script must run in the provided conda environment, `permutation_test.yaml`:

```bash
conda env create -f permutation_test.yaml
conda activate permutation_test
```

You may use only the following:
- **numpy**, for numerical operations and random number generation
- **pandas**, for data manipulation and file output
- Python standard library modules (`random`, `sys`, `argparse`, `bisect`, `collections`, etc.)

**Do not** install or use bioinformatics or interval libraries (e.g. `pysam`, `pybedtools`, `pyranges`, `intervaltree`). The point of the assignment is to implement the overlap calculation yourself.

## Data Files Provided

All input files are in the `data/` folder:

| File | Contents |
|------|----------|
| `data/SetA.bed` | Transcription factor binding sites (3-column BED) |
| `data/SetB.bed` | Active chromatin regions (3-column BED) |
| `data/genome.fa.fai` | Chromosome names and lengths (FASTA index) |
| `data/toy/` | The small worked example from [Key Definitions](#key-definitions), for testing |

Things to know about the provided data:
- The genome has 5 chromosomes (`chr1`–`chr5`). Every chromosome that appears in a BED file is also listed in `genome.fa.fai`.
- Each BED file has 5,000 regions. There is no header line.
- **The BED files are not sorted.** Sort the regions yourself if your algorithm needs sorted input.
- **Regions within the same file can overlap each other.** You must handle this (see [Key Definitions](#key-definitions)).
- All regions lie inside their chromosome (`0 <= start < end <= chromosome length`).

Your script is graded on this data and on the toy example in `data/toy/`, which has only 2 chromosomes, different chromosome lengths, and far fewer regions. Do not hard-code chromosome names, region counts, or chromosome lengths.

### BED File Format
- **3-column format:** chromosome, start, end (tab-separated)
- **0-based, half-open coordinates:** start is inclusive, end is exclusive
- **Length** of a region is `end - start`
- **Example:** `chr1	1000	2000` covers bases 1000–1999 on chr1 (1,000 bases)
- Two regions `[s1, e1)` and `[s2, e2)` on the same chromosome share `max(0, min(e1, e2) - max(s1, s2))` bases

### FASTA Index Format
- **5-column format:** chromosome, length, offset, bases_per_line, bytes_per_line
- **You only need** the chromosome name (column 1) and length (column 2)
- **Example:** `chr1	1876424	0	50	50` means chr1 is 1,876,424 bases long, so valid positions are 0–1,876,423

## Key Definitions

Read this section carefully. Your numbers will only match the expected answers if you use these definitions.

**Merging (union) of a region set.** Within one set, on each chromosome, combine regions that overlap into a single region. For example, `chr1 10 30` and `chr1 20 40` merge into `chr1 10 40`. You may also merge regions that touch end-to-start (`chr1 10 20` and `chr1 20 30` → `chr1 10 30`). Either choice gives the same base counts. After merging, no base is covered twice within a set.

**Overlap.** The overlap between Set A and Set B is the **number of bases covered by at least one region in Set A and at least one region in Set B**. Each base counts at most once, no matter how many regions cover it. The simplest correct way to compute this is to merge each set first, then add up the lengths of the intersections between the merged A regions and the merged B regions on each chromosome.

**Region count vs. total bases.**
- `setX_regions` is the number of regions in the BED file as given (number of lines, **before** merging).
- `setX_total_bases` is the number of bases covered by the set **after** merging.

### Worked toy example

`genome.fa.fai`:
```
chr1	100	0	50	50
chr2	50	0	50	50
```

`SetA.bed`:
```
chr1	10	30
chr1	20	40
chr1	60	70
chr2	0	10
```

`SetB.bed`:
```
chr1	25	65
chr2	5	20
chr1	90	100
```

1. Merge Set A: `chr1 10 40`, `chr1 60 70`, `chr2 0 10`, which covers 30 + 10 + 10 = **50 bases**.
2. Merge Set B: nothing overlaps, so it stays `chr1 25 65`, `chr1 90 100`, `chr2 5 20`, which covers 40 + 10 + 15 = **65 bases**.
3. Intersections: `chr1 25 40` (15) + `chr1 60 65` (5) + `chr2 5 10` (5) = **25 bases of overlap**.

The expected values for the deterministic metrics are:

```
observed_overlap	25
setA_regions	4
setB_regions	3
setA_total_bases	50
setB_total_bases	65
```

> **Common mistake:** If you skip merging and add up the overlap of every A–B pair, you get 30 instead of 25, because bases 25–29 on chr1 are covered by two Set A regions and get counted twice.

These files are provided in `data/toy/`. Test your script on them before running it on the real data (the arguments are explained under [Command Line Interface](#command-line-interface)):

```bash
python firstname_lastname_permutation_test.py data/toy/SetA.bed data/toy/SetB.bed data/toy/genome.fa.fai toy_results.tsv
```

## Assignment Tasks

### Part 1: Overlap Calculation (40 points)

1. **Parse** the two BED files and the genome index. Use only the first two columns of the `.fai` file (chromosome name and length).
2. **Merge** the overlapping regions within each set, separately for each chromosome.
3. **Compute the observed overlap** in base pairs between merged Set A and merged Set B, summed over all chromosomes.
4. **Compute** the region counts and total bases for each set, as defined above.

**Requirements:**
- Count overlapping **base pairs**, not the number of overlapping regions.
- Only regions on the same chromosome can overlap.
- Use an efficient algorithm. Comparing every A region with every B region (a nested loop) is far too slow once it runs inside 10,000 permutations. Instead, sort the merged regions on each chromosome and walk through both lists together in one pass (the "two-pointer" or "sweep" approach).

### Part 2: Permutation Testing (60 points)

**Null hypothesis:** Set B's regions are placed along their chromosomes independently of Set A. Any overlap with Set A happens by chance.

Build a null distribution by repeatedly shuffling Set B and recomputing the overlap:

1. **Keep Set A fixed.** Merge it once and reuse it in every permutation.
2. **Shuffle Set B.** For each **original, unmerged** Set B region with length `L = end - start` on chromosome `c`:
   - Keep it on chromosome `c`.
   - Choose a new start uniformly at random from the integers `0, 1, ..., chrom_len - L` (both ends included), where `chrom_len` is the length of chromosome `c`.
   - Set the new end to `new_start + L`.

   This keeps the number of regions, the size of each region, and the number of regions per chromosome the same as in the real data. The shuffled region always fits inside its chromosome.
3. **Allow shuffled regions to overlap each other.** Do not reject or redraw placements. Instead, **merge the shuffled Set B** before computing its overlap with Set A, exactly as you did for the observed data.
4. **Record the overlap** between merged Set A and the merged, shuffled Set B.
5. **Repeat** steps 2–4 `N` times. `N` is the number of permutations, 10,000 by default.

**p-value.** We are testing for **enrichment** (more overlap than chance), so use a one-sided, upper-tail p-value with a +1 correction:

```
p_value = (1 + number of permutations with overlap >= observed_overlap) / (1 + N)
```

The +1 counts the observed arrangement as one possible arrangement under the null. This means the p-value can never be exactly 0. With `N` permutations the smallest possible p-value is `1 / (N + 1)`. **Do not report `p_value` as 0.**

**Watch out for off-by-one errors in your random number generator:**
- Python's `random.randint(a, b)` **includes** `b`: use `random.randint(0, chrom_len - L)`.
- NumPy's `rng.integers(low, high)` **excludes** `high`: use `rng.integers(0, chrom_len - L + 1)`.
- If a region is exactly as long as its chromosome, the only valid start is 0. An upper bound that is one too small gives an empty range and crashes. One that is one too large lets regions run past the end of the chromosome.

## Technical Requirements

### Command Line Interface

Your script must accept these positional command-line arguments, in this order:

```bash
python firstname_lastname_permutation_test.py <setA_bed> <setB_bed> <genome_fai> <output_file> [num_permutations]
```

- `num_permutations` is optional. If it is not given, use **10000**.
- Do not require any other arguments, and do not prompt for input.
- You may set a fixed random seed inside your script for reproducibility, but it is not required.

Example:
```bash
python firstname_lastname_permutation_test.py data/SetA.bed data/SetB.bed data/genome.fa.fai results.tsv 10000
```

### Output File

Your script must write a **tab-separated** file to `<output_file>` with:
- a header line `metric	value`, then
- exactly these 7 rows, in this order:

| metric | meaning |
|--------|---------|
| `observed_overlap` | Bases covered by both merged Set A and merged Set B (integer) |
| `p_value` | One-sided permutation p-value, `(1 + #{null >= observed}) / (1 + N)` |
| `num_permutations` | `N`, the number of permutations performed (integer) |
| `setA_regions` | Number of regions in `SetA.bed` before merging (integer) |
| `setB_regions` | Number of regions in `SetB.bed` before merging (integer) |
| `setA_total_bases` | Bases covered by Set A after merging (integer) |
| `setB_total_bases` | Bases covered by Set B after merging (integer) |

Integers must be written without decimals (`25`, not `25.0`). Any standard decimal or scientific notation is fine for `p_value` (e.g. `0.0001` or `9.999e-05`).

Example of the format, using the toy example above with `N = 1000`. Your `p_value` will differ because the permutations are random.

```
metric	value
observed_overlap	25
p_value	0.4116
num_permutations	1000
setA_regions	4
setB_regions	3
setA_total_bases	50
setB_total_bases	65
```

Only the output file is graded. Results printed to the screen are ignored, although printing progress messages is fine.

### Performance Requirements

- **Runtime limit:** With the provided data and 10,000 permutations, your script must finish within **5 minutes**.
- **Feasibility:** A straightforward pure-Python solution that sorts, merges, and sweeps each chromosome finishes in under a minute on a laptop, so NumPy is not required. If your script is much slower, check for nested loops over regions.
- **Memory:** Do not build per-base arrays (e.g. one entry per base for the whole genome) in every permutation. Work with region coordinates.

## Grading

Your script is run on the provided data (with no `num_permutations` argument, so the default is used) and on the toy example in `data/toy/` (with 10,000 permutations). There are no hidden datasets.

| Part | What is checked | Points |
|------|-----------------|--------|
| 1 | On both datasets: `observed_overlap`, `setA_total_bases`, `setB_total_bases`, the two region counts, and the output format must match **exactly** (20 points per dataset) | 40 |
| 2 | Provided data: `p_value` equals `1 / (N + 1)` (the observed overlap is larger than every permutation) | 20 |
| 2 | Provided data: `num_permutations` defaults to 10,000 | 10 |
| 2 | Provided data: finishes within 5 minutes | 10 |
| 2 | Toy example: `p_value` is within random-sampling error of the correct one-sided p-value (about 0.40) | 20 |

Because the permutations are random, your toy p-value will not exactly match ours. A correct implementation with 10,000 permutations passes comfortably. A two-sided p-value will not, and neither will shuffling Set A instead of Set B. Leaving out the +1 gives `p_value 0` on the provided data, which loses points.

Your results file is required but not scored (see [Submission Requirements](#submission-requirements)).

## Submission Requirements

Submit these two files:

| File | Contents |
|------|----------|
| `firstname_lastname_permutation_test.py` | Your script |
| `firstname_lastname_results.tsv` | Your script's output on the provided data, using the default 10,000 permutations |

Create the results file by running your final script with exactly this command, and do not edit the file by hand:

```bash
python firstname_lastname_permutation_test.py data/SetA.bed data/SetB.bed data/genome.fa.fai firstname_lastname_results.tsv
```

The grader re-runs your script and compares its output with your results file. The results file is **required for your submission to count, but it is not scored**. Your grade comes from running your script (see [Grading](#grading)).

Before submitting, check that:
- [ ] Your script runs in the `permutation_test` conda environment using only the allowed libraries
- [ ] It accepts the arguments above and defaults to 10,000 permutations
- [ ] It reproduces the toy example values exactly (`25`, `4`, `3`, `50`, `65`)
- [ ] It writes the 7 metrics, tab-separated, with the header line, in the order shown
- [ ] It never reports a `p_value` of 0
- [ ] It finishes on the provided data with 10,000 permutations in under 5 minutes
- [ ] `firstname_lastname_results.tsv` was produced by your final script and not edited
