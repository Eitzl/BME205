#!/usr/bin/env python
"""RNA-seq read counting and differential expression for E. coli BAM files.

Usage:
    python evan_fitzhugh_rnaseq.py <sample_info_tsv> <annotation_gff> <output_directory>

Outputs (in <output_directory>):
    gene_counts.tsv            raw integer fragment (read-pair) counts per gene and sample
                               (all GFF genes, including zeros)
    differential_expression.tsv  TPM-based control vs treatment statistics for genes with
                               >= 10 total reads (means, medians, log2 FC, Mann-Whitney p)

Background for readers new to RNA-seq
-------------------------------------
RNA-seq measures how active each gene is.  RNA is extracted from cells, converted to
cDNA, chopped into short fragments, and sequenced from both ends ("paired-end"
sequencing), giving two reads (mates) per fragment.  Each read is then aligned
(mapped) to the reference genome and the alignments are stored in a BAM file.
The more RNA a gene produces, the more fragments align to it, so counting
fragments per gene gives a rough measure of expression.

Raw counts cannot be compared directly, though:
  * Longer genes produce more fragments than short genes at the same expression.
  * A sample sequenced more deeply has more reads of everything.
TPM (Transcripts Per Million) corrects for both, so TPM values can be compared
across genes and across samples.

Pipeline implemented here
-------------------------
  Part 1 (counting):  for each BAM file, pair up mates, decide which single gene
                      each read pair belongs to, and tally an integer count per gene.
  Part 2 (statistics): drop genes with too little data, convert counts to TPM,
                      then compare control vs treatment samples per gene with the
                      mean/median TPM, a log2 fold change and a Mann-Whitney U test.

Coordinate convention used throughout (0-based, half-open)
----------------------------------------------------------
A region [start, end) includes position `start` but NOT position `end`, and the
first base of a chromosome is position 0.  This is the convention pysam uses for
BAM alignments (reference_start / reference_end).  GFF files instead use 1-based,
inclusive coordinates (first base is 1, and `end` is included).  Everything is
converted to 0-based half-open once, at GFF-parsing time, so that BAM and GFF
coordinates are directly comparable.  A big practical benefit of half-open
intervals is that lengths are simply end - start and overlaps are simply
min(ends) - max(starts), with no "+1" corrections to get wrong.
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd
import pysam
from scipy import stats

# Minimum MAPQ (mapping quality) for a read to be trusted.  MAPQ is a Phred-scaled
# confidence that the read was placed at the right genomic position; 10 means about
# a 10% chance of being misplaced.  Reads below this are often multi-mappers
# (identical sequence in several places), so we discard them to avoid miscounting.
MIN_MAPQ = 10
# Genes whose counts summed over ALL samples are below this are dropped before the
# statistics: with so few reads, TPM estimates and rank tests are mostly noise.
MIN_TOTAL_COUNT = 10
MIN_OVERLAP_FRAC = 0.5  # min fraction of the template that must lie in the assigned gene
PSEUDOCOUNT = 0.01  # added to group mean TPM before the log2 ratio


def parse_gff(path):
    """Return a DataFrame of `gene` features: gene_id, start (0-based), end, length (bp).

    GFF3 is a tab-separated annotation format.  Each non-comment line has 9 columns:
        0 seqid (chromosome)   1 source   2 type (gene, CDS, ...)   3 start
        4 end                  5 score    6 strand                  7 phase
        8 attributes (semicolon-separated key=value pairs, e.g. "ID=gene123;Name=abc")
    We keep only rows of type "gene" and record their ID, chromosome and interval.

    Coordinates: GFF start/end are 1-based and inclusive on both ends.  We convert to
    0-based half-open [start, end) by subtracting 1 from the start and leaving the end
    unchanged (the 1-based inclusive end equals the 0-based exclusive end).  This
    matches pysam's BAM coordinates, so reads and genes can be compared directly.
    The resulting columns are: gene_id, chrom, start, end, length.
    """
    # One tuple per gene; turned into a DataFrame once at the end (much faster than
    # appending to a DataFrame row by row).
    rows = []
    with open(path) as fh:
        for line in fh:
            # Lines beginning with '#' are headers/directives (e.g. ##gff-version 3,
            # ##sequence-region), not annotation records.
            if line.startswith("#"):
                continue
            # GFF columns are separated by tabs; strip the trailing newline first so
            # it does not stick to the last (attributes) column.
            f = line.rstrip("\n").split("\t")
            # Skip malformed/short lines and any feature that is not a gene
            # (exons, CDS, rRNA, etc. would otherwise double-count the same locus).
            if len(f) < 9 or f[2] != "gene":
                continue
            # Parse "k1=v1;k2=v2" into a dict.  split("=", 1) splits only on the first
            # '=' so values that themselves contain '=' stay intact; the `if "=" in kv`
            # guard ignores empty pieces (e.g. from a trailing ';').
            attrs = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
            start, end = int(f[3]) - 1, int(f[4])  # GFF is 1-based inclusive -> 0-based half-open
            # f[0] is the chromosome/contig name; attrs["ID"] is the unique gene id.
            rows.append((attrs["ID"], f[0], start, end))
    genes = pd.DataFrame(rows, columns=["gene_id", "chrom", "start", "end"])
    # Because the interval is half-open, length is just end - start (no +1 needed).
    # The length is later used to normalise counts for gene size when computing TPM.
    genes["length"] = genes["end"] - genes["start"]
    return genes


def read_samples(path):
    """Return (sample names, bam paths, sample types) in file order.

    The sample-info file is a tab-separated table with (at least) the columns
    `bam_path` (location of each sample's alignment file) and `sample_type`
    (the experimental group, "control" or "treatment").  File order is preserved so
    that names, paths and types stay aligned index-for-index, which also fixes the
    column order of the output tables.
    """
    info = pd.read_csv(path, sep="\t")
    # Derive a short, readable sample name from each BAM filename by dropping the
    # directory and the ".sorted.bam" / ".bam" suffix.  The longer suffix is removed
    # first; if it is absent that replace is a no-op and the plain ".bam" one applies.
    names = [
        os.path.basename(p).replace(".sorted.bam", "").replace(".bam", "")
        for p in info["bam_path"]
    ]
    # Normalise group labels (trim whitespace, lower-case) so that " Control" and
    # "control" are treated the same later when selecting columns by group.
    return names, list(info["bam_path"]), list(info["sample_type"].str.strip().str.lower())


class GeneIndex:
    """Interval lookup over genes sorted by start with a running max of ends.

    Problem: for every read pair we need all genes whose interval overlaps the pair's
    template.  Testing every gene for every pair is O(genes) per pair, far too slow
    for millions of pairs.  This class answers each query with two binary searches.

    Data structure:
      * Genes are sorted by `start` (self.starts, with self.ends reordered to match).
      * self.max_end[i] = max(ends[0..i]) -- the running (prefix) maximum of the ends.
        It is non-decreasing by construction, so it can be binary-searched.

    Why the search window is correct for a query [start, end):
      * Upper bound `hi`: a gene can only overlap the query if it begins before the
        query ends (gene.start < end).  Because starts are sorted, all such genes form
        a prefix [0, hi) of the array; hi comes from searchsorted(starts, end, "left")
        which counts the starts strictly less than `end`.
      * Lower bound `lo`: a gene can only overlap if it finishes after the query begins
        (gene.end > start).  Ends are NOT sorted (a long gene can end after a later,
        shorter one), so we use the prefix max.  max_end[i] <= start means that every
        gene at positions 0..i ends at or before `start`, so none of them can overlap.
        Because max_end is non-decreasing, the first index where max_end > start is
        found with searchsorted(max_end, start, "right"); all genes before it are
        guaranteed non-overlapping, so it is safe to skip them.
      * Genes inside [lo, hi) are only CANDIDATES (some may still end before `start`
        because the running max can be driven by a different, longer gene), so each is
        verified with an exact overlap calculation.  No true overlap is ever missed,
        and the window is tiny for typical (non-nested) bacterial genes.
    """

    def __init__(self, genes):
        # Indices that would sort genes by start.  kind="stable" keeps ties (equal
        # starts) in their original file order so results are deterministic.
        self.order = np.argsort(genes["start"].to_numpy(), kind="stable")
        # Starts and ends rearranged into sorted order (same permutation for both so
        # that starts[i] and ends[i] still describe the same gene).
        self.starts = genes["start"].to_numpy()[self.order]
        self.ends = genes["end"].to_numpy()[self.order]
        # Prefix maximum of ends: max_end[i] = largest end among the first i+1 genes.
        self.max_end = np.maximum.accumulate(self.ends)

    def overlaps(self, start, end):
        """Yield (original gene index, overlap bp) for genes overlapping [start, end).

        The yielded index refers to the row in the ORIGINAL `genes` DataFrame (not the
        sorted position), so it can be used directly to index the counts array.
        """
        hi = np.searchsorted(self.starts, end, side="left")  # genes starting before `end`
        lo = np.searchsorted(self.max_end, start, side="right")  # first gene that can reach `start`
        # Examine only the candidate window [lo, hi) found above.
        for i in range(lo, hi):
            # Overlap length of two half-open intervals: the shared stretch runs from
            # the larger start to the smaller end.  Negative/zero means no overlap
            # (the candidate window can contain a few genes that end too early).
            ov = min(end, self.ends[i]) - max(start, self.starts[i])
            if ov > 0:
                # Map the sorted position back to the gene's original row number.
                yield self.order[i], ov


def count_bam(bam_path, index, n_genes):
    """Count read pairs per gene in one BAM (integer counts). A pair's template is the
    span from the leftmost aligned base to the rightmost aligned base of both mates; the
    pair is assigned wholly to the gene with the greatest template overlap, provided that
    overlap is >= 50% of the template length. Ties are ambiguous and not counted.

    Why count pairs, not reads: both mates come from one RNA fragment, so counting each
    mate would double-count the fragment.  Each qualifying pair adds exactly 1.

    Why the TEMPLATE span (not just the bases covered by the reads): the template is the
    whole stretch of the original fragment, including the unsequenced gap between the
    mates.  Using the full template span (leftmost to rightmost aligned base) follows
    the TA hints and treats the fragment as one interval.  The gene is therefore judged
    on how much of the fragment's extent it contains.

    Parameters
      bam_path : path to a coordinate- or name-sorted BAM file.
      index    : GeneIndex built from the annotation.
      n_genes  : number of genes (length of the returned counts vector).
    Returns a numpy int64 array of length n_genes, in the original GFF gene order.
    """
    counts = np.zeros(n_genes, dtype=np.int64)
    # Holds the first-seen mate of each pair until its partner shows up.  Keyed by
    # query name (the read name), which both mates of a fragment share.
    pending = {}  # query_name -> first mate seen
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        # until_eof=True streams every record from start to end of the file without
        # needing a BAM index (.bai) and including unmapped reads; we filter ourselves.
        for read in bam.fetch(until_eof=True):
            # Quality filters: skip anything that is not a clean, confidently mapped,
            # properly paired primary alignment.
            if (
                read.is_unmapped  # read did not align anywhere
                or read.is_secondary  # alternative alignment of a multi-mapper
                or read.is_supplementary  # chimeric piece of a split alignment
                or read.is_qcfail  # failed the sequencer's quality check
                or read.is_duplicate  # PCR/optical duplicate (would inflate counts)
                or not read.is_proper_pair  # mates not aligned in the expected orientation/distance
                or read.mapping_quality < MIN_MAPQ  # position not trusted (see MIN_MAPQ)
            ):
                continue
            # Mate pairing via query name: if the partner was already seen, pop it out
            # of the dictionary (completing the pair); otherwise this is the first mate.
            mate = pending.pop(read.query_name, None)
            if mate is None:
                # First mate of the pair: park it and wait for its partner.  If the
                # partner was filtered out above, this entry just stays unused, so
                # half-filtered pairs are never counted.
                pending[read.query_name] = read
                continue
            # Template span as a 0-based half-open interval: leftmost aligned base of
            # either mate (min of starts) to one past the rightmost aligned base
            # (max of reference_end, which is already exclusive).
            t_start = min(mate.reference_start, read.reference_start)
            t_end = max(mate.reference_end, read.reference_end)
            # All genes touching the template as (gene index, overlap bp), best
            # (largest overlap) first.
            hits = sorted(index.overlaps(t_start, t_end), key=lambda h: h[1], reverse=True)
            # The fragment lies in no annotated gene (e.g. intergenic): not counted.
            if not hits:
                continue
            # Tie handling: if the top two genes overlap by exactly the same amount
            # (e.g. overlapping genes on opposite strands) we cannot say which one the
            # fragment came from, so we drop it rather than guess.
            if len(hits) > 1 and hits[0][1] == hits[1][1]:
                continue  # ambiguous: equal overlap with two genes
            # 50% overlap gate: the winning gene must contain at least half of the
            # template, otherwise the fragment mostly lies elsewhere (e.g. spans a
            # gene boundary or is mostly intergenic) and is not assigned to the gene.
            if hits[0][1] < MIN_OVERLAP_FRAC * (t_end - t_start):
                continue  # best gene covers < 50% of the template
            # Passed every test: credit exactly one fragment to the winning gene
            # (hits[0][0] is its original row index).
            counts[hits[0][0]] += 1
    return counts


def compute_tpm(counts, lengths):
    """TPM per sample (column). counts: DataFrame genes x samples; lengths: bp per gene.

    TPM = Transcripts Per Million.  Steps (done independently for each sample):
      1. RPK (reads per kilobase) = count / (gene length in kb).  This removes the
         gene-length bias: a gene twice as long is expected to yield twice the reads.
      2. Per-sample scaling factor = (sum of all RPK values in that sample) / 1,000,000.
         This is the sample's total "expression mass" in millions, capturing
         sequencing depth.
      3. TPM = RPK / scaling factor.  Dividing by the scaling factor rescales every
         sample so its TPM values sum to 1,000,000, making samples comparable (a TPM
         of 50 means 50 out of every million transcripts in that sample).
    """
    # Step 1: divide each row (gene) by its length in kilobases.  axis=0 aligns the
    # `lengths` Series with the row index (gene ids) rather than with the columns.
    rpk = counts.div(lengths / 1000.0, axis=0)
    # Step 2: column sums give each sample's total RPK; /1e6 puts it in "millions".
    scaling = rpk.sum(axis=0) / 1e6
    # Step 3: divide each column (sample) by its own scaling factor (axis=1 aligns
    # on column names).  A sample with zero total would divide by zero, so its factor
    # is replaced with NaN first (result NaN rather than inf) and the NaNs are then
    # set to 0.0, giving an all-zero TPM column for an empty sample.
    return rpk.div(scaling.replace(0, np.nan), axis=1).fillna(0.0)


def differential_expression(tpm, types):
    """Compare control vs treatment TPM for every gene and return a results DataFrame.

    Parameters
      tpm   : DataFrame (genes x samples) of TPM values.
      types : list of group labels ("control"/"treatment"), one per sample column, in
              the same order as the columns of `tpm`.

    For each gene the output has: mean and median TPM of each group, the log2 fold
    change (treatment relative to control; positive = higher in treatment), and a
    two-sided Mann-Whitney U test p-value.

    Why Mann-Whitney U: with only a few replicates per group, the distribution of
    expression values cannot be assumed normal, so a rank-based non-parametric test
    (which asks whether values from one group tend to be larger than the other) is
    safer than a t-test.
    """
    # Array of labels so that boolean masks can be built with elementwise ==.
    types = np.array(types)
    # Select the columns belonging to each group; .to_numpy() gives plain
    # (genes x replicates) arrays, which are faster to iterate and aggregate.
    ctrl = tpm.loc[:, types == "control"].to_numpy()
    trt = tpm.loc[:, types == "treatment"].to_numpy()
    # Per-gene mean across replicates (axis=1 averages along the sample dimension).
    mean_c, mean_t = ctrl.mean(axis=1), trt.mean(axis=1)
    pvals = []
    # Walk the genes in lockstep: each c / t is the vector of replicate TPMs of one gene.
    for c, t in zip(ctrl, trt):
        try:
            p = stats.mannwhitneyu(t, c, alternative="two-sided").pvalue
        except ValueError:  # e.g. all values identical
            # Mann-Whitney fallback: scipy raises ValueError when every value is
            # identical (no rank information), so there is no evidence of a
            # difference; report the least significant p-value, 1.0.
            p = 1.0
        # Also guard against a NaN p-value (e.g. degenerate input) by mapping it to
        # 1.0 so the output never contains missing p-values.
        pvals.append(1.0 if np.isnan(p) else p)
    return pd.DataFrame(
        {
            "gene_id": tpm.index,
            "mean_control": mean_c,
            "median_control": np.median(ctrl, axis=1),
            "mean_treatment": mean_t,
            "median_treatment": np.median(trt, axis=1),
            # log2 fold change = log2(treatment / control) on the group MEAN TPMs.
            # The pseudocount (a small constant) is added to both means so that a
            # gene with zero expression in one group gives a finite value instead of
            # log2(0) = -inf or a division by zero.  It is small relative to typical
            # TPM values, so it barely affects well-expressed genes.
            # A value of +1 means 2x higher in treatment, -1 means 2x lower, 0 = equal.
            "log2_fold_change": np.log2((mean_t + PSEUDOCOUNT) / (mean_c + PSEUDOCOUNT)),
            "p_value": pvals,
        }
    )


def main():
    """Command-line entry point: parse arguments, count reads, then run the statistics."""
    # Three positional command-line arguments; see the module docstring for usage.
    ap = argparse.ArgumentParser(description="RNA-seq counting and differential expression")
    ap.add_argument("sample_info_tsv")
    ap.add_argument("annotation_gff")
    ap.add_argument("output_directory")
    args = ap.parse_args()

    # Create the output directory if needed; exist_ok avoids an error on re-runs.
    os.makedirs(args.output_directory, exist_ok=True)
    # Load annotation and sample table, and build the fast interval index once so
    # that all BAM files can reuse it.
    genes = parse_gff(args.annotation_gff)
    names, bams, types = read_samples(args.sample_info_tsv)
    index = GeneIndex(genes)

    # Part 1: read counting
    # Start with an empty table whose rows are gene ids (every GFF gene, even those
    # that end up with zero reads); one column of counts is added per sample.  The
    # count vectors are in GFF order, matching this index.
    counts = pd.DataFrame(index=genes["gene_id"])
    for name, bam in zip(names, bams):
        # Progress goes to stderr so it does not mix with any stdout output.
        print(f"counting {name} ...", file=sys.stderr)
        counts[name] = count_bam(bam, index, len(genes))
    counts.index.name = "gene_id"
    # Raw integer counts for all genes (unfiltered), tab-separated.
    counts.to_csv(os.path.join(args.output_directory, "gene_counts.tsv"), sep="\t")

    # Part 2: filter, TPM, statistics
    # Low-count filter: keep genes whose total over all samples reaches the minimum.
    # This is a boolean Series aligned with `counts`.
    keep = counts.sum(axis=1) >= MIN_TOTAL_COUNT
    # Gene lengths keyed by gene_id so they can be looked up for the kept genes.
    lengths = genes.set_index("gene_id")["length"]
    kept = counts[keep]
    # TPM is computed AFTER filtering, so the per-sample normalisation (the "million"
    # total) is over the retained genes only.  lengths[kept.index] picks the lengths
    # in the same gene order as the filtered counts.
    tpm = compute_tpm(kept, lengths[kept.index])
    de = differential_expression(tpm, types)
    de.to_csv(
        os.path.join(args.output_directory, "differential_expression.tsv"),
        sep="\t",
        index=False,  # gene_id is already a regular column, so skip the row numbers
        float_format="%.10g",  # 10 significant digits: compact but precise
    )
    # Report how many genes survived the filter, e.g. "4100 / 4300 genes ...".
    print(f"{keep.sum()} / {len(genes)} genes passed the count filter", file=sys.stderr)


if __name__ == "__main__":
    main()
