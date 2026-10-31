"""Permutation test for the overlap between two sets of genomic regions.

This script calculates the number of bases that merged Set A and merged
Set B have in common. It then moves each Set B region to a random start on
its chromosome, and calculates the overlap again for each permutation.
The p-value is one-sided (upper tail) with a +1 correction.

Usage:
    python Evan_Fitzhugh_permutation_test.py <setA_bed> <setB_bed> \
        <genome_fai> <output_file> [num_permutations]
"""

import random
import sys
from collections import defaultdict

DEFAULT_PERMUTATIONS = 10000
SEED = 205


def read_fai(path):
    """Read the FASTA index. Give a dict of chromosome name to length."""
    lengths = {}
    with open(path) as fh:
        for line in fh:
            fields = line.split("\t")
            if len(fields) >= 2:
                lengths[fields[0]] = int(fields[1])
    return lengths


def read_bed(path):
    """Read a 3-column BED file.

    The function gives a dict of chromosome name to a list of (start, end)
    tuples, and the number of regions in the file before the merge.
    """
    regions = defaultdict(list)
    count = 0
    with open(path) as fh:
        for line in fh:
            fields = line.split()
            if len(fields) < 3:
                continue
            regions[fields[0]].append((int(fields[1]), int(fields[2])))
            count += 1
    return regions, count


def merge(intervals):
    """Sort and merge one list of (start, end) tuples on one chromosome.

    Regions that overlap or touch become one region. After the merge, no
    base is in more than one region.
    """
    if not intervals:
        return []
    intervals = sorted(intervals)
    merged = []
    cur_start, cur_end = intervals[0]
    for start, end in intervals[1:]:
        if start <= cur_end:
            if end > cur_end:
                cur_end = end
        else:
            merged.append((cur_start, cur_end))
            cur_start, cur_end = start, end
    merged.append((cur_start, cur_end))
    return merged


def merge_all(regions):
    """Merge the regions on each chromosome."""
    return {chrom: merge(ivs) for chrom, ivs in regions.items()}


def total_bases(merged):
    """Count the bases in a set of merged regions."""
    return sum(end - start for ivs in merged.values() for start, end in ivs)


def overlap_sorted(a, b):
    """Count the bases that two sorted, merged lists have in common.

    The function uses two pointers and moves through each list one time.
    At each step, the pointer for the region that ends first moves forward.
    """
    i = j = 0
    n_a, n_b = len(a), len(b)
    total = 0
    while i < n_a and j < n_b:
        a_start, a_end = a[i]
        b_start, b_end = b[j]
        lo = a_start if a_start > b_start else b_start
        hi = a_end if a_end < b_end else b_end
        if hi > lo:
            total += hi - lo
        if a_end < b_end:
            i += 1
        else:
            j += 1
    return total


def overlap(merged_a, merged_b):
    """Add the overlap for all chromosomes that are in the two sets."""
    return sum(
        overlap_sorted(merged_a[chrom], ivs)
        for chrom, ivs in merged_b.items()
        if chrom in merged_a
    )


def shuffle_lengths(regions_b, chrom_lengths):
    """Get the data that each permutation uses for Set B.

    For each chromosome, the function gives the region lengths and the
    largest permitted start for each region (chrom_len - L).
    """
    plan = {}
    for chrom, ivs in regions_b.items():
        chrom_len = chrom_lengths[chrom]
        lengths = [end - start for start, end in ivs]
        plan[chrom] = [(chrom_len - length, length) for length in lengths]
    return plan


def shuffled_merged_b(plan, rng):
    """Put each Set B region at a random start on its chromosome, then merge.

    Each new start is a random integer from 0 to chrom_len - L. The two
    limits are included. Shuffled regions can overlap, and the merge
    removes the bases that they have in common.
    """
    randint = rng.randint
    merged = {}
    for chrom, items in plan.items():
        new_ivs = []
        for max_start, length in items:
            start = randint(0, max_start)
            new_ivs.append((start, start + length))
        merged[chrom] = merge(new_ivs)
    return merged


def write_results(path, results):
    """Write the metrics to a tab-separated file with a header line."""
    with open(path, "w") as fh:
        fh.write("metric\tvalue\n")
        for name, value in results:
            fh.write(f"{name}\t{value}\n")


def main(argv):
    if len(argv) not in (5, 6):
        sys.exit(
            "usage: python Evan_Fitzhugh_permutation_test.py <setA_bed> "
            "<setB_bed> <genome_fai> <output_file> [num_permutations]"
        )
    set_a_path, set_b_path, fai_path, out_path = argv[1:5]
    n_perm = int(argv[5]) if len(argv) == 6 else DEFAULT_PERMUTATIONS
    if n_perm < 1:
        sys.exit("num_permutations must be 1 or more")

    chrom_lengths = read_fai(fai_path)
    regions_a, count_a = read_bed(set_a_path)
    regions_b, count_b = read_bed(set_b_path)

    # Set A stays the same. Merge it one time only.
    merged_a = merge_all(regions_a)
    merged_b = merge_all(regions_b)

    observed = overlap(merged_a, merged_b)

    rng = random.Random(SEED)
    plan = shuffle_lengths(regions_b, chrom_lengths)
    n_ge = 0
    for k in range(n_perm):
        if overlap(merged_a, shuffled_merged_b(plan, rng)) >= observed:
            n_ge += 1
        if (k + 1) % 1000 == 0:
            print(f"{k + 1}/{n_perm} permutations done", file=sys.stderr)

    p_value = (1 + n_ge) / (1 + n_perm)

    write_results(
        out_path,
        [
            ("observed_overlap", observed),
            ("p_value", p_value),
            ("num_permutations", n_perm),
            ("setA_regions", count_a),
            ("setB_regions", count_b),
            ("setA_total_bases", total_bases(merged_a)),
            ("setB_total_bases", total_bases(merged_b)),
        ],
    )
    print(f"observed_overlap={observed} p_value={p_value}", file=sys.stderr)


if __name__ == "__main__":
    main(sys.argv)
