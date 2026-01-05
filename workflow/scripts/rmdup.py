#!/usr/bin/env python3
"""
Cell-barcode aware unbiased duplicate removal for paired-end reads.
Removes duplicate FRAGMENTS (read pairs) at random within each cell.
"""

import pysam
import sys
import random
from collections import defaultdict

def get_fragment_key(read, use_cb=True):
    """
    Generate a key for identifying duplicate fragments (read pairs).
    Key includes: CB, chromosome, fragment start/end positions, strand.
    """
    if read.is_unmapped:
        return None, None
    
    # Get cell barcode
    cb = ""
    if use_cb and read.has_tag("CB"):
        cb = read.get_tag("CB")
    
    # For paired-end reads, use the leftmost position of the pair as fragment start
    if read.is_paired and not read.mate_is_unmapped:
        # Determine fragment boundaries
        if read.reference_start < read.next_reference_start:
            # Read1 is leftmost
            frag_start = read.reference_start
            frag_end = read.next_reference_start + abs(read.template_length)
            strand1 = read.is_reverse
            strand2 = read.mate_is_reverse
        else:
            # Read2 is leftmost
            frag_start = read.next_reference_start
            frag_end = read.reference_start + abs(read.template_length)
            strand1 = read.mate_is_reverse
            strand2 = read.is_reverse
        
        # Create fragment key (same for both read1 and read2)
        key = (
            cb,
            read.reference_name,
            frag_start,
            frag_end,
            strand1,
            strand2
        )
    else:
        # Single-end or unpaired
        if read.is_reverse:
            pos = read.reference_end
        else:
            pos = read.reference_start
        key = (cb, read.reference_name, pos, read.is_reverse)
    
    return key, cb

def rmdup_per_cell_fragments(input_bam, output_bam, seed=None):
    """
    Remove duplicate fragments (read pairs) at random per cell barcode.
    Keeps both read1 and read2 for selected fragments.
    """
    if seed is not None:
        random.seed(seed)
    
    # First pass: group reads by fragment key and read name
    print("Reading BAM and grouping fragments...", file=sys.stderr)
    fragment_groups = defaultdict(lambda: defaultdict(list))  # key -> qname -> [reads]
    cb_total_fragments = defaultdict(int)
    
    with pysam.AlignmentFile(input_bam, "rb") as inbam:
        for read in inbam:
            key, cb = get_fragment_key(read, use_cb=True)
            if key is not None:
                qname = read.query_name
                fragment_groups[key][qname].append(read)
                if cb and read.is_read1:  # Count fragments once per pair
                    cb_total_fragments[cb] += 1
    
    # Calculate statistics
    cb_kept_fragments = defaultdict(int)
    cb_removed_fragments = defaultdict(int)
    
    total_fragments = 0
    kept_fragments = 0
    removed_fragments = 0
    
    selected_qnames = set()  # Track which fragments to keep
    
    for key, qname_dict in fragment_groups.items():
        cb = key[0]
        n_fragments = len(qname_dict)
        total_fragments += n_fragments
        
        # Randomly select ONE fragment (qname) from this duplicate group
        selected_qname = random.choice(list(qname_dict.keys()))
        selected_qnames.add(selected_qname)
        
        kept_fragments += 1
        removed_fragments += n_fragments - 1
        
        if cb:
            cb_kept_fragments[cb] += 1
            cb_removed_fragments[cb] += n_fragments - 1
    
    # Overall statistics
    print(f"Found {len(fragment_groups)} unique fragment positions", file=sys.stderr)
    print(f"Total fragments: {total_fragments}", file=sys.stderr)
    print(f"Kept fragments: {kept_fragments}", file=sys.stderr)
    print(f"Removed fragments: {removed_fragments} ({100*removed_fragments/total_fragments:.2f}%)", file=sys.stderr)
    
    # Per-cell statistics
    print("\nPer-cell barcode statistics:", file=sys.stderr)
    print(f"{'Cell Barcode':<25} {'Total':<10} {'Kept':<10} {'Removed':<10} {'% Dup':<8}", file=sys.stderr)
    print("-" * 70, file=sys.stderr)
    
    for cb in sorted(cb_total_fragments.keys()):
        total = cb_total_fragments[cb]
        kept = cb_kept_fragments[cb]
        removed = cb_removed_fragments[cb]
        pct = 100 * removed / total if total > 0 else 0
        print(f"{cb:<25} {total:<10} {kept:<10} {removed:<10} {pct:<8.2f}", file=sys.stderr)
    
    # Second pass: write only selected fragments (both read1 and read2)
    print("\nWriting output...", file=sys.stderr)
    reads_written = 0
    with pysam.AlignmentFile(input_bam, "rb") as inbam:
        with pysam.AlignmentFile(output_bam, "wb", template=inbam) as outbam:
            for read in inbam:
                if read.query_name in selected_qnames:
                    outbam.write(read)
                    reads_written += 1
    
    print(f"Wrote {reads_written} reads ({reads_written//2} fragments)", file=sys.stderr)
    print("Done!", file=sys.stderr)

def main():
    if len(sys.argv) < 3:
        print("Usage: python rmdup.py <input.bam> <output.bam> [random_seed]")
        print("\nPerforms unbiased fragment-level duplicate removal per cell barcode (CB tag).")
        print("Randomly selects one fragment (read pair) from each duplicate set within each cell.")
        sys.exit(1)
    
    input_bam = sys.argv[1]
    output_bam = sys.argv[2]
    seed = int(sys.argv[3]) if len(sys.argv) > 3 else None
    
    rmdup_per_cell_fragments(input_bam, output_bam, seed)

if __name__ == "__main__":
    main()

