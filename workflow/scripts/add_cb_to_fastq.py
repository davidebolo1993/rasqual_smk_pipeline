#!/usr/bin/env python3
"""
Add CB (cell barcode) tags from BAM to WASP-formatted FASTQ headers.
The script matches reads by parsing the original read name from the WASP format.
"""

import gzip
import pysam
import sys
from collections import defaultdict

def parse_wasp_readname(fastq_name):
    """
    Parse WASP fastq read name format:
    <orig_name>.<coordinate>.<read_number>.<total_read_number>
    Returns the original read name.
    """
    parts = fastq_name.split('.')
    if len(parts) >= 4:
        return parts[0]
    return fastq_name

def build_cb_dict(bam_file):
    """
    Build a dictionary mapping read names to CB tags from BAM file.
    """
    cb_dict = {}
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for read in bam:
            if read.has_tag("CB"):
                cb_tag = read.get_tag("CB")
                cb_dict[read.query_name] = cb_tag
    return cb_dict

def process_fastq(input_fq, output_fq, cb_dict, is_gzipped=True):
    """
    Process FASTQ file and add CB tag to headers.
    """
    open_func = gzip.open if is_gzipped else open
    mode_in = 'rt' if is_gzipped else 'r'
    mode_out = 'wt' if is_gzipped else 'w'
    
    with open_func(input_fq, mode_in) as fin, open_func(output_fq, mode_out) as fout:
        while True:
            # Read FASTQ record (4 lines)
            header = fin.readline()
            if not header:
                break
            seq = fin.readline()
            plus = fin.readline()
            qual = fin.readline()
            
            # Parse original read name from WASP format
            fastq_name = header.strip()[1:]  # Remove '@'
            orig_name = parse_wasp_readname(fastq_name)
            
            # Add CB tag if available
            if orig_name in cb_dict:
                new_header = f"@{fastq_name} CB:Z:{cb_dict[orig_name]}\n"
            else:
                new_header = header
            
            # Write modified record
            fout.write(new_header)
            fout.write(seq)
            fout.write(plus)
            fout.write(qual)

def main():
    if len(sys.argv) != 5:
        print("Usage: python add_cb_to_fastq.py <input.bam> <input.fq1.gz> <input.fq2.gz> <output_prefix>")
        print("Example: python add_cb_to_fastq.py test.to.remap.bam test.remap.fq1.gz test.remap.fq2.gz test.remap.cb")
        sys.exit(1)
    
    bam_file = sys.argv[1]
    fq1_in = sys.argv[2]
    fq2_in = sys.argv[3]
    output_prefix = sys.argv[4]
    
    fq1_out = f"{output_prefix}.fq1.gz"
    fq2_out = f"{output_prefix}.fq2.gz"
    
    print(f"Building CB dictionary from {bam_file}...")
    cb_dict = build_cb_dict(bam_file)
    print(f"Found {len(cb_dict)} reads with CB tags")
    
    print(f"Processing {fq1_in}...")
    process_fastq(fq1_in, fq1_out, cb_dict, is_gzipped=True)
    
    print(f"Processing {fq2_in}...")
    process_fastq(fq2_in, fq2_out, cb_dict, is_gzipped=True)
    
    print(f"Done! Output files: {fq1_out}, {fq2_out}")

if __name__ == "__main__":
    main()

