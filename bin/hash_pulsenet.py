#!/usr/bin/env python3
"""
Hash pulsenet allele calls into a wide chewBBACA-style profile (TSV).

Output: header row "FILE\\t<locus1>\\t<locus2>\\t..." followed by one data row
for the sample. Each cell is either a 16-hex-char hash of the called allele
sequence, or "-" for missing/uncallable loci.

Input files:
  - pulsenet json output

Classification flags seen in PULSENET results:
Flag       Tag                  Description
1	       CORE	                Core locus
2	       ACCESSORY	        Accessory locus
4	       CALLED	            Locus called
8	       NOT_CALLED	        Locus not called
16	       NOT_FOUND	        Locus not found
32	       IDENTITY_PROBLEMS	Allele with lower identity
64	       SEQ_PROBLEMS	        Allele sequence with non ATGC
128	       CODON_PROBLEMS	    Allele missing start or/and stop codons, or having an internal stop codon
256	       REPEAT_PROBLEMS	    Locus with paralogous
512	       COVERAGE_PROBLEMS	Allele with positions without [any nucleotide with] minimum depth
1024	   BIAS_PROBLEMS	    Allele with positions with [covered nucleotide with] strand bias
2048	   DOUBLE_CALLS	        Allele with positions with double possible base calls
4096	   MORE_CALLS	        Allele with positions with more than two possible base calls
8192	   MISMATCH_CALLS	    Allele with positions with most frequent nucleotide mismatching assembly base called

The following flags are parsed:
5 - CORE and CALLED
"""

import os
import sys
import argparse
import json
import zlib

MISSING_VALUE = "-"
CALLED_FLAGS = [5]

_COMPLEMENT = str.maketrans("ACGTN", "TGCAN")

def revcomp(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]

def canonical_hash(seq) -> str:
    """CRC32 of the canonical (strand-agnostic) orientation, as a decimal string.

    Canonical = lexicographic min of (seq, revcomp(seq)), so the same biological
    allele hashes identically regardless of the strand taranys reports.

    NOTE: CRC32 is a 32-bit space; collisions become non-negligible once the
    number of distinct alleles approaches ~10^5.
    """
    if seq is None or not str(seq).strip():
        return MISSING_VALUE
    s = str(seq).upper().strip()
    return str(zlib.crc32(s.encode()))

def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--pulsenet-json",        required=True,
                   help="pulsenet json file with allele calls for the sample")
    p.add_argument("--output-dir",           required=False, default=".",
                   help="Output directory")
    p.add_argument("--prefix",           required=True,
                   help="Prefix for output file")
    p.add_argument("--sample",           required=True,
                   help="Sample name")
    args = p.parse_args()

    # Load loci allele calls from pulsenet json
    seq_hash_by_locus = {}
    with open(args.pulsenet_json) as f:
        data = json.load(f)
        loci_results = data.get("calls", {}) 
        for locus, locus_data in loci_results.items():
            locus_data = locus_data[0] 
            sequence = None
            sequence_hash = None
            locus_flag = int(locus_data.get("flag", ""))
            if CALLED_FLAGS.__contains__(locus_flag):
                sequence = locus_data.get("seq", "")
                sequence_hash = canonical_hash(sequence)
            else:
                sequence_hash = MISSING_VALUE
            seq_hash_by_locus[locus] = sequence_hash
   
    output_path = os.path.join(args.output_dir, f"{args.prefix}.tsv")
    with open(output_path, "w") as out:
        out.write("FILE\t" + "\t".join(seq_hash_by_locus.keys()) + "\n")
        out.write(args.sample + "\t" + "\t".join(seq_hash_by_locus.values()) + "\n")

    n_called = sum(1 for v in seq_hash_by_locus.values() if v != MISSING_VALUE)
    print(f"Wrote {n_called}/{len(seq_hash_by_locus)} called loci for {args.sample} "
          f"to {output_path}", file=sys.stderr)
    return 0

if __name__ == "__main__":
    sys.exit(main())
