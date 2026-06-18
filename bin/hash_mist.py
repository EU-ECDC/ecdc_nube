#!/usr/bin/env python3
"""
Hash mist allele calls into a wide chewBBACA-style profile (TSV).

Output: header row "FILE\\t<locus1>\\t<locus2>\\t..." followed by one data row
for the sample. Each cell is either a 16-hex-char hash of the called allele
sequence, or "-" for missing/uncallable loci.

Input files:
  - mist json output
  - schema-dir

Classification codes seen in MIST results:
  - EXACT            exact match to database     -> hash sequence
  - NOVEL            novel allele                -> hash sequence
  - MULTI            multiple exact matches      -> "-"
  - EDGE             contig edge (incomplete)    -> "-"
  - INDEL            length anomalies             -> "-"
  - ABSENT           locus not found             -> "-"
  
"""

import argparse
import json
import os
import sys
import zlib
from glob import glob
from Bio import SeqIO

MISSING_TAGS = ("MULTI", "EDGE", "INDEL", "ABSENT")
MISSING_VALUE = "-"

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
    canonical = min(s, revcomp(s))
    return str(zlib.crc32(canonical.encode()))

def list_schema_loci(schema_dir: str) -> list:
    """One locus per FASTA file in the schema, sorted for stable column order."""
    extensions = ("*.fasta", "*.fa", "*.fna")
    locus_fasta_files = {}

    fastas = sorted(
        f
        for ext in extensions
        for f in glob(os.path.join(schema_dir, "*", ext))
        if os.path.basename(os.path.dirname(f)) ==
        os.path.splitext(os.path.basename(f))[0]
    )

    if not fastas:
        raise SystemExit(f"ERROR: no FASTA files in schema dir {schema_dir}")

    for f in fastas:
        locus = os.path.basename(os.path.dirname(f))
        locus_fasta_files[locus] = f

    return locus_fasta_files

def is_missing(codes: list) -> bool:
    """A code is missing if it matches a missing category"""
    if not codes:
        return True
    for code in codes:
        if not code or code.lower() == "nan":
            return True
        if code in MISSING_TAGS:
            return True
    return False

def get_sequence_from_schema(locus_fasta_files: dict, locus: str, allele: str) -> str:
    """Extract the allele sequence from the locus fasta file in the schema"""
    fasta_file = locus_fasta_files.get(locus)
    try:
        with open(fasta_file, "r") as f:
            for record in SeqIO.parse(f, "fasta"):
                if record.id == f"{locus}_{allele}":
                    return str(record.seq)
    except FileNotFoundError:
        print(f"Error: File '{fasta_file}' not found.")
    except Exception as e:
        print(f"An error occurred: {e}")

def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--schema-dir",       required=True,
                   help="Directory of locus FASTAs (defines locus column universe)")
    p.add_argument("--mist-json",        required=True,
                   help="mist json file with allele calls for the sample")
    p.add_argument("--output",           required=True,
                   help="Output tsv file")
    p.add_argument("--prefix",           required=True,
                   help="Prefix for output file")
    p.add_argument("--sample",           required=True,
                   help="Sample name")
    args = p.parse_args()

    locus_fasta_files = list_schema_loci(args.schema_dir)
    loci = sorted(locus_fasta_files.keys())

    # Load loci allele calls from mist json
    seq_hash_by_locus = {}
    with open(args.mist_json) as f:
        data = json.load(f)
        loci_results = data.get("alleles", {})
        for locus, locus_data in loci_results.items():
            sequence = None
            sequence_hash = None
            locus_allele = locus_data.get("allele_str", "")
            locus_tags = locus_data.get("tags", [])
            if is_missing(locus_tags):
                sequence_hash = MISSING_VALUE
            elif "NOVEL" in locus_tags:
                locus_results = locus_data.get("allele_results", {})
                sequence = locus_results[0].get("sequence", "")
            elif "EXACT" in locus_tags:
                print(locus, locus_allele)
                sequence = get_sequence_from_schema(locus_fasta_files, locus, locus_allele)
            else:
                print("ERROR: Something went wrong.")
                sys.exit(1)
            if sequence is not None:
                sequence_hash = canonical_hash(sequence)
            seq_hash_by_locus[locus] = sequence_hash

    # Build the row in schema-defined locus order.
    row_values = []
    for locus in loci:
        row_values.append(seq_hash_by_locus[locus])
   
    # Write TSV (chewBBACA-style: FILE column + locus columns).
    with open(os.path.join(args.output, f"{args.prefix}.tsv"), "w") as out:
        out.write("FILE\t" + "\t".join(loci) + "\n")
        out.write(args.sample + "\t" + "\t".join(row_values) + "\n")

    n_called = sum(1 for v in row_values if v != MISSING_VALUE)
    print(f"Wrote {n_called}/{len(loci)} called loci for {args.sample} "
          f"to {args.output}", file=sys.stderr)
    return 0

if __name__ == "__main__":
    sys.exit(main())
