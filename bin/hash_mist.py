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
    return str(zlib.crc32(s.encode()))

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
    loci = []
    for f in fastas:
        loci.append(os.path.basename(os.path.dirname(f)))
    return loci

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

def get_sequence_from_fasta(fasta_file: str, seq_id: str) -> str:
    """Extract a fasta sequence based on sequence id"""
    try:
        with open(fasta_file, "r") as f:
            for record in SeqIO.parse(f, "fasta"):
                if record.id == seq_id:
                    return str(record.seq)
    except FileNotFoundError:
        print(f"Error: File '{fasta_file}' not found.")
    except Exception as e:
        print(f"An error occurred: {e}")

def get_sequence_from_assembly(assembly: str, align_start: int, align_end: int, seq_id: str, strand: str):
    """Extract the allele sequence from the assembly fasta file"""
    contig_sequence = get_sequence_from_fasta(assembly, seq_id)
    allele_sequence = contig_sequence[align_start - 1: align_end]
    if strand == "-":
        return revcomp(allele_sequence)
    return allele_sequence


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--schema-dir",       required=True,
                   help="Directory of locus FASTAs (defines locus column universe)")
    p.add_argument("--mist-json",        required=True,
                   help="mist json file with allele calls for the sample")
    p.add_argument("--output-dir",           required=False, default=".",
                   help="Output directory")
    p.add_argument("--prefix",           required=True,
                   help="Prefix for output file")
    p.add_argument("--sample",           required=True,
                   help="Sample name")
    p.add_argument("--assembly",         required=True,
                   help="Assembly fasta")
    args = p.parse_args()

    loci_sequential = list_schema_loci(args.schema_dir)
    loci = sorted(loci_sequential)

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
            elif (len(locus_tags) == 1) and (locus_tags[0] == "NOVEL"):
                locus_results = locus_data.get("allele_results", {})
                sequence = locus_results[0].get("sequence", "")
            elif (len(locus_tags) == 1) and (locus_tags[0] == "EXACT"):
                locus_results = locus_data.get("allele_results", {})[0]
                alignment_results = locus_results.get("alignment", {})
                align_start = alignment_results.get("start", int)
                align_end = alignment_results.get("end", int)
                seq_id = alignment_results.get("seq_id", str)
                strand = alignment_results.get("strand", str)
                sequence = get_sequence_from_assembly(args.assembly, align_start, align_end, seq_id, strand)
            else:
                print(f"ERROR: Parsing is not supported for specified locus tags: {locus_tags}")
                sys.exit(1)
            if sequence is not None:
                sequence_hash = canonical_hash(sequence)
            seq_hash_by_locus[locus] = sequence_hash

    # Build the row in schema-defined locus order.
    row_values = []
    for locus in loci:
        row_values.append(seq_hash_by_locus[locus])
   
    # Write TSV (chewBBACA-style: FILE column + locus columns).
    output_path = os.path.join(args.output_dir, f"{args.prefix}.tsv")
    with open(output_path, "w") as out:
        out.write("FILE\t" + "\t".join(loci) + "\n")
        out.write(args.sample + "\t" + "\t".join(row_values) + "\n")

    n_called = sum(1 for v in row_values if v != MISSING_VALUE)
    print(f"Wrote {n_called}/{len(loci)} called loci for {args.sample} "
          f"to {output_path}", file=sys.stderr)
    return 0

if __name__ == "__main__":
    sys.exit(main())
