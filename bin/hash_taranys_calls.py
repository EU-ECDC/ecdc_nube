#!/usr/bin/env python3
"""
Hash taranys allele calls into a wide chewBBACA-style profile (TSV).

Output: header row "FILE\\t<locus1>\\t<locus2>\\t..." followed by one data row
for the sample. Each cell is either a 16-hex-char hash of the called allele
sequence, or "-" for missing/uncallable loci.

Input files (taranys v3.0.1 output):
  - contig_alignment_info.csv  (per-(sample,locus,hit) rows with sample_allele_seq)
  - allele_calling_match.csv   (wide: rows=samples, cols=loci, values=codes)
  - schema-dir                 (defines the locus column universe)

Classification codes seen in allele_calling_match.csv:
  - EXC_<locus>_<n>  exact match           -> hash sample_allele_seq
  - INF_<n>          novel/inferred allele -> hash sample_allele_seq
  - TPR_<n>          truncated / partial   -> "-"
  - NIPH_<...>       paralog               -> "-"
  - PLOT             contig edge           -> "-"
  - ASM, ALM         length anomalies      -> "-"
  - LNF              locus not found       -> "-"
"""

import argparse
import csv
import zlib
import os
import sys
from glob import glob


# Code prefixes that count as missing / uncallable.
MISSING_PREFIXES = ("LNF", "PLOT", "ASM", "ALM", "NIPH", "NIPHEM", "TPR")
MISSING_VALUE = "-"

# DNA reverse complement (with N).
_COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


def revcomp(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def get_hash(seq) -> str:
    """CRC32 hash of the sequence, as a decimal string.

    NOTE: CRC32 is a 32-bit space; collisions become non-negligible once the
    number of distinct alleles approaches ~10^5.
    """
    if seq is None or not str(seq).strip():
        return MISSING_VALUE
    s = str(seq).upper().strip()
    return str(zlib.crc32(s.encode()))


def list_schema_loci(schema_dir: str) -> list:
    """One locus per FASTA file in the schema, sorted for stable column order."""
    fastas = sorted(
        glob(os.path.join(schema_dir, "*.fasta"))
        + glob(os.path.join(schema_dir, "*.fa"))
        + glob(os.path.join(schema_dir, "*.fna"))
    )
    if not fastas:
        raise SystemExit(f"ERROR: no FASTA files in schema dir {schema_dir}")
    return [os.path.splitext(os.path.basename(f))[0] for f in fastas]


def is_missing(code: str) -> bool:
    """A code is missing if its prefix (before '_') matches a missing category,
    or it's empty / blank."""
    if not code or code.lower() == "nan":
        return True
    prefix = code.split("_", 1)[0]
    return prefix in MISSING_PREFIXES


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--schema-dir",       required=True,
                   help="Directory of locus FASTAs (defines locus column universe)")
    p.add_argument("--contig-alignment", required=True,
                   help="taranys contig_alignment_info.csv")
    p.add_argument("--match-csv",        required=True,
                   help="taranys allele_calling_match.csv")
    p.add_argument("--sample",           required=True,
                   help="Sample identifier (must match Sample column in match-csv)")
    p.add_argument("--output",           required=True)
    args = p.parse_args()

    loci = list_schema_loci(args.schema_dir)

    # contig_alignment_info.csv: load sample_allele_seq per locus.
    # Paralogs (NIPH) have multiple rows; they'll be overridden to "-"
    # via the match-csv code anyway, so last-write-wins is harmless.
    seq_by_locus = {}
    with open(args.contig_alignment, newline="") as fh:
        reader = csv.DictReader(fh)
        for col in ("locus_name", "sample_allele_seq"):
            if col not in reader.fieldnames:
                print(f"ERROR: column '{col}' missing from {args.contig_alignment}; "
                      f"got {reader.fieldnames}", file=sys.stderr)
                return 1
        for row in reader:
            seq_by_locus[row["locus_name"]] = row["sample_allele_seq"]

    # allele_calling_match.csv: pull codes for our sample.
    with open(args.match_csv, newline="") as fh:
        reader = csv.DictReader(fh)
        codes = None
        sample_col = reader.fieldnames[0]  # "Sample"
        for row in reader:
            if row[sample_col] == args.sample:
                codes = {k: row[k] for k in reader.fieldnames if k != sample_col}
                break
        if codes is None:
            # Fall back: if there's only one sample in the file, use it.
            fh.seek(0)
            reader = csv.DictReader(fh)
            rows = list(reader)
            if len(rows) == 1:
                only = rows[0][sample_col]
                print(f"NOTE: '{args.sample}' not found in {args.match_csv}; "
                      f"using single available row '{only}'", file=sys.stderr)
                codes = {k: rows[0][k] for k in reader.fieldnames if k != sample_col}
            else:
                print(f"ERROR: sample '{args.sample}' not in {args.match_csv} "
                      f"(found {len(rows)} samples); cannot proceed",
                      file=sys.stderr)
                return 1

    # Build the row in schema-defined locus order.
    row_values = []
    for locus in loci:
        code = codes.get(locus, "")
        if is_missing(code):
            row_values.append(MISSING_VALUE)
        else:
            row_values.append(get_hash(seq_by_locus.get(locus)))

    # Write TSV (chewBBACA-style: FILE column + locus columns).
    with open(args.output, "w") as out:
        out.write("FILE\t" + "\t".join(loci) + "\n")
        out.write(args.sample + "\t" + "\t".join(row_values) + "\n")

    n_called = sum(1 for v in row_values if v != MISSING_VALUE)
    print(f"Wrote {n_called}/{len(loci)} called loci for {args.sample} "
          f"to {args.output}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())