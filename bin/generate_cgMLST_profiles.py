#!/usr/bin/env python3

import argparse
import random, string
import shutil
import subprocess as sp
import sys
import glob
import re

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input_fasta", type=str,
                    help="Input fasta", required=True)
parser.add_argument("-a", "--accession", type=str,
                    help="Accession", required=True)
parser.add_argument("-g", "--schema_directory", type=str,
                    help="Path to the schema directory", required=True)
parser.add_argument("-p", "--ptf", type=str,
                    help="Path to the Prodigal training file", required=True)
parser.add_argument("-gl", "--gene_list", type=str,
                    help="Path to a file with the list of genes/loci to perform allele calling", required=True)
parser.add_argument("-ao", "--adv_options", type=str,
                    help="Advanced options", required=True)
parser.add_argument("-f", "--prefix", type=str,
                    help="Output file prefix", required=False)
args = parser.parse_args()

random_tag = ''.join(random.choices(string.ascii_letters + string.digits, k=8))
print(f"- I generate a random tag: {random_tag}")

print(f"- making a copy of {args.input_fasta} [{args.input_fasta} -> {random_tag}.fasta]")
shutil.copy(args.input_fasta, f"{random_tag}.fasta")

with open("assembly_to_process.txt","w") as inpf:
    inpf.write(f"{random_tag}.fasta\n")

cmd_chewbbaca = [
    "chewBBACA.py",
    "AlleleCall",
    "-i", "assembly_to_process.txt",
    "-g", args.schema_directory,
    "-o", ".",
    "--no-inferred",
    "--ptf", args.ptf,
    "--genes-list", args.gene_list,
    "--hash-profiles", "crc32"
    ]
if args.adv_options != "None":
    cmd_chewbbaca.extend(args.adv_options.split())

print("- running chewBBACA")
print(" ".join(cmd_chewbbaca))
sp.run(cmd_chewbbaca,text=True, stdout=sys.stdout, stderr=sys.stderr)

files_hashed_profiles = glob.glob("results_*/results_alleles_hashed.tsv")

if files_hashed_profiles and len(files_hashed_profiles) == 1:
    target_file = glob.glob("results_*/results_alleles_hashed.tsv")[0]
    with open(target_file,"r") as tf:
        if args.prefix:
            output_prefix = args.prefix
        else:
            output_prefix =  f"{args.accession}_cgMLST"
        with open(f"{output_prefix}.tsv","w") as outf:
            for line in tf:
                if line.startswith(f"{random_tag}\t"):
                    outf.write(re.sub(rf"^{random_tag}", f"{args.accession}", line))
                else:
                    outf.write(line)
else:
    sys.exit(80)