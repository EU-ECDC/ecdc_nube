#!/usr/bin/env python

import subprocess as sp
import argparse
import sys
import os
import logging
from tqdm import tqdm
from Bio import SeqIO
import glob

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S"
)

parser = argparse.ArgumentParser()
parser.add_argument("-a", "--alignments_dir", type=str,
                    help="path to the directory with the existing alignments", required=True)
parser.add_argument("-d", "--distmat", type=str, 
                    help="directory with current distance matrices")
parser.add_argument("-o", "--outdir", type=str,
                    help="path to the directory with the output distance matrices", required=True)
args = parser.parse_args()

def hamming_distance_ignore_gaps(seq1, seq2):
    mismatches = 0
    comparable_positions = 0
    for c1, c2 in zip(seq1, seq2):
        if c1 != '-' and c2 != '-':
            comparable_positions += 1
            if c1 != c2:
                mismatches += 1
    proportion_mismatches = mismatches / comparable_positions if comparable_positions > 0 else 0
    return mismatches, comparable_positions, proportion_mismatches

def find_sequence_offsets(fasta_file):
    offsets = []
    with open(fasta_file, "rb") as f:
        while True:
            pos = f.tell()
            line = f.readline()
            if not line:
                break
            if line.startswith(b">"):
                offsets.append(pos)
    return offsets

def get_sequence_ids(fasta_file):
    return [record.id for record in SeqIO.parse(fasta_file, "fasta")]

def stream_incremental_distances(fasta_file, output_file, new_ids):
    print(f"Finding sequence positions in {fasta_file}...")
    sequence_offsets = find_sequence_offsets(fasta_file)
    n_seq = len(sequence_offsets)
    print(f"Found {n_seq} sequences.")

    all_ids = []
    new_ids_offsets=[]

    with open(fasta_file, "r") as f:
        for offset in sequence_offsets:
            f.seek(offset)
            seq = next(SeqIO.parse(f, "fasta"))
            all_ids.append((offset, seq.id))
            if seq.id in new_ids:
                new_ids_offsets.append((offset, seq.id))

    if not new_ids_offsets:
        print("No new sequences found.")
        return

    with open(fasta_file, "r") as f_in, open(output_file, "a") as f_out:
        print(f"Computing distances for {len(new_ids_offsets)} new sequences and appending to file...")
        processed_pairs=set()
        for offset_i, id_i in tqdm(new_ids_offsets):
            f_in.seek(offset_i)
            seq1 = next(SeqIO.parse(f_in, "fasta"))
            seq1_str = str(seq1.seq).upper()

            for offset_j, id_j in all_ids:
                if id_i == id_j:
                    continue

                pair = tuple(sorted((id_i, id_j)))
                if pair in processed_pairs:
                    continue
                else:
                    processed_pairs.add(pair)

                f_in.seek(offset_j)
                seq2 = next(SeqIO.parse(f_in, "fasta"))
                seq2_str = str(seq2.seq).upper()

                mm, c_pos, pdist,  = hamming_distance_ignore_gaps(seq1_str, seq2_str)
                f_out.write(f"{id_i}\t{id_j}\t{mm}\t{c_pos}\t{pdist}\n")


    print("Update complete.")


alignments_to_consider = glob.glob(os.path.join(args.alignments_dir, '**', '*_aligned.fasta'), recursive=True)

os.makedirs(f"{args.outdir}", exist_ok=True)

for ali in alignments_to_consider:    
    # I get the txt file with the new sequences
    bname = os.path.basename(ali).replace("_aligned.fasta", "")
    dir_ali = os.path.dirname(ali)

    ids_new_seqs = []
    with open(f"{dir_ali}/inp_seqs_matching_{bname}.txt", "r") as fin:
        ids_new_seqs = [line.strip() for line in fin]


    # I need to remove the rows where there are the new ids
    if args.distmat and os.path.exists(f"{args.distmat}/{bname}_dist.tsv"):        
        with open(f"{args.distmat}/{bname}_dist.tsv", "r") as fin, open(f"{args.outdir}/{bname}_dist.tsv", "w") as fout:
            for line in fin:
                fields = line.strip().split('\t')
                if not any(id_ in ids_new_seqs for id_ in fields):
                    fout.write(line)
    else:
        with open(f"{args.outdir}/{bname}_dist.tsv", "w") as fout:
            pass

    stream_incremental_distances(ali, f"{args.outdir}/{bname}_dist.tsv", ids_new_seqs)



