#!/usr/bin/env python

import subprocess as sp
import argparse
import sys
import os
import logging
import glob
import io

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S"
)

parser = argparse.ArgumentParser()
parser.add_argument("-in", "--inp", type=str, 
                    help="fasta file with the sequences to add to the alignment", required=True)
parser.add_argument("-r", "--ref", type=str,
                    help="path to fasta file of reference sequences", required=True)
parser.add_argument("-a", "--alignments_dir", type=str,
                    help="path to the directory with the existing alignments", required=True)
parser.add_argument("-o", "--outdir", type=str,
                    help="path to the directory with the output alignments", required=True)
parser.add_argument("-f", "--flags", type=str,
                    help="Flags [force-align, regenerate]")
args = parser.parse_args()

list_flags = []
if args.flags:
    list_flags = args.flags.rstrip("\n").split(",")


def extract_ids_from_fasta(fasta_path):
    ids = []
    with io.open(fasta_path, 'r', buffering=1024*1024) as file:  # 1MB buffer
        for line in file:
            if line.startswith('>'):
                ids.append(line[1:].strip().split()[0])
    return ids

def extract_sequences_by_ids(fasta_path, id_list, output_path):
    wanted = set(id_list)  # For O(1) lookups
    write = False

    with open(fasta_path, 'r') as infile, open(output_path, 'w') as outfile:
        for line in infile:
            if line.startswith('>'):
                seq_id = line[1:].strip().split()[0]  # Take ID before first space
                write = seq_id in wanted
            if write:
                outfile.write(line)

# Creates a db for the references if it does not exist
rname = os.path.basename(args.ref).split(".")[0]
name_db = rname + "_db"
if not os.path.exists("/".join([args.alignments_dir, name_db + ".nsq"])):
    logging.info("Creating blast db")
    cmd_makeblastdb = ["makeblastdb", 
        "-in", args.ref,
        "-dbtype", "nucl",
        "-out", "/".join([args.alignments_dir, name_db])
        ]
    print(" ".join(cmd_makeblastdb))
    sp.run(cmd_makeblastdb, check=True, capture_output=False, text=True)

# Running blast to find the best reference
logging.info("Running blast")
cmd_blast = ["blastn", 
        "-query", args.inp,
        "-db", "/".join([args.alignments_dir, name_db]),
        "-out", "results_blast.tsv",
        "-outfmt", "6",
        "-max_target_seqs", "1"
        ] 
print(" ".join(cmd_blast))
sp.run(cmd_blast, check=True, capture_output=False, text=True)

# Getting the list of the ids of the reference sequences
#logging.info("Getting the ids of the reference sequences")
#cdata = sp.run(["seqkit", 
#        "seq","-n", args.ref
#        ], check=True, capture_output=True, text=True)
#l_ids = cdata.stdout.strip().split("\n")
#print(l_ids)

# I split the input sequences based on the best hit in the 
# pool of reference sequences
logging.info("Analyzing blast results")
d_split_inp_seqs = {}
with open("results_blast.tsv","r") as inp:
    for row in inp:
        explode_row = row.rstrip().split("\t")
        id_ref_seq = explode_row[1]
        id_query_seq = explode_row[0]
        if id_ref_seq not in d_split_inp_seqs:
            d_split_inp_seqs[id_ref_seq] = []    
        d_split_inp_seqs[id_ref_seq].append(id_query_seq)
print(d_split_inp_seqs)

for ref_seq_id in d_split_inp_seqs.keys():
    logging.info(f"Filtering input sequences based on the best-matching reference (ref: {ref_seq_id})")
    os.makedirs(f"{args.outdir}/", exist_ok=True)
    
    ids_sequences_inp_fasta_matching_ref = d_split_inp_seqs[ref_seq_id]
    
    with open(f"{args.outdir}/inp_seqs_matching_{ref_seq_id}.txt","w") as tout:
        tout.write("\n".join(ids_sequences_inp_fasta_matching_ref))

    extract_sequences_by_ids(args.inp, 
                             ids_sequences_inp_fasta_matching_ref,
                             f"{args.outdir}/inp_seqs_matching_{ref_seq_id}.fasta")
    # if an alignment exists
    if os.path.exists(f"{args.alignments_dir}/{ref_seq_id}_aligned.fasta") and ("regenerate" not in list_flags):
        if "force-align" in list_flags:
            # Removing sequences from existing alignment (due to force-align option)
            seqs_to_remove = d_split_inp_seqs[ref_seq_id]
            seqs_in_existing_alignment = extract_ids_from_fasta(f"{args.alignments_dir}/{ref_seq_id}_aligned.fasta")

            seqs_to_extract = list(set(seqs_in_existing_alignment) - set(seqs_to_remove))
            if len(seqs_to_extract) == 0:
                raise Exception("ERROR: No sequences in the existing alignment")
            
            extract_sequences_by_ids(f"{args.alignments_dir}/{ref_seq_id}_aligned.fasta", 
                             seqs_to_extract,
                             f"{args.outdir}/base_ali_seqs_{ref_seq_id}.fasta")
            
            # I perform the alignment
            logging.info(f"Generating the sequence alignment")
            cmd_alignment = ["augur", "align", 
                "--sequences", f"{args.outdir}/inp_seqs_matching_{ref_seq_id}.fasta",
                "--existing-alignment", f"{args.outdir}/base_ali_seqs_{ref_seq_id}.fasta",
                "--debug",
                "--output", f"{args.outdir}/{ref_seq_id}_aligned.fasta"
                ]
            sp.run(cmd_alignment, check=True, capture_output=False, text=True)
            
        else:
            # I perform the alignment
            logging.info(f"Generating the sequence alignment")
            cmd_alignment = ["augur", "align", 
                "--sequences", f"{args.outdir}/inp_seqs_matching_{ref_seq_id}.fasta",
                "--existing-alignment", f"{args.alignments_dir}/{ref_seq_id}_aligned.fasta",
                "--debug",
                "--output", f"{args.outdir}/{ref_seq_id}_aligned.fasta"
                ]
            print(" ".join(cmd_alignment))
            sp.run(cmd_alignment, check=True, capture_output=False, text=True)

    # if there is no alignment or there is the flag "regenerate": 
    else:
        # The reference sequence constitutes the existing alignment
        extract_sequences_by_ids(args.ref, 
                                 [ref_seq_id],
                                 f"{args.outdir}/{ref_seq_id}.fasta")
        
        # I perform the alignment
        logging.info(f"Generating the sequence alignment")
        cmd_alignment = ["augur", "align", 
            "--sequences", f"{args.outdir}/inp_seqs_matching_{ref_seq_id}.fasta",
            "--existing-alignment", f"{args.outdir}/{ref_seq_id}.fasta",
            "--debug",
            "--output", f"{args.outdir}/{ref_seq_id}_aligned.fasta"
            ] 
        print(" ".join(cmd_alignment))
        sp.run(cmd_alignment, check=True, capture_output=False, text=True)
