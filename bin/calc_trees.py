#!/usr/bin/env python

import subprocess as sp
import argparse
import sys
import os
import logging
import glob


logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S"
)

parser = argparse.ArgumentParser()
parser.add_argument("-a", "--alignments_dir", type=str,
                    help="path to the directory with the existing alignments", required=True)
parser.add_argument("-o", "--outdir", type=str,
                    help="path to the directory with the output alignments", required=True)
args = parser.parse_args()

alignments_to_consider = glob.glob(os.path.join(args.alignments_dir, '**', '*_aligned.fasta'), recursive=True)

os.makedirs(f"{args.outdir}", exist_ok=True)

for ali in alignments_to_consider:    
    logging.info(f"Generating tree for {ali}")
    bname = os.path.basename(ali).replace("_aligned.fasta", "")
    cmd_alignment = ["iqtree",
        "-m", "GTR",  
        "-s", ali,
        "--prefix", f"{args.outdir}/{bname}"
        ]
    try:
        cmd_out = sp.run(cmd_alignment, check=True, capture_output=True, text=True)
    except sp.CalledProcessError as e:
        if e.returncode != 0:
            logging.error(f"iqtree failed with exit code {e.returncode}")
            #print(f"[STDOUT]\n{e.stdout}")
            print(f"[STDERR]\n{e.stderr}")
 